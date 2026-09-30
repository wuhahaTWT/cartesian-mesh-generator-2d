#include "cartmesh2d/fv/ReactingImplicit2D.hpp"
#include "cartmesh2d/chemistry/SpeciesMassClosure.hpp"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_ls.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_config.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <utility>

#if SUNDIALS_VERSION_MAJOR != 5
#error "Coupled implicit development backend currently requires the verified SUNDIALS 5 API"
#endif

namespace cartmesh2d::fv {
namespace {
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error("implicit reacting flow: " + message);
}
void check(int flag, const char* name) {
    require(flag >= 0, std::string(name) + " failed with flag " + std::to_string(flag));
}
struct EvaluationBudgetExceeded : std::runtime_error {
    EvaluationBudgetExceeded() : std::runtime_error("residual-evaluation budget exhausted before endpoint") {}
};
struct EvaluationCanceled : std::runtime_error {
    EvaluationCanceled() : std::runtime_error("canceled during a trial solve; restart from the last accepted state with fresh BDF history") {}
};
struct SundialsResources {
    void* ode = nullptr;
    N_Vector y = nullptr, absolute = nullptr, quadrature = nullptr, quadAbsolute = nullptr;
    N_Vector localError = nullptr, errorWeights = nullptr;
    SUNMatrix matrix = nullptr;
    SUNLinearSolver linear = nullptr;
    SUNNonlinearSolver nonlinear = nullptr;
    ~SundialsResources() {
        if (ode) CVodeFree(&ode);
        if (nonlinear) SUNNonlinSolFree(nonlinear);
        if (linear) SUNLinSolFree(linear);
        if (matrix) SUNMatDestroy(matrix);
        for (auto v : {y, absolute, quadrature, quadAbsolute, localError, errorWeights}) if (v) N_VDestroy(v);
    }
};
std::vector<std::vector<std::size_t>> adjacency(const FvMesh2D& mesh) {
    std::vector<std::vector<std::size_t>> graph(mesh.cells.size());
    for (const auto& f : mesh.faces) if (f.neighbour) {
        graph[f.owner].push_back(*f.neighbour); graph[*f.neighbour].push_back(f.owner);
    }
    for (auto& row : graph) { std::sort(row.begin(), row.end()); row.erase(std::unique(row.begin(), row.end()), row.end()); }
    return graph;
}
std::vector<std::size_t> reverseCuthillMcKee(const std::vector<std::vector<std::size_t>>& graph) {
    std::vector<std::size_t> permutation;
    std::vector<bool> seen(graph.size());
    const auto less = [&](std::size_t a, std::size_t b) {
        return std::pair{graph[a].size(), a} < std::pair{graph[b].size(), b};
    };
    while (permutation.size() < graph.size()) {
        auto first = graph.size();
        for (std::size_t i = 0; i < graph.size(); ++i) if (!seen[i] && (first == graph.size() || less(i, first))) first = i;
        std::queue<std::size_t> todo; todo.push(first); seen[first] = true;
        while (!todo.empty()) {
            const auto i = todo.front(); todo.pop(); permutation.push_back(i);
            auto neighbours = graph[i]; std::sort(neighbours.begin(), neighbours.end(), less);
            for (auto j : neighbours) if (!seen[j]) { seen[j] = true; todo.push(j); }
        }
    }
    std::reverse(permutation.begin(), permutation.end());
    return permutation;
}
}

struct ReactingImplicitIntegrator2D::Impl {
    // Decorate the public SUNDIALS Newton solver interface. Whole Newton
    // updates are backtracked to an admissible iterate; no component is
    // projected/clipped, and the original correction controls convergence.
    struct GuardedNewton {
        Impl& owner;
        SUNNonlinearSolver inner = nullptr;
        SUNNonlinSolSysFn system = nullptr;
        void* integrator = nullptr;
        N_Vector prediction = nullptr, base = nullptr, candidate = nullptr, full = nullptr, residual = nullptr;
        static inline thread_local GuardedNewton* active = nullptr;
        explicit GuardedNewton(Impl& i) : owner(i) {}
        ~GuardedNewton() {
            if (inner) SUNNonlinSolFree(inner);
            for (auto v : {base, candidate, full, residual}) if (v) N_VDestroy(v);
        }
        static GuardedNewton& get(SUNNonlinearSolver n) { return *static_cast<GuardedNewton*>(n->content); }
        static int systemCallback(N_Vector correction, N_Vector result, void* memory) {
            if (!active || active->integrator != memory) return SUN_NLS_ILL_INPUT;
            const int status = active->system(correction, result, memory);
            if (status == 0) N_VScale(1, correction, active->base);
            return status;
        }
        bool admissible(N_Vector correction) {
            N_VLinearSum(1, prediction, 1, correction, full);
            try {
                owner.decode(owner.progress.lastAccepted.time, NV_DATA_S(full), owner.scratch);
                owner.solver.validate(owner.scratch);
                return true;
            } catch (const std::exception& error) { owner.lastRhsError = error.what(); return false; }
        }
        static int test(SUNNonlinearSolver, N_Vector correction, N_Vector delta,
                        realtype tolerance, N_Vector weights, void* context) {
            auto& g = *static_cast<GuardedNewton*>(context);
            try {
                double fraction = 1;
                if (!g.admissible(correction)) {
                    g.owner.progress.lastDampedTrialFailure = g.owner.lastRhsError;
                    bool found = false;
                    for (unsigned attempt = 0; attempt < 64; ++attempt) {
                        fraction *= .5;
                        N_VLinearSum(1, g.base, fraction, delta, g.candidate);
                        if (g.admissible(g.candidate)) { found = true; break; }
                    }
                    // A feasible base iterate may already solve the nonlinear
                    // system within its unchanged local tolerance. Accepting
                    // that iterate requires BOTH the full Newton correction
                    // and the actual nonlinear residual to meet that tolerance.
                    if (!found) {
                        fraction = 0; N_VScale(1, g.base, g.candidate);
                        if (!g.admissible(g.candidate)) return SUN_NLS_CONV_RECVR;
                    }
                    N_VScale(1, g.candidate, correction);
                    ++g.owner.progress.dampedNewtonUpdates;
                    g.owner.progress.minimumNewtonFraction = std::min(g.owner.progress.minimumNewtonFraction, fraction);
                }
                // Evaluate the final iterate even when the stock correction
                // test would stop without evaluating it. This also synchronizes
                // CVODES' current state before the quadrature callback.
                const int status = g.system(correction, g.residual, g.integrator);
                if (status != 0) return status;
                auto& p = g.owner.progress;
                p.lastNewtonCorrectionNorm = N_VWrmsNorm(delta, weights);
                p.lastNewtonResidualNorm = N_VWrmsNorm(g.residual, weights);
                p.lastNewtonTolerance = tolerance;
                p.lastPredictorCorrectionNorm = N_VWrmsNorm(correction, weights);
                // The stock CVODES convergence test caches |delta| as |acor|
                // on the first Newton iteration. Our feasible initial guess
                // and backtracking invalidate delta == acor. Never invoke
                // that shortcut: leaving CVODES' acnrmcur flag false makes
                // it compute the actual total predictor correction for LTE.
                // Both the original Newton update and the actual residual
                // must meet the supplied nonlinear tolerance.
                if (p.lastNewtonCorrectionNorm <= tolerance && p.lastNewtonResidualNorm <= tolerance)
                    return SUN_NLS_SUCCESS;
                if (fraction < 1) {
                    if (!g.owner.controls.continueDampedNewton || fraction == 0) return SUN_NLS_CONV_RECVR;
                    ++p.continuedDampedNewtonUpdates;
                }
                return SUN_NLS_CONTINUE;
            } catch (const std::exception& error) { g.owner.lastRhsError = error.what(); return SUN_NLS_EXT_FAIL; }
            catch (...) { g.owner.lastRhsError = "unknown nonlinear convergence exception"; return SUN_NLS_EXT_FAIL; }
        }
        static SUNNonlinearSolver create(Impl& owner, N_Vector model) {
            std::unique_ptr<GuardedNewton> g = std::make_unique<GuardedNewton>(owner);
            g->inner = SUNNonlinSol_Newton(model);
            g->base = N_VClone(model); g->candidate = N_VClone(model); g->full = N_VClone(model); g->residual = N_VClone(model);
            require(g->inner && g->base && g->candidate && g->full && g->residual, "nonlinear solver allocation failed");
            auto n = SUNNonlinSolNewEmpty(); require(n, "nonlinear interface allocation failed");
            n->content = g.release();
            n->ops->gettype = [](auto) { return SUNNONLINEARSOLVER_ROOTFIND; };
            n->ops->initialize = [](auto s) { return SUNNonlinSolInitialize(get(s).inner); };
            n->ops->free = [](auto s) { delete &get(s); SUNNonlinSolFreeEmpty(s); return SUN_NLS_SUCCESS; };
            n->ops->setsysfn = [](auto s, SUNNonlinSolSysFn f) {
                get(s).system = f; return SUNNonlinSolSetSysFn(get(s).inner, systemCallback);
            };
            n->ops->setlsetupfn = [](auto s, SUNNonlinSolLSetupFn f) { return SUNNonlinSolSetLSetupFn(get(s).inner, f); };
            n->ops->setlsolvefn = [](auto s, SUNNonlinSolLSolveFn f) { return SUNNonlinSolSetLSolveFn(get(s).inner, f); };
            n->ops->setctestfn = [](auto s, SUNNonlinSolConvTestFn, void*) {
                auto& a = get(s);
                return SUNNonlinSolSetConvTestFn(a.inner, test, &a);
            };
            n->ops->setmaxiters = [](auto s, int count) { return SUNNonlinSolSetMaxIters(get(s).inner, count); };
            n->ops->getnumiters = [](auto s, long* count) { return SUNNonlinSolGetNumIters(get(s).inner, count); };
            n->ops->getcuriter = [](auto s, int* count) { return SUNNonlinSolGetCurIter(get(s).inner, count); };
            n->ops->getnumconvfails = [](auto s, long* count) { return SUNNonlinSolGetNumConvFails(get(s).inner, count); };
            n->ops->solve = [](auto s, N_Vector prediction, N_Vector correction, N_Vector weights,
                               realtype tolerance, booleantype setup, void* memory) {
                auto& a = get(s); a.prediction = prediction; a.integrator = memory;
                // The C system callback has no solver-context argument. Scope
                // its adapter per thread and restore nested calls on exit.
                struct ActiveScope {
                    GuardedNewton* previous = active;
                    explicit ActiveScope(GuardedNewton& value) { active = &value; }
                    ~ActiveScope() { active = previous; }
                } scope(a);
                try {
                    if (!a.admissible(correction)) {
                        // Choose a feasible initial Newton iterate when a BDF
                        // predictor crosses a bound. This changes the nonlinear
                        // starting guess only, never the accepted solution.
                        for (std::size_t i = 0; i < a.owner.mapping.size(); ++i) {
                            const auto [cell, component] = a.owner.mapping[i];
                            NV_Ith_S(correction, i) = a.owner.progress.lastAccepted.cells[cell][component]
                                / a.owner.scale[i] - NV_Ith_S(prediction, i);
                        }
                        if (!a.admissible(correction)) return SUN_NLS_CONV_RECVR;
                    }
                    return SUNNonlinSolSolve(a.inner, prediction, correction, weights, tolerance, setup, memory);
                } catch (const std::exception& error) { a.owner.lastRhsError = error.what(); return SUN_NLS_EXT_FAIL; }
                catch (...) { a.owner.lastRhsError = "unknown nonlinear solve exception"; return SUN_NLS_EXT_FAIL; }
            };
            return n;
        }
    };
    chemistry::DetailedGas& gas;
    ReactingFlowStepper2D& solver;
    ReactingImplicitControls2D controls;
    ReactingState2D initial, scratch;
    ReactingImplicitProgress2D progress;
    SundialsResources resources;
    std::vector<std::pair<std::size_t,std::size_t>> mapping;
    std::vector<std::size_t> dependent;
    std::vector<double> scale, integralScale, cachedY;
    ReactingResidual2D cachedResidual;
    double cachedTime = -1;
    unsigned cachedOrder = 0;
    std::string lastRhsError, integratorError;
    const std::function<void(const ReactingImplicitProgress2D&, double)>* evaluationObserver = nullptr;
    const std::function<bool()>* evaluationCancel = nullptr;
    bool terminalFailure = false;
    std::size_t nv = 0;

    Impl(chemistry::DetailedGas& chemistry, ReactingFlowStepper2D& residual,
         const ReactingState2D& input, const ReactingImplicitControls2D& c)
        : gas(chemistry), solver(residual), controls(c), initial(input), scratch(input) {
        require(&gas == &solver.gas_, "different chemistry contexts for residual and integration");
        solver.validate(input);
        require(std::isfinite(c.relativeTolerance) && c.relativeTolerance > 0 && c.relativeTolerance < 1
            && std::isfinite(c.absoluteConservedTolerance) && c.absoluteConservedTolerance > 0
            && std::isfinite(c.absoluteSpeciesFraction) && c.absoluteSpeciesFraction > 0, "invalid local error weights");
        require(std::isfinite(c.maximumStep) && c.maximumStep > 0 && std::isfinite(c.initialStep) && c.initialStep >= 0
            && c.initialStep <= c.maximumStep && c.maximumAcceptedSteps > 0 && c.maximumResidualEvaluations > 0, "invalid step controls");
        require((c.spatialOrder == 1 || c.spatialOrder == 2) && c.maximumBdfOrder >= 1 && c.maximumBdfOrder <= 5,
            "invalid spatial/BDF order");
        require(c.jacobianAdvectionOrder == 1 || c.jacobianAdvectionOrder == 2, "invalid Newton advection linearization order");
        require(c.maximumNonlinearIterations > 0 && c.maximumNonlinearIterations <= static_cast<unsigned>(std::numeric_limits<int>::max()),
            "invalid nonlinear iteration budget");
        const auto& mesh = solver.mesh(); const auto nc = mesh.cells.size();
        nv = gas.mechanism().species.size() + 4;
        const auto independent = nv - 1;
        const auto graph = adjacency(mesh);
        const auto permutation = reverseCuthillMcKee(graph);
        std::vector<std::size_t> position(nc);
        for (std::size_t i = 0; i < nc; ++i) position[permutation[i]] = i;
        std::size_t cellBand = 0;
        // MUSCL and nonorthogonal viscous/diffusive gradients couple at most
        // neighbours of neighbours. Chemistry remains cell-local.
        for (std::size_t i = 0; i < nc; ++i) {
            const auto include = [&](std::size_t j) { cellBand = std::max(cellBand,
                position[i] > position[j] ? position[i] - position[j] : position[j] - position[i]); };
            for (auto j : graph[i]) { include(j); for (auto k : graph[j]) include(k); }
        }
        require(nc <= static_cast<std::size_t>(std::numeric_limits<sunindextype>::max()) / independent, "equation count overflow");
        const auto neq = nc * independent;
        progress.bandHalfWidth = std::min(neq - 1, (cellBand + 1) * independent - 1);
        const auto storedWidth = std::min(neq - 1, 2 * progress.bandHalfWidth) + progress.bandHalfWidth + 1;
        require(storedWidth <= std::numeric_limits<std::size_t>::max() / sizeof(double) / neq, "band storage overflow");
        progress.bandBytes = storedWidth * neq * sizeof(double);
        require(progress.bandBytes <= c.maximumBandBytes, "direct band matrix exceeds configured storage budget; no reduced-physics fallback");
        dependent.resize(nc); integralScale.assign(nv, 0);
        std::vector<double> initialConstraint(nv);
        for (auto i : permutation) {
            auto represented = input.cells[i];
            const auto closure = chemistry::closeSpeciesMassRoundoff(represented[0], represented, 4);
            dependent[i] = closure.species + 4;
            initialConstraint[dependent[i]] += closure.change * mesh.cells[i].area;
            const auto p = reactingPrimitive2D(gas, represented);
            const double density = represented[0], energy = std::max(std::abs(represented[3]), density * p.properties.cv * p.properties.temperature);
            for (std::size_t k = 0; k < nv; ++k) {
                const double unit = k == 3 ? energy : (k == 1 || k == 2) ? density * p.soundSpeed : density;
                integralScale[k] += unit * mesh.cells[i].area;
                if (k != dependent[i]) { mapping.emplace_back(i, k); scale.push_back(unit); }
            }
        }
        require(mapping.size() == neq, "invalid independent-state layout");
        auto& r = resources;
        r.y = N_VNew_Serial(static_cast<sunindextype>(neq)); r.absolute = N_VNew_Serial(static_cast<sunindextype>(neq));
        require(r.y && r.absolute, "state vector allocation failed");
        r.localError = N_VClone(r.y); r.errorWeights = N_VClone(r.y);
        r.quadrature = N_VNew_Serial(static_cast<sunindextype>(3 * nv)); r.quadAbsolute = N_VNew_Serial(static_cast<sunindextype>(3 * nv));
        require(r.y && r.absolute && r.quadrature && r.quadAbsolute && r.localError && r.errorWeights, "vector allocation failed");
        for (std::size_t a = 0; a < neq; ++a) {
            const auto [i, k] = mapping[a];
            NV_Ith_S(r.y, a) = input.cells[i][k] / scale[a];
            NV_Ith_S(r.absolute, a) = k >= 4 ? c.absoluteSpeciesFraction : c.absoluteConservedTolerance;
        }
        for (std::size_t a = 0; a < 3 * nv; ++a) {
            NV_Ith_S(r.quadrature, a) = a >= 2 * nv ? initialConstraint[a % nv] / integralScale[a % nv] : 0;
            NV_Ith_S(r.quadAbsolute, a) = c.absoluteConservedTolerance;
        }
        r.ode = CVodeCreate(CV_BDF); require(r.ode, "CVODES allocation failed");
        check(CVodeSetUserData(r.ode, this), "CVodeSetUserData");
        check(CVodeSetErrHandlerFn(r.ode, errorCallback, this), "CVodeSetErrHandlerFn");
        check(CVodeInit(r.ode, rhsCallback, input.time, r.y), "CVodeInit");
        check(CVodeSVtolerances(r.ode, c.relativeTolerance, r.absolute), "CVodeSVtolerances");
        // Do not enable CVODES inequality constraints: cvCheckConstraints can
        // project small negative components to zero. Physical admissibility
        // is checked on unmodified trial and accepted states instead.
        check(CVodeSetMaxOrd(r.ode, static_cast<int>(c.maximumBdfOrder)), "CVodeSetMaxOrd");
        check(CVodeSetMaxStep(r.ode, c.maximumStep), "CVodeSetMaxStep");
        if (c.initialStep > 0) check(CVodeSetInitStep(r.ode, c.initialStep), "CVodeSetInitStep");
        r.matrix = SUNBandMatrix(static_cast<sunindextype>(neq), static_cast<sunindextype>(progress.bandHalfWidth),
            static_cast<sunindextype>(progress.bandHalfWidth));
        require(r.matrix, "band allocation failed"); r.linear = SUNLinSol_Band(r.y, r.matrix); require(r.linear, "band solver allocation failed");
        check(CVodeSetLinearSolver(r.ode, r.linear, r.matrix), "CVodeSetLinearSolver");
        check(CVodeSetJacFn(r.ode, jacobianCallback), "CVodeSetJacFn");
        r.nonlinear = GuardedNewton::create(*this, r.y);
        check(CVodeSetNonlinearSolver(r.ode, r.nonlinear), "CVodeSetNonlinearSolver");
        check(CVodeSetMaxNonlinIters(r.ode, static_cast<int>(c.maximumNonlinearIterations)), "CVodeSetMaxNonlinIters");
        check(CVodeQuadInit(r.ode, quadCallback, r.quadrature), "CVodeQuadInit");
        check(CVodeQuadSVtolerances(r.ode, c.relativeTolerance, r.quadAbsolute), "CVodeQuadSVtolerances");
        check(CVodeSetQuadErrCon(r.ode, 1), "CVodeSetQuadErrCon");
        progress.lastAccepted = input;
        progress.boundaryImpulse.assign(nv, 0); progress.chemistryChange.assign(nv, 0); progress.constraintChange.assign(nv, 0);
    }
    void decode(double t, const double* values, ReactingState2D& state) const {
        state.time = t;
        for (std::size_t a = 0; a < mapping.size(); ++a) {
            const auto [i, k] = mapping[a]; state.cells[i][k] = values[a] * scale[a];
        }
        for (std::size_t i = 0; i < state.cells.size(); ++i) {
            // Compensated positive sum; negative trial states remain errors.
            double sum = 0, correction = 0;
            for (std::size_t k = 4; k < nv; ++k) if (k != dependent[i]) {
                const double v = state.cells[i][k];
                if (!(std::isfinite(v) && v >= 0)) {
                    std::ostringstream message; message << std::setprecision(17) << "trial cell " << i
                        << " species " << gas.mechanism().species[k - 4] << " density=" << v;
                    require(false, message.str());
                }
                const double next = sum + v; correction += sum >= v ? (sum - next) + v : (v - next) + sum; sum = next;
            }
            state.cells[i][dependent[i]] = state.cells[i][0] - (sum + correction);
        }
    }
    const ReactingResidual2D& evaluate(double t, N_Vector y, unsigned order) {
        if (evaluationCancel && *evaluationCancel && (*evaluationCancel)()) throw EvaluationCanceled();
        const auto* values = NV_DATA_S(y);
        if (cachedTime == t && cachedOrder == order && cachedY.size() == mapping.size() && std::equal(cachedY.begin(), cachedY.end(), values)) return cachedResidual;
        if (progress.rhsCalls >= controls.maximumResidualEvaluations) throw EvaluationBudgetExceeded();
        ++progress.rhsCalls;
        if (evaluationObserver && *evaluationObserver) (*evaluationObserver)(progress, t);
        decode(t, values, scratch);
        auto residual = solver.evaluateResidual(scratch, order, false);
        cachedY.assign(values, values + mapping.size()); cachedTime = t; cachedOrder = order; cachedResidual = std::move(residual);
        return cachedResidual;
    }
    static int rhsCallback(realtype t, N_Vector y, N_Vector derivative, void* context) {
        return rhsForOrder(t, y, derivative, context, static_cast<Impl*>(context)->controls.spatialOrder);
    }
    static int rhsForOrder(realtype t, N_Vector y, N_Vector derivative, void* context, unsigned order) {
        auto& self = *static_cast<Impl*>(context);
        try {
            const auto& f = self.evaluate(t, y, order);
            for (std::size_t a = 0; a < self.mapping.size(); ++a) {
                const auto [i, k] = self.mapping[a]; NV_Ith_S(derivative, a) = f.derivative[i][k] / self.scale[a];
            }
            return 0;
        } catch (const EvaluationCanceled& error) { self.progress.canceled = true; self.lastRhsError = error.what(); return -1; }
        catch (const EvaluationBudgetExceeded& error) { self.lastRhsError = error.what(); return -1; }
        catch (const std::exception& error) { ++self.progress.rejectedRhsCalls; self.lastRhsError = error.what(); return 1; }
        catch (...) { self.lastRhsError = "unknown RHS exception"; return -1; }
    }
    static int jacobianCallback(realtype t, N_Vector y, N_Vector baseDerivative, SUNMatrix matrix,
                                void* context, N_Vector trial, N_Vector derivative, N_Vector linearizationBase) {
        auto& self = *static_cast<Impl*>(context);
        try {
            check(SUNMatZero(matrix), "SUNMatZero");
            const auto count = self.mapping.size(), band = self.progress.bandHalfWidth;
            if (self.controls.jacobianAdvectionOrder != self.controls.spatialOrder) {
                const int status = rhsForOrder(t, y, linearizationBase, context, self.controls.jacobianAdvectionOrder);
                if (status != 0) return status;
                baseDerivative = linearizationBase;
            }
            const auto groups = std::min(count, 2 * band + 1);
            std::vector<double> increments(count);
            // Difference the dimensionless conservative variables using their
            // physical unit scales. The stock h-dependent increment can fall
            // far below bulk-state roundoff for absent/trace species, yielding
            // enormous spurious Jacobian entries as h decreases. This affects
            // only the Newton matrix; all actual RHS states and error weights
            // retain the full mechanism and unmodified trace concentrations.
            const double relativeIncrement = std::sqrt(std::numeric_limits<double>::epsilon());
            for (std::size_t group = 0; group < groups; ++group) {
                int status = 1;
                for (unsigned attempt = 0; attempt < 24; ++attempt) {
                    N_VScale(1, y, trial);
                    for (std::size_t column = group; column < count; column += groups) {
                        const double value = NV_Ith_S(y, column);
                        NV_Ith_S(trial, column) = value + std::ldexp(relativeIncrement * std::max(1., std::abs(value)), -static_cast<int>(attempt));
                        increments[column] = NV_Ith_S(trial, column) - value;
                        require(increments[column] > 0 && std::isfinite(increments[column]), "Jacobian perturbation is not representable");
                    }
                    status = rhsForOrder(t, trial, derivative, context, self.controls.jacobianAdvectionOrder);
                    if (status <= 0) break;
                }
                if (status != 0) return status;
                for (std::size_t column = group; column < count; column += groups) {
                    const auto first = column > band ? column - band : 0;
                    const auto last = std::min(count - 1, column + band);
                    for (auto row = first; row <= last; ++row) {
                        const double value = (NV_Ith_S(derivative, row) - NV_Ith_S(baseDerivative, row)) / increments[column];
                        require(std::isfinite(value), "nonfinite band Jacobian entry");
                        SM_ELEMENT_B(matrix, row, column) = value;
                    }
                }
            }
            return 0;
        } catch (const std::exception& error) { self.lastRhsError = error.what(); return -1; }
        catch (...) { self.lastRhsError = "unknown Jacobian exception"; return -1; }
    }
    static int quadCallback(realtype t, N_Vector y, N_Vector derivative, void* context) {
        auto& self = *static_cast<Impl*>(context);
        try {
            const auto& f = self.evaluate(t, y, self.controls.spatialOrder); const auto nv = self.nv;
            for (std::size_t k = 0; k < nv; ++k) {
                NV_Ith_S(derivative, k) = f.boundaryFlux[k] / self.integralScale[k];
                NV_Ith_S(derivative, nv + k) = f.chemistryIntegral[k] / self.integralScale[k];
                NV_Ith_S(derivative, 2 * nv + k) = 0;
            }
            for (std::size_t i = 0; i < self.scratch.cells.size(); ++i) {
                long double defect = f.derivative[i][0];
                for (std::size_t k = 4; k < nv; ++k) defect -= f.derivative[i][k];
                const auto k = self.dependent[i];
                NV_Ith_S(derivative, 2 * nv + k) += static_cast<double>(defect) * self.solver.mesh().cells[i].area / self.integralScale[k];
            }
            return 0;
        } catch (const EvaluationCanceled& error) { self.progress.canceled = true; self.lastRhsError = error.what(); return -1; }
        catch (const EvaluationBudgetExceeded& error) { self.lastRhsError = error.what(); return -1; }
        // A failed final-state evaluation after nonlinear convergence must
        // stop explicitly. SUNDIALS 5 resets its retry counter before every
        // quadrature evaluation, so recoverable failure here can otherwise
        // repeat indefinitely without an accepted step.
        catch (const std::exception& error) { ++self.progress.rejectedRhsCalls; self.lastRhsError = error.what(); return -1; }
        catch (...) { self.lastRhsError = "unknown quadrature exception"; return -1; }
    }
    static void errorCallback(int, const char*, const char*, char* message, void* context) {
        static_cast<Impl*>(context)->integratorError = message;
    }
    void statistics() {
        CVodeGetNumErrTestFails(resources.ode, &progress.errorTestFailures);
        CVodeGetNumLinSolvSetups(resources.ode, &progress.linearSetups);
        CVodeGetNumJacEvals(resources.ode, &progress.jacobianEvaluations);
        CVodeGetNumNonlinSolvIters(resources.ode, &progress.nonlinearIterations);
        CVodeGetNumNonlinSolvConvFails(resources.ode, &progress.nonlinearConvergenceFailures);
    }
    ReactingImplicitProgress2D advance(double endTime, const std::function<bool()>& cancel,
        const std::function<void(const ReactingImplicitProgress2D&)>& onAccepted,
        const std::function<void(const ReactingImplicitProgress2D&, double)>& onEvaluation) {
        progress.reachedEnd = false; progress.canceled = false;
        if (terminalFailure) return progress;
        evaluationObserver = &onEvaluation;
        evaluationCancel = &cancel;
        struct ObserverScope {
            decltype(evaluationObserver)& observer;
            decltype(evaluationCancel)& cancel;
            ~ObserverScope() { observer = nullptr; cancel = nullptr; }
        } observerScope{evaluationObserver, evaluationCancel};
        try {
            require(std::isfinite(endTime) && endTime > progress.lastAccepted.time, "endpoint is not ahead of accepted clock");
            check(CVodeSetStopTime(resources.ode, endTime), "CVodeSetStopTime");
            while (progress.lastAccepted.time < endTime) {
                if (cancel && cancel()) { progress.canceled = true; return progress; }
                require(static_cast<std::size_t>(progress.internalSteps) < controls.maximumAcceptedSteps, "accepted-step budget exhausted before endpoint");
                double reached = progress.lastAccepted.time;
                const int flag = CVode(resources.ode, endTime, resources.y, &reached, CV_ONE_STEP);
                require(flag >= 0, "CVODES flag " + std::to_string(flag) + ": " + integratorError + "; last RHS error: " + lastRhsError);
                require(std::isfinite(reached) && reached > progress.lastAccepted.time && reached <= endTime, "invalid accepted clock");
                ReactingState2D candidate = progress.lastAccepted; decode(reached, NV_DATA_S(resources.y), candidate);
                check(CVodeGetNumSteps(resources.ode, &progress.internalSteps), "CVodeGetNumSteps");
                require(static_cast<std::size_t>(progress.internalSteps) < std::numeric_limits<std::size_t>::max() - initial.steps, "step counter overflow");
                candidate.steps = initial.steps + static_cast<std::size_t>(progress.internalSteps); solver.validate(candidate);
                double quadratureTime = 0; check(CVodeGetQuad(resources.ode, &quadratureTime, resources.quadrature), "CVodeGetQuad");
                require(quadratureTime == reached, "quadrature clock differs from accepted state");
                ReactingConservative2D boundary(nv), chemistry(nv), constraint(nv);
                for (std::size_t k = 0; k < nv; ++k) {
                    boundary[k] = NV_Ith_S(resources.quadrature, k) * integralScale[k];
                    chemistry[k] = NV_Ith_S(resources.quadrature, nv + k) * integralScale[k];
                    constraint[k] = NV_Ith_S(resources.quadrature, 2 * nv + k) * integralScale[k];
                    require(std::isfinite(boundary[k]) && std::isfinite(chemistry[k]) && std::isfinite(constraint[k]), "nonfinite integral budget");
                }
                const double dt = reached - progress.lastAccepted.time;
                check(CVodeGetEstLocalErrors(resources.ode, resources.localError), "CVodeGetEstLocalErrors");
                check(CVodeGetErrWeights(resources.ode, resources.errorWeights), "CVodeGetErrWeights");
                const double localErrorNorm = N_VWrmsNorm(resources.localError, resources.errorWeights);
                // This is CVODES' own dimensionless weighted LTE acceptance
                // bound, not a flame-accuracy threshold. Allow accumulated
                // binary64 norm/rescaling roundoff only.
                const double normRoundoff = 64 * std::numeric_limits<double>::epsilon() * static_cast<double>(mapping.size());
                require(std::isfinite(localErrorNorm) && localErrorNorm <= 1 + normRoundoff,
                    "accepted state exceeds weighted local-error bound");
                progress.minimumAcceptedStep = progress.minimumAcceptedStep > 0 ? std::min(progress.minimumAcceptedStep, dt) : dt;
                progress.maximumAcceptedStep = std::max(progress.maximumAcceptedStep, dt);
                progress.lastAccepted = std::move(candidate); progress.boundaryImpulse = std::move(boundary);
                progress.chemistryChange = std::move(chemistry); progress.constraintChange = std::move(constraint);
                check(CVodeGetLastOrder(resources.ode, &progress.lastBdfOrder), "CVodeGetLastOrder");
                require(progress.lastBdfOrder >= 1 && progress.lastBdfOrder <= static_cast<int>(controls.maximumBdfOrder), "invalid accepted BDF order");
                ++progress.acceptedByBdfOrder[static_cast<std::size_t>(progress.lastBdfOrder)];
                progress.lastLocalErrorNorm = localErrorNorm;
                progress.maximumLocalErrorNorm = std::max(progress.maximumLocalErrorNorm, progress.lastLocalErrorNorm);
                statistics();
                if (onAccepted) onAccepted(progress);
            }
            progress.reachedEnd = progress.lastAccepted.time == endTime;
        } catch (const std::exception& error) { terminalFailure = true; progress.failure = error.what(); }
        // Preserve backend statistics on failure as well as successful exit.
        // A failed candidate never replaces lastAccepted or its quadratures.
        statistics();
        return progress;
    }
};

ReactingImplicitIntegrator2D::ReactingImplicitIntegrator2D(chemistry::DetailedGas& gas, ReactingFlowStepper2D& solver,
    const ReactingState2D& initial, const ReactingImplicitControls2D& controls)
    : impl_(std::make_unique<Impl>(gas, solver, initial, controls)) {}
ReactingImplicitIntegrator2D::~ReactingImplicitIntegrator2D() = default;
ReactingImplicitProgress2D ReactingImplicitIntegrator2D::advance(double endTime, const std::function<bool()>& cancel,
    const std::function<void(const ReactingImplicitProgress2D&)>& onAccepted,
    const std::function<void(const ReactingImplicitProgress2D&, double)>& onEvaluation) {
    return impl_->advance(endTime, cancel, onAccepted, onEvaluation);
}

} // namespace cartmesh2d::fv
