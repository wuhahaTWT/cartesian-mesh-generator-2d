#include "cartmesh2d/fv/ReactingFlow2D.hpp"
#ifdef CARTMESH2D_IMPLICIT_FLAME_PROBE
#include "cartmesh2d/fv/ReactingImplicit2D.hpp"
#endif
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::chemistry;
namespace {
void require(bool value, const char* text) { if (!value) throw std::runtime_error(text); }
template<class V> void numbers(std::ostream& out, const V& v) {
    out << '['; bool first = true; for (auto x : v) { if (!first) out << ','; first = false; out << x; } out << ']';
}
template<class V> void matrix(std::ostream& out, const V& v) {
    out << '['; bool first = true; for (const auto& x : v) { if (!first) out << ','; first = false; numbers(out, x); } out << ']';
}
void string(std::ostream& out, const std::string& value) {
    constexpr char hex[] = "0123456789abcdef"; out << '"';
    for (char character : value) {
        const auto c = static_cast<unsigned char>(character);
        if (c == '"' || c == '\\') out << '\\' << static_cast<char>(c);
        else if (c < 32) out << "\\u00" << hex[c / 16] << hex[c % 16];
        else out << static_cast<char>(c);
    }
    out << '"';
}
FvMesh2D strip(const std::vector<double>& x, std::size_t ny, double height) {
    require(x.size() > 2 && ny > 0 && height > 0 && std::isfinite(height), "invalid strip geometry");
    TopologyMesh2D t; const auto nx = x.size() - 1;
    for (std::size_t j = 0; j <= ny; ++j) for (double a : x)
        t.vertices.push_back({t.vertices.size(), {a, height * static_cast<double>(j) / static_cast<double>(ny)}});
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> edges;
    for (std::size_t j = 0; j < ny; ++j) for (std::size_t i = 0; i < nx; ++i) {
        require(std::isfinite(x[i]) && std::isfinite(x[i + 1]) && x[i + 1] > x[i], "non-increasing strip grid");
        TopologyCell2D c; c.id = t.cells.size(); c.geometryArea = (x[i + 1] - x[i]) * height / static_cast<double>(ny);
        const auto a = j * (nx + 1) + i; c.vertices = {a, a + 1, a + nx + 2, a + nx + 1};
        for (std::size_t k = 0; k < 4; ++k) {
            const auto u = c.vertices[k], v = c.vertices[(k + 1) % 4];
            const auto [it, inserted] = edges.emplace(std::minmax(u, v), t.edges.size());
            if (inserted) t.edges.push_back({it->second, u, v, c.id, {}, BoundaryPatch2D::DomainBoundary});
            else { t.edges[it->second].neighbour = c.id; t.edges[it->second].patch = BoundaryPatch2D::None; }
            c.edges.push_back(it->second);
        }
        t.cells.push_back(c);
    }
    return makeFvMesh2D(t);
}
void state(std::ostream& out, DetailedGas& gas, const ReactingState2D& s, const ReactingResidual2D& r) {
    out << "{\"time\":" << s.time << ",\"steps\":" << s.steps << ",\"U\":"; matrix(out, s.cells);
    out << ",\"temperature\":[";
    for (std::size_t i = 0; i < s.cells.size(); ++i) { if (i) out << ','; out << reactingPrimitive2D(gas, s.cells[i]).properties.temperature; }
    out << "],\"pressure\":[";
    for (std::size_t i = 0; i < s.cells.size(); ++i) { if (i) out << ','; out << reactingPrimitive2D(gas, s.cells[i]).properties.pressure; }
    out << "],\"transportDerivative\":"; matrix(out, r.transportDerivative);
    out << ",\"chemistryDerivative\":"; matrix(out, r.chemistryDerivative);
    out << ",\"derivative\":"; matrix(out, r.derivative);
    out << ",\"faceFlux\":"; matrix(out, r.faceFlux);
    out << ",\"transportRate\":"; numbers(out, r.transportRate);
    out << ",\"boundaryFlux\":"; numbers(out, r.boundaryFlux);
    out << ",\"chemistryIntegral\":"; numbers(out, r.chemistryIntegral);
    out << ",\"hlleFallbacks\":" << r.hlleFallbacks << '}';
}
#ifdef CARTMESH2D_IMPLICIT_FLAME_PROBE
void integration(std::ostream& out, const ReactingImplicitControls2D& c,
                 const ReactingImplicitProgress2D& r, std::size_t sampleEvery, std::size_t sampleEvaluations,
                 std::size_t maximumSamples, const std::string& samplingFailure) {
    out << "{\"method\":\"CVODES-BDF\",\"errorControlCorrection\":\"fullPredictorCorrection\",\"maximumOrder\":" << c.maximumBdfOrder
        << ",\"canceled\":" << (r.canceled ? "true" : "false")
        << ",\"maximumNonlinearIterations\":" << c.maximumNonlinearIterations
        << ",\"continueDampedNewton\":" << (c.continueDampedNewton ? "true" : "false")
        << ",\"reflectSpeciesNewton\":" << (c.reflectSpeciesNewton ? "true" : "false")
        << ",\"relativeTolerance\":" << c.relativeTolerance
        << ",\"spatialOrder\":" << c.spatialOrder
        << ",\"jacobianAdvectionOrder\":" << c.jacobianAdvectionOrder
        << ",\"absoluteConservedTolerance\":" << c.absoluteConservedTolerance
        << ",\"absoluteSpeciesFraction\":" << c.absoluteSpeciesFraction
        << ",\"maximumStep\":" << c.maximumStep
        << ",\"maximumResidualEvaluations\":" << c.maximumResidualEvaluations
        << ",\"rhsCalls\":" << r.rhsCalls << ",\"rejectedRhsCalls\":" << r.rejectedRhsCalls
        << ",\"errorTestFailures\":" << r.errorTestFailures << ",\"linearSetups\":" << r.linearSetups
        << ",\"jacobianEvaluations\":" << r.jacobianEvaluations
        << ",\"nonlinearIterations\":" << r.nonlinearIterations
        << ",\"nonlinearConvergenceFailures\":" << r.nonlinearConvergenceFailures
        << ",\"dampedNewtonUpdates\":" << r.dampedNewtonUpdates
        << ",\"continuedDampedNewtonUpdates\":" << r.continuedDampedNewtonUpdates
        << ",\"reflectedNewtonUpdates\":" << r.reflectedNewtonUpdates
        << ",\"reflectionAttempts\":" << r.reflectionAttempts
        << ",\"maximumReflectedSpeciesScaledChange\":" << r.maximumReflectedSpeciesScaledChange
        << ",\"maximumLocalErrorNorm\":" << r.maximumLocalErrorNorm
        << ",\"minimumNewtonFraction\":" << r.minimumNewtonFraction
        << ",\"bandHalfWidth\":" << r.bandHalfWidth << ",\"bandBytes\":" << r.bandBytes
        << ",\"sampleEveryAcceptedSteps\":" << sampleEvery
        << ",\"sampleResidualEvaluations\":" << sampleEvaluations
        << ",\"maximumSamples\":" << maximumSamples << ",\"samplingFailure\":"; string(out, samplingFailure);
    out << ",\"acceptedByBdfOrder\":"; numbers(out, r.acceptedByBdfOrder);
    out << ",\"lastDampedTrialFailure\":"; string(out, r.lastDampedTrialFailure); out << '}';
}
#endif
}
int main(int argc, char** argv) {
    try {
        const bool regression = argc == 4 && std::string(argv[3]) == "--regression";
#ifdef CARTMESH2D_IMPLICIT_FLAME_PROBE
        require(regression || (argc >= 5 && (argc - 5) % 2 == 0),
            "expected mechanism fixture output duration [--restart checkpoint] [--rtol value] [--conserved-atol value] [--species-atol value] [--bdf-order 1..5] [--jacobian-order 1|2] [--newton-iterations count] [--continue-damped 0|1] [--reflect-species 0|1] [--sample-every acceptedSteps] [--max-samples count], or mechanism fixture --regression");
#else
        require(argc == 5 || (argc == 7 && std::string(argv[5]) == "--restart") || regression,
            "expected mechanism fixture NEW-output duration [--restart checkpoint], or mechanism fixture --regression");
#endif
        const auto start = std::chrono::steady_clock::now();
        std::string restartPath, restartData;
        for (int argument = 5; argument < argc; argument += 2) if (std::string(argv[argument]) == "--restart") {
            require(restartPath.empty() && std::string(argv[argument + 1]).size() > 0, "invalid or duplicate restart checkpoint");
            restartPath = std::filesystem::absolute(argv[argument + 1]).string();
        }
        DetailedGas gas(argv[1]); std::ifstream input(argv[2]);
        std::string token; unsigned version = 0; std::size_t nx = 0, ny = 0, ns = 0; double height = 0;
        require(static_cast<bool>(input >> token >> version) && token == "CM2D_FLAME_FIXTURE" && version == 1, "invalid fixture version");
        require(static_cast<bool>(input >> nx >> ny >> ns >> height) && nx > 1 && nx <= 20000 && ny > 0 && ny <= 8
                && ns == gas.mechanism().species.size(), "invalid fixture sizes");
        for (const auto& species : gas.mechanism().species)
            require(static_cast<bool>(input >> token) && token == species, "fixture species order differs from mechanism");
        std::vector<double> x(nx + 1); for (double& v : x) require(static_cast<bool>(input >> v), "truncated fixture grid");
        std::vector<ReactingConservative2D> columns(nx, ReactingConservative2D(ns + 4));
        for (auto& u : columns) for (double& v : u) require(static_cast<bool>(input >> v), "truncated fixture cells");
        ReactingConservative2D inlet(ns + 4); for (double& v : inlet) require(static_cast<bool>(input >> v), "truncated inlet state");
        require(static_cast<bool>(input >> token) && token == "END" && !(input >> token), "invalid fixture ending");
        const auto mesh = strip(x, ny, height);
        std::vector<ReactingBoundary2D> bc; const auto in = reactingPrimitive2D(gas, inlet);
        for (std::size_t id = 0; id < mesh.faces.size(); ++id) if (!mesh.faces[id].neighbour) {
            ReactingBoundary2D b; b.face = id; const auto& f = mesh.faces[id];
            if (f.areaVector.x < 0) { b.kind = ReactingBoundaryKind2D::Reservoir; b.reservoir = in.gas; b.velocity = in.velocity; }
            else if (f.areaVector.x > 0) b.kind = ReactingBoundaryKind2D::ExtrapolatedOutflow;
            bc.push_back(b);
        }
        ReactingFlowStepper2D solver(gas, mesh, bc);
        std::vector<ReactingConservative2D> cells;
        for (std::size_t j = 0; j < ny; ++j) cells.insert(cells.end(), columns.begin(), columns.end());
        auto initial = solver.initialState(cells);
        if (!restartPath.empty()) {
            std::ifstream checkpoint(restartPath, std::ios::binary);
            require(static_cast<bool>(checkpoint), "cannot open restart checkpoint");
            std::ostringstream saved; saved << checkpoint.rdbuf();
            require(!checkpoint.bad() && static_cast<bool>(saved), "cannot read restart checkpoint");
            restartData = saved.str(); std::istringstream reader(restartData);
            initial = solver.readCheckpoint(reader);
        }
        auto current = initial;
        const auto initialResidual = solver.evaluateResidual(initial);
        if (regression) {
#ifdef CARTMESH2D_IMPLICIT_FLAME_PROBE
            ReactingImplicitControls2D controls; controls.maximumBdfOrder = 5;
            controls.maximumResidualEvaluations = 3000;
            ReactingImplicitIntegrator2D integrator(gas, solver, initial, controls);
            const auto result = integrator.advance(1e-8);
            std::cout << "implicit trace: time=" << result.lastAccepted.time << " steps=" << result.internalSteps
                << " RHS=" << result.rhsCalls << " damping=" << result.dampedNewtonUpdates
                << " continued=" << result.continuedDampedNewtonUpdates
                << " delta=" << result.lastNewtonCorrectionNorm << " residual=" << result.lastNewtonResidualNorm
                << " tolerance=" << result.lastNewtonTolerance << '\n';
            if (!result.reachedEnd) throw std::runtime_error(result.failure);
            require(result.lastAccepted.time == 1e-8 && result.lastAccepted.steps > 0, "implicit trace endpoint mismatch");
            std::vector<double> before(ns + 4), after(ns + 4);
            double energyScale = 0, sound = 0;
            for (std::size_t i = 0; i < cells.size(); ++i) {
                const auto p = reactingPrimitive2D(gas, initial.cells[i]);
                energyScale += mesh.cells[i].area * std::max(std::abs(initial.cells[i][3]),
                    initial.cells[i][0] * p.properties.cv * p.properties.temperature);
                sound = std::max(sound, p.soundSpeed);
                for (std::size_t k = 0; k < ns + 4; ++k) {
                    before[k] += mesh.cells[i].area * initial.cells[i][k];
                    after[k] += mesh.cells[i].area * result.lastAccepted.cells[i][k];
                }
            }
            std::vector<double> physical(ns + 4);
            for (std::size_t k = 0; k < ns + 4; ++k) physical[k] = after[k] - before[k] + result.boundaryImpulse[k];
            require(std::abs(physical[0]) / before[0] < 1e-8 && std::abs(physical[3]) / energyScale < 1e-8,
                    "implicit trace mass/energy budget failed");
            require(std::max(std::abs(physical[1]), std::abs(physical[2])) / (before[0] * sound) < 1e-8,
                    "implicit trace momentum budget failed");
            for (std::size_t k = 4; k < ns + 4; ++k)
                require(std::abs(physical[k] - result.chemistryChange[k] - result.constraintChange[k]) / before[0] < 1e-8,
                        "implicit trace species quadrature budget failed");
            const auto& m = gas.mechanism();
            for (std::size_t e = 0; e < m.elements.size(); ++e) {
                double element = 0;
                for (std::size_t k = 0; k < ns; ++k)
                    element += physical[k + 4] * m.atomCounts[k * m.elements.size() + e] * m.atomicWeights[e] / m.molecularWeights[k];
                require(std::abs(element) / before[0] < 1e-8, "implicit trace physical element budget failed");
            }
            (void)solver.evaluateResidual(result.lastAccepted);
            require(initial.cells == cells && initial.steps == 0 && initial.time == 0, "implicit trace modified input");
#else
            ReactingStepControls2D controls; controls.endTime = 1e-8;
            while (current.time < *controls.endTime && current.steps < 128) {
                const auto r = solver.advance(current, controls);
                if (!r.accepted) throw std::runtime_error(r.failure);
                require(r.rejectedReasons.empty(), "trace flame fixture requires rejected transport candidates");
                for (std::size_t k = 0; k < ns + 4; ++k) {
                    double scale = std::abs(r.beforeIntegral[k]) + std::abs(r.afterIntegral[k]) + std::abs(r.chemistryChange[k]);
                    for (const auto& f : r.faceFlux) scale += r.step * std::abs(f[k]);
                    require(std::abs(r.balanceError[k]) <= 3e-11 * scale, "trace regression conservation failed");
                }
                for (double v : r.elementalBalanceError) require(std::abs(v) < 2e-8, "trace regression elemental drift");
                current = *r.accepted;
            }
            require(current.time == *controls.endTime, "trace regression did not reach physical endpoint");
            require(initial.cells == cells && initial.steps == 0 && initial.time == 0, "trace regression modified initial state");
            (void)solver.evaluateResidual(current);
            std::cout << "Trace flame regression: " << mesh.cells.size() << " cells, " << current.steps << " accepted steps, no rejections\n";
#endif
            return 0;
        }
        std::size_t used = 0; const std::string durationText(argv[4]); const double duration = std::stod(durationText, &used);
        require(used == durationText.size() && duration >= 0 && std::isfinite(duration), "invalid physical duration");
        const double endTime = initial.time + duration;
        require(std::isfinite(endTime) && (duration == 0 || endTime > initial.time), "physical duration cannot advance checkpoint clock");
        const std::filesystem::path directory(argv[3]); require(!std::filesystem::exists(directory), "output directory already exists");
        std::filesystem::create_directories(directory);
        if (!restartPath.empty()) {
            std::ofstream checkpoint(directory / "restart.checkpoint", std::ios::binary); checkpoint << restartData;
            checkpoint.close(); require(static_cast<bool>(checkpoint), "cannot preserve restart checkpoint");
        }
        { std::ofstream definition(directory / "resolved-mechanism.yaml"); definition << gas.mechanism().resolvedDefinition;
          require(static_cast<bool>(definition), "cannot write resolved mechanism"); }
        std::ofstream log(directory / "steps.jsonl"); log << std::setprecision(17);
        std::vector<double> boundaryImpulse(ns + 4), chemistryChange(ns + 4);
        std::string failure;
#ifdef CARTMESH2D_IMPLICIT_FLAME_PROBE
        ReactingImplicitControls2D implicitControls;
        implicitControls.maximumBdfOrder = 5; implicitControls.maximumStep = std::max(duration, 1e-15);
        std::size_t sampleEvery = 0, sampleEvaluations = 0, lastSampleSteps = initial.steps, maximumSamples = 256;
        std::string samplingFailure;
        std::set<std::string> suppliedOptions;
        for (int argument = 5; argument < argc; argument += 2) {
            const std::string option(argv[argument]), value(argv[argument + 1]);
            require(suppliedOptions.insert(option).second, "duplicate implicit control option");
            if (option == "--restart") continue;
            std::size_t read = 0;
            if (option == "--rtol" || option == "--conserved-atol" || option == "--species-atol") {
                const double tolerance = std::stod(value, &read);
                require(read == value.size() && std::isfinite(tolerance) && tolerance > 0,
                    "invalid local integration tolerance");
                if (option == "--rtol") {
                    require(tolerance < 1, "relative local tolerance must be less than one");
                    implicitControls.relativeTolerance = tolerance;
                } else if (option == "--conserved-atol") implicitControls.absoluteConservedTolerance = tolerance;
                else implicitControls.absoluteSpeciesFraction = tolerance;
            } else if (option == "--continue-damped") {
                require(value == "0" || value == "1", "invalid damped continuation control");
                implicitControls.continueDampedNewton = value == "1";
            } else if (option == "--reflect-species") {
                require(value == "0" || value == "1", "invalid species reflection control");
                implicitControls.reflectSpeciesNewton = value == "1";
            } else {
                require(option == "--bdf-order" || option == "--jacobian-order" || option == "--newton-iterations"
                    || option == "--sample-every" || option == "--max-samples", "unknown implicit control option");
                const auto count = std::stoul(value, &read);
                require(read == value.size() && count > 0 && count <= static_cast<unsigned long>(std::numeric_limits<int>::max()), "invalid implicit integer control");
                if (option == "--bdf-order") {
                    require(count <= 5, "invalid maximum BDF order"); implicitControls.maximumBdfOrder = static_cast<unsigned>(count);
                } else if (option == "--jacobian-order") {
                    require(count <= 2, "invalid Jacobian advection order"); implicitControls.jacobianAdvectionOrder = static_cast<unsigned>(count);
                } else if (option == "--sample-every") sampleEvery = static_cast<std::size_t>(count);
                else if (option == "--max-samples") maximumSamples = static_cast<std::size_t>(count);
                else implicitControls.maximumNonlinearIterations = static_cast<unsigned>(count);
            }
        }
        ReactingImplicitProgress2D implicitResult;
        implicitResult.constraintChange.assign(ns + 4, 0);
        implicitResult.boundaryImpulse = boundaryImpulse; implicitResult.chemistryChange = chemistryChange;
        std::unique_ptr<DetailedGas> sampleGas;
        std::unique_ptr<ReactingFlowStepper2D> sampleSolver;
        std::ofstream samples;
        if (sampleEvery) {
            // Diagnostic queries have their own mutable thermodynamic and
            // transport context. Sampling never changes the integrator's
            // state, BDF history, requested endpoint, or gas-query sequence.
            sampleGas = std::make_unique<DetailedGas>(argv[1]);
            sampleSolver = std::make_unique<ReactingFlowStepper2D>(*sampleGas, mesh, bc);
            samples.open(directory / "samples.jsonl"); samples << std::setprecision(17);
            require(static_cast<bool>(samples), "cannot open accepted-state samples");
        }
        const auto writeSample = [&](const ReactingState2D& saved, const ReactingImplicitProgress2D& progress) {
            if (!samplingFailure.empty()) return;
            try {
                require(sampleEvaluations < maximumSamples, "accepted-state sample budget exhausted");
                const auto residual = sampleSolver->evaluateResidual(saved); ++sampleEvaluations;
                samples << "{\"kind\":\"" << (saved.steps == initial.steps ? "initial" : "accepted") << "\",\"state\":";
                state(samples, *sampleGas, saved, residual);
                samples << ",\"boundaryImpulse\":"; numbers(samples, progress.boundaryImpulse);
                samples << ",\"chemistryChange\":"; numbers(samples, progress.chemistryChange);
                samples << ",\"constraintChange\":"; numbers(samples, progress.constraintChange);
                samples << ",\"integration\":"; integration(samples, implicitControls, progress, sampleEvery, sampleEvaluations, maximumSamples, samplingFailure);
                samples << "}\n"; samples.flush(); require(static_cast<bool>(samples), "accepted-state sample write failed");
                lastSampleSteps = saved.steps;
            } catch (const std::exception& error) {
                // Stop at an accepted state and let the normal checkpoint and
                // final-state writers run even when optional diagnostics fail.
                samplingFailure = std::string("accepted-state sampling failed: ") + error.what();
            }
        };
        if (sampleEvery) writeSample(initial, implicitResult);
        std::ofstream evaluations(directory / "evaluation-progress.jsonl"); evaluations << std::setprecision(17);
        if (duration > 0) {
            ReactingImplicitIntegrator2D implicit(gas, solver, initial, implicitControls);
            implicitResult = implicit.advance(endTime, [&] {
                return !samplingFailure.empty() || std::filesystem::exists(directory / "cancel.request");
            }, [&](const auto& r) {
                const double dt = r.lastAccepted.time - current.time;
                log << "{\"time\":" << current.time << ",\"dt\":" << dt << ",\"accepted\":true,\"internalSteps\":"
                    << r.internalSteps << ",\"rhsCalls\":" << r.rhsCalls << ",\"rejectedRhsCalls\":" << r.rejectedRhsCalls
                    << ",\"bdfOrder\":" << r.lastBdfOrder << ",\"localErrorNorm\":" << r.lastLocalErrorNorm
                    << ",\"predictorCorrectionNorm\":" << r.lastPredictorCorrectionNorm << ",\"failure\":\"\"}\n";
                log.flush(); current = r.lastAccepted;
                if (sampleEvery && (current.steps - initial.steps) % sampleEvery == 0) writeSample(current, r);
                if (current.steps % 10 == 0) std::cerr << "implicit accepted " << current.steps << " t=" << current.time << '\n';
            }, [&](const auto& r, double trialTime) {
                if (r.rhsCalls == 1 || r.rhsCalls % 100 == 0) {
                    std::cerr << "implicit residuals=" << r.rhsCalls << " rejected=" << r.rejectedRhsCalls
                        << " accepted=" << r.lastAccepted.steps << " trial_t=" << trialTime
                        << " band=" << r.bandHalfWidth << '\n';
                    evaluations << "{\"rhsCalls\":" << r.rhsCalls << ",\"rejectedRhsCalls\":" << r.rejectedRhsCalls
                        << ",\"acceptedSteps\":" << r.lastAccepted.steps << ",\"acceptedTime\":" << r.lastAccepted.time
                        << ",\"trialTime\":" << trialTime << ",\"bandHalfWidth\":" << r.bandHalfWidth
                        << ",\"dampedNewtonUpdates\":" << r.dampedNewtonUpdates
                        << ",\"continuedDampedNewtonUpdates\":" << r.continuedDampedNewtonUpdates
                        << ",\"reflectedNewtonUpdates\":" << r.reflectedNewtonUpdates
                        << ",\"reflectionAttempts\":" << r.reflectionAttempts
                        << ",\"maximumReflectedSpeciesScaledChange\":" << r.maximumReflectedSpeciesScaledChange
                        << ",\"jacobianEvaluations\":" << r.jacobianEvaluations
                        << ",\"nonlinearIterations\":" << r.nonlinearIterations
                        << ",\"nonlinearConvergenceFailures\":" << r.nonlinearConvergenceFailures
                        << ",\"lastBdfOrder\":" << r.lastBdfOrder
                        << ",\"newtonCorrectionNorm\":" << r.lastNewtonCorrectionNorm
                        << ",\"newtonResidualNorm\":" << r.lastNewtonResidualNorm
                        << ",\"newtonTolerance\":" << r.lastNewtonTolerance << ",\"lastDampedTrialFailure\":";
                    string(evaluations, r.lastDampedTrialFailure); evaluations << "}\n";
                    evaluations.flush(); require(static_cast<bool>(evaluations), "evaluation log write failed");
                }
            });
            current = implicitResult.lastAccepted; failure = implicitResult.failure;
            if (!samplingFailure.empty()) failure = samplingFailure;
            else if (implicitResult.canceled && failure.empty()) failure = "canceled between accepted steps";
            boundaryImpulse = implicitResult.boundaryImpulse; chemistryChange = implicitResult.chemistryChange;
            if (!implicitResult.reachedEnd) {
                log << "{\"time\":" << current.time << ",\"dt\":0,\"accepted\":false,\"failure\":";
                string(log, failure); log << "}\n"; log.flush();
            }
        }
        if (sampleEvery && current.steps > lastSampleSteps) writeSample(current, implicitResult);
        if (!samplingFailure.empty()) failure = samplingFailure;
#else
        ReactingStepControls2D controls; controls.endTime = endTime;
        std::vector<double> transportClosure(ns + 4), chemistryClosure(ns + 4);
        double maximumClosure = 0, absoluteClosure = 0;
        std::size_t rejected = 0, sourceCalls = 0;
        while (current.time < endTime && current.steps - initial.steps < 20000) {
            const auto r = solver.advance(current, controls); sourceCalls += r.sourceCalls; rejected += r.rejectedReasons.size();
            log << "{\"time\":" << current.time << ",\"dt\":" << r.step << ",\"accepted\":" << (r.accepted ? "true" : "false")
                << ",\"courant\":" << r.combinedCourant << ",\"rejectedReasons\":[";
            for (std::size_t j = 0; j < r.rejectedReasons.size(); ++j) { if (j) log << ','; string(log, r.rejectedReasons[j]); }
            log << "],\"failure\":"; string(log, r.failure);
            log << ",\"maximumMassClosureFraction\":" << r.maximumMassClosureFraction
                << ",\"absoluteMassClosureIntegral\":" << r.absoluteMassClosureIntegral
                << ",\"transportMassClosureChange\":"; numbers(log, r.transportMassClosureChange);
            log << ",\"chemistryMassClosureChange\":"; numbers(log, r.chemistryMassClosureChange);
            log << "}\n"; log.flush();
            if (!r.accepted) { failure = r.failure; break; }
            for (std::size_t k = 0; k < ns + 4; ++k) {
                boundaryImpulse[k] += r.step * r.boundaryFlux[k]; chemistryChange[k] += r.chemistryChange[k];
                transportClosure[k] += r.transportMassClosureChange[k]; chemistryClosure[k] += r.chemistryMassClosureChange[k];
            }
            maximumClosure = std::max(maximumClosure, r.maximumMassClosureFraction);
            absoluteClosure += r.absoluteMassClosureIntegral;
            current = *r.accepted;
            if (current.steps % 50 == 0) std::cerr << "accepted " << current.steps << " t=" << current.time << '\n';
        }
#endif
        if (current.time < endTime && failure.empty()) failure = "step budget exhausted before physical endpoint";
        if (current.steps > 0) { std::ofstream checkpoint(directory / "accepted.checkpoint"); solver.writeCheckpoint(checkpoint, current); }
        // Preserve raw final accepted values before optional residual
        // diagnostics, which can themselves expose an interpolation failure.
        { std::ofstream saved(directory / "accepted-state.json"); saved << std::setprecision(17)
              << "{\"time\":" << current.time << ",\"steps\":" << current.steps << ",\"U\":";
          matrix(saved, current.cells); saved << ",\"failure\":"; string(saved, failure); saved << "}\n";
          require(static_cast<bool>(saved), "cannot save final accepted state"); }
        const auto finalResidual = solver.evaluateResidual(current);
        std::ofstream out(directory / "field.json"); out << std::setprecision(17);
        out << "{\"complete\":" << (current.time == endTime ? "true" : "false") << ",\"qualifiedFlame\":false,\"duration\":" << duration
            << ",\"endTime\":" << endTime << ",\"restart\":";
        if (restartPath.empty()) out << "null";
        else {
            out << "{\"sourcePath\":"; string(out, restartPath);
            out << ",\"checkpointFile\":\"restart.checkpoint\",\"restoresIntegratorHistory\":false,\"budgetOrigin\":\"initial-state\"}";
        }
        out << ",\"cells\":" << mesh.cells.size() << ",\"nx\":" << nx << ",\"ny\":" << ny << ",\"height\":" << height << ",\"xEdges\":";
        numbers(out, x); out << ",\"areas\":[";
        for (std::size_t i = 0; i < mesh.cells.size(); ++i) { if (i) out << ','; out << mesh.cells[i].area; }
        out << "],\"faces\":[";
        for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
            if (i) out << ','; const auto& f = mesh.faces[i]; out << "{\"owner\":" << f.owner << ",\"neighbour\":";
            if (f.neighbour) out << *f.neighbour; else out << "null";
            out << ",\"centre\":[" << f.centre.x << ',' << f.centre.y << "],\"S\":[" << f.areaVector.x << ',' << f.areaVector.y << "]}";
        }
        out << "],\"boundaryImpulse\":"; numbers(out, boundaryImpulse); out << ",\"chemistryChange\":"; numbers(out, chemistryChange);
#ifdef CARTMESH2D_IMPLICIT_FLAME_PROBE
        out << ",\"constraintChange\":"; numbers(out, implicitResult.constraintChange);
        out << ",\"integration\":"; integration(out, implicitControls, implicitResult, sampleEvery, sampleEvaluations, maximumSamples, samplingFailure);
#else
        out << ",\"massClosure\":{\"maximumFraction\":" << maximumClosure << ",\"absoluteIntegral\":" << absoluteClosure
            << ",\"transport\":"; numbers(out, transportClosure); out << ",\"chemistry\":"; numbers(out, chemistryClosure); out << '}';
#endif
        out << ",\"initial\":"; state(out, gas, initial, initialResidual);
        out << ",\"final\":"; state(out, gas, current, finalResidual);
#ifdef CARTMESH2D_IMPLICIT_FLAME_PROBE
        out << ",\"rejections\":null,\"sourceCalls\":null,\"failure\":";
#else
        out << ",\"rejections\":" << rejected << ",\"sourceCalls\":" << sourceCalls << ",\"failure\":";
#endif
        string(out, failure);
        out << ",\"elapsedSecondsBeforeFinalFlush\":" << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << "}\n";
        out.close(); require(static_cast<bool>(out) && static_cast<bool>(log), "verification output write failed");
        if (!failure.empty()) std::cerr << failure << '\n';
        return current.time == endTime && failure.empty() ? 0 : 1;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
