#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <optional>

#include "cartmesh2d/fv/detail/FlowAggregation2D.hpp"

namespace cartmesh2d::fv::detail {

using LinearVector2D = std::vector<double>;
enum class LinearPressureMethod2D { Jacobi, IC0, Aggregation };
enum class LinearSolveMethod2D { Jacobi, ILU0 };

inline void linearEnsure(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

inline double linearFinite(double value) {
    linearEnsure(std::isfinite(value), "Flow numerical range exceeded");
    return value;
}

// Explicit two-component expansion. No implicit conversion to double: callers
// must choose when to round and then audit the field they actually retain.
struct LinearTwofold2D {
    double high=0, low=0;
    void add(double value) {
        const double sum=linearFinite(high+value);
        const double split=sum-high;
        const double error=(high-(sum-split))+(value-split);
        const double tail=linearFinite(low+error);
        const double next=linearFinite(sum+tail);
        const double nextSplit=next-sum;
        low=linearFinite((sum-(next-nextSplit))+(tail-nextSplit));
        high=next;
    }
    void scale(double factor) {
        const double product=linearFinite(high*factor);
        const double tail=linearFinite(std::fma(high,factor,-product)+low*factor);
        high=product;low=0;add(tail);
    }
};

struct LinearCandidate2D {
    LinearVector2D high, low; // empty low means ordinary double arithmetic
    std::size_t iterations=0;
    double relaxedDouble(std::size_t i,double previous,double relaxation) const {
        if(low.empty() || low[i]==0)
            return linearFinite(previous+relaxation*(high[i]-previous));
        LinearTwofold2D value{high[i],low[i]};
        value.add(-previous);value.scale(relaxation);value.add(previous);
        return value.high;
    }
};

inline double linearProduct(const LinearVector2D& a, const LinearVector2D& b) {
    long double sum = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        sum += static_cast<long double>(a[i]) * b[i];
    return linearFinite(static_cast<double>(sum));
}

inline double linearNorm(const LinearVector2D& a) {
    // Scaled sum of squares avoids overflow/underflow without calling hypot
    // for every component. This is the same Euclidean norm and stopping target,
    // not a cheaper approximate residual. Standard scaled-norm principle:
    // https://www.netlib.org/lapack/explore-html/d8/d76/group__lassq.html
    // Independently implemented here; no LAPACK dependency. Accumulate in long double
    // (which may equal double on a given platform).
    long double scale = 0, sumSquares = 1;
    for (double value : a) {
        const long double magnitude = std::abs(static_cast<long double>(linearFinite(value)));
        if (magnitude == 0) continue;
        if (scale < magnitude) {
            const long double ratio = scale / magnitude;
            sumSquares = 1 + sumSquares * ratio * ratio;
            scale = magnitude;
        } else {
            const long double ratio = magnitude / scale;
            sumSquares += ratio * ratio;
        }
    }
    return linearFinite(static_cast<double>(scale * std::sqrt(sumSquares)));
}

// Fixed, sorted off-diagonal CSR graph. Duplicate faces between a cell pair
// share one coefficient, matching additive finite-volume assembly.
struct SparsePattern2D {
    std::vector<std::size_t> rows, columns, transpose, lowerEnd;

    SparsePattern2D(std::size_t n,
                    const std::vector<std::pair<std::size_t, std::size_t>>& connections) {
        linearEnsure(n > 0, "Flow sparse pattern is empty");
        std::vector<std::vector<std::size_t>> neighbours(n);
        for (auto [i, j] : connections) {
            linearEnsure(i < n && j < n && i != j, "Flow sparse connection invalid");
            neighbours[i].push_back(j);
            neighbours[j].push_back(i);
        }
        rows.push_back(0);
        for (std::size_t i = 0; i < n; ++i) {
            auto& row = neighbours[i];
            std::sort(row.begin(), row.end());
            row.erase(std::unique(row.begin(), row.end()), row.end());
            lowerEnd.push_back(columns.size() + static_cast<std::size_t>(
                std::lower_bound(row.begin(), row.end(), i) - row.begin()));
            columns.insert(columns.end(), row.begin(), row.end());
            rows.push_back(columns.size());
        }
        transpose.resize(columns.size());
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = rows[i]; k < rows[i + 1]; ++k)
                transpose[k] = slot(columns[k], i);
    }

    std::size_t slot(std::size_t row, std::size_t column) const {
        linearEnsure(row < rows.size() - 1, "Flow sparse row invalid");
        const auto first = columns.begin() + static_cast<std::ptrdiff_t>(rows[row]);
        const auto last = columns.begin() + static_cast<std::ptrdiff_t>(rows[row + 1]);
        const auto found = std::lower_bound(first, last, column);
        linearEnsure(found != last && *found == column, "Flow sparse connection missing");
        return static_cast<std::size_t>(found - columns.begin());
    }
};

// Shared by sequential U, V and pressure solves; no vector allocation per
// Krylov step. The original matrix is retained for true-residual checks.
struct LinearWorkspace2D {
    LinearVector2D r, r0, p, v, s, t, z, zs, ax;
    explicit LinearWorkspace2D(std::size_t n)
        : r(n), r0(n), p(n), v(n), s(n), t(n), z(n), zs(n), ax(n) {}
};

struct SparseSystem2D {
    const SparsePattern2D& pattern;
    LinearVector2D diag, rhs, off;
    explicit SparseSystem2D(const SparsePattern2D& graph)
        : pattern(graph), diag(graph.rows.size() - 1), rhs(diag.size()), off(graph.columns.size()) {}

    std::size_t ic0Builds() const { return factorBuilds_; }
    std::size_t ic0Reuses() const { return factorReuses_; }
    std::size_t hierarchyBuilds() const { return hierarchyBuilds_; }
    std::size_t hierarchyReuses() const { return hierarchyReuses_; }
    std::size_t hierarchyRefreshes() const { return hierarchyRefreshes_; }
    std::size_t hierarchyLevels() const { return hierarchy_ ? hierarchy_->levels() : 0; }
    std::size_t hierarchyCoarseCells() const { return hierarchy_ ? hierarchy_->coarseCells() : 0; }

    std::size_t ilu0Builds() const {return iluBuilds_;}
    std::size_t ilu0Reuses() const {return iluReuses_;}

    // Natural-order nonsymmetric ILU(0): eliminate only entries in the original
    // CSR graph. Independent implementation of the standard zero-fill method;
    // see https://www.netlib.org/templates/templates.pdf, section 3.4.
    // No diagonal shift, pivot replacement, symmetrization or fallback.
    void factorILU0() const {
        linearEnsure(diag.size()+1==pattern.rows.size() && off.size()==pattern.columns.size(),
                     "Flow ILU0 dimensions invalid");
        if(iluReady_ && iluMatrixDiagonal_==diag && iluMatrixOff_==off) {++iluReuses_;return;}
        iluReady_=false;
        for(double d:diag)linearEnsure(std::isfinite(d)&&d>0,"Flow ILU0 nonpositive/nonfinite diagonal");
        for(double v:off)linearEnsure(std::isfinite(v),"Flow ILU0 nonfinite coefficient");
        iluDiagonal_=diag;iluOff_=off;
        for(std::size_t i=0;i<diag.size();++i) {
            for(auto ij=pattern.rows[i];ij<pattern.lowerEnd[i];++ij) {
                const auto j=pattern.columns[ij];
                const double factor=linearFinite(iluOff_[ij]/iluDiagonal_[j]);
                iluOff_[ij]=factor;
                auto ik=ij+1;
                for(auto jk=pattern.lowerEnd[j];jk<pattern.rows[j+1];++jk) {
                    const auto k=pattern.columns[jk];
                    if(k==i)iluDiagonal_[i]=linearFinite(iluDiagonal_[i]-factor*iluOff_[jk]);
                    else {
                        while(ik<pattern.rows[i+1] && pattern.columns[ik]<k)++ik;
                        if(ik<pattern.rows[i+1] && pattern.columns[ik]==k)
                            iluOff_[ik]=linearFinite(iluOff_[ik]-factor*iluOff_[jk]);
                    }
                }
            }
            linearEnsure(std::isfinite(iluDiagonal_[i])&&iluDiagonal_[i]>0,
                         "Flow ILU0 nonpositive/nonfinite pivot");
        }
        iluMatrixDiagonal_=diag;iluMatrixOff_=off;iluReady_=true;++iluBuilds_;
    }
    // Public diagnostic application rejects stale factors after direct coefficient
    // mutation. Krylov uses the private application after factorILU0 once per
    // solve; coefficients remain unchanged throughout that solve.
    void preconditionILU0(const LinearVector2D& input,LinearVector2D& output) const {
        linearEnsure(iluMatrixDiagonal_==diag && iluMatrixOff_==off,
                     "Flow ILU0 factors are stale");
        applyILU0(input,output);
    }

    void reset() {
        factorReady_ = false;iluReady_=false;
        hierarchyReady_ = false;
        std::fill(diag.begin(), diag.end(), 0.);
        std::fill(rhs.begin(), rhs.end(), 0.);
        std::fill(off.begin(), off.end(), 0.);
    }
    void add(std::size_t row, std::size_t column, double value) {
        factorReady_ = false;iluReady_=false;
        hierarchyReady_ = false;
        off[pattern.slot(row, column)] += value;
    }
    void apply(const LinearVector2D& x, LinearVector2D& y) const {
        linearEnsure(x.size() == diag.size() && &x != &y, "Flow matrix vector size or alias invalid");
        y.resize(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) {
            y[i] = diag[i] * x[i];
            for (std::size_t k = pattern.rows[i]; k < pattern.rows[i + 1]; ++k)
                y[i] += off[k] * x[pattern.columns[k]];
            linearFinite(y[i]);
        }
    }
    // Compensated b_i - A_i*x: preserve both multiplication and addition
    // roundoff. In particular, do not first round A_i*x and then subtract b_i.
    // Error-free TwoSum and FMA-product principles (independent implementation):
    // Ogita, Rump, Oishi, SIAM SISC 26(6), doi:10.1137/030601818.
    // Finite normal arithmetic is required; this is not arbitrary precision.
    double compensatedResidualRow(std::size_t i, const LinearVector2D& x,
                                  const LinearVector2D* low=nullptr) const {
        double sum=rhs[i], correction=0;
        const auto subtractProduct=[&](double coefficient,double value) {
            const double product=linearFinite(-coefficient*value);
            const double productError=std::fma(-coefficient,value,-product);
            const double next=linearFinite(sum+product);
            const double split=next-sum;
            const double additionError=(sum-(next-split))+(product-split);
            correction=linearFinite(correction+(additionError+productError));
            sum=next;
        };
        subtractProduct(diag[i],x[i]);
        if(low)subtractProduct(diag[i],(*low)[i]);
        for(auto k=pattern.rows[i];k<pattern.rows[i+1];++k)
        {
            subtractProduct(off[k],x[pattern.columns[k]]);
            if(low)subtractProduct(off[k],(*low)[pattern.columns[k]]);
        }
        return linearFinite(sum+correction);
    }
    void pin(std::size_t id) {
        linearEnsure(id < diag.size(), "Flow pressure gauge invalid");
        factorReady_ = false;iluReady_=false;
        hierarchyReady_ = false;
        for (std::size_t k = pattern.rows[id]; k < pattern.rows[id + 1]; ++k)
            off[k] = off[pattern.transpose[k]] = 0.;
        rhs[id] = 0.;
    }

    // Natural-order incomplete LDL^T with zero fill (IC(0)). No diagonal shift,
    // threshold relaxation or silent fallback. A failed pivot is explicit.
    void factorIC0() const {
        // Coefficients remain public for assembly. Exact snapshots guard direct
        // writes as well as RHS-only reuse; no hash collision or approximate
        // matrix comparison may authorize a stale preconditioner. The pattern
        // is immutable for the lifetime of this system, as for CSR assembly.
        if (factorReady_ && factorMatrixDiagonal_ == diag && factorMatrixOff_ == off) {
            ++factorReuses_;
            return;
        }
        factorReady_ = false; // A failed factorization must not leave a usable cache.
        factorLower_.assign(off.size(), 0.);
        factorDiagonal_.resize(diag.size());
        for (std::size_t k = 0; k < off.size(); ++k)
            linearEnsure(std::isfinite(off[k]) && off[k] == off[pattern.transpose[k]],
                         "Flow IC0 requires a finite symmetric matrix");
        for (std::size_t i = 0; i < diag.size(); ++i) {
            double pivot = diag[i];
            for (std::size_t ij = pattern.rows[i]; ij < pattern.lowerEnd[i]; ++ij) {
                const auto j = pattern.columns[ij];
                double value = off[ij];
                auto ik = pattern.rows[i];
                auto jk = pattern.rows[j];
                while (ik < ij && jk < pattern.lowerEnd[j]) {
                    const auto ci = pattern.columns[ik], cj = pattern.columns[jk];
                    if (ci == cj) {
                        value -= factorLower_[ik] * factorDiagonal_[ci] * factorLower_[jk];
                        ++ik;
                        ++jk;
                    } else if (ci < cj) ++ik;
                    else ++jk;
                }
                const double factor = linearFinite(value / factorDiagonal_[j]);
                factorLower_[ij] = factor;
                pivot -= factor * factor * factorDiagonal_[j];
            }
            linearEnsure(std::isfinite(pivot) && pivot > 0, "Flow IC0 nonpositive/nonfinite pivot");
            factorDiagonal_[i] = pivot;
        }
        factorMatrixDiagonal_ = diag;
        factorMatrixOff_ = off;
        factorReady_ = true;
        ++factorBuilds_;
    }

    void precondition(LinearWorkspace2D& w, LinearPressureMethod2D method) const {
        if (method == LinearPressureMethod2D::Aggregation) {
            linearEnsure(hierarchyReady_ && hierarchy_.has_value(), "Flow aggregation hierarchy unavailable");
            hierarchy_->apply(w.r,w.z);
            return;
        }
        if (method == LinearPressureMethod2D::Jacobi) {
            for (std::size_t i = 0; i < diag.size(); ++i) w.z[i] = w.r[i] / diag[i];
            return;
        }
        for (std::size_t i = 0; i < diag.size(); ++i) {
            double value = w.r[i];
            for (std::size_t k = pattern.rows[i]; k < pattern.lowerEnd[i]; ++k)
                value -= factorLower_[k] * w.z[pattern.columns[k]];
            w.z[i] = linearFinite(value);
        }
        for (std::size_t i = 0; i < diag.size(); ++i)
            w.z[i] = linearFinite(w.z[i] / factorDiagonal_[i]);
        for (std::size_t i = diag.size(); i-- > 0;)
            for (std::size_t k = pattern.rows[i]; k < pattern.lowerEnd[i]; ++k)
                w.z[pattern.columns[k]] -= factorLower_[k] * w.z[i];
    }

    std::size_t solvePressure(LinearVector2D& x, LinearWorkspace2D& w, bool ic0 = true) const {
        return solvePressure(x,w,ic0?LinearPressureMethod2D::IC0:LinearPressureMethod2D::Jacobi);
    }
    std::size_t solvePressure(LinearVector2D& x, LinearWorkspace2D& w, LinearPressureMethod2D method) const {
        linearEnsure(method==LinearPressureMethod2D::Jacobi || method==LinearPressureMethod2D::IC0 ||
                     method==LinearPressureMethod2D::Aggregation, "Flow pressure preconditioner invalid");
        linearEnsure(x.size() == diag.size() && w.r.size() == diag.size(), "Flow linear workspace size invalid");
        for (double d : diag)
            linearEnsure(d > 0 && std::isfinite(d), "Flow pressure matrix diagonal invalid");
        auto& residual = w.r;
        auto& z = w.z;
        auto& direction = w.p;
        auto& ax = w.ax;
        apply(x, ax);
        for (std::size_t i = 0; i < x.size(); ++i) residual[i] = rhs[i] - ax[i];
        const double stop = linearFinite(1e-13 + 1e-11 * linearNorm(rhs));
        if (linearNorm(residual) <= stop) return 0;
        if (method == LinearPressureMethod2D::IC0) factorIC0();
        else if (method == LinearPressureMethod2D::Aggregation) {
            hierarchyReady_=false;
            try {
                if (hierarchy_ && hierarchy_->matches(diag,off)) ++hierarchyReuses_;
                else if (hierarchy_ && refreshesSinceBuild_<8 && hierarchy_->refresh(diag,off)) {
                    ++hierarchyRefreshes_; ++refreshesSinceBuild_;
                } else {
                    // Periodic regrouping bounds how long the interpolation
                    // topology lags changing strengths. This is a setup policy,
                    // not a relaxation of the true-residual stopping criterion.
                    hierarchy_.reset();
                    hierarchy_.emplace(pattern.rows,pattern.columns,diag,off);
                    refreshesSinceBuild_=0; ++hierarchyBuilds_;
                }
                hierarchyReady_=true;
            } catch (...) {
                // A failed numerical refresh must never remain a usable cache.
                hierarchy_.reset(); refreshesSinceBuild_=0; throw;
            }
        }
        precondition(w, method);
        direction = z;
        double rz = linearProduct(residual, z);
        for (std::size_t iteration = 0; iteration < 3000; ++iteration) {
            auto& ad = w.v;
            apply(direction, ad);
            const double denominator = linearProduct(direction, ad);
            linearEnsure(denominator > 0 && rz > 0, "Flow pressure PCG lost positive definiteness");
            const double alpha = rz / denominator;
            for (std::size_t i = 0; i < x.size(); ++i) {
                x[i] += alpha * direction[i];
                residual[i] -= alpha * ad[i];
            }
            if (linearNorm(residual) <= stop) {
                apply(x, ax);
                for (std::size_t i = 0; i < x.size(); ++i) residual[i] = rhs[i] - ax[i];
                if (linearNorm(residual) <= stop) return iteration + 1;
                precondition(w, method);
                direction = z;
                rz = linearProduct(residual, z);
                continue;
            }
            precondition(w, method);
            const double next = linearProduct(residual, z);
            const double beta = next / rz;
            for (std::size_t i = 0; i < x.size(); ++i)
                direction[i] = z[i] + beta * direction[i];
            rz = next;
        }
        apply(x, ax);
        for (std::size_t i = 0; i < x.size(); ++i) residual[i] = rhs[i] - ax[i];
        std::ostringstream message;
        message << "Flow pressure PCG iteration limit reached: true residual=" << linearNorm(residual)
                << ", target=" << stop << ", rhs=" << linearNorm(rhs);
        throw std::runtime_error(message.str());
    }

    std::size_t solve(LinearVector2D& x, LinearWorkspace2D& w,
                      double diagonalScaledStop = std::numeric_limits<double>::infinity(),
                      double residualNormStop = std::numeric_limits<double>::infinity(),
                      LinearSolveMethod2D method=LinearSolveMethod2D::Jacobi) const {
        return solveImpl(x,w,diagonalScaledStop,residualNormStop,nullptr,method);
    }

    // The original gates apply to high+low, not to high alone. A transport
    // caller must round after relaxation and recheck its original equations.
    LinearCandidate2D solveCandidate(LinearVector2D initial,LinearWorkspace2D& w,
                      double diagonalScaledStop,double residualNormStop,
                      LinearSolveMethod2D method=LinearSolveMethod2D::Jacobi) const {
        LinearCandidate2D candidate{std::move(initial),{},0};
        candidate.iterations=solveImpl(candidate.high,w,diagonalScaledStop,
                                      residualNormStop,&candidate.low,method);
        return candidate;
    }
private:
    void applyILU0(const LinearVector2D& input,LinearVector2D& output) const {
        linearEnsure(iluReady_ && input.size()==diag.size() && output.size()==diag.size() && &input!=&output,
                     "Flow ILU0 factors/workspace unavailable");
        for(std::size_t i=0;i<diag.size();++i) {
            double value=input[i];
            for(auto k=pattern.rows[i];k<pattern.lowerEnd[i];++k)value-=iluOff_[k]*output[pattern.columns[k]];
            output[i]=linearFinite(value);
        }
        for(std::size_t i=diag.size();i-- >0;) {
            double value=output[i];
            for(auto k=pattern.lowerEnd[i];k<pattern.rows[i+1];++k)value-=iluOff_[k]*output[pattern.columns[k]];
            output[i]=linearFinite(value/iluDiagonal_[i]);
        }
    }

    std::size_t solveImpl(LinearVector2D& x,LinearWorkspace2D& w,
                         double diagonalScaledStop,double residualNormStop,
                         LinearVector2D* tail,LinearSolveMethod2D method) const {
        linearEnsure(method==LinearSolveMethod2D::Jacobi || method==LinearSolveMethod2D::ILU0,
                     "Flow linear preconditioner invalid");
        linearEnsure(diagonalScaledStop > 0 && !std::isnan(diagonalScaledStop), "Invalid momentum diagonal-scaled residual tolerance");
        linearEnsure(residualNormStop > 0 && !std::isnan(residualNormStop), "Invalid linear residual norm tolerance");
        linearEnsure(x.size() == diag.size() && w.r.size() == diag.size(), "Flow linear workspace size invalid");
        for (double d : diag)
            linearEnsure(d > 0 && std::isfinite(d), "Flow singular/nonpositive matrix diagonal");
        auto& ax = w.ax; auto& r = w.r; auto& r0 = w.r0;
        auto& p = w.p; auto& v = w.v; auto& s = w.s; auto& t = w.t;
        auto& z = w.z; auto& zs = w.zs;
        apply(x, ax);
        for (std::size_t i = 0; i < x.size(); ++i) r[i] = rhs[i] - ax[i];
        r0 = r;
        std::fill(p.begin(), p.end(), 0.);
        std::fill(v.begin(), v.end(), 0.);
        const double stop = std::min(linearFinite(1e-13 + 1e-11 * linearNorm(rhs)),residualNormStop);
        const auto smallResidual = [&](const LinearVector2D& residual) {
            if (linearNorm(residual)>stop) return false;
            // In addition to the existing global norm, constrain every row in
            // solution units. A large far-field RHS must not mask a tiny cell.
            for (std::size_t i=0;i<diag.size();++i)
                if (std::abs(residual[i])/diag[i]>diagonalScaledStop) return false;
            return true;
        };
        const auto nearRoundoff=[&] {
            double magnitude=0;
            for(double value:x)magnitude=std::max(magnitude,std::abs(value));
            return diagonalScaledStop<16*std::numeric_limits<double>::epsilon()*magnitude;
        };
        bool compensated=nearRoundoff(), refined=false;
        if(tail) {
            tail->clear();
            if(compensated)tail->assign(x.size(),0.);
        }
        const auto rowResidual=[&](std::size_t i) {
            return compensatedResidualRow(i,x,tail&&!tail->empty()?tail:nullptr);
        };
        const auto addUpdate=[&](std::size_t i,double update) {
            if(tail&&!tail->empty()) {
                LinearTwofold2D value{x[i],(*tail)[i]};value.add(update);
                x[i]=value.high;(*tail)[i]=value.low;
            } else x[i]=linearFinite(x[i]+update);
        };
        if(compensated)
            for(std::size_t i=0;i<x.size();++i)r[i]=rowResidual(i);
        r0=r;
        if (smallResidual(r)) return 0;
        if(method==LinearSolveMethod2D::ILU0)factorILU0();
        double rhoOld = 1, alpha = 1, omega = 1;
        const auto restart = [&] {
            // r is the freshly recomputed b-Ax, not a recursively updated
            // estimate. A new shadow residual can recover exact biorthogonal
            // breakdown without changing the matrix, iterate or stopping gates.
            r0=r;
            std::fill(p.begin(),p.end(),0.);
            std::fill(v.begin(),v.end(),0.);
            rhoOld=alpha=omega=1;
        };
        for (std::size_t step = 0; step < 3000; ++step) {
            double rho = linearProduct(r0, r);
            if(rho==0 || omega==0) {
                restart();rho=linearProduct(r0,r);
                linearEnsure(rho>0,"Flow BiCGStab residual product underflow");
            }
            const double beta = (rho / rhoOld) * (alpha / omega);
            if(method==LinearSolveMethod2D::Jacobi) {
                for (std::size_t i = 0; i < x.size(); ++i) {
                    p[i] = r[i] + beta * (p[i] - omega * v[i]);
                    z[i] = p[i] / diag[i];
                }
            } else {
                for(std::size_t i=0;i<x.size();++i)p[i]=r[i]+beta*(p[i]-omega*v[i]);
                applyILU0(p,z);
            }
            apply(z, v);
            const double rv = linearProduct(r0, v);
            linearEnsure(rv != 0, "Flow BiCGStab singular projection");
            alpha = rho / rv;
            for (std::size_t i = 0; i < x.size(); ++i) s[i] = r[i] - alpha * v[i];
            if (smallResidual(s)) {
                for (std::size_t i = 0; i < x.size(); ++i) addUpdate(i,alpha*z[i]);
            } else {
                if(method==LinearSolveMethod2D::ILU0)applyILU0(s,zs);
                else for (std::size_t i = 0; i < x.size(); ++i) zs[i] = s[i] / diag[i];
                apply(zs, t);
                const double tt = linearProduct(t, t);
                linearEnsure(tt > 0, "Flow BiCGStab null update");
                omega = linearProduct(t, s) / tt;
                for (std::size_t i = 0; i < x.size(); ++i)
                    addUpdate(i,alpha*z[i]+omega*zs[i]);
            }
            if(!compensated && nearRoundoff()) {
                compensated=true;
                if(tail)tail->assign(x.size(),0.);
            }
            if(compensated) {
                for(std::size_t i=0;i<x.size();++i)r[i]=rowResidual(i);
            } else {
                apply(x, ax);
                for (std::size_t i = 0; i < x.size(); ++i) r[i] = rhs[i] - ax[i];
            }
            if (smallResidual(r)) return step + 1;
            // Near the floating-point floor, use accurately recomputed b-Ax
            // for every acceptance. Once only, try bounded coordinate refinement
            // for a diagonally dominant matrix. Gates and matrix stay unchanged;
            // unrepresentable accuracy requests still fail. Iteration counts
            // remain Krylov steps, not a count of these additional row sweeps.
            if(compensated && !refined && step%40==39 && linearNorm(r)<=stop) {
                refined=true;
                bool diagonallyDominant=true;
                for(std::size_t i=0;i<x.size();++i) {
                    double rowSum=0;
                    for(auto k=pattern.rows[i];k<pattern.rows[i+1];++k)rowSum+=std::abs(off[k]);
                    if(rowSum>diag[i])diagonallyDominant=false;
                }
                for(std::size_t sweep=0;diagonallyDominant && sweep<8;++sweep) {
                    for(std::size_t i=0;i<x.size();++i)
                        addUpdate(i,rowResidual(i)/diag[i]);
                    for(std::size_t i=0;i<x.size();++i)r[i]=rowResidual(i);
                    if(smallResidual(r))return step+1;
                }
            }
            if (step % 40 == 39) {
                r0 = r;
                std::fill(p.begin(), p.end(), 0.);
                std::fill(v.begin(), v.end(), 0.);
                rhoOld = alpha = omega = 1;
            } else rhoOld = rho;
        }
        double scaled=0;
        for (std::size_t i=0;i<diag.size();++i) scaled=std::max(scaled,std::abs(r[i])/diag[i]);
        std::ostringstream message;
        message << "Flow linear solver iteration limit reached: true residual=" << linearNorm(r)
                << ", target=" << stop << ", diagonal-scaled residual=" << scaled
                << ", diagonal-scaled target=" << diagonalScaledStop << ", rhs=" << linearNorm(rhs);
        throw std::runtime_error(message.str());
    }
private:
    // Owned by this matrix, so another system using the same Krylov workspace
    // cannot overwrite its factors. Storage is allocated lazily for pressure.
    mutable LinearVector2D iluDiagonal_,iluOff_,iluMatrixDiagonal_,iluMatrixOff_;
    mutable bool iluReady_=false;
    mutable std::size_t iluBuilds_=0,iluReuses_=0;
    mutable LinearVector2D factorDiagonal_, factorLower_;
    mutable LinearVector2D factorMatrixDiagonal_, factorMatrixOff_;
    mutable std::optional<AggregationHierarchy2D> hierarchy_;
    mutable bool hierarchyReady_ = false;
    mutable std::size_t hierarchyBuilds_ = 0, hierarchyReuses_ = 0;
    mutable std::size_t hierarchyRefreshes_ = 0, refreshesSinceBuild_ = 0;
    mutable bool factorReady_ = false;
    mutable std::size_t factorBuilds_ = 0, factorReuses_ = 0;

};

} // namespace cartmesh2d::fv::detail
