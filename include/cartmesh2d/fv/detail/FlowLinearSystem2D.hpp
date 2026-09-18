#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cartmesh2d::fv::detail {

using LinearVector2D = std::vector<double>;

inline void linearEnsure(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

inline double linearFinite(double value) {
    linearEnsure(std::isfinite(value), "Flow numerical range exceeded");
    return value;
}

inline double linearProduct(const LinearVector2D& a, const LinearVector2D& b) {
    long double sum = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        sum += static_cast<long double>(a[i]) * b[i];
    return linearFinite(static_cast<double>(sum));
}

inline double linearNorm(const LinearVector2D& a) {
    long double result = 0;
    for (double value : a)
        result = std::hypot(result, static_cast<long double>(linearFinite(value)));
    return linearFinite(static_cast<double>(result));
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
    LinearVector2D r, r0, p, v, s, t, z, zs, ax, factorDiagonal, factorLower;
    explicit LinearWorkspace2D(std::size_t n)
        : r(n), r0(n), p(n), v(n), s(n), t(n), z(n), zs(n), ax(n), factorDiagonal(n) {}
};

struct SparseSystem2D {
    const SparsePattern2D& pattern;
    LinearVector2D diag, rhs, off;
    explicit SparseSystem2D(const SparsePattern2D& graph)
        : pattern(graph), diag(graph.rows.size() - 1), rhs(diag.size()), off(graph.columns.size()) {}

    void reset() {
        std::fill(diag.begin(), diag.end(), 0.);
        std::fill(rhs.begin(), rhs.end(), 0.);
        std::fill(off.begin(), off.end(), 0.);
    }
    void add(std::size_t row, std::size_t column, double value) {
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
    void pin(std::size_t id) {
        linearEnsure(id < diag.size(), "Flow pressure gauge invalid");
        for (std::size_t k = pattern.rows[id]; k < pattern.rows[id + 1]; ++k)
            off[k] = off[pattern.transpose[k]] = 0.;
        rhs[id] = 0.;
    }

    // Natural-order incomplete LDL^T with zero fill (IC(0)). No diagonal shift,
    // threshold relaxation or silent fallback. A failed pivot is explicit.
    void factorIC0(LinearWorkspace2D& w) const {
        w.factorLower.assign(off.size(), 0.);
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
                        value -= w.factorLower[ik] * w.factorDiagonal[ci] * w.factorLower[jk];
                        ++ik;
                        ++jk;
                    } else if (ci < cj) ++ik;
                    else ++jk;
                }
                const double factor = linearFinite(value / w.factorDiagonal[j]);
                w.factorLower[ij] = factor;
                pivot -= factor * factor * w.factorDiagonal[j];
            }
            linearEnsure(std::isfinite(pivot) && pivot > 0, "Flow IC0 nonpositive/nonfinite pivot");
            w.factorDiagonal[i] = pivot;
        }
    }

    void precondition(LinearWorkspace2D& w, bool ic0) const {
        if (!ic0) {
            for (std::size_t i = 0; i < diag.size(); ++i) w.z[i] = w.r[i] / diag[i];
            return;
        }
        for (std::size_t i = 0; i < diag.size(); ++i) {
            double value = w.r[i];
            for (std::size_t k = pattern.rows[i]; k < pattern.lowerEnd[i]; ++k)
                value -= w.factorLower[k] * w.z[pattern.columns[k]];
            w.z[i] = linearFinite(value);
        }
        for (std::size_t i = 0; i < diag.size(); ++i)
            w.z[i] = linearFinite(w.z[i] / w.factorDiagonal[i]);
        for (std::size_t i = diag.size(); i-- > 0;)
            for (std::size_t k = pattern.rows[i]; k < pattern.lowerEnd[i]; ++k)
                w.z[pattern.columns[k]] -= w.factorLower[k] * w.z[i];
    }

    std::size_t solvePressure(LinearVector2D& x, LinearWorkspace2D& w, bool ic0 = true) const {
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
        if (ic0) factorIC0(w);
        precondition(w, ic0);
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
                precondition(w, ic0);
                direction = z;
                rz = linearProduct(residual, z);
                continue;
            }
            precondition(w, ic0);
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
                      double diagonalScaledStop = std::numeric_limits<double>::infinity()) const {
        linearEnsure(diagonalScaledStop > 0 && !std::isnan(diagonalScaledStop), "Invalid momentum diagonal-scaled residual tolerance");
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
        const double stop = linearFinite(1e-13 + 1e-11 * linearNorm(rhs));
        const auto smallResidual = [&](const LinearVector2D& residual) {
            if (linearNorm(residual)>stop) return false;
            // In addition to the existing global norm, constrain every row in
            // solution units. A large far-field RHS must not mask a tiny cell.
            for (std::size_t i=0;i<diag.size();++i)
                if (std::abs(residual[i])/diag[i]>diagonalScaledStop) return false;
            return true;
        };
        if (smallResidual(r)) return 0;
        double rhoOld = 1, alpha = 1, omega = 1;
        for (std::size_t step = 0; step < 3000; ++step) {
            const double rho = linearProduct(r0, r);
            linearEnsure(rho != 0 && omega != 0, "Flow BiCGStab breakdown");
            const double beta = (rho / rhoOld) * (alpha / omega);
            for (std::size_t i = 0; i < x.size(); ++i) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
                z[i] = p[i] / diag[i];
            }
            apply(z, v);
            const double rv = linearProduct(r0, v);
            linearEnsure(rv != 0, "Flow BiCGStab singular projection");
            alpha = rho / rv;
            for (std::size_t i = 0; i < x.size(); ++i) s[i] = r[i] - alpha * v[i];
            if (smallResidual(s)) {
                for (std::size_t i = 0; i < x.size(); ++i) x[i] += alpha * z[i];
            } else {
                for (std::size_t i = 0; i < x.size(); ++i) zs[i] = s[i] / diag[i];
                apply(zs, t);
                const double tt = linearProduct(t, t);
                linearEnsure(tt > 0, "Flow BiCGStab null update");
                omega = linearProduct(t, s) / tt;
                for (std::size_t i = 0; i < x.size(); ++i)
                    x[i] += alpha * z[i] + omega * zs[i];
            }
            apply(x, ax);
            for (std::size_t i = 0; i < x.size(); ++i) r[i] = rhs[i] - ax[i];
            if (smallResidual(r)) return step + 1;
            if (step % 40 == 39) {
                r0 = r;
                std::fill(p.begin(), p.end(), 0.);
                std::fill(v.begin(), v.end(), 0.);
                rhoOld = alpha = omega = 1;
            } else rhoOld = rho;
        }
        throw std::runtime_error("Flow linear solver iteration limit reached");
    }
};

} // namespace cartmesh2d::fv::detail
