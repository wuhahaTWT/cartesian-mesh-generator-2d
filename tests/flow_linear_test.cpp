#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using cartmesh2d::fv::detail::LinearWorkspace2D;
using cartmesh2d::fv::detail::LinearPressureMethod2D;
using cartmesh2d::fv::detail::linearNorm;
using cartmesh2d::fv::detail::SparsePattern2D;
using cartmesh2d::fv::detail::SparseSystem2D;

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

template<class F>
void rejects(F&& function, const std::string& message) {
    try {
        function();
        check(false, message);
    } catch (const std::exception&) {
        check(true, message);
    }
}

using Dense = std::vector<std::vector<double>>;

Dense zeroMatrix(std::size_t n) {
    return Dense(n, std::vector<double>(n, 0.0));
}

std::vector<double> denseApply(const Dense& matrix, const std::vector<double>& x) {
    std::vector<double> result(matrix.size(), 0.0);
    for (std::size_t row = 0; row < matrix.size(); ++row)
        for (std::size_t col = 0; col < matrix.size(); ++col)
            result[row] += matrix[row][col] * x[col];
    return result;
}

double norm2(const std::vector<double>& values) {
    double sum = 0.0;
    for (const double value : values) sum += value * value;
    return std::sqrt(sum);
}

double residualNorm(const Dense& matrix, const std::vector<double>& x,
                    const std::vector<double>& rhs) {
    const auto actual = denseApply(matrix, x);
    double sum = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const double error = actual[i] - rhs[i];
        sum += error * error;
    }
    return std::sqrt(sum);
}

double errorNorm(const std::vector<double>& actual,
                 const std::vector<double>& expected) {
    check(actual.size() == expected.size(), "known-solution vector sizes agree");
    double sum = 0.0;
    for (std::size_t i = 0; i < actual.size() && i < expected.size(); ++i) {
        const double error = actual[i] - expected[i];
        sum += error * error;
    }
    return std::sqrt(sum);
}

double solverTolerance(const std::vector<double>& rhs) {
    // The implementation's true-residual target is 1e-13 + 1e-11 ||rhs||.
    // Keep only a 5% margin for an independently evaluated dense product.
    return 1.05 * (1e-13 + 1e-11 * norm2(rhs));
}

double solutionTolerance(const std::vector<double>& exact) {
    return 5e-11 * (1.0 + norm2(exact));
}

void linearNormRegression() {
    check(linearNorm({}) == 0.0, "linear norm of an empty vector is zero");
    check(linearNorm({0.0, -0.0, 0.0}) == 0.0,
          "linear norm of a zero vector is zero");

    const auto checkScaled345 = [](double scale, const std::string& label) {
        const double actual = linearNorm({3.0 * scale, 4.0 * scale});
        const double expected = 5.0 * scale;
        check(std::abs(actual / expected - 1.) <= 8*std::numeric_limits<double>::epsilon(),
              label + " preserves the analytic 3-4-5 norm within floating-point rounding");
    };
    checkScaled345(1e200, "large-scale");
    checkScaled345(1e-200, "tiny-scale");

    check(linearNorm({std::numeric_limits<double>::denorm_min(), 0.0}) ==
              std::numeric_limits<double>::denorm_min(),
          "representable subnormal norm remains representable");
    check(linearNorm({1e200, 1e-200, std::numeric_limits<double>::denorm_min()}) ==
              1e200,
          "mixed dynamic range is dominated by the large finite component");

    const double nearLimitComponent =
        (std::numeric_limits<double>::max() * 0.999) / std::sqrt(2.0);
    const double nearLimitExpected = std::sqrt(2.0) * nearLimitComponent;
    check(std::isfinite(nearLimitExpected) &&
              std::abs(linearNorm({nearLimitComponent, nearLimitComponent}) / nearLimitExpected - 1.)
                  <= 8*std::numeric_limits<double>::epsilon(),
          "finite norm near the double limit remains finite and accurate");

    rejects([] {
        (void)linearNorm({std::numeric_limits<double>::max(),
                          std::numeric_limits<double>::max()});
    }, "mathematically overflowing norm is rejected");
    rejects([] { (void)linearNorm({std::numeric_limits<double>::quiet_NaN()}); },
            "NaN norm input is rejected");
    rejects([] { (void)linearNorm({std::numeric_limits<double>::infinity()}); },
            "infinite norm input is rejected");
}

void checkKnownSolution(const Dense& dense, const std::vector<double>& exact,
                        const std::vector<double>& rhs,
                        const std::vector<double>& actual,
                        const std::string& label) {
    check(residualNorm(dense, actual, rhs) <= solverTolerance(rhs),
          label + " reaches the independently computed true-residual target");
    check(errorNorm(actual, exact) <= solutionTolerance(exact),
          label + " matches the independently supplied known solution");
}

void addEntry(SparseSystem2D& system, Dense& dense, std::size_t row,
              std::size_t col, double value) {
    if (row == col)
        system.diag[row] += value;
    else
        system.add(row, col, value);
    dense[row][col] += value;
}

void addSymmetricEdge(SparseSystem2D& system, Dense& dense, std::size_t a,
                      std::size_t b, double value) {
    addEntry(system, dense, a, b, value);
    addEntry(system, dense, b, a, value);
}

void setRhs(SparseSystem2D& system, const Dense& dense,
            const std::vector<double>& exact) {
    const auto rhs = denseApply(dense, exact);
    system.rhs = rhs;
}

void tridiagonalPressureRegression() {
    constexpr std::size_t n = 9;
    std::vector<std::pair<std::size_t, std::size_t>> connections;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        connections.push_back({i, i + 1});
        connections.push_back({i + 1, i});
        connections.push_back({i, i + 1});
    }
    SparsePattern2D pattern(n, connections);
    SparseSystem2D system(pattern);
    Dense dense = zeroMatrix(n);
    for (std::size_t i = 0; i < n; ++i) addEntry(system, dense, i, i, 4.0);
    for (std::size_t i = 0; i + 1 < n; ++i)
        addSymmetricEdge(system, dense, i, i + 1, -1.0);

    const std::vector<double> exact{0.2, -0.4, 0.9, 1.3, -0.7, 0.1, 0.8, -1.1, 0.5};
    setRhs(system, dense, exact);
    LinearWorkspace2D workspace(n);

    std::vector<double> x(n, 0.0);
    const auto ic0Steps = system.solvePressure(x, workspace, true);
    check(ic0Steps > 0, "IC(0) PCG reports actual iterations on nonzero RHS");
    checkKnownSolution(dense, exact, system.rhs, x, "IC(0) tridiagonal PCG");

    x.assign(n, 0.0);
    const auto jacobiSteps = system.solvePressure(x, workspace, false);
    check(jacobiSteps > 0, "Jacobi PCG reports actual iterations on nonzero RHS");
    checkKnownSolution(dense, exact, system.rhs, x, "Jacobi tridiagonal PCG");
}

void irregularGraphRegression() {
    constexpr std::size_t n = 8;
    const std::vector<std::pair<std::size_t, std::size_t>> edges{
        {0, 1}, {1, 0}, {0, 2}, {2, 0}, {0, 2}, {1, 2}, {2, 1},
        {1, 4}, {4, 1}, {2, 3}, {3, 2}, {2, 5}, {5, 2}, {3, 5},
        {5, 3}, {3, 6}, {6, 3}, {4, 5}, {5, 4}, {4, 7}, {7, 4},
        {5, 6}, {6, 5}, {6, 7}, {7, 6}, {5, 6}, {6, 7}};
    SparsePattern2D pattern(n, edges);
    SparseSystem2D system(pattern);
    Dense dense = zeroMatrix(n);
    const auto edge = [&](std::size_t a, std::size_t b, double weight) {
        addEntry(system, dense, a, a, weight);
        addEntry(system, dense, b, b, weight);
        addSymmetricEdge(system, dense, a, b, -weight);
    };
    edge(0, 1, 1.0);
    edge(0, 2, 0.7);
    edge(1, 2, 1.3);
    edge(1, 4, 0.8);
    edge(2, 3, 1.1);
    edge(2, 5, 0.6);
    edge(3, 5, 1.4);
    edge(3, 6, 0.9);
    edge(4, 5, 1.2);
    edge(4, 7, 0.5);
    edge(5, 6, 0.75);
    edge(6, 7, 1.05);
    for (std::size_t i = 0; i < n; ++i) addEntry(system, dense, i, i, 0.35);

    const std::vector<double> exact{0.5, -0.1, 0.8, -0.4, 1.2, 0.3, -0.7, 0.6};
    setRhs(system, dense, exact);
    LinearWorkspace2D workspace(n);
    std::vector<double> x(n, 0.0);
    const auto ic0Steps = system.solvePressure(x, workspace, true);
    check(ic0Steps > 0, "IC(0) reports iterations on the irregular graph");
    checkKnownSolution(dense, exact, system.rhs, x,
                       "IC(0) irregular Laplacian PCG");
    const auto ic0Solution = x;
    x.assign(n, 0.0);
    const auto jacobiSteps = system.solvePressure(x, workspace, false);
    check(jacobiSteps > 0, "Jacobi reports iterations on the irregular graph");
    checkKnownSolution(dense, exact, system.rhs, x,
                       "Jacobi irregular Laplacian PCG");
    check(errorNorm(x, ic0Solution) <= 2.0 * solutionTolerance(exact),
          "Jacobi and IC(0) agree on the irregular graph solution");
}

void nonsymmetricBiCGRegression() {
    constexpr std::size_t n = 7;
    std::vector<std::pair<std::size_t, std::size_t>> connections;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        connections.push_back({i, i + 1});
        connections.push_back({i + 1, i});
        connections.push_back({i, i + 1});
    }
    SparsePattern2D pattern(n, connections);
    SparseSystem2D system(pattern);
    Dense dense = zeroMatrix(n);
    for (std::size_t i = 0; i < n; ++i) addEntry(system, dense, i, i, 3.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        addEntry(system, dense, i, i + 1, -1.2);
        addEntry(system, dense, i + 1, i, -0.7);
    }
    const std::vector<double> exact{1.1, -0.8, 0.4, 1.7, -0.2, 0.9, -1.3};
    setRhs(system, dense, exact);
    LinearWorkspace2D workspace(n);
    std::vector<double> x(n, 0.0);
    const auto steps = system.solve(x, workspace);
    check(steps > 0, "BiCGStab reports actual iterations for nonsymmetric system");
    checkKnownSolution(dense, exact, system.rhs, x, "BiCGStab convection matrix");
    x.assign(n, 0.0);
    rejects([&] { (void)system.solvePressure(x, workspace, true); },
            "IC(0) rejects the nonsymmetric convection matrix");
}

void compensatedResidualRegression() {
    SparsePattern2D single(1,{});SparseSystem2D product(single);
    const double e=std::ldexp(1.,-27);
    product.diag={1+e};product.rhs={1};
    // Exact binary arithmetic: 1-(1+2^-27)(1-2^-27)=2^-54.
    check(product.compensatedResidualRow(0,{1-e})==std::ldexp(1.,-54),
          "compensated residual retains rounded product error");
    rejects([&] {
        LinearWorkspace2D workspace(1);std::vector<double> x={1-e};
        (void)product.solve(x,workspace,1e-18);
    },"unrepresentable row accuracy must not pass through a rounded zero");
    SparsePattern2D pair(2,{{0,1}});SparseSystem2D cancellation(pair);
    cancellation.diag={1e16,1};cancellation.rhs={1e16,1};cancellation.add(0,1,1);
    check(cancellation.compensatedResidualRow(0,{1,1})==-1,
          "compensated residual does not turn cancellation into false zero");
    check(cancellation.compensatedResidualRow(1,{1,1})==0,
          "compensated residual retains exact zero");
}

void biorthogonalBreakdownRegression() {
    // Three one-way coupled control volumes: A is nonsingular and positive
    // diagonal. After the first update, r=(0,.5,.5) is exactly orthogonal to
    // the original shadow (1,0,0); the old solver aborted despite x != (1,1,1).
    SparsePattern2D pattern(3,{{0,1},{1,2}});
    SparseSystem2D system(pattern);Dense dense=zeroMatrix(3);
    for(std::size_t i=0;i<3;++i)addEntry(system,dense,i,i,1.);
    addEntry(system,dense,1,0,-1.);addEntry(system,dense,2,1,-1.);
    const std::vector<double> exact{1,1,1};setRhs(system,dense,exact);
    std::vector<double> x(3,0.);LinearWorkspace2D workspace(3);
    const auto steps=system.solve(x,workspace);
    check(steps>1,"breakdown recovery reports actual work");
    checkKnownSolution(dense,exact,system.rhs,x,"BiCGStab orthogonal shadow restart");
}

void structureAndPinRegression() {
    const std::vector<std::pair<std::size_t, std::size_t>> connections{
        {0, 1}, {1, 0}, {0, 1}, {0, 1}, {1, 2}, {2, 1}, {2, 1}};
    SparsePattern2D pattern(3, connections);
    SparseSystem2D system(pattern);
    Dense dense = zeroMatrix(3);
    check(system.off.size() == 4,
          "duplicate undirected connections occupy one off entry per directed edge");
    for (std::size_t i = 0; i < 3; ++i) addEntry(system, dense, i, i, 2.5);
    addSymmetricEdge(system, dense, 0, 1, -0.5);
    addSymmetricEdge(system, dense, 1, 2, -0.8);
    addEntry(system, dense, 0, 1, -0.25);
    addEntry(system, dense, 1, 0, -0.25);
    check(std::abs(system.diag[0] - 2.5) < 1e-14 &&
              std::abs(system.diag[1] - 2.5) < 1e-14,
          "diagonal assembly remains separate from off-diagonal storage");
    std::vector<double> assembledApply(3, 0.0);
    system.apply({1.0, 2.0, 3.0}, assembledApply);
    check(residualNorm(dense, {1.0, 2.0, 3.0}, assembledApply) <= 1e-12,
          "repeated off-diagonal add calls accumulate in the unique storage slot");

    const std::vector<double> exact{0.0, 1.25, -0.75};
    setRhs(system, dense, exact);
    system.pin(0);
    dense[0][0] = 2.5;
    for (std::size_t i = 1; i < dense.size(); ++i) {
        dense[i][0] = 0.0;
        dense[0][i] = 0.0;
    }
    check(std::abs(system.rhs[0]) <= 1e-15,
          "pin clears the constrained RHS entry");
    LinearWorkspace2D workspace(3);
    std::vector<double> x(3, 0.0);
    const auto steps = system.solvePressure(x, workspace, true);
    check(steps > 0, "pinned PCG reports actual iterations");
    checkKnownSolution(dense, exact, system.rhs, x, "symmetric pinned PCG");
    check(std::abs(x[0]) <= 1e-12, "pinned unknown is zero in the solved state");

    std::vector<double> duplicateApply(3, 0.0);
    system.apply({1.0, 2.0, 3.0}, duplicateApply);
    check(duplicateApply.size() == 3 &&
              residualNorm(dense, {1.0, 2.0, 3.0}, duplicateApply) <= 1e-12,
          "apply agrees with the independently assembled pinned matrix");
}

void zeroResetAndPivotRegression() {
    SparsePattern2D pattern(2, {{0, 1}, {1, 0}, {0, 1}});
    SparseSystem2D system(pattern);
    Dense dense = zeroMatrix(2);
    addEntry(system, dense, 0, 0, 2.0);
    addEntry(system, dense, 1, 1, 3.0);
    addSymmetricEdge(system, dense, 0, 1, -0.4);
    system.rhs = {0.0, 0.0};
    LinearWorkspace2D workspace(2);
    std::vector<double> x{0.0, 0.0};
    const auto zeroSteps = system.solvePressure(x, workspace, true);
    check(zeroSteps == 0 && x == std::vector<double>({0.0, 0.0}),
          "zero RHS converges in zero iterations from the exact zero guess");

    const std::vector<double> exactA{1.0, -2.0};
    setRhs(system, dense, exactA);
    x.assign(2, 0.0);
    (void)system.solvePressure(x, workspace, true);
    checkKnownSolution(dense, exactA, system.rhs, x, "workspace first solve");
    system.reset();
    dense = zeroMatrix(2);
    addEntry(system, dense, 0, 0, 4.0);
    addEntry(system, dense, 1, 1, 5.0);
    addSymmetricEdge(system, dense, 0, 1, -0.2);
    const std::vector<double> exactB{-0.7, 0.6};
    setRhs(system, dense, exactB);
    x.assign(2, 0.0);
    (void)system.solvePressure(x, workspace, true);
    checkKnownSolution(dense, exactB, system.rhs, x,
                       "reset workspace second solve");

    SparsePattern2D badPattern(2, {{0, 1}});
    SparseSystem2D bad(badPattern);
    bad.diag[0] = 1.0;
    bad.diag[1] = 0.5;
    bad.add(0, 1, 1.0);
    bad.add(1, 0, 1.0);
    bad.rhs[0] = 1.0;
    std::vector<double> badX(2, 0.0);
    LinearWorkspace2D badWorkspace(2);
    rejects([&] { (void)bad.solvePressure(badX, badWorkspace, true); },
            "IC(0) rejects a non-positive pivot explicitly");
}

void ic0CacheRegression() {
    const std::vector<std::pair<std::size_t, std::size_t>> connections{
        {0, 1}, {1, 0}, {1, 2}, {2, 1}};
    SparsePattern2D pattern(3, connections);
    SparseSystem2D system(pattern);
    Dense dense = zeroMatrix(3);
    addEntry(system, dense, 0, 0, 3.0);
    addEntry(system, dense, 1, 1, 4.0);
    addEntry(system, dense, 2, 2, 5.0);
    addSymmetricEdge(system, dense, 0, 1, -0.5);
    addSymmetricEdge(system, dense, 1, 2, -0.7);
    LinearWorkspace2D workspace(3);
    std::vector<double> x(3, 0.0);

    const std::vector<double> exactA{1.2, -0.4, 0.8};
    setRhs(system, dense, exactA);
    const auto firstSteps = system.solvePressure(x, workspace, true);
    check(firstSteps > 0, "IC(0) cache first solve iterates");
    check(system.ic0Builds() == 1 && system.ic0Reuses() == 0,
          "first IC(0) solve builds exactly one factorization");
    checkKnownSolution(dense, exactA, system.rhs, x,
                       "IC(0) cache first known solution");

    // A new RHS must reuse the matrix factorization and still solve the new
    // independent system, rather than returning a cached solution.
    const std::vector<double> exactB{-0.3, 1.1, -0.6};
    setRhs(system, dense, exactB);
    x.assign(3, 0.0);
    const auto secondSteps = system.solvePressure(x, workspace, true);
    check(secondSteps > 0, "IC(0) cache changed-RHS solve iterates");
    check(system.ic0Builds() == 1 && system.ic0Reuses() == 1,
          "same matrix changed RHS reuses the cached factorization");
    checkKnownSolution(dense, exactB, system.rhs, x,
                       "IC(0) cache changed-RHS known solution");

    // Coefficients are public for assembly, so direct writes must be covered
    // by exact snapshots as well as the explicit add/reset invalidators.
    system.diag[1] += 0.6;
    dense[1][1] += 0.6;
    const std::vector<double> exactC{0.4, -0.9, 1.3};
    setRhs(system, dense, exactC);
    x.assign(3, 0.0);
    (void)system.solvePressure(x, workspace, true);
    check(system.ic0Builds() == 2,
          "direct diagonal mutation rebuilds IC(0) exactly once");
    checkKnownSolution(dense, exactC, system.rhs, x,
                       "IC(0) direct diagonal mutation solution");

    const auto edge01 = pattern.slot(0, 1);
    const auto edge10 = pattern.slot(1, 0);
    system.off[edge01] -= 0.1;
    system.off[edge10] -= 0.1;
    dense[0][1] -= 0.1;
    dense[1][0] -= 0.1;
    const std::vector<double> exactD{-0.8, 0.2, 0.7};
    setRhs(system, dense, exactD);
    x.assign(3, 0.0);
    (void)system.solvePressure(x, workspace, true);
    check(system.ic0Builds() == 3,
          "direct off-diagonal mutation rebuilds IC(0) exactly once");
    checkKnownSolution(dense, exactD, system.rhs, x,
                       "IC(0) direct off-diagonal mutation solution");

    // A different system may share the Krylov workspace without sharing
    // factors. Solving it must not corrupt the original system's cache.
    SparseSystem2D other(pattern);
    Dense otherDense = zeroMatrix(3);
    addEntry(other, otherDense, 0, 0, 6.0);
    addEntry(other, otherDense, 1, 1, 7.0);
    addEntry(other, otherDense, 2, 2, 8.0);
    addSymmetricEdge(other, otherDense, 0, 1, -0.4);
    addSymmetricEdge(other, otherDense, 1, 2, -0.9);
    const std::vector<double> otherExact{0.6, -1.4, 0.3};
    setRhs(other, otherDense, otherExact);
    x.assign(3, 0.0);
    (void)other.solvePressure(x, workspace, true);
    check(other.ic0Builds() == 1 && other.ic0Reuses() == 0,
          "independent system owns an independent IC(0) factorization");
    const std::vector<double> exactAfterOther{0.9, 0.5, -1.1};
    setRhs(system, dense, exactAfterOther);
    x.assign(3, 0.0);
    (void)system.solvePressure(x, workspace, true);
    check(system.ic0Builds() == 3 && system.ic0Reuses() == 2,
          "shared workspace does not replace the original cached factors");
    checkKnownSolution(dense, exactAfterOther, system.rhs, x,
                       "original system survives another-system solve");

    // reset and pin both invalidate a previously built factorization.
    system.reset();
    dense = zeroMatrix(3);
    addEntry(system, dense, 0, 0, 4.0);
    addEntry(system, dense, 1, 1, 5.0);
    addEntry(system, dense, 2, 2, 6.0);
    addSymmetricEdge(system, dense, 0, 1, -0.3);
    addSymmetricEdge(system, dense, 1, 2, -0.6);
    const std::vector<double> exactReset{0.7, -0.2, 0.9};
    setRhs(system, dense, exactReset);
    x.assign(3, 0.0);
    (void)system.solvePressure(x, workspace, true);
    check(system.ic0Builds() == 4,
          "reset invalidates and rebuilds the IC(0) factorization");
    checkKnownSolution(dense, exactReset, system.rhs, x,
                       "reset matrix known solution");

    system.pin(0);
    dense[0][1] = dense[1][0] = 0.0;
    const std::vector<double> exactPin{0.0, -0.4, 1.2};
    setRhs(system, dense, exactPin);
    system.rhs[0] = 0.0;
    x.assign(3, 0.0);
    (void)system.solvePressure(x, workspace, true);
    check(system.ic0Builds() == 5,
          "pin invalidates and rebuilds the IC(0) factorization");
    checkKnownSolution(dense, exactPin, system.rhs, x,
                       "pinned matrix known solution after invalidation");

    // Fail after entering factorization (asymmetric direct off write), then
    // restore the matrix. The failed attempt must not leave stale factors
    // usable for the repaired matrix.
    const auto edge12 = pattern.slot(1, 2);
    system.off[edge12] += 0.05;
    const auto buildsBeforeFailure = system.ic0Builds();
    rejects([&] {
        std::vector<double> failedX(3, 0.0);
        (void)system.solvePressure(failedX, workspace, true);
    }, "failed IC(0) factorization is reported");
    check(system.ic0Builds() == buildsBeforeFailure,
          "failed IC(0) factorization is not counted as a build");
    system.off[edge12] -= 0.05;
    x.assign(3, 0.0);
    (void)system.solvePressure(x, workspace, true);
    check(system.ic0Builds() == buildsBeforeFailure + 1,
          "repaired matrix rebuilds after failed IC(0) factorization");
    checkKnownSolution(dense, exactPin, system.rhs, x,
                       "repaired matrix does not use stale failed factors");

    // Neither a Jacobi solve nor an exact zero-residual solve should create or
    // consume IC(0) factorization statistics.
    const auto buildsBeforeNoIc0 = system.ic0Builds();
    const auto reusesBeforeNoIc0 = system.ic0Reuses();
    x.assign(3, 0.0);
    (void)system.solvePressure(x, workspace, false);
    check(system.ic0Builds() == buildsBeforeNoIc0 &&
              system.ic0Reuses() == reusesBeforeNoIc0,
          "Jacobi pressure solve does not affect IC(0) counters");
    system.rhs.assign(3, 0.0);
    x.assign(3, 0.0);
    check(system.solvePressure(x, workspace, true) == 0,
          "zero-residual pressure solve skips Krylov iterations");
    check(system.ic0Builds() == buildsBeforeNoIc0 &&
              system.ic0Reuses() == reusesBeforeNoIc0,
          "zero-residual pressure solve does not fake IC(0) use");

    // A zero residual is a completed solve for every pressure preconditioner:
    // the caller's already-valid guess must remain byte-identical and no
    // preconditioner work should be inferred from the zero iteration count.
    for (const auto method : {LinearPressureMethod2D::Jacobi,
                              LinearPressureMethod2D::IC0,
                              LinearPressureMethod2D::Aggregation}) {
        SparsePattern2D diagonalPattern(2, {});
        SparseSystem2D diagonalSystem(diagonalPattern);
        diagonalSystem.diag={2.,3.}; diagonalSystem.rhs={.6,-.9};
        LinearWorkspace2D diagonalWorkspace(2);
        std::vector<double> guess{.3,-.3};
        const auto before=guess;
        check(diagonalSystem.solvePressure(guess,diagonalWorkspace,method)==0,
              "zero-residual pressure solve reports zero iterations for each method");
        check(guess==before,
              "zero-residual pressure solve preserves the caller guess for each method");
    }
}

void localResidualScaleRegression() {
    SparsePattern2D pattern(2, {});
    SparseSystem2D system(pattern);
    LinearWorkspace2D workspace(2);
    system.diag={1e6,1e-4}; system.rhs={1e6,1e-10};
    std::vector<double> x{1.,0.};
    check(system.solve(x,workspace)==0,"global RHS tolerance can mask a small-cell residual");
    check(system.solve(x,workspace,1e-10)>0,"diagonal-scaled residual check triggers a real correction");
    check(std::abs(x[1]-1e-6)<1e-15,"small-cell solution recovered without relaxing nonlinear tolerance");
    rejects([&] {system.solve(x,workspace,0);},"nonpositive scaled residual target rejected");
}

void tightResidualNormRegression() {
    SparsePattern2D pattern(4, {{0,1},{1,0},{1,2},{2,1},{2,3},{3,2}});
    SparseSystem2D system(pattern);
    Dense dense = zeroMatrix(4);
    for (std::size_t i = 0; i < 4; ++i) addEntry(system, dense, i, i, 3.0);
    for (std::size_t i = 0; i + 1 < 4; ++i) addSymmetricEdge(system, dense, i, i + 1, -.9);
    const std::vector<double> exact{.7, -1.1, .4, 1.3};
    setRhs(system, dense, exact);
    LinearWorkspace2D workspace(4);
    std::vector<double> defaultGuess(4, 0.0);
    const auto defaultSteps = system.solve(defaultGuess, workspace);
    std::vector<double> tightGuess(4, 0.0);
    const auto tightSteps = system.solve(tightGuess, workspace, std::numeric_limits<double>::infinity(), 1e-15);
    check(tightSteps >= defaultSteps, "explicit tight residual norm does not stop earlier than default");
    check(residualNorm(dense, tightGuess, system.rhs) <= 1.1e-15,
          "explicit tight residual norm is independently reached");
    check(errorNorm(tightGuess, exact) <= 2e-14,
          "tight residual norm retains the correct linear solution");
    auto nearGuess=exact; nearGuess[0]+=1e-12;
    check(system.solve(nearGuess,workspace)==0,
          "legacy norm accepts an already close guess");
    check(system.solve(nearGuess,workspace,std::numeric_limits<double>::infinity(),1e-15)>0,
          "tight norm corrects a guess accepted by the legacy norm");
    check(residualNorm(dense,nearGuess,system.rhs)<=1.1e-15,
          "warm-start true residual meets the explicitly tight norm");
    rejects([&] { std::vector<double> guess(4); (void)system.solve(guess, workspace, std::numeric_limits<double>::infinity(), 0.); },
            "invalid linear residual norm tolerance");
    rejects([&] { std::vector<double> guess(4); (void)system.solve(guess, workspace, std::numeric_limits<double>::infinity(), -.1); },
            "invalid linear residual norm tolerance");
    rejects([&] { std::vector<double> guess(4); (void)system.solve(guess, workspace, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()); },
            "invalid linear residual norm tolerance");
}

} // namespace

int main() {
    try {
        linearNormRegression();
        tridiagonalPressureRegression();
        irregularGraphRegression();
        nonsymmetricBiCGRegression();
        biorthogonalBreakdownRegression();
        compensatedResidualRegression();
        structureAndPinRegression();
        zeroResetAndPivotRegression();
        ic0CacheRegression();
        localResidualScaleRegression();
        tightResidualNormRegression();
    } catch (const std::exception& error) {
        std::cerr << "UNEXPECTED EXCEPTION: " << error.what() << '\n';
        return 1;
    }
    std::cout << "flow linear regression failures: " << failures << '\n';
    return failures == 0 ? 0 : 1;
}
