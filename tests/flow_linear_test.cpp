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
using cartmesh2d::fv::detail::LinearSolveMethod2D;
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

std::vector<double> denseMaskedILU0Solve(
    const Dense& matrix, const std::vector<std::pair<std::size_t, std::size_t>>& pattern,
    const std::vector<double>& rhs) {
    const std::size_t n = matrix.size();
    Dense lower = zeroMatrix(n);
    Dense upper = zeroMatrix(n);
    std::vector<std::vector<bool>> present(n, std::vector<bool>(n, false));
    for (const auto [row, col] : pattern) present[row][col] = true;
    for (std::size_t i = 0; i < n; ++i) {
        lower[i][i] = 1.0;
        for (std::size_t j = 0; j <= i; ++j) {
            if (j < i && !present[i][j]) continue;
            double value = matrix[i][j];
            for (std::size_t k = 0; k < j; ++k)
                value -= lower[i][k] * upper[k][j];
            if (j == i) upper[i][i] = value;
            else lower[i][j] = value / upper[j][j];
        }
        for (std::size_t j = i + 1; j < n; ++j) {
            if (!present[i][j]) continue;
            double value = matrix[i][j];
            for (std::size_t k = 0; k < i; ++k)
                value -= lower[i][k] * upper[k][j];
            upper[i][j] = value;
        }
    }
    std::vector<double> y(n, 0.0), result(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        y[i] = rhs[i];
        for (std::size_t j = 0; j < i; ++j) y[i] -= lower[i][j] * y[j];
    }
    for (std::size_t i = n; i-- > 0;) {
        result[i] = y[i];
        for (std::size_t j = i + 1; j < n; ++j) result[i] -= upper[i][j] * result[j];
        result[i] /= upper[i][i];
    }
    return result;
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
    LinearWorkspace2D candidateWorkspace(1);
    const auto candidate=product.solveCandidate({1-e},candidateWorkspace,1e-25,1e-25);
    check(candidate.high[0]==1-e && candidate.low.size()==1 && candidate.low[0]>0,
          "accurate candidate exposes low bits instead of silently discarding them");
    // Independent closed-form residual: d*high=1-2^-54 exactly for this input.
    const double expandedResidual=std::ldexp(1.,-54)-(1+e)*candidate.low[0];
    check(std::abs(expandedResidual)<1e-25,
          "high+low candidate satisfies the original equation at the requested gate");
    const auto late=product.solveCandidate({0},candidateWorkspace,1e-25,1e-25);
    check(late.high[0]==1-e && late.low.size()==1 && late.low[0]>0 &&
          std::abs(std::ldexp(1.,-54)-(1+e)*late.low[0])<1e-25,
          "late expansion activation recovers lost bits from the recomputed residual");
    const double rounded=candidate.relaxedDouble(0,0,1);
    check(rounded==1-e && std::ldexp(1.,-54)/(1+e)>1e-18,
          "rounded field demonstrably fails the gate satisfied by the candidate");
    // Exact binary tie: premature rounding to 1 loses the positive low part.
    cartmesh2d::fv::detail::LinearCandidate2D tie{{1},{std::ldexp(1.,-54)},0};
    const double previous=1+std::ldexp(1.,-52);
    check(tie.relaxedDouble(0,previous,.5)==previous,
          "relaxation retains low bits until the final field rounding");
    const auto ordinary=product.solveCandidate({1-e},candidateWorkspace,1e-10,1e-10);
    check(ordinary.low.empty(),"ordinary accuracy does not allocate an expansion");
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

void ilu0SelectionAndCacheRegression() {
    // This deliberately has unequal upper/lower coefficients.  It must be
    // exercised as a genuinely nonsymmetric factorization, rather than being
    // accepted by the symmetric IC(0) path after silently symmetrizing it.
    constexpr std::size_t n = 4;
    SparsePattern2D pattern(n, {{0,1},{1,0},{1,2},{2,1},{2,3},{3,2}});
    SparseSystem2D system(pattern);
    Dense dense = zeroMatrix(n);
    for (std::size_t i = 0; i < n; ++i) addEntry(system, dense, i, i, 4.0);
    addEntry(system, dense, 0, 1, -1.25);
    addEntry(system, dense, 1, 0, -0.35);
    addEntry(system, dense, 1, 2, -1.10);
    addEntry(system, dense, 2, 1, -0.55);
    addEntry(system, dense, 2, 3, -0.90);
    addEntry(system, dense, 3, 2, -0.25);

    const std::vector<double> exactA{1.2, -0.7, 0.35, 1.1};
    setRhs(system, dense, exactA);
    LinearWorkspace2D workspace(n);
    std::vector<double> x(n, 0.0);
    const auto steps = system.solve(x, workspace,
                                   std::numeric_limits<double>::infinity(),
                                   std::numeric_limits<double>::infinity(),
                                   LinearSolveMethod2D::ILU0);
    check(steps > 0, "ILU(0) reports work on a nonsymmetric known system");
    checkKnownSolution(dense, exactA, system.rhs, x,
                       "ILU(0) nonsymmetric known solution");
    check(system.ilu0Builds() == 1 && system.ilu0Reuses() == 0,
          "first ILU(0) solve builds one factorization");

    // The direct preconditioner API is independently checked on a vector;
    // this catches a solver that happens to converge while applying the wrong
    // triangular factors.
    std::vector<double> preconditioned(n, 0.0);
    const std::vector<double> probe{1.0, -2.0, 0.5, 3.0};
    system.preconditionILU0(probe, preconditioned);
    for (double value : preconditioned)
        check(std::isfinite(value), "ILU(0) preconditioner returns finite values");
    check(residualNorm(dense, preconditioned, probe) <= 2e-14,
          "tridiagonal ILU(0) preconditioner is an exact dense solve");

    // The next pattern contains a ring and a triangle but omits the fill
    // edge (1,3).  Compare against an independently coded dense masked
    // Doolittle factorization, so a preconditioner that accidentally keeps
    // fill or uses the wrong triangular order is detected.
    const std::vector<std::pair<std::size_t, std::size_t>> maskedConnections{
        {0,1},{1,0},{1,2},{2,1},{2,3},{3,2},{0,2},{2,0},{0,3},{3,0}};
    SparsePattern2D maskedPattern(n, maskedConnections);
    SparseSystem2D masked(maskedPattern);
    Dense maskedDense = zeroMatrix(n);
    const std::vector<double> maskedDiagonal{5.,4.,3.,6.};
    for (std::size_t i = 0; i < n; ++i) addEntry(masked, maskedDense, i, i,maskedDiagonal[i]);
    const auto addMasked = [&](std::size_t row, std::size_t col, double value) {
        addEntry(masked, maskedDense, row, col, value);
    };
    addMasked(0,1,1.0); addMasked(1,0,-0.5);
    addMasked(1,2,1.2); addMasked(2,1,-0.4);
    addMasked(2,3,0.8); addMasked(3,2,-0.3);
    addMasked(0,2,-0.6); addMasked(2,0,0.25);
    addMasked(0,3,0.7); addMasked(3,0,-0.2);
    LinearWorkspace2D maskedWorkspace(n);
    masked.factorILU0();
    const std::vector<double> maskedProbe{0.8, -1.1, 2.0, 0.35};
    std::vector<double> maskedActual(n, 0.0);
    masked.preconditionILU0(maskedProbe, maskedActual);
    const auto maskedExpected = denseMaskedILU0Solve(maskedDense, maskedConnections,
                                                     maskedProbe);
    check(errorNorm(maskedActual, maskedExpected) <= 2e-13,
          "masked nonsymmetric ILU(0) matches independent Doolittle factors");
    check(residualNorm(maskedDense, maskedActual, maskedProbe) > 1e-8,
          "masked ILU(0) exposes the intentionally dropped fill in dense residual");

    const std::vector<double> exactB{-0.4, 0.8, -1.3, 0.6};
    setRhs(system, dense, exactB);
    x.assign(n, 0.0);
    (void)system.solve(x, workspace,
                       std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(),
                       LinearSolveMethod2D::ILU0);
    check(system.ilu0Builds() == 1 && system.ilu0Reuses() == 1,
          "changed RHS reuses the exact ILU(0) coefficient snapshot");
    checkKnownSolution(dense, exactB, system.rhs, x,
                       "ILU(0) changed-RHS known solution");

    // Public assembly arrays are intentionally mutable.  A changed diagonal
    // or off coefficient must invalidate the factorization independently of
    // the RHS and solve the new matrix.
    system.diag[1] += 0.7;
    dense[1][1] += 0.7;
    const auto edge = pattern.slot(1, 2);
    system.off[edge] += 0.13;
    dense[1][2] += 0.13;
    const std::vector<double> exactC{0.3, -0.9, 1.4, -0.2};
    setRhs(system, dense, exactC);
    x.assign(n, 0.0);
    (void)system.solve(x, workspace,
                       std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(),
                       LinearSolveMethod2D::ILU0);
    check(system.ilu0Builds() == 2 && system.ilu0Reuses() == 1,
          "ILU(0) coefficient mutation rebuilds instead of reusing stale factors");
    checkKnownSolution(dense, exactC, system.rhs, x,
                       "ILU(0) changed-coefficient known solution");

    // A bad pivot must fail closed.  Repairing it must force a fresh build;
    // no partial factor from the failed attempt may survive.
    system.diag[2] = 0.0;
    rejects([&] { system.factorILU0(); },
            "ILU(0) rejects a zero pivot without shifting or fallback");
    check(system.ilu0Builds() == 2,
          "failed ILU(0) factorization is not counted as a successful build");
    system.diag[2] = dense[2][2];
    setRhs(system, dense, exactA);
    x.assign(n, 0.0);
    (void)system.solve(x, workspace,
                       std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(),
                       LinearSolveMethod2D::ILU0);
    check(system.ilu0Builds() == 3,
          "repaired ILU(0) matrix rebuilds after failed factorization");

    SparsePattern2D eliminationPattern(2, {{0,1},{1,0}});
    system.diag[0] += .25;
    rejects([&] { system.preconditionILU0(probe, preconditioned); },
            "direct ILU(0) application rejects stale factors after coefficient mutation");

    SparseSystem2D elimination(eliminationPattern);
    elimination.diag = {1.0, 1.0};
    elimination.add(0, 1, 2.0);
    elimination.add(1, 0, 1.0);
    rejects([&] { elimination.factorILU0(); },
            "ILU(0) rejects a pivot made nonpositive by elimination");
    SparseSystem2D nanSystem(eliminationPattern);
    nanSystem.diag = {1.0, std::numeric_limits<double>::quiet_NaN()};
    rejects([&] { nanSystem.factorILU0(); }, "ILU(0) rejects a NaN coefficient");
    SparseSystem2D infSystem(eliminationPattern);
    infSystem.diag = {1.0, std::numeric_limits<double>::infinity()};
    rejects([&] { infSystem.factorILU0(); }, "ILU(0) rejects an infinite coefficient");

    // Default Jacobi remains a separate method and must not create or consume
    // ILU(0) factors.  The invalid enum check also prevents accidental
    // fallback to a different preconditioner.
    const auto builds = system.ilu0Builds();
    const auto reuses = system.ilu0Reuses();
    system.rhs = denseApply(dense, exactB);
    x.assign(n, 0.0);
    (void)system.solve(x, workspace);
    check(system.ilu0Builds() == builds && system.ilu0Reuses() == reuses,
          "default Jacobi solve does not touch ILU(0) cache counters");
    rejects([&] {
        (void)system.solve(x, workspace,
                           std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity(),
                           static_cast<LinearSolveMethod2D>(99));
    }, "invalid linear solve method is rejected");
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
        ilu0SelectionAndCacheRegression();
    } catch (const std::exception& error) {
        std::cerr << "UNEXPECTED EXCEPTION: " << error.what() << '\n';
        return 1;
    }
    std::cout << "flow linear regression failures: " << failures << '\n';
    return failures == 0 ? 0 : 1;
}
