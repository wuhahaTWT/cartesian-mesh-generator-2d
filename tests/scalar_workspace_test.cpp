#include "cartmesh2d/fv/ScalarTransport2D.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {
int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

TopologyMesh2D fromPolygons(const std::vector<Polygon2D>& polygons) {
    TopologyMesh2D result;
    std::map<std::pair<double, double>, std::size_t> vertices;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edges;
    for (const auto& polygon : polygons) {
        TopologyCell2D cell;
        cell.id = result.cells.size();
        for (const auto& point : polygon.vertices) {
            const auto key = std::make_pair(point.x, point.y);
            auto [it, inserted] = vertices.emplace(key, result.vertices.size());
            if (inserted) result.vertices.push_back({it->second, point});
            cell.vertices.push_back(it->second);
        }
        cell.geometryArea = polygon.area();
        for (std::size_t i = 0; i < cell.vertices.size(); ++i) {
            const auto a = cell.vertices[i];
            const auto b = cell.vertices[(i + 1) % cell.vertices.size()];
            auto [it, inserted] = edges.emplace(std::minmax(a, b), result.edges.size());
            if (inserted) {
                result.edges.push_back({it->second, a, b, cell.id, {}, BoundaryPatch2D::DomainBoundary});
            } else {
                auto& edge = result.edges[it->second];
                edge.neighbour = cell.id;
                edge.patch = BoundaryPatch2D::None;
            }
            cell.edges.push_back(it->second);
        }
        result.cells.push_back(cell);
    }
    return result;
}

FvMesh2D grid(int nx, int ny, double shear = 0.) {
    const auto point = [=](int i, int j) {
        const double x = static_cast<double>(i) / nx;
        const double y = static_cast<double>(j) / ny;
        return Point2D{x + shear * y, y};
    };
    std::vector<Polygon2D> polygons;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            polygons.push_back({{point(i, j), point(i + 1, j), point(i + 1, j + 1), point(i, j + 1)}});
    return makeFvMesh2D(fromPolygons(polygons));
}

ScalarBoundary2D valueBoundary(double value) {
    return {ScalarBoundaryKind2D::Value, value, std::nullopt};
}

ScalarTransportControls2D tightControls(ScalarPreconditioner2D method = ScalarPreconditioner2D::Jacobi) {
    ScalarTransportControls2D controls;
    controls.maxCorrections = 1200;
    controls.relaxation = 1.;
    controls.relativeTolerance = 1e-11;
    controls.absoluteTolerance = 1e-13;
    controls.cellTolerance = 1e-11;
    controls.preconditioner = method;
    controls.profile = true;
    return controls;
}

ScalarTransportProblem2D uniformProblem(const FvMesh2D& mesh, double source = 1., double sink = .5,
                                        double boundary = 2.) {
    ScalarTransportProblem2D problem;
    problem.diffusivity = .17;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.sourceDensity.assign(mesh.cells.size(), source);
    problem.sinkRate.assign(mesh.cells.size(), sink);
    problem.boundaryData.assign(mesh.faces.size(), valueBoundary(boundary));
    return problem;
}

ScalarTransportProblem2D varyingProblem(const FvMesh2D& mesh, double boundaryShift = 0.) {
    ScalarTransportProblem2D problem;
    problem.diffusivity = .11;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.sourceDensity.resize(mesh.cells.size());
    problem.sinkRate.resize(mesh.cells.size());
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        problem.sourceDensity[i] = .35 + .27 * mesh.cells[i].centre.x + .13 * mesh.cells[i].centre.y;
        problem.sinkRate[i] = .25 + .08 * mesh.cells[i].centre.x;
    }
    problem.faceDiffusivity.resize(mesh.faces.size());
    problem.boundaryData.resize(mesh.faces.size());
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        problem.faceDiffusivity[id] = .08 + .04 * (1. + face.centre.x);
        problem.boundaryData[id] = valueBoundary(1.1 + boundaryShift + .2 * face.centre.x);
    }
    return problem;
}

void compareResults(const ScalarTransportResult2D& expected, const ScalarTransportResult2D& actual,
                    const char* message) {
    check(expected.converged && actual.converged, message);
    check(expected.values == actual.values, message);
    check(expected.advectiveFlux == actual.advectiveFlux, message);
    check(expected.diffusiveFlux == actual.diffusiveFlux, message);
    check(expected.sourceIntegrals == actual.sourceIntegrals, message);
    check(expected.temporalIntegrals == actual.temporalIntegrals, message);
    check(expected.sinkIntegrals == actual.sinkIntegrals, message);
    check(expected.history.size() == actual.history.size(), message);
    if (expected.history.size() == actual.history.size())
        for (std::size_t i = 0; i < expected.history.size(); ++i) {
            check(expected.history[i].iteration == actual.history[i].iteration, message);
            check(expected.history[i].linearIterations == actual.history[i].linearIterations, message);
            check(expected.history[i].residualNorm == actual.history[i].residualNorm, message);
            check(expected.history[i].maxCellImbalance == actual.history[i].maxCellImbalance, message);
        }
    check(expected.boundaryFlux == actual.boundaryFlux && expected.sourceIntegral == actual.sourceIntegral &&
              expected.temporalIntegral == actual.temporalIntegral && expected.globalBalance == actual.globalBalance,
          message);
}

void exactUniformSteadyAndTransient() {
    const auto mesh = grid(3, 2, .15);
    const auto problem = uniformProblem(mesh);
    auto controls = tightControls();
    ScalarTransportWorkspace2D workspace;

    const auto fresh = solveScalarTransport2D(mesh, problem, controls);
    const auto reused = solveScalarTransport2D(mesh, problem, controls, {}, 0., &workspace);
    compareResults(fresh, reused, "workspace steady solve matches fresh solve exactly");
    for(double value:reused.values)check(std::abs(value-2.)<1e-10,
          "steady source/sink problem has the analytical uniform solution");
    check(reused.performance.patternBuilds == 1, "first workspace solve builds one pattern");

    const std::vector<double> previous(mesh.cells.size(), .5);
    const double dt = .25;
    const double exact = (previous.front() / dt + 1.) / (1. / dt + .5);
    auto transientProblem=problem;
    transientProblem.boundaryData.assign(mesh.faces.size(),valueBoundary(exact));
    const auto transientFresh = solveScalarTransport2D(mesh, transientProblem, controls, previous, dt);
    const auto transientReused = solveScalarTransport2D(mesh, transientProblem, controls, previous, dt, &workspace);
    compareResults(transientFresh, transientReused, "workspace transient solve matches fresh solve exactly");
    for(double value:transientReused.values)check(std::abs(value-exact)<1e-10,
          "transient source/sink problem has the analytical backward-Euler solution");
    check(transientReused.performance.patternReuses >= 1, "same graph reuses the workspace pattern");
}

void coefficientsAndPreconditionerChangesResetAllValues() {
    const auto mesh = grid(4, 3, .1);
    auto problem = varyingProblem(mesh);
    auto controls = tightControls(ScalarPreconditioner2D::Jacobi);
    ScalarTransportWorkspace2D workspace;

    const auto jacobiFresh = solveScalarTransport2D(mesh, problem, controls);
    const auto jacobiReused = solveScalarTransport2D(mesh, problem, controls, {}, 0., &workspace);
    compareResults(jacobiFresh, jacobiReused, "Jacobi workspace solve matches fresh variable-coefficient solve");

    problem = varyingProblem(mesh, .37);
    for(double& value:problem.sourceDensity)value*=1.7;
    for(double& value:problem.sinkRate)value*=.6;
    for(double& value:problem.faceDiffusivity)value*=2.3;
    controls.preconditioner = ScalarPreconditioner2D::ILU0;
    const auto iluFresh = solveScalarTransport2D(mesh, problem, controls);
    const auto iluReused = solveScalarTransport2D(mesh, problem, controls, {}, 0., &workspace);
    compareResults(iluFresh, iluReused, "ILU0 after changed source/BC matches fresh solve");
    check(iluReused.performance.patternReuses >= 1, "coefficient changes retain only the graph pattern");

    controls.preconditioner = ScalarPreconditioner2D::Jacobi;
    const auto backJacobiFresh = solveScalarTransport2D(mesh, problem, controls);
    const auto backToJacobi = solveScalarTransport2D(mesh, problem, controls, {}, 0., &workspace);
    compareResults(backJacobiFresh, backToJacobi, "preconditioner switch does not reuse stale scalar fields");
}

void geometryAndGraphChangesDoNotReuseStaleMatrix() {
    const auto sameGraph = grid(2, 2);
    const auto changedGeometry = grid(2, 2, .23);
    const auto changedGraph = grid(4, 1);
    auto controls = tightControls(ScalarPreconditioner2D::ILU0);
    ScalarTransportWorkspace2D workspace;

    const auto first = solveScalarTransport2D(sameGraph, varyingProblem(sameGraph), controls, {}, 0., &workspace);
    check(first.performance.patternBuilds >= 1, "first graph workspace call builds a pattern");
    const auto geometryFresh = solveScalarTransport2D(changedGeometry, varyingProblem(changedGeometry), controls);
    const auto geometryReused = solveScalarTransport2D(changedGeometry, varyingProblem(changedGeometry), controls, {}, 0., &workspace);
    compareResults(geometryFresh, geometryReused, "same graph with changed geometry uses current coefficients");
    check(geometryReused.performance.patternReuses >= 1, "same ordered graph reuses only the pattern after geometry change");

    const auto graphFresh = solveScalarTransport2D(changedGraph, varyingProblem(changedGraph), controls);
    const auto graphReused = solveScalarTransport2D(changedGraph, varyingProblem(changedGraph), controls, {}, 0., &workspace);
    compareResults(graphFresh, graphReused, "different ordered graph matches a fresh solve");
    check(graphReused.performance.patternBuilds >= 1, "different ordered graph rebuilds the pattern");
}

void recursiveCallbackFailureReleasesWorkspace() {
    const auto mesh = grid(2, 2);
    auto controls = tightControls();
    ScalarTransportWorkspace2D workspace;
    ScalarTransportProblem2D recursiveProblem = varyingProblem(mesh);
    recursiveProblem.sourceDensity.clear();
    recursiveProblem.source = [&](Point2D) {
        (void)solveScalarTransport2D(mesh, recursiveProblem, controls, {}, 0., &workspace);
        return 0.;
    };

    bool rejected = false;
    try {
        (void)solveScalarTransport2D(mesh, recursiveProblem, controls, {}, 0., &workspace);
    } catch (const std::runtime_error& error) {
        rejected = std::string(error.what()).find("already in use") != std::string::npos;
    }
    check(rejected, "recursive source callback use of one workspace is explicitly rejected");

    auto validProblem = varyingProblem(mesh);
    validProblem.sourceDensity.clear();
    validProblem.source = [](Point2D point) {
        return .35 + .27 * point.x + .13 * point.y;
    };
    const auto fresh = solveScalarTransport2D(mesh, validProblem, controls);
    const auto recovered = solveScalarTransport2D(mesh, validProblem, controls, {}, 0., &workspace);
    compareResults(fresh, recovered, "workspace is reusable after recursive callback failure");
}

void evaluationFailureAndMoveLifecycle() {
    const auto mesh = grid(3, 2);
    const auto problem = varyingProblem(mesh);
    auto controls = tightControls();
    controls.maxCorrections = 1;
    controls.relaxation = .5;
    ScalarTransportWorkspace2D workspace;
    const auto failed = solveSteadyScalarTransportFromInitial2D(
        mesh, problem, std::vector<double>(mesh.cells.size(), 0.), controls, &workspace);
    check(!failed.converged, "one-correction workspace solve remains honestly unconverged");

    controls.maxCorrections = 1200;
    const auto fresh = solveSteadyScalarTransportFromInitial2D(
        mesh, problem, std::vector<double>(mesh.cells.size(), 0.), controls);
    const auto recovered = solveSteadyScalarTransportFromInitial2D(
        mesh, problem, std::vector<double>(mesh.cells.size(), 0.), controls, &workspace);
    compareResults(fresh, recovered, "workspace recovers exactly after a failed solve");

    const auto evaluatedFresh = evaluateScalarTransport2D(mesh, problem, recovered.values, controls);
    const auto evaluatedWorkspace = evaluateScalarTransport2D(mesh, problem, recovered.values, controls, {}, 0., &workspace);
    compareResults(evaluatedFresh, evaluatedWorkspace, "workspace evaluation matches fresh evaluation");
    check(evaluatedWorkspace.values == recovered.values && !evaluatedWorkspace.history.empty() &&
              evaluatedWorkspace.history.back().linearIterations == 0,
          "evaluation neither mutates nor solves the supplied field");

    static_assert(std::is_move_constructible_v<ScalarTransportWorkspace2D>);
    static_assert(std::is_move_assignable_v<ScalarTransportWorkspace2D>);
    static_assert(!std::is_copy_constructible_v<ScalarTransportWorkspace2D>);
    static_assert(!std::is_copy_assignable_v<ScalarTransportWorkspace2D>);
    ScalarTransportWorkspace2D movedTo(std::move(workspace));
    const auto movedResult = solveScalarTransport2D(mesh, problem, controls, {}, 0., &movedTo);
    compareResults(fresh, movedResult, "moved workspace remains reusable");
    const auto movedFromResult = solveScalarTransport2D(mesh, problem, controls, {}, 0., &workspace);
    compareResults(fresh, movedFromResult, "moved-from workspace is safely reusable");

    ScalarTransportWorkspace2D assigned;
    assigned = std::move(movedTo);
    const auto assignedResult = solveScalarTransport2D(mesh, problem, controls, {}, 0., &assigned);
    compareResults(fresh, assignedResult, "move-assigned workspace remains reusable");
}
}

int main() {
    try {
        exactUniformSteadyAndTransient();
        coefficientsAndPreconditionerChangesResetAllValues();
        geometryAndGraphChangesDoNotReuseStaleMatrix();
        recursiveCallbackFailureReleasesWorkspace();
        evaluationFailureAndMoveLifecycle();
    } catch (const std::exception& error) {
        std::cerr << "EXCEPTION: " << error.what() << '\n';
        return 1;
    }
    if (failures) {
        std::cerr << failures << " scalar workspace checks failed\n";
        return 1;
    }
    std::cout << "Scalar transport workspace checks passed\n";
    return 0;
}
