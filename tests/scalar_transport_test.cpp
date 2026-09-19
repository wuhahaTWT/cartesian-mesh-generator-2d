#include "cartmesh2d/fv/ScalarTransport2D.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {
int failures = 0;

void check(bool value, const std::string& message) {
    if (!value) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

template<class F>
void rejects(F&& f, const std::string& expected, const std::string& message) {
    try {
        f();
        check(false, message);
    } catch (const std::exception& error) {
        check(std::string(error.what()).find(expected) != std::string::npos,
              message + " (wrong reason: " + error.what() + ")");
    }
}

TopologyMesh2D fromPolygons(const std::vector<Polygon2D>& polygons) {
    TopologyMesh2D result;
    std::map<std::pair<double, double>, std::size_t> vertices;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edges;
    for (const auto& polygon : polygons) {
        TopologyCell2D cell;
        cell.id = result.cells.size();
        cell.geometryArea = polygon.area();
        for (const auto& point : polygon.vertices) {
            const auto key = std::make_pair(point.x, point.y);
            auto [it, inserted] = vertices.emplace(key, result.vertices.size());
            if (inserted) result.vertices.push_back({it->second, point});
            cell.vertices.push_back(it->second);
        }
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

TopologyMesh2D grid(int nx, int ny, bool skew = false) {
    auto point = [=](int i, int j) {
        const double x = static_cast<double>(i) / nx;
        const double y = static_cast<double>(j) / ny;
        return Point2D{x + (skew ? .2 * y : 0.), y};
    };
    std::vector<Polygon2D> polygons;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            polygons.push_back({{point(i, j), point(i + 1, j), point(i + 1, j + 1), point(i, j + 1)}});
        }
    }
    return fromPolygons(polygons);
}

FvMesh2D fvGrid(int nx, int ny, bool skew = false) {
    return makeFvMesh2D(grid(nx, ny, skew));
}

double faceNormalFlux(const Face& face, const Vector2D& velocity) {
    return velocity.x * face.areaVector.x + velocity.y * face.areaVector.y;
}

ScalarBoundary2D valueBoundary(double value) {
    return {ScalarBoundaryKind2D::Value, value, std::nullopt};
}

void affineMixedDiffusion() {
    const auto mesh = fvGrid(5, 4, true);
    constexpr double diffusivity = .7;
    const auto exact = [](Point2D p) { return .8 + 1.3 * p.x - .6 * p.y; };
    ScalarTransportProblem2D problem;
    problem.diffusivity = diffusivity;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.source = [](Point2D) { return 0.; };
    problem.boundary = [&](std::size_t, const Face& face) {
        if (face.centre.x < .3 || face.centre.y > .8) return valueBoundary(exact(face.centre));
        const double length = std::hypot(face.areaVector.x, face.areaVector.y);
        const double outwardFluxPerLength = -diffusivity *
            (1.3 * face.areaVector.x - .6 * face.areaVector.y) / length;
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux, outwardFluxPerLength, std::nullopt};
    };
    ScalarTransportControls2D controls;
    controls.maxCorrections = 2000;
    controls.relativeTolerance = 1e-11;
    controls.absoluteTolerance = 1e-13;
    controls.cellTolerance = 1e-11;
    const auto result = solveScalarTransport2D(mesh, problem, controls);
    check(result.converged, "mixed value/flux affine diffusion converges");
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        check(std::abs(result.values[i] - exact(mesh.cells[i].centre)) < 5e-8,
              "mixed affine diffusion recovers the exact scalar");
    }
    check(std::abs(result.globalBalance) < 1e-8, "affine diffusion has zero balance");
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        double balance = -result.sourceIntegrals[i];
        for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
            const auto& face = mesh.faces[id];
            const double flux = result.advectiveFlux[id] + result.diffusiveFlux[id];
            if (face.owner == i) balance += flux;
            if (face.neighbour && *face.neighbour == i) balance -= flux;
        }
        check(std::abs(balance) < 2e-8, "independent affine cell flux balance closes");
    }
}

void insulatingTransientSource() {
    const auto mesh = fvGrid(4, 3);
    ScalarTransportProblem2D problem;
    problem.diffusivity = .2;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.source = [](Point2D) { return 2.; };
    problem.boundary = [](std::size_t, const Face&) {
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux, 0., std::nullopt};
    };
    std::vector<double> previous(mesh.cells.size(), 0.);
    const auto result = solveScalarTransport2D(mesh, problem, {}, previous, .125);
    check(result.converged, "insulating source transient converges");
    for (double value : result.values) check(std::abs(value - .25) < 2e-9, "uniform source gives exact constant rise");
    check(std::abs(result.sourceIntegral - 2.) < 1e-9,
          "source integral is reported independently");
    check(std::abs(result.temporalIntegral - result.sourceIntegral) < 1e-8,
          "transient source and temporal balance closes");
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        double balance = result.temporalIntegrals[i] - result.sourceIntegrals[i];
        for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
            const auto& face = mesh.faces[id];
            if (face.owner == i) balance += result.advectiveFlux[id] + result.diffusiveFlux[id];
            if (face.neighbour && *face.neighbour == i) balance -= result.advectiveFlux[id] + result.diffusiveFlux[id];
        }
        check(std::abs(balance) < 1e-8, "independent transient cell balance closes");
    }
}

void constantCarrierPreservesConstant() {
    const auto mesh = fvGrid(6, 5, true);
    const Vector2D velocity{.73, -.41};
    ScalarTransportProblem2D problem;
    problem.diffusivity = .03;
    problem.volumeFlux.resize(mesh.faces.size());
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) problem.volumeFlux[i] = faceNormalFlux(mesh.faces[i], velocity);
    problem.source = [](Point2D) { return 0.; };
    problem.boundary = [](std::size_t, const Face&) { return valueBoundary(3.25); };
    const auto result = solveScalarTransport2D(mesh, problem);
    check(result.converged, "divergence-free constant-carrier solve converges");
    for (double value : result.values) check(std::abs(value - 3.25) < 2e-8, "constant scalar is preserved by carrier");
    check(result.maxCarrierImbalance < 1e-10, "constant carrier is independently divergence-free");
}

void boundedUpwindFront() {
    const auto mesh = fvGrid(24, 8);
    const Vector2D velocity{1., 0.};
    ScalarTransportProblem2D problem;
    problem.diffusivity = 1e-5;
    problem.volumeFlux.resize(mesh.faces.size());
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) problem.volumeFlux[i] = faceNormalFlux(mesh.faces[i], velocity);
    problem.source = [](Point2D) { return 0.; };
    problem.boundary = [](std::size_t, const Face& face) {
        if (face.centre.x < .01) return valueBoundary(1.);
        if (face.centre.x > .99) return valueBoundary(0.);
        return valueBoundary(.5);
    };
    ScalarTransportControls2D controls;
    controls.convection = ConvectionScheme2D::Upwind;
    const auto result = solveScalarTransport2D(mesh, problem, controls);
    check(result.converged, "upwind front converges");
    check(result.minValue > -1e-10 && result.maxValue < 1. + 1e-10, "upwind front remains bounded");
}

double sineError(int n, ConvectionScheme2D scheme) {
    const auto mesh = fvGrid(n, n, false);
    constexpr double pi = std::numbers::pi;
    constexpr double diffusivity = .08;
    const Vector2D velocity{.35, -.17};
    const auto exact = [](Point2D p) { return std::sin(std::numbers::pi * p.x) * std::sin(std::numbers::pi * p.y); };
    ScalarTransportProblem2D problem;
    problem.diffusivity = diffusivity;
    problem.volumeFlux.resize(mesh.faces.size());
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) problem.volumeFlux[i] = faceNormalFlux(mesh.faces[i], velocity);
    problem.source = [&](Point2D p) {
        return velocity.x * pi * std::cos(pi * p.x) * std::sin(pi * p.y)
            + velocity.y * pi * std::sin(pi * p.x) * std::cos(pi * p.y)
            + 2. * diffusivity * pi * pi * exact(p);
    };
    problem.boundary = [&](std::size_t, const Face& face) { return valueBoundary(exact(face.centre)); };
    ScalarTransportControls2D controls;
    controls.convection = scheme;
    const auto result = solveScalarTransport2D(mesh, problem, controls);
    check(result.converged, "manufactured sine transport converges");
    double squared = 0.;
    double area = 0.;
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        squared += mesh.cells[i].area * std::pow(result.values[i] - exact(mesh.cells[i].centre), 2);
        area += mesh.cells[i].area;
    }
    return std::sqrt(squared / area);
}

void manufacturedRefinement() {
    for (const auto scheme : {ConvectionScheme2D::Upwind, ConvectionScheme2D::LimitedLinearUpwind}) {
        const double e8 = sineError(8, scheme);
        const double e16 = sineError(16, scheme);
        const double e32 = sineError(32, scheme);
        check(std::isfinite(e8) && std::isfinite(e16) && std::isfinite(e32), "manufactured errors are finite");
        check(e16 < e8 && e32 < e16, "manufactured sine error decreases under refinement");
        check(e32 < .025, "manufactured sine fine-grid error is small");
        std::cout << (scheme == ConvectionScheme2D::Upwind ? "upwind" : "limited-linear")
                  << " sine L2: " << e8 << " -> " << e16 << " -> " << e32 << '\n';
    }
}

void transientSineTimeRefinement() {
    const auto mesh = fvGrid(16, 16);
    constexpr double diffusivity = .2;
    const auto exact = [](Point2D p, double time) {
        return std::sin(std::numbers::pi * p.x) * std::sin(std::numbers::pi * p.y)
            * std::exp(-2. * std::numbers::pi * std::numbers::pi * diffusivity * time);
    };
    ScalarTransportProblem2D problem;
    problem.diffusivity = diffusivity;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.source = [](Point2D) { return 0.; };
    problem.boundary = [](std::size_t, const Face&) { return valueBoundary(0.); };
    std::vector<double> previous;
    for (const auto& cell : mesh.cells) previous.push_back(exact(cell.centre, 0.));
    auto errorAt = [&](int steps) {
        const double finalTime = .5;
        const double dt = finalTime / static_cast<double>(steps);
        auto state = previous;
        ScalarTransportResult2D result;
        for (int step = 0; step < steps; ++step) {
            result = solveScalarTransport2D(mesh, problem, {}, state, dt);
            check(result.converged, "transient sine step converges");
            state = result.values;
        }
        // Compare every run at the same physical time against the continuous
        // sine decay; the spatial grid is held fixed, so dt is the only
        // refinement variable.
        const double expectedAmplitude = std::exp(-2. * std::numbers::pi * std::numbers::pi * diffusivity * finalTime);
        double e = 0., area = 0.;
        for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
            e += mesh.cells[i].area * std::pow(state[i] - exact(mesh.cells[i].centre, 0.) * expectedAmplitude, 2);
            area += mesh.cells[i].area;
        }
        return std::sqrt(e / area);
    };
    const double coarse = errorAt(5);
    const double middle = errorAt(10);
    const double fine = errorAt(20);
    std::cout << "transient sine time-refinement errors: " << coarse << " -> " << middle << " -> " << fine << '\n';
    check(middle < coarse && fine < middle, "transient sine error decreases with time refinement at common final time");
}

void invalidInputs() {
    const auto mesh = fvGrid(3, 3);
    ScalarTransportProblem2D valid;
    valid.diffusivity = .1;
    valid.volumeFlux.assign(mesh.faces.size(), 0.);
    valid.source = [](Point2D) { return 0.; };
    valid.boundary = [](std::size_t, const Face&) { return valueBoundary(0.); };

    auto nonfinite = valid;
    nonfinite.diffusivity = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { (void)solveScalarTransport2D(mesh, nonfinite); }, "diffusivity must be finite positive", "nonfinite diffusivity is rejected");

    for (const auto bad : {std::vector<double>(mesh.cells.size() - 1, .1),
                           std::vector<double>(mesh.cells.size(), -.1),
                           std::vector<double>(mesh.cells.size(), std::numeric_limits<double>::quiet_NaN()),
                           std::vector<double>(mesh.cells.size(), std::numeric_limits<double>::infinity())}) {
        auto invalidSink = valid;
        invalidSink.sinkRate = bad;
        rejects([&] { (void)solveScalarTransport2D(mesh, invalidSink); }, "sink",
                "invalid sink field is rejected");
    }

    auto wrongFaceCount = valid;
    wrongFaceCount.faceDiffusivity.assign(mesh.faces.size() - 1, .1);
    rejects([&] { (void)solveScalarTransport2D(mesh, wrongFaceCount); },
            "face diffusivity dimensions mismatch", "face diffusivity size is rejected");
    for (const double value : {0., -.1, std::numeric_limits<double>::quiet_NaN()}) {
        auto badFace = valid;
        badFace.faceDiffusivity.assign(mesh.faces.size(), .1);
        badFace.faceDiffusivity[0] = value;
        rejects([&] { (void)solveScalarTransport2D(mesh, badFace); },
                "face diffusivity must be finite positive", "invalid face diffusivity is rejected");
    }

    auto divergent = valid;
    divergent.volumeFlux[0] = 1.;
    rejects([&] { (void)solveScalarTransport2D(mesh, divergent); }, "carrier flux violates cell continuity", "non-divergence-free carrier is rejected");

    auto unanchored = valid;
    unanchored.boundary = [](std::size_t, const Face&) {
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux, 0., std::nullopt};
    };
    rejects([&] { (void)solveScalarTransport2D(mesh, unanchored); }, "unanchored steady connected component", "unanchored steady Neumann problem is rejected");

    auto missingInflowValue = valid;
    const Vector2D carrier{1., 0.};
    for (std::size_t i = 0; i < mesh.faces.size(); ++i)
        missingInflowValue.volumeFlux[i] = faceNormalFlux(mesh.faces[i], carrier);
    missingInflowValue.boundary = [](std::size_t, const Face&) {
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux, 0., std::nullopt};
    };
    rejects([&] { (void)solveScalarTransport2D(mesh, missingInflowValue); }, "inflow requires a prescribed scalar value", "inflow flux without inflow value is rejected");

    auto oneCorrection = valid;
    oneCorrection.boundary = [](std::size_t, const Face& face) {
        return valueBoundary(face.centre.x < .5 ? 1. : 0.);
    };
    ScalarTransportControls2D controls;
    controls.maxCorrections = 1;
    const auto result = solveScalarTransport2D(mesh, oneCorrection, controls);
    check(!result.converged, "one correction reports nonconvergence");
}

void transientUniformReaction() {
    const auto mesh = fvGrid(4, 3);
    ScalarTransportProblem2D problem;
    problem.diffusivity = .2;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.sinkRate.assign(mesh.cells.size(), .7);
    problem.source = [](Point2D) { return 2.; };
    problem.boundary = [](std::size_t, const Face&) {
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux, 0., std::nullopt};
    };
    const std::vector<double> previous(mesh.cells.size(), 1.);
    const double dt = 10.;
    const auto result = solveScalarTransport2D(mesh, problem, {}, previous, dt);
    const double expected = (1. + dt * 2.) / (1. + dt * .7);
    check(result.converged, "implicit uniform reaction converges for large dt");
    for (double value : result.values) {
        check(std::abs(value - expected) < 2e-9, "implicit reaction matches backward-Euler exact value");
        check(value > 0., "implicit reaction remains nonnegative for large dt");
    }
    check(std::abs(result.sinkIntegral - .7 * expected) < 2e-9,
          "reaction sink integral is reported");
    check(std::abs(result.globalBalance) < 2e-8,
          "reaction transient global balance closes");
}

void steadyNeumannReaction() {
    const auto mesh = fvGrid(3, 2);
    ScalarTransportProblem2D problem;
    problem.diffusivity = .3;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.sinkRate.assign(mesh.cells.size(), .4);
    problem.source = [](Point2D) { return 2.; };
    problem.boundary = [](std::size_t, const Face&) {
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux, 0., std::nullopt};
    };
    const auto result = solveScalarTransport2D(mesh, problem);
    check(result.converged, "positive reaction anchors steady Neumann problem");
    for (double value : result.values)
        check(std::abs(value - 5.) < 2e-8, "steady Neumann reaction has source-over-sink constant solution");
    check(std::abs(result.sinkIntegral - 2.) < 2e-8, "steady reaction sink balances source");
    check(std::abs(result.globalBalance) < 2e-8, "steady reaction global balance closes");
}

void spatialReactionVariableDiffusionBalance() {
    const auto mesh = fvGrid(5, 4, true);
    const auto exact = [](Point2D p) { return 1.7 + p.x + .3 * p.y; };
    const Vector2D gradient{1., .3};
    ScalarTransportProblem2D problem;
    problem.diffusivity = .1;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.faceDiffusivity.resize(mesh.faces.size());
    problem.sinkRate.resize(mesh.cells.size());
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        problem.faceDiffusivity[id] = .1 * (1. + mesh.faces[id].centre.x);
    for (std::size_t i = 0; i < mesh.cells.size(); ++i)
        problem.sinkRate[i] = .2 * (1. + .3 * mesh.cells[i].centre.x);
    problem.source = [&](Point2D p) { return -.1 + .2 * (1. + .3 * p.x) * exact(p); };
    problem.boundary = [&](std::size_t, const Face& face) {
        const double d = .1 * (1. + face.centre.x);
        if (face.centre.x - .2 * face.centre.y < 1e-10)
            return valueBoundary(exact(face.centre));
        const double length = std::hypot(face.areaVector.x, face.areaVector.y);
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux,
            -d * (gradient.x * face.areaVector.x + gradient.y * face.areaVector.y) / length,
            std::nullopt};
    };
    ScalarTransportControls2D controls;
    controls.maxCorrections = 3000;
    controls.relativeTolerance = 1e-11;
    controls.absoluteTolerance = 1e-13;
    controls.cellTolerance = 1e-11;
    const auto result = solveScalarTransport2D(mesh, problem, controls);
    check(result.converged, "spatial reaction variable-D affine solve converges");
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        check(std::abs(result.values[i] - exact(mesh.cells[i].centre)) < 3e-7,
              "spatial reaction variable-D affine solution is recovered");
        double balance = result.sinkIntegrals[i] - result.sourceIntegrals[i];
        for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
            const auto& face = mesh.faces[id];
            const double flux = result.advectiveFlux[id] + result.diffusiveFlux[id];
            if (face.owner == i) balance += flux;
            if (face.neighbour && *face.neighbour == i) balance -= flux;
        }
        check(std::abs(balance) < 3e-8, "spatial reaction cell balance includes sink");
    }
}

void zeroSinkCompatibilityAndLocalAnchor() {
    const auto mesh = fvGrid(3, 2);
    ScalarTransportProblem2D base;
    base.diffusivity = .2;
    base.volumeFlux.assign(mesh.faces.size(), 0.);
    base.source = [](Point2D) { return 0.; };
    base.boundary = [](std::size_t, const Face&) {
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux, 0., std::nullopt};
    };
    auto previous = std::vector<double>(mesh.cells.size(), 1.25);
    const auto oldResult = solveScalarTransport2D(mesh, base, {}, previous, .1);
    auto zero = base;
    zero.sinkRate.assign(mesh.cells.size(), 0.);
    const auto zeroResult = solveScalarTransport2D(mesh, zero, {}, previous, .1);
    check(oldResult.values == zeroResult.values && oldResult.diffusiveFlux == zeroResult.diffusiveFlux,
          "empty and explicit zero sink preserve old outputs");
    auto local = base;
    local.sinkRate.assign(mesh.cells.size(), 0.);
    local.sinkRate[0] = .5;
    const auto anchored = solveScalarTransport2D(mesh, local);
    check(anchored.converged, "one positive sink anchors a connected Neumann component");
    check(std::abs(anchored.globalBalance) < 2e-8, "local sink anchor balance closes");
}

void evaluateScalarTransportRegression() {
    const auto mesh = fvGrid(4, 3, true);
    const auto exact = [](Point2D p) { return 1.4 + p.x + .25 * p.y; };
    ScalarTransportProblem2D problem;
    problem.diffusivity = .1;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.faceDiffusivity.resize(mesh.faces.size());
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        problem.faceDiffusivity[id] = .1 * (1. + mesh.faces[id].centre.x);
    problem.sinkRate.resize(mesh.cells.size());
    for (std::size_t i = 0; i < mesh.cells.size(); ++i)
        problem.sinkRate[i] = .2 * (1. + .1 * mesh.cells[i].centre.x);
    problem.source = [&](Point2D p) { return -.1 + .2 * (1. + .1 * p.x) * exact(p); };
    problem.boundary = [&](std::size_t, const Face& face) {
        const double d = .1 * (1. + face.centre.x);
        if (face.centre.x - .2 * face.centre.y < 1e-10) return valueBoundary(exact(face.centre));
        const double length = std::hypot(face.areaVector.x, face.areaVector.y);
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux,
            -d * (face.areaVector.x + .25 * face.areaVector.y) / length, std::nullopt};
    };
    ScalarTransportControls2D controls;
    controls.maxCorrections = 2500;
    controls.relativeTolerance = 1e-11;
    controls.absoluteTolerance = 1e-13;
    controls.cellTolerance = 1e-11;
    const auto solved = solveScalarTransport2D(mesh, problem, controls);
    check(solved.converged, "nonlinear re-evaluation fixture converges");
    const auto evaluated = evaluateScalarTransport2D(mesh, problem, solved.values, controls);
    check(evaluated.values == solved.values, "evaluation preserves supplied solved values");
    check(evaluated.advectiveFlux == solved.advectiveFlux && evaluated.diffusiveFlux == solved.diffusiveFlux,
          "evaluation reproduces solved face fluxes");
    check(evaluated.sinkIntegrals == solved.sinkIntegrals && evaluated.sourceIntegrals == solved.sourceIntegrals,
          "evaluation reproduces solved source and sink integrals");
    check(evaluated.globalBalance == solved.globalBalance && evaluated.converged == solved.converged,
          "evaluation reproduces solved balance and convergence");
    check(evaluated.history.size() == 1 && evaluated.history[0].iteration == 0 &&
          evaluated.history[0].linearIterations == 0,
          "evaluation records one zero-linear-iteration history entry");

    auto perturbedValues = solved.values;
    perturbedValues[0] += .5;
    const auto perturbed = evaluateScalarTransport2D(mesh, problem, perturbedValues, controls);
    check(perturbed.values == perturbedValues, "evaluation does not silently solve a perturbed field");
    check(!perturbed.converged && perturbed.history.size() == 1 &&
          perturbed.history[0].linearIterations == 0,
          "perturbed evaluation reports the nonlinear residual without solving");
    double maxResidual = 0., residualNormSquared = 0.;
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        double residual = perturbed.sinkIntegrals[i] - perturbed.sourceIntegrals[i];
        for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
            const auto& face = mesh.faces[id];
            const double flux = perturbed.advectiveFlux[id] + perturbed.diffusiveFlux[id];
            if (face.owner == i) residual += flux;
            if (face.neighbour && *face.neighbour == i) residual -= flux;
        }
        maxResidual = std::max(maxResidual, std::abs(residual));
        residualNormSquared += residual * residual;
    }
    check(std::abs(perturbed.history[0].maxCellImbalance - maxResidual) < 1e-12,
          "evaluation history reports independently reconstructed cell residual");
    check(std::abs(perturbed.history[0].residualNorm - std::sqrt(residualNormSquared)) < 1e-12,
          "evaluation history reports independently reconstructed residual norm");

    auto changedSink = problem;
    changedSink.sinkRate[0] *= 20.;
    const auto sinkChanged = evaluateScalarTransport2D(mesh, changedSink, solved.values, controls);
    check(!sinkChanged.converged, "changed sink invalidates a previously converged field");
    auto changedFace = problem;
    changedFace.faceDiffusivity[0] *= 4.;
    const auto faceChanged = evaluateScalarTransport2D(mesh, changedFace, solved.values, controls);
    check(!faceChanged.converged, "changed face diffusivity invalidates a previously converged field");

    rejects([&] { (void)evaluateScalarTransport2D(mesh, problem, std::vector<double>(mesh.cells.size() - 1), controls); },
            "field dimensions", "evaluation rejects wrong value size");
    auto nonfinite = solved.values;
    nonfinite[0] = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { (void)evaluateScalarTransport2D(mesh, problem, nonfinite, controls); },
            "numerical range", "evaluation rejects nonfinite values");
}

void uniformFaceDiffusivityIsIdentical() {
    const auto mesh = fvGrid(5, 4, true);
    ScalarTransportProblem2D scalar;
    scalar.diffusivity = .17;
    scalar.volumeFlux.assign(mesh.faces.size(), 0.);
    scalar.source = [](Point2D p) { return .3 + p.x - .4 * p.y; };
    scalar.boundary = [](std::size_t, const Face& face) { return valueBoundary(.2 + face.centre.x); };
    auto faceValues = scalar;
    faceValues.faceDiffusivity.assign(mesh.faces.size(), scalar.diffusivity);
    const auto a = solveScalarTransport2D(mesh, scalar);
    const auto b = solveScalarTransport2D(mesh, faceValues);
    check(a.converged && b.converged, "uniform and face-valued scalar solves converge");
    check(a.values == b.values, "uniform face diffusivity preserves scalar values exactly");
    check(a.advectiveFlux == b.advectiveFlux && a.diffusiveFlux == b.diffusiveFlux,
          "uniform face diffusivity preserves fluxes exactly");
    check(a.globalBalance == b.globalBalance, "uniform face diffusivity preserves balance exactly");
}

void variableFaceDiffusivityConservesHarmonicInterface() {
    const auto mesh = fvGrid(2, 1);
    constexpr double leftD = .1, rightD = .4;
    constexpr double interfaceD = 2. * leftD * rightD / (leftD + rightD);
    constexpr double flux = .16;
    const std::vector<double> exact{.6, .1};
    ScalarTransportProblem2D problem;
    problem.diffusivity = 1.;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.faceDiffusivity.resize(mesh.faces.size(), 1.);
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        if (face.neighbour) problem.faceDiffusivity[id] = interfaceD;
        else if (face.centre.x < 1e-12) problem.faceDiffusivity[id] = leftD;
        else if (face.centre.x > 1. - 1e-12) problem.faceDiffusivity[id] = rightD;
    }
    problem.source = [](Point2D) { return 0.; };
    problem.boundary = [&](std::size_t id, const Face& face) {
        if (face.neighbour) {
            return valueBoundary(0.);
        }
        if (face.centre.x < 1e-12) {
            return valueBoundary(1.);
        }
        if (face.centre.x > 1. - 1e-12) {
            return valueBoundary(0.);
        }
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux, 0., std::nullopt};
    };
    // Tight algebraic controls make this exact two-cell resistance check
    // independent of the production default stopping error.
    ScalarTransportControls2D controls;
    controls.relativeTolerance = 1e-12;
    controls.absoluteTolerance = 1e-14;
    controls.cellTolerance = 1e-12;
    const auto result = solveScalarTransport2D(mesh, problem, controls);
    check(result.converged, "two-material harmonic-interface solve converges");
    for (std::size_t i = 0; i < exact.size(); ++i)
        check(std::abs(result.values[i] - exact[i]) < 2e-10,
              "two-material interface recovers finite-volume analytic values");
    std::size_t interface = mesh.faces.size();
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        if (mesh.faces[id].neighbour) interface = id;
    check(interface < mesh.faces.size(), "two-material mesh has an internal interface");
    if (interface < mesh.faces.size())
        check(std::abs(result.diffusiveFlux[interface] - flux) < 2e-10,
              "harmonic interface carries the analytic conserved flux");
    check(std::abs(result.globalBalance) < 2e-10, "variable face diffusivity closes global balance");
}

void variableFaceDiffusivityManufacturedSkewAffine() {
    const auto mesh = fvGrid(6, 5, true);
    const auto exact = [](Point2D p) { return 1.7 + p.x + .3 * p.y; };
    const Vector2D gradient{1., .3};
    ScalarTransportProblem2D problem;
    problem.diffusivity = .1;
    problem.volumeFlux.assign(mesh.faces.size(), 0.);
    problem.faceDiffusivity.resize(mesh.faces.size());
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        problem.faceDiffusivity[id] = .1 * (1. + mesh.faces[id].centre.x);
    problem.source = [](Point2D) { return -.1; };
    problem.boundary = [&](std::size_t, const Face& face) {
        const double d = .1 * (1. + face.centre.x);
        const double length = std::hypot(face.areaVector.x, face.areaVector.y);
        // The left boundary is the only value anchor; all other boundaries
        // prescribe the exact variable-D normal diffusive flux.
        if (face.centre.x - .2 * face.centre.y < 1e-10)
            return valueBoundary(exact(face.centre));
        return ScalarBoundary2D{ScalarBoundaryKind2D::DiffusiveFlux,
                                -d * (gradient.x * face.areaVector.x + gradient.y * face.areaVector.y) / length,
                                std::nullopt};
    };
    ScalarTransportControls2D controls;
    controls.maxCorrections = 3000;
    controls.relativeTolerance = 1e-11;
    controls.absoluteTolerance = 1e-13;
    controls.cellTolerance = 1e-11;
    const auto result = solveScalarTransport2D(mesh, problem, controls);
    check(result.converged, "variable-D skew affine solve converges");
    for (std::size_t i = 0; i < mesh.cells.size(); ++i)
        check(std::abs(result.values[i] - exact(mesh.cells[i].centre)) < 2e-7,
              "variable-D skew affine manufactured solution is recovered");
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        const double expectedFlux = -problem.faceDiffusivity[id] * dot(gradient, face.areaVector);
        check(std::abs(result.diffusiveFlux[id] - expectedFlux) < 2e-8,
              "variable-D skew affine shared-face flux agrees with exact field");
    }
    check(std::abs(result.globalBalance) < 2e-8,
          "variable-D skew affine flux balance closes");
}
}

int main() {
    try {
        affineMixedDiffusion();
        insulatingTransientSource();
        constantCarrierPreservesConstant();
        boundedUpwindFront();
        manufacturedRefinement();
        transientSineTimeRefinement();
        invalidInputs();
        uniformFaceDiffusivityIsIdentical();
        variableFaceDiffusivityConservesHarmonicInterface();
        variableFaceDiffusivityManufacturedSkewAffine();
        transientUniformReaction();
        steadyNeumannReaction();
        spatialReactionVariableDiffusionBalance();
        zeroSinkCompatibilityAndLocalAnchor();
        evaluateScalarTransportRegression();
    } catch (const std::exception& error) {
        std::cerr << "unexpected scalar transport exception: " << error.what() << '\n';
        return 1;
    }
    std::cout << "scalar transport failures: " << failures << '\n';
    return failures == 0 ? 0 : 1;
}
