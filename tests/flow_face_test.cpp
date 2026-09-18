#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cartmesh2d;
using cartmesh2d::fv::FvMesh2D;
using cartmesh2d::fv::detail::conservativePressureGradient;
using cartmesh2d::fv::detail::faceReconstructionLimiter;
using cartmesh2d::fv::detail::flowGradient;
using cartmesh2d::fv::detail::pressureFaceValues;
using cartmesh2d::fv::detail::upwindFaceValue;
using cartmesh2d::fv::detail::symmetricViscousCorrection;

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

void near(double actual, double expected, double tolerance, const std::string& message) {
    check(std::abs(actual - expected) <= tolerance, message);
}

TopologyMesh2D fromPolygons(const std::vector<Polygon2D>& polygons) {
    TopologyMesh2D topology;
    std::map<std::pair<double, double>, std::size_t> vertices;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edges;
    for (const auto& polygon : polygons) {
        TopologyCell2D cell;
        cell.id = topology.cells.size();
        cell.geometryArea = polygon.area();
        for (const auto point : polygon.vertices) {
            const auto key = std::make_pair(point.x, point.y);
            const auto [it, inserted] = vertices.emplace(key, topology.vertices.size());
            if (inserted) topology.vertices.push_back({it->second, point});
            cell.vertices.push_back(it->second);
        }
        for (std::size_t i = 0; i < cell.vertices.size(); ++i) {
            const auto a = cell.vertices[i];
            const auto b = cell.vertices[(i + 1) % cell.vertices.size()];
            const auto key = std::minmax(a, b);
            const auto [it, inserted] = edges.emplace(key, topology.edges.size());
            if (inserted) {
                topology.edges.push_back({it->second, a, b, cell.id, {},
                                          BoundaryPatch2D::DomainBoundary});
            } else {
                auto& edge = topology.edges[it->second];
                edge.neighbour = cell.id;
                edge.patch = BoundaryPatch2D::None;
            }
            cell.edges.push_back(it->second);
        }
        topology.cells.push_back(std::move(cell));
    }
    return topology;
}

FvMesh2D rectangularMesh(std::size_t nx, std::size_t ny, double shear = 0.0) {
    std::vector<Polygon2D> polygons;
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const auto point = [shear](std::size_t x, std::size_t y) {
                const double yy = static_cast<double>(y);
                return Point2D{static_cast<double>(x) + shear * yy, yy};
            };
            polygons.push_back({{point(i, j), point(i + 1, j), point(i + 1, j + 1),
                                 point(i, j + 1)}});
        }
    }
    return cartmesh2d::fv::makeFvMesh2D(fromPolygons(polygons));
}

FvMesh2D normalizedSkewMesh(std::size_t n) {
    std::vector<Polygon2D> polygons;
    const double inverse = 1.0 / static_cast<double>(n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const auto point = [inverse](std::size_t x, std::size_t y) {
                return Point2D{inverse * (static_cast<double>(x) +
                                          0.3 * static_cast<double>(y)),
                               inverse * static_cast<double>(y)};
            };
            polygons.push_back({{point(i, j), point(i + 1, j), point(i + 1, j + 1),
                                 point(i, j + 1)}});
        }
    }
    return cartmesh2d::fv::makeFvMesh2D(fromPolygons(polygons));
}

FvMesh2D splitFaceMesh() {
    // The left cell has a collinear vertex on its right side.  Its one coarse
    // side is therefore represented by two actual shared faces with the two
    // right-hand fine cells.
    return cartmesh2d::fv::makeFvMesh2D(fromPolygons({
        {{{0, 0}, {1, 0}, {1, 1}, {1, 2}, {0, 2}}},
        {{{1, 0}, {2, 0}, {2, 1}, {1, 1}}},
        {{{1, 1}, {2, 1}, {2, 2}, {1, 2}}},
    }));
}

FvMesh2D triangularTipMesh() {
    // Cell 0 has one face-neighbour (the diagonal triangle).  The two
    // quadrilaterals are second-ring cells reached through that neighbour,
    // and provide non-collinear centre offsets for a linear reconstruction.
    return cartmesh2d::fv::makeFvMesh2D(fromPolygons({
        {{{0, 0}, {1, 0}, {0, 1}}},
        {{{1, 0}, {1, 1}, {0, 1}}},
        {{{1, 0}, {2, 0}, {2, 1}, {1, 1}}},
        {{{0, 1}, {1, 1}, {1, 2}, {0, 2}}},
    }));
}

double linear(Point2D point) {
    return 1.25 + 2.0 * point.x - 3.0 * point.y;
}

Vector2D linearGradient() {
    return {2.0, -3.0};
}

std::vector<double> cellValues(const FvMesh2D& mesh, double (*function)(Point2D)) {
    std::vector<double> values;
    values.reserve(mesh.cells.size());
    for (const auto& cell : mesh.cells) values.push_back(function(cell.centre));
    return values;
}

void constantPressureAndInternalCancellation(const FvMesh2D& mesh) {
    const std::vector<double> pressure(mesh.cells.size(), 7.5);
    const std::vector<Vector2D> gradients(mesh.cells.size(), {0.0, 0.0});
    const std::vector<double> boundary(mesh.faces.size(), 7.5);
    const std::vector<bool> fixed(mesh.faces.size(), true);
    const auto facePressure = pressureFaceValues(mesh, pressure, gradients, boundary, fixed);
    const auto result = conservativePressureGradient(mesh, facePressure);
    for (std::size_t i = 0; i < result.size(); ++i) {
        near(result[i].x, 0.0, 1e-12, "constant pressure has zero x force in every cell");
        near(result[i].y, 0.0, 1e-12, "constant pressure has zero y force in every cell");
    }
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        if (!face.neighbour) continue;
        near(facePressure[id], 7.5, 1e-12, "one constant pressure value is shared by an internal face");
        // The two cell contributions of an internal face are opposite because
        // the final mesh stores one owner-oriented normal for that face.
        const double ownerForce = facePressure[id] * face.areaVector.x;
        const double neighbourForce = -facePressure[id] * face.areaVector.x;
        near(ownerForce + neighbourForce, 0.0, 1e-12,
             "constant internal-face pressure forces cancel");
    }
}

void linearPressureIsExact(const FvMesh2D& mesh, const std::string& label) {
    const auto pressure = cellValues(mesh, linear);
    const auto gradients = std::vector<Vector2D>(mesh.cells.size(), linearGradient());
    std::vector<double> boundary(mesh.faces.size());
    std::vector<bool> fixed(mesh.faces.size(), true);
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        boundary[id] = linear(mesh.faces[id].centre);
    }
    const auto facePressure = pressureFaceValues(mesh, pressure, gradients, boundary, fixed);
    const auto result = conservativePressureGradient(mesh, facePressure);
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        near(facePressure[id], linear(mesh.faces[id].centre), 1e-12,
             label + ": reconstructed pressure is exact at every face centre");
    }
    for (const auto value : result) {
        near(value.x, 2.0, 1e-11, label + ": Gauss pressure gradient has exact x component");
        near(value.y, -3.0, 1e-11, label + ": Gauss pressure gradient has exact y component");
    }
}

void reconstructedLinearPressureIsExact(const FvMesh2D& mesh, const std::string& label) {
    const auto pressure = cellValues(mesh, linear);
    std::vector<double> boundary(mesh.faces.size());
    const std::vector<bool> fixed(mesh.faces.size(), true);
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        boundary[id] = linear(mesh.faces[id].centre);

    const auto gradients = flowGradient(mesh, pressure, boundary, fixed, true);
    for (const auto& gradient : gradients) {
        near(gradient.x, linearGradient().x, 1e-11,
             label + ": extrapolated linear gradient has exact x component");
        near(gradient.y, linearGradient().y, 1e-11,
             label + ": extrapolated linear gradient has exact y component");
    }
    const auto facePressure = pressureFaceValues(mesh, pressure, gradients, boundary, fixed);
    const auto result = conservativePressureGradient(mesh, facePressure);
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        near(facePressure[id], linear(mesh.faces[id].centre), 1e-11,
             label + ": reconstructed pressure is exact at every fixed face");
    }
    for (const auto& gradient : result) {
        near(gradient.x, linearGradient().x, 1e-11,
             label + ": Gauss recovery has exact x component");
        near(gradient.y, linearGradient().y, 1e-11,
             label + ": Gauss recovery has exact y component");
    }
}

void unknownBoundaryAndRankChecks(const FvMesh2D& mesh) {
    const auto pressure = cellValues(mesh, linear);
    const std::vector<double> unknownBoundary(mesh.faces.size(), 0.0);
    const std::vector<bool> unknown(mesh.faces.size(), false);
    const auto constrained = flowGradient(mesh, pressure, unknownBoundary, unknown, false);
    bool changed = false;
    for (const auto& gradient : constrained) {
        changed = changed || std::abs(gradient.x - linearGradient().x) > 1e-8
                          || std::abs(gradient.y - linearGradient().y) > 1e-8;
    }
    check(changed, "unknown wall rows with extrapolateUnknown=false do not preserve arbitrary pressure gradient");

    bool threw = false;
    try {
        const auto rankDeficient = rectangularMesh(1, 3);
        const auto values = cellValues(rankDeficient, linear);
        (void)flowGradient(rankDeficient, values,
                           std::vector<double>(rankDeficient.faces.size(), 0.0),
                           std::vector<bool>(rankDeficient.faces.size(), false), true);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    check(threw, "1xN unknown-boundary gradient with extrapolation explicitly reports rank deficiency");
}

void triangularTipUsesSecondRing() {
    const auto mesh = triangularTipMesh();
    std::size_t directNeighbours = 0;
    for (const auto id : mesh.cells[0].faces)
        if (mesh.faces[id].neighbour) ++directNeighbours;
    check(directNeighbours == 1, "triangular tip has exactly one face-neighbour");

    const auto pressure = cellValues(mesh, linear);
    const std::vector<double> unknownBoundary(mesh.faces.size(), 0.0);
    const std::vector<bool> unknown(mesh.faces.size(), false);
    const auto gradients = flowGradient(mesh, pressure, unknownBoundary, unknown, true);
    for (const auto& gradient : gradients) {
        near(gradient.x, linearGradient().x, 1e-11,
             "second-ring tip reconstruction has exact x gradient");
        near(gradient.y, linearGradient().y, 1e-11,
             "second-ring tip reconstruction has exact y gradient");
    }

    const auto facePressure = pressureFaceValues(mesh, pressure, gradients,
                                                  unknownBoundary, unknown);
    const auto result = conservativePressureGradient(mesh, facePressure);
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        near(facePressure[id], linear(mesh.faces[id].centre), 1e-11,
             "second-ring tip pressure is exact on internal and extrapolated faces");
    for (const auto& gradient : result) {
        near(gradient.x, linearGradient().x, 1e-11,
             "second-ring tip Gauss recovery has exact x component");
        near(gradient.y, linearGradient().y, 1e-11,
             "second-ring tip Gauss recovery has exact y component");
    }
}

void checkerboardPressureKeepsDirectDifference(const FvMesh2D& mesh) {
    check(mesh.cells.size() == 2, "checkerboard fixture has two cells");
    const std::vector<double> pressure{1.0, -1.0};
    const std::vector<Vector2D> gradients(mesh.cells.size(), {0.0, 0.0});
    const std::vector<double> boundary(mesh.faces.size(), 0.0);
    const std::vector<bool> fixed(mesh.faces.size(), false);
    const auto facePressure = pressureFaceValues(mesh, pressure, gradients, boundary, fixed);
    std::size_t internal = mesh.faces.size();
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        if (mesh.faces[id].neighbour) {
            internal = id;
            break;
        }
    }
    check(internal < mesh.faces.size(), "checkerboard fixture has an internal face");
    if (internal == mesh.faces.size()) return;
    near(facePressure[internal], 0.0, 1e-12,
         "equal-centre interpolation makes checkerboard internal face average vanish");
    const auto& face = mesh.faces[internal];
    const double directDifference = pressure[*face.neighbour] - pressure[face.owner];
    check(std::abs(directDifference) > 1.0,
          "checkerboard direct pressure difference remains nonzero");
}

void limiterAndUpwindSelection() {
    const auto mesh = rectangularMesh(3, 3);
    const auto constantField = std::vector<double>(mesh.cells.size(), 4.0);
    const auto zeroGradient = std::vector<Vector2D>(mesh.cells.size(), {0.0, 0.0});
    const auto constantBoundary = std::vector<double>(mesh.faces.size(), 4.0);
    const auto fixed = std::vector<bool>(mesh.faces.size(), true);
    const auto constantLimiter = faceReconstructionLimiter(
        mesh, constantField, zeroGradient, constantBoundary, fixed);
    for (const auto value : constantLimiter) near(value, 1.0, 1e-12,
        "limiter preserves a constant field");

    const auto smooth = cellValues(mesh, linear);
    const auto smoothGradient = std::vector<Vector2D>(mesh.cells.size(), linearGradient());
    const auto smoothBoundary = [&]() {
        std::vector<double> values(mesh.faces.size());
        for (std::size_t id = 0; id < mesh.faces.size(); ++id)
            values[id] = linear(mesh.faces[id].centre);
        return values;
    }();
    const auto smoothLimiter = faceReconstructionLimiter(
        mesh, smooth, smoothGradient, smoothBoundary, fixed);
    near(smoothLimiter[4], 1.0, 1e-12,
         "smooth linear reconstruction is unlimited in the interior cell");
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        if (!face.neighbour || (face.owner != 4 && *face.neighbour != 4)) continue;
        const bool cellIsOwner = face.owner == 4;
        const double flux = cellIsOwner ? 1.0 : -1.0;
        near(upwindFaceValue(mesh, id, flux, smooth, smoothGradient, smoothLimiter),
             linear(face.centre), 1e-12,
             "positive/negative flux selects the correct linear upwind state");
    }
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        near(upwindFaceValue(mesh, id, 1.0, constantField, zeroGradient, constantLimiter),
             4.0, 1e-12, "positive flux preserves constant face value");
        near(upwindFaceValue(mesh, id, -1.0, constantField, zeroGradient, constantLimiter),
             4.0, 1e-12, "negative flux preserves constant face value");
    }

    auto step = std::vector<double>(mesh.cells.size(), 0.0);
    for (std::size_t i = 0; i < mesh.cells.size(); ++i)
        if (mesh.cells[i].centre.x > 1.5) step[i] = 1.0;
    auto steep = std::vector<Vector2D>(mesh.cells.size(), {0.0, 0.0});
    steep[4] = {10.0, 0.0};
    const auto noBoundary = std::vector<double>(mesh.faces.size(), 0.0);
    const auto noFixed = std::vector<bool>(mesh.faces.size(), false);
    const auto stepLimiter = faceReconstructionLimiter(mesh, step, steep, noBoundary, noFixed);
    check(stepLimiter[4] < 1.0 && stepLimiter[4] >= 0.0,
          "step limiter reduces an overshooting interior reconstruction");
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        if (!face.neighbour || (face.owner != 4 && *face.neighbour != 4)) continue;
        const bool cellIsOwner = face.owner == 4;
        const double flux = cellIsOwner ? 1.0 : -1.0;
        const double reconstructed = upwindFaceValue(mesh, id, flux, step, steep, stepLimiter);
        const double neighbour = cellIsOwner ? step[*face.neighbour] : step[face.owner];
        const double lo = std::min(step[4], neighbour);
        const double hi = std::max(step[4], neighbour);
        check(reconstructed >= lo - 1e-12 && reconstructed <= hi + 1e-12,
              "limited step reconstruction stays within local stencil bounds");
    }
}

double exponential(Point2D point) {
    return std::exp(point.x + 0.3 * point.y);
}

Vector2D exponentialGradient(Point2D point) {
    const double value = exponential(point);
    return {value, 0.3 * value};
}

struct FaceRefinementError {
    double upwind = 0.0;
    double limited = 0.0;
};

FaceRefinementError exponentialFaceError(std::size_t n) {
    const auto mesh = normalizedSkewMesh(n);
    const auto field = cellValues(mesh, exponential);
    std::vector<Vector2D> gradients;
    gradients.reserve(mesh.cells.size());
    for (const auto& cell : mesh.cells) gradients.push_back(exponentialGradient(cell.centre));
    std::vector<double> boundary(mesh.faces.size());
    const std::vector<bool> fixed(mesh.faces.size(), true);
    std::vector<bool> interior(mesh.cells.size(), true);
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        for (const auto id : mesh.cells[i].faces)
            if (!mesh.faces[id].neighbour) interior[i] = false;
    }
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        boundary[id] = exponential(mesh.faces[id].centre);
    const auto limiter = faceReconstructionLimiter(mesh, field, gradients, boundary, fixed);
    double upwindSquared = 0.0;
    double limitedSquared = 0.0;
    std::size_t samples = 0;
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        if (!face.neighbour || !interior[face.owner] || !interior[*face.neighbour]) continue;
        const double geometricFlux = face.areaVector.x; // U=(1,0), fixed for both runs below.
        if (std::abs(geometricFlux) < 1e-14) continue;
        const double exact = exponential(face.centre);
        for (const double flux : {geometricFlux, -geometricFlux}) {
            const double upwind = upwindFaceValue(mesh, id, flux, field, gradients, {});
            const double limited = upwindFaceValue(mesh, id, flux, field, gradients, limiter);
            upwindSquared += std::pow(upwind - exact, 2);
            limitedSquared += std::pow(limited - exact, 2);
            ++samples;
        }
    }
    check(samples > 0, "exponential refinement has interior positive and negative flux samples");
    return {std::sqrt(upwindSquared / static_cast<double>(samples)),
            std::sqrt(limitedSquared / static_cast<double>(samples))};
}

void exponentialFaceRefinement() {
    const auto coarse = exponentialFaceError(8);
    const auto medium = exponentialFaceError(16);
    const auto fine = exponentialFaceError(32);
    const double upwindOrder = std::log(coarse.upwind / medium.upwind) / std::log(2.0);
    const double limitedOrder = std::log(coarse.limited / medium.limited) / std::log(2.0);
    const double upwindOrderFine = std::log(medium.upwind / fine.upwind) / std::log(2.0);
    const double limitedOrderFine = std::log(medium.limited / fine.limited) / std::log(2.0);
    std::cout << "exp face refinement: upwind " << coarse.upwind << " -> " << medium.upwind
              << " -> " << fine.upwind << " (orders " << upwindOrder << ", "
              << upwindOrderFine << "); limited " << coarse.limited << " -> "
              << medium.limited << " -> " << fine.limited << " (orders " << limitedOrder
              << ", " << limitedOrderFine << ")\n";
    check(medium.limited < medium.upwind && fine.limited < fine.upwind,
          "limited exponential face values improve over upwind values");
    check(limitedOrder > 1.5 && limitedOrderFine > 1.5,
          "limited exponential face reconstruction is approximately second order");
    check(upwindOrder < 1.4 && upwindOrderFine < 1.4,
          "upwind exponential face reconstruction remains approximately first order");
}

void splitFacesAreBothUsed(const FvMesh2D& mesh) {
    const auto pressure = cellValues(mesh, linear);
    const auto gradients = std::vector<Vector2D>(mesh.cells.size(), linearGradient());
    std::vector<double> boundary(mesh.faces.size());
    const std::vector<bool> fixed(mesh.faces.size(), true);
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        boundary[id] = linear(mesh.faces[id].centre);
    const auto facePressure = pressureFaceValues(mesh, pressure, gradients, boundary, fixed);
    std::size_t splitCount = 0;
    const std::vector<double> cellField{10.0, 20.0, 30.0};
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        if (!face.neighbour || std::abs(face.centre.x - 1.0) > 1e-12) continue;
        ++splitCount;
        near(facePressure[id], linear(face.centre), 1e-12,
             "each coarse/fine split face gets its own exact pressure value");
        near(upwindFaceValue(mesh, id, 1.0, cellField, {}, {}),
             cellField[face.owner], 1e-12, "positive flux uses split-face owner cell");
        near(upwindFaceValue(mesh, id, -1.0, cellField, {}, {}),
             cellField[*face.neighbour], 1e-12, "negative flux uses split-face neighbour cell");
    }
    check(splitCount == 2, "coarse side is represented by two actual shared faces");
}

void affineSymmetricStress(const FvMesh2D& mesh, const std::string& label) {
    constexpr double nu = 0.37;
    const std::vector<double> boundary(mesh.faces.size(), 0.0);
    const std::vector<bool> fixed(mesh.faces.size(), true);
    const std::vector<bool> constant(mesh.faces.size(), false);
    const auto cellField = [&](double ax, double ay, double bx, double by) {
        std::vector<double> value;
        for (const auto& c : mesh.cells) value.push_back(ax*c.centre.x + ay*c.centre.y + bx);
        return value;
    };
    const auto faceValue = [&](const cartmesh2d::fv::Face& f, double ax, double ay, double b) {
        return ax*f.centre.x + ay*f.centre.y + b;
    };
    const auto checkCase = [&](const std::string& name, Vector2D gu, Vector2D gv,
                               double bu0, double bv0) {
        const auto u = cellField(gu.x, gu.y, bu0, 0.0);
        const auto v = cellField(gv.x, gv.y, bv0, 0.0);
        std::vector<double> bu(mesh.faces.size()), bv(mesh.faces.size());
        for (std::size_t id=0; id<mesh.faces.size(); ++id) {
            bu[id]=faceValue(mesh.faces[id],gu.x,gu.y,bu0);
            bv[id]=faceValue(mesh.faces[id],gv.x,gv.y,bv0);
        }
        const std::vector<Vector2D> guCell(mesh.cells.size(),gu), gvCell(mesh.cells.size(),gv);
        const auto correction=symmetricViscousCorrection(mesh,u,v,guCell,gvCell,bu,bv,
                                                         fixed,fixed,constant,constant,nu);
        for (std::size_t id=0; id<mesh.faces.size(); ++id) {
            const auto& f=mesh.faces[id];
            if (!f.neighbour) continue;
            const auto j=*f.neighbour;
            const auto gfu=gu;
            const auto gfv=gv;
            const double baseU=-nu*(f.transmissibility*(u[j]-u[f.owner])+
                                    dot(gfu,f.correction));
            const double baseV=-nu*(f.transmissibility*(v[j]-v[f.owner])+
                                    dot(gfv,f.correction));
            const double expectedU=-nu*(2*gu.x*f.areaVector.x+(gu.y+gv.x)*f.areaVector.y);
            const double expectedV=-nu*((gu.y+gv.x)*f.areaVector.x+2*gv.y*f.areaVector.y);
            near(baseU+correction[id].x, expectedU, 1e-11,
                 label+" "+name+": x shared-face symmetric stress");
            near(baseV+correction[id].y, expectedV, 1e-11,
                 label+" "+name+": y shared-face symmetric stress");
        }
    };
    // Rigid rotation has an antisymmetric gradient and therefore zero stress.
    checkCase("rotation", {0.0,-1.3}, {1.3,0.0}, 0.0, 0.0);
    checkCase("extension", {1.2,0.0}, {0.0,-1.2}, 0.0, 0.0);
    checkCase("shear", {0.0,1.7}, {0.0,0.0}, 0.0, 0.0);
}

void wallTraceAndSlipStress(const FvMesh2D& mesh) {
    constexpr double nu=0.5;
    const auto u=cellValues(mesh, [](Point2D p) { return p.x; });
    const std::vector<double> v(mesh.cells.size(),0.0), bu(mesh.faces.size(),0.0), bv(mesh.faces.size(),0.0);
    const std::vector<Vector2D> gu(mesh.cells.size(),{1.0,0.0}), gv(mesh.cells.size(),{0.0,0.0});
    std::vector<bool> fu(mesh.faces.size(),true), fv(mesh.faces.size(),true);
    std::vector<bool> constantU(mesh.faces.size(),false), constantV(mesh.faces.size(),false);
    bool foundNormal=false;
    for(std::size_t id=0;id<mesh.faces.size();++id) if(!mesh.faces[id].neighbour &&
        std::abs(mesh.faces[id].areaVector.x)>0.9) {
        // A prescribed constant wall trace still retains the normal stress from
        // the prescribed owner-to-wall difference; its tangential trace is zero.
        constantU[id]=true;
        const auto& f=mesh.faces[id];
        const auto d=f.centre-mesh.cells[f.owner].centre;
        const double area=std::hypot(f.areaVector.x,f.areaVector.y);
        const auto n=f.areaVector*(1/area);
        const std::vector<Vector2D> misleadingGradient(mesh.cells.size(),{3.7,-2.4});
        const auto faceG=cartmesh2d::fv::detail::viscousFaceGradient(mesh,id,u,misleadingGradient,bu,fu,constantU);
        near(faceG.x*(-n.y)+faceG.y*n.x,0,1e-12,"constant wall removes cell tangential gradient");
        near(faceG.x*n.x+faceG.y*n.y,(bu[id]-u[f.owner])/(d.x*n.x+d.y*n.y),
             1e-12,"constant wall normal derivative matches prescribed value");
        const auto one=symmetricViscousCorrection(mesh,u,v,gu,gv,bu,bv,fu,fv,constantU,constantV,nu);
        near(one[id].y,0.0,1e-12,"constant wall trace has zero tangential shear");
        check(std::abs(one[id].x)>1e-12,"constant wall trace retains prescribed normal stress");
        foundNormal=true;
    }
    check(foundNormal,"wall trace test found an axis-aligned vertical wall");

    // At a horizontal slip wall, v is fixed to zero while u varies tangentially.
    const auto tangential=cellValues(mesh, [](Point2D p) { return p.x; });
    std::fill(fu.begin(),fu.end(),false); std::fill(fv.begin(),fv.end(),true);
    std::fill(constantU.begin(),constantU.end(),false); std::fill(constantV.begin(),constantV.end(),false);
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        if(!mesh.faces[id].neighbour && std::abs(mesh.faces[id].areaVector.x)<1e-12) {
            constantV[id]=true;
            const auto slip=symmetricViscousCorrection(mesh,tangential,v,gu,gv,bu,bv,fu,fv,constantU,constantV,nu);
            near(slip[id].x,0.0,1e-12,"axis-aligned slip has zero tangential shear");
        }
    }
}

} // namespace

int main() {
    try {
        const auto skew = rectangularMesh(2, 2, 0.3);
        constantPressureAndInternalCancellation(skew);
        linearPressureIsExact(skew, "skew mesh");
        reconstructedLinearPressureIsExact(skew, "skew mesh");
        unknownBoundaryAndRankChecks(skew);

        const auto split = splitFaceMesh();
        constantPressureAndInternalCancellation(split);
        linearPressureIsExact(split, "coarse/fine mesh");
        reconstructedLinearPressureIsExact(split, "coarse/fine mesh");
        splitFacesAreBothUsed(split);

        triangularTipUsesSecondRing();

        checkerboardPressureKeepsDirectDifference(rectangularMesh(2, 1));
        limiterAndUpwindSelection();
        exponentialFaceRefinement();
        affineSymmetricStress(skew, "skew mesh");
        wallTraceAndSlipStress(rectangularMesh(3, 3));
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "flow face operators: failures=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
