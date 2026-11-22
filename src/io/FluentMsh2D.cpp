#include "cartmesh2d/io/FluentMsh2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <utility>

namespace cartmesh2d {
namespace {

struct FluentFaceRecord2D {
    std::size_t n0 = 0;
    std::size_t n1 = 0;
    std::size_t c0 = 0;
    std::size_t c1 = 0;
};

struct FluentPatchBucket2D {
    std::string name;
    std::string bcName;
    int bcTypeCode = 3;
    BoundaryConditionRole2D role = BoundaryConditionRole2D::Wall;
    std::string sourceLayer;
    std::vector<FluentFaceRecord2D> faces;
};

void setError(std::string* error, std::string message) {
    if (error != nullptr) *error = std::move(message);
}

std::string toHex(std::size_t value) {
    std::ostringstream out;
    out << std::hex << value;
    return out.str();
}

std::pair<int, const char*> fluentBcForRole(BoundaryConditionRole2D role) noexcept {
    switch (role) {
        case BoundaryConditionRole2D::Wall: return {3, "wall"};
        case BoundaryConditionRole2D::Inlet: return {10, "velocity-inlet"};
        case BoundaryConditionRole2D::Outlet: return {5, "pressure-outlet"};
        case BoundaryConditionRole2D::Slip:
        case BoundaryConditionRole2D::Symmetry: return {7, "symmetry"};
    }
    return {3, "wall"};
}

[[nodiscard]] bool edgeDirectionForCell(const TopologyCell2D& cell,
                                        std::size_t v0, std::size_t v1,
                                        std::size_t& from, std::size_t& to) noexcept {
    const std::size_t count = cell.vertices.size();
    for (std::size_t i = 0; i < count; ++i) {
        const std::size_t a = cell.vertices[i];
        const std::size_t b = cell.vertices[(i + 1U) % count];
        if ((a == v0 && b == v1) || (a == v1 && b == v0)) {
            from = a;
            to = b;
            return true;
        }
    }
    return false;
}

[[nodiscard]] std::string domainPatchName(const Edge2D& edge,
                                          const TopologyMesh2D& topology,
                                          const Domain2D& domain,
                                          const TolerancePolicy& tol) {
    const auto& a = topology.vertices[edge.v0].point;
    const auto& b = topology.vertices[edge.v1].point;
    const auto near = [&](double x, double y) {
        return std::abs(x - y) <= tol.scale(std::max(std::abs(x), std::abs(y)));
    };
    if (near(a.x, domain.bounds.min.x) && near(b.x, domain.bounds.min.x)) return "left";
    if (near(a.x, domain.bounds.max.x) && near(b.x, domain.bounds.max.x)) return "right";
    if (near(a.y, domain.bounds.min.y) && near(b.y, domain.bounds.min.y)) return "bottom";
    if (near(a.y, domain.bounds.max.y) && near(b.y, domain.bounds.max.y)) return "top";
    return {};
}

[[nodiscard]] FluentPatchBucket2D* findBucket(std::vector<FluentPatchBucket2D>& buckets,
                                              const std::string& name) {
    for (auto& bucket : buckets) if (bucket.name == name) return &bucket;
    return nullptr;
}

[[nodiscard]] double pointSegmentDistance(const Point2D& point,
                                          const Segment2D& segment) noexcept {
    const Vector2D direction = segment.b - segment.a;
    const double lengthSquared = squaredNorm(direction);
    if (!(lengthSquared > 0.0)) return std::sqrt(squaredNorm(point - segment.a));
    const double parameter = std::clamp(
        dot(point - segment.a, direction) / lengthSquared, 0.0, 1.0);
    return std::sqrt(squaredNorm(point - (segment.a + direction * parameter)));
}

[[nodiscard]] std::optional<std::size_t> embeddedLoopId(
    const Edge2D& edge, const TopologyMesh2D& topology,
    const BoundaryRegion2D& boundary, const TolerancePolicy& tol) {
    const Segment2D face{topology.vertices[edge.v0].point,
                         topology.vertices[edge.v1].point};
    const Vector2D direction = face.b - face.a;
    const double lengthSquared = squaredNorm(direction);
    if (!(lengthSquared > 0.0)) return std::nullopt;
    const double parameterTolerance = tol.scale(std::sqrt(lengthSquared)) / std::sqrt(lengthSquared);
    std::optional<std::size_t> exactLoop;
    for (std::size_t loopId = 0; loopId < boundary.loops().size(); ++loopId) {
        const auto& vertices = boundary.loops()[loopId].vertices();
        std::vector<std::pair<double, double>> intervals;
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            const auto intersection = intersectSegments(
                face, {vertices[i], vertices[(i + 1U) % vertices.size()]}, tol);
            if (!intersection.overlap) continue;
            double a = dot(intersection.overlap->a - face.a, direction) / lengthSquared;
            double b = dot(intersection.overlap->b - face.a, direction) / lengthSquared;
            if (a > b) std::swap(a, b);
            intervals.emplace_back(a, b);
        }
        std::sort(intervals.begin(), intervals.end());
        double covered = 0.0;
        for (const auto& [a, b] : intervals) {
            if (a > covered + parameterTolerance) break;
            covered = std::max(covered, b);
        }
        if (covered < 1.0 - parameterTolerance) continue;
        if (exactLoop) return std::nullopt;
        exactLoop = loopId;
    }
    if (exactLoop) return exactLoop;

    if (!topology.constructionRegistry ||
        topology.canonicalVertexIds.size() != topology.vertices.size()) {
        return std::nullopt;
    }
    const auto& registry = *topology.constructionRegistry;
    const auto handle0 = topology.canonicalVertexIds[edge.v0];
    const auto handle1 = topology.canonicalVertexIds[edge.v1];
    if (handle0 >= registry.vertices().size() || handle1 >= registry.vertices().size()) {
        return std::nullopt;
    }
    const double arithmeticBudget =
        tol.scale(std::max({1.0, std::abs(face.a.x), std::abs(face.a.y),
                            std::abs(face.b.x), std::abs(face.b.y)}));
    double distanceBudgetA = arithmeticBudget;
    double distanceBudgetB = arithmeticBudget;
    for (const auto& event : registry.events()) {
        if (event.canonicalVertex == handle0)
            distanceBudgetA = std::max(distanceBudgetA, event.displacement + arithmeticBudget);
        if (event.canonicalVertex == handle1)
            distanceBudgetB = std::max(distanceBudgetB, event.displacement + arithmeticBudget);
    }
    std::optional<std::size_t> matchedLoop;
    for (std::size_t loopId = 0; loopId < boundary.loops().size(); ++loopId) {
        const auto& vertices = boundary.loops()[loopId].vertices();
        bool matches = false;
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            const Segment2D segment{vertices[i], vertices[(i + 1U) % vertices.size()]};
            if (pointSegmentDistance(face.a, segment) <= distanceBudgetA &&
                pointSegmentDistance(face.b, segment) <= distanceBudgetB) {
                matches = true;
                break;
            }
        }
        if (!matches) continue;
        if (matchedLoop && *matchedLoop != loopId) return std::nullopt;
        matchedLoop = loopId;
    }
    return matchedLoop;
}

} // namespace

FluentMshWriteReport2D writeFluentMsh2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const std::filesystem::path& mshFile,
    std::string* error,
    const TolerancePolicy& tol,
    const std::vector<EmbeddedBoundaryPatch2D>& embeddedPatches) {
    FluentMshWriteReport2D report;
    if (!topology.valid() || topology.cells.empty() || topology.vertices.empty() ||
        topology.edges.empty() || !domain.valid(tol) || !boundary.diagnose(tol).valid()) {
        setError(error, "Fluent 2D export requires valid non-empty topology, domain, and boundary");
        return report;
    }
    if (!embeddedPatches.empty() && embeddedPatches.size() != boundary.loops().size()) {
        setError(error, "embedded boundary patch metadata count does not match loops");
        return report;
    }

    // Validate finite vertex coordinates and positive cell polygon areas.
    for (const auto& vertex : topology.vertices) {
        if (!std::isfinite(vertex.point.x) || !std::isfinite(vertex.point.y)) {
            setError(error, "Fluent 2D export encountered non-finite vertex coordinate");
            return {};
        }
    }

    std::vector<int> cellTypes;
    cellTypes.reserve(topology.cells.size());
    std::size_t triangles = 0U, quads = 0U, polygons = 0U;
    for (const auto& cell : topology.cells) {
        if (cell.vertices.size() < 3U) {
            setError(error, "Fluent 2D export encountered cell with fewer than 3 vertices");
            return {};
        }
        Polygon2D poly;
        poly.vertices.reserve(cell.vertices.size());
        for (const auto vid : cell.vertices) {
            if (vid >= topology.vertices.size()) {
                setError(error, "Fluent 2D export cell vertex ID out of range");
                return {};
            }
            poly.vertices.push_back(topology.vertices[vid].point);
        }
        if (!(poly.signedArea() > 0.0)) {
            setError(error, "Fluent 2D export requires positive CCW cell orientation");
            return {};
        }
        if (cell.vertices.size() == 3U) { cellTypes.push_back(1); ++triangles; }
        else if (cell.vertices.size() == 4U) { cellTypes.push_back(3); ++quads; }
        else { cellTypes.push_back(7); ++polygons; }
    }

    // Build boundary patch buckets matching OpenFoam2D patch naming & roles.
    std::vector<FluentPatchBucket2D> boundaryBuckets;
    std::vector<std::string> loopPatchNames;
    loopPatchNames.reserve(boundary.loops().size());
    for (std::size_t loopId = 0; loopId < boundary.loops().size(); ++loopId) {
        EmbeddedBoundaryPatch2D spec;
        if (embeddedPatches.empty()) {
            spec = {"wall_" + std::to_string(loopId), "wall", BoundaryConditionRole2D::Wall, ""};
        } else {
            spec = embeddedPatches[loopId];
        }
        if (!spec.valid() || spec.name == "interior" || spec.name == "fluid" ||
            spec.name == "left" || spec.name == "right" ||
            spec.name == "bottom" || spec.name == "top") {
            setError(error, "invalid or reserved embedded boundary patch metadata");
            return {};
        }
        const auto [bcCode, bcName] = fluentBcForRole(spec.role);
        if (const auto* existing = findBucket(boundaryBuckets, spec.name); existing != nullptr) {
            if (existing->role != spec.role || existing->sourceLayer != spec.sourceLayer) {
                setError(error, "duplicate embedded patch name has conflicting metadata");
                return {};
            }
        } else {
            boundaryBuckets.push_back({spec.name, bcName, bcCode, spec.role, spec.sourceLayer, {}});
        }
        loopPatchNames.push_back(spec.name);
    }

    const std::vector<FluentPatchBucket2D> domainBuckets{
        {"left", "velocity-inlet", 10, BoundaryConditionRole2D::Inlet, "", {}},
        {"right", "pressure-outlet", 5, BoundaryConditionRole2D::Outlet, "", {}},
        {"bottom", "symmetry", 7, BoundaryConditionRole2D::Slip, "", {}},
        {"top", "symmetry", 7, BoundaryConditionRole2D::Slip, "", {}}
    };
    boundaryBuckets.insert(boundaryBuckets.end(), domainBuckets.begin(), domainBuckets.end());

    std::vector<FluentFaceRecord2D> internalFaces;
    internalFaces.reserve(topology.edges.size());

    for (const auto& edge : topology.edges) {
        if (edge.owner >= topology.cells.size()) {
            setError(error, "Fluent 2D edge owner out of range");
            return {};
        }
        std::size_t from = 0U, to = 0U;
        if (!edgeDirectionForCell(topology.cells[edge.owner], edge.v0, edge.v1, from, to)) {
            setError(error, "Fluent 2D edge absent from owner cell loop");
            return {};
        }
        // In Fluent's 2D face-based format, walking n0 -> n1 keeps c0 on the
        // left (matching our CCW owner loop, whose outward normal is (dy,-dx)).
        const auto& p0 = topology.vertices[from].point;
        const auto& p1 = topology.vertices[to].point;
        if (!(squaredNorm(p1 - p0) > 0.0)) {
            setError(error, "Fluent 2D export encountered zero-length face");
            return {};
        }
        if (edge.neighbour) {
            if (*edge.neighbour >= topology.cells.size() || *edge.neighbour == edge.owner) {
                setError(error, "Fluent 2D edge neighbour invalid");
                return {};
            }
            internalFaces.push_back({from + 1U, to + 1U, edge.owner + 1U, *edge.neighbour + 1U});
            continue;
        }

        std::string patchName;
        if (edge.patch == BoundaryPatch2D::EmbeddedBoundary) {
            const auto loopId = embeddedLoopId(edge, topology, boundary, tol);
            if (!loopId) {
                setError(error, "Fluent 2D embedded boundary edge could not be matched to a unique loop");
                return {};
            }
            patchName = loopPatchNames[*loopId];
        } else if (edge.patch == BoundaryPatch2D::DomainBoundary) {
            patchName = domainPatchName(edge, topology, domain, tol);
            if (patchName.empty()) {
                setError(error, "Fluent 2D domain boundary edge does not lie on domain rectangle");
                return {};
            }
        } else {
            setError(error, "Fluent 2D boundary edge has unclassified patch");
            return {};
        }
        auto* bucket = findBucket(boundaryBuckets, patchName);
        if (bucket == nullptr) {
            setError(error, "Fluent 2D boundary patch bucket not found: " + patchName);
            return {};
        }
        bucket->faces.push_back({from + 1U, to + 1U, edge.owner + 1U, 0U});
    }

    const std::filesystem::path parent = mshFile.parent_path();
    if (!parent.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(parent, ec);
        if (ec) {
            setError(error, "failed to create output directory for Fluent .msh: " + ec.message());
            return {};
        }
    }

    std::ofstream out(mshFile);
    if (!out) {
        setError(error, "failed to open Fluent .msh file for writing: " + mshFile.string());
        return {};
    }

    out << "(0 \"CartMesh2D native 2D face-based Fluent ASCII mesh\")\n";
    out << "(2 2)\n";

    // Nodes: zone 1, 2D coordinates at exact 17-digit precision.
    const std::size_t nNodes = topology.vertices.size();
    out << "(10 (0 1 " << toHex(nNodes) << " 0))\n";
    out << "(10 (1 1 " << toHex(nNodes) << " 1 2)(\n";
    out << std::setprecision(17);
    for (const auto& vertex : topology.vertices) {
        out << vertex.point.x << ' ' << vertex.point.y << '\n';
    }
    out << "))\n";

    // Cells: zone 2 (fluid), element-type 0 (mixed list: 1=tri, 3=quad, 7=poly).
    const std::size_t nCells = topology.cells.size();
    out << "(12 (0 1 " << toHex(nCells) << " 0))\n";
    out << "(12 (2 1 " << toHex(nCells) << " 1 0)(\n";
    for (std::size_t i = 0; i < nCells; ++i) {
        out << cellTypes[i] << ((i + 1U) % 32U == 0U || i + 1U == nCells ? '\n' : ' ');
    }
    out << "))\n";

    // Faces: header + interior zone 3 + non-empty boundary zones 4, 5, ...
    const std::size_t nFaces = topology.edges.size();
    out << "(13 (0 1 " << toHex(nFaces) << " 0))\n";

    std::size_t nextFace = 1U;
    std::size_t nextZoneId = 3U;
    std::vector<FluentZoneSummary2D> zones;

    if (!internalFaces.empty()) {
        const std::size_t start = nextFace;
        const std::size_t end = start + internalFaces.size() - 1U;
        const std::size_t zoneId = nextZoneId++;
        out << "(13 (" << toHex(zoneId) << ' ' << toHex(start) << ' ' << toHex(end) << " 2 2)(\n";
        for (const auto& f : internalFaces) {
            out << toHex(f.n0) << ' ' << toHex(f.n1) << ' '
                << toHex(f.c0) << ' ' << toHex(f.c1) << '\n';
        }
        out << "))\n";
        zones.push_back({zoneId, "interior", "interior", BoundaryConditionRole2D::Wall, internalFaces.size(), start});
        nextFace = end + 1U;
    }

    for (const auto& bucket : boundaryBuckets) {
        if (bucket.faces.empty()) continue;
        const std::size_t start = nextFace;
        const std::size_t end = start + bucket.faces.size() - 1U;
        const std::size_t zoneId = nextZoneId++;
        out << "(13 (" << toHex(zoneId) << ' ' << toHex(start) << ' ' << toHex(end)
            << ' ' << toHex(static_cast<std::size_t>(bucket.bcTypeCode)) << " 2)(\n";
        for (const auto& f : bucket.faces) {
            out << toHex(f.n0) << ' ' << toHex(f.n1) << ' '
                << toHex(f.c0) << " 0\n";
        }
        out << "))\n";
        zones.push_back({zoneId, bucket.name, bucket.bcName, bucket.role, bucket.faces.size(), start});
        nextFace = end + 1U;
    }

    // Zone section (45 ...): zone IDs are DECIMAL per Fluent ASCII mesh spec.
    out << "(45 (2 fluid fluid)())\n";
    for (const auto& zone : zones) {
        out << "(45 (" << zone.zoneId << ' ' << zone.bcType << ' ' << zone.name << ")())\n";
    }
    out.flush();
    if (!out.good()) {
        setError(error, "failed while writing Fluent .msh file");
        return {};
    }

    report.nodeCount = nNodes;
    report.cellCount = nCells;
    report.faceCount = nFaces;
    report.internalFaceCount = internalFaces.size();
    report.triangleCount = triangles;
    report.quadCount = quads;
    report.polygonCount = polygons;
    report.zones = std::move(zones);
    return report;
}

} // namespace cartmesh2d
