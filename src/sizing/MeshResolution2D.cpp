#include "cartmesh2d/sizing/MeshResolution2D.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace cartmesh2d {
namespace {
void optionalNumber(std::ostream& out, std::optional<double> value) {
    if (value) out << *value;
    else out << "null";
}

void distribution(std::ostream& out, std::vector<double> values) {
    if (values.empty()) { out << "null"; return; }
    std::sort(values.begin(), values.end());
    const auto quantile = [&values](double p) {
        const double index = p * static_cast<double>(values.size() - 1U);
        const auto lower = static_cast<std::size_t>(std::floor(index));
        const auto upper = static_cast<std::size_t>(std::ceil(index));
        return values[lower] + (index - static_cast<double>(lower)) * (values[upper] - values[lower]);
    };
    out << "{\"min\":" << values.front() << ",\"p50\":" << quantile(0.5)
        << ",\"p95\":" << quantile(0.95) << ",\"max\":" << values.back() << '}';
}
} // namespace

std::string meshResolutionReportToJson2D(const TopologyMesh2D& mesh,
                                         const MeshResolutionTargets2D& targets,
                                         const TolerancePolicy& tol) {
    const double reference = targets.referenceLength;
    if (!(reference > 0.0) || !std::isfinite(reference)) {
        throw std::invalid_argument("resolution reference length must be finite and positive");
    }
    for (const auto value : {targets.wallSize, targets.backgroundSize, targets.firstLayerHeight}) {
        if (value && (!std::isfinite(*value) || !(*value > 0.0))) {
            throw std::invalid_argument("resolution target must be finite and positive");
        }
    }
    std::vector<double> areaSize, wallEdgeSize, wallTangent, wallNormal;
    std::vector<std::size_t> exceeded;
    double wallLength = 0.0, exceededLength = 0.0;
    areaSize.reserve(mesh.cells.size());
    for (const auto& cell : mesh.cells) {
        if (!(cell.geometryArea > 0.0) || !std::isfinite(cell.geometryArea)) {
            throw std::invalid_argument("invalid cell area in resolution measurement");
        }
        areaSize.push_back(std::sqrt(cell.geometryArea) / reference);
    }
    for (const auto& edge : mesh.edges) {
        if (edge.neighbour || edge.patch != BoundaryPatch2D::EmbeddedBoundary) continue;
        const auto a = mesh.vertices.at(edge.v0).point;
        const auto b = mesh.vertices.at(edge.v1).point;
        const double length = std::hypot(b.x-a.x, b.y-a.y);
        if (!(length > 0.0) || !std::isfinite(length)) {
            throw std::invalid_argument("invalid wall edge in resolution measurement");
        }
        const double tx = (b.x-a.x)/length, ty = (b.y-a.y)/length;
        double tmin = std::numeric_limits<double>::infinity(), tmax = -tmin;
        double nmin = tmin, nmax = -tmin;
        const auto& cell = mesh.cells.at(edge.owner);
        for (const auto vertex : cell.vertices) {
            const auto p = mesh.vertices.at(vertex).point;
            const double x = p.x-a.x, y = p.y-a.y;
            const double t = x*tx+y*ty, n = -x*ty+y*tx;
            tmin = std::min(tmin,t); tmax = std::max(tmax,t);
            nmin = std::min(nmin,n); nmax = std::max(nmax,n);
        }
        wallLength += length;
        wallEdgeSize.push_back(length/reference);
        wallTangent.push_back((tmax-tmin)/reference);
        wallNormal.push_back((nmax-nmin)/reference);
        // Compare nondimensional lengths so the tolerance cannot change merely
        // because a user supplies metres instead of millimetres.
        if (targets.wallSize && (tmax-tmin)/reference > *targets.wallSize/reference &&
            !tol.nearlyEqual((tmax-tmin)/reference, *targets.wallSize/reference)) {
            exceeded.push_back(edge.id);
            exceededLength += length;
        }
    }
    std::ostringstream out;
    out << std::setprecision(17)
        << "{\n  \"format\":\"cartmesh2d-mesh-resolution-v1\",\n"
        << "  \"basis\":\"final_solver_topology\",\n"
        << "  \"reference_length\":" << reference << ",\n"
        << "  \"reference_source\":\"" << (targets.explicitReferenceLength ? "explicit" : "bounding_box_span") << "\",\n"
        << "  \"cell_count\":" << mesh.cells.size() << ",\n"
        << "  \"requested\":{\"wall_h_over_reference\":";
    optionalNumber(out, targets.wallSize ? std::optional<double>(*targets.wallSize/reference) : std::nullopt);
    out << ",\"background_h_over_reference\":";
    optionalNumber(out, targets.backgroundSize ? std::optional<double>(*targets.backgroundSize/reference) : std::nullopt);
    out << ",\"first_layer_h_over_reference\":";
    optionalNumber(out, targets.firstLayerHeight ? std::optional<double>(*targets.firstLayerHeight/reference) : std::nullopt);
    out << "},\n  \"actual\":{\"sqrt_area_over_reference\":";
    distribution(out, std::move(areaSize));
    out << ",\"wall_edge_length_over_reference\":";
    distribution(out, std::move(wallEdgeSize));
    out << ",\"wall_owner_tangential_extent_over_reference\":";
    distribution(out, std::move(wallTangent));
    out << ",\"wall_owner_normal_extent_over_reference\":";
    distribution(out, std::move(wallNormal));
    out << "},\n  \"wall_length_over_reference\":" << wallLength/reference
        << ",\n  \"wall_owner_tangential_exceedance_length_fraction\":";
    optionalNumber(out, targets.wallSize && wallLength > 0 ? std::optional<double>(exceededLength/wallLength) : std::nullopt);
    out << ",\n  \"wall_owner_tangential_exceedance_edge_ids\":[";
    for (std::size_t i = 0; i < exceeded.size(); ++i) { if (i) out << ','; out << exceeded[i]; }
    out << "],\n  \"boundary_layer_coverage_status\":\"not_evaluated\",\n"
        << "  \"notes\":[\"Wall statistics sample final embedded edges; split edges do not replace owner extents.\","
        << "\"Normal owner extent is not a certified first-layer height or y-plus.\","
        << "\"Sizing exceedance is a diagnostic independent of topology, Solver and Q1 gates.\"]\n}\n";
    return out.str();
}
} // namespace cartmesh2d
