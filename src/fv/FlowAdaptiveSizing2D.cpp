#include "cartmesh2d/fv/FlowAdaptiveSizing2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>

namespace cartmesh2d::fv {
namespace {

bool finiteValues(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); });
}

AABB2D unionBox(const AABB2D& a, const AABB2D& b) noexcept {
    return {
        {std::min(a.min.x, b.min.x), std::min(a.min.y, b.min.y)},
        {std::max(a.max.x, b.max.x), std::max(a.max.y, b.max.y)}
    };
}

double boxArea(const AABB2D& box) noexcept {
    const double dx = std::max(0.0, box.max.x - box.min.x);
    const double dy = std::max(0.0, box.max.y - box.min.y);
    return dx * dy;
}

double overlapRatio(const AABB2D& a, const AABB2D& b) noexcept {
    const double ix = std::min(a.max.x, b.max.x) - std::max(a.min.x, b.min.x);
    const double iy = std::min(a.max.y, b.max.y) - std::max(a.min.y, b.min.y);
    if (!(ix > 0.0) || !(iy > 0.0)) return 0.0;
    const double minArea = std::min(boxArea(a), boxArea(b));
    return minArea > 0.0 ? (ix * iy) / minArea : 0.0;
}

AABB2D clampBox(const AABB2D& box, const AABB2D& limit) noexcept {
    const double xmin = std::clamp(box.min.x, limit.min.x, limit.max.x);
    const double ymin = std::clamp(box.min.y, limit.min.y, limit.max.y);
    const double xmax = std::clamp(box.max.x, limit.min.x, limit.max.x);
    const double ymax = std::clamp(box.max.y, limit.min.y, limit.max.y);
    return {{std::min(xmin, xmax), std::min(ymin, ymax)},
            {std::max(xmin, xmax), std::max(ymin, ymax)}};
}

double percentileSorted(const std::vector<double>& sorted, double fraction) noexcept {
    if (sorted.empty()) return 0.0;
    const double idx = std::clamp(fraction, 0.0, 1.0) * static_cast<double>(sorted.size() - 1U);
    const auto lo = static_cast<std::size_t>(std::floor(idx));
    const auto hi = static_cast<std::size_t>(std::ceil(idx));
    const double t = idx - static_cast<double>(lo);
    return sorted[lo] * (1.0 - t) + sorted[hi] * t;
}

AABB2D cellPolygonBounds(const TopologyMesh2D& topology, std::size_t cellId) noexcept {
    AABB2D box{
        {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}
    };
    for (const auto vid : topology.cells[cellId].vertices) {
        const auto p = topology.vertices[vid].point;
        box.min.x = std::min(box.min.x, p.x);
        box.min.y = std::min(box.min.y, p.y);
        box.max.x = std::max(box.max.x, p.x);
        box.max.y = std::max(box.max.y, p.y);
    }
    return box;
}

struct Cluster {
    std::vector<std::size_t> cells;
    AABB2D tightBounds;
    double areaSum = 0.0;
};

Cluster makeCluster(const std::vector<std::size_t>& cellIds,
                    const std::vector<AABB2D>& cellBounds,
                    const FvMesh2D& fvMesh) {
    Cluster cluster;
    cluster.cells = cellIds;
    cluster.tightBounds = cellBounds[cellIds.front()];
    for (const auto id : cellIds) {
        cluster.tightBounds = unionBox(cluster.tightBounds, cellBounds[id]);
        cluster.areaSum += fvMesh.cells[id].area;
    }
    return cluster;
}

void splitSparseCluster(const Cluster& cluster,
                        const std::vector<AABB2D>& cellBounds,
                        const FvMesh2D& fvMesh,
                        double fillThreshold,
                        std::size_t depth,
                        std::vector<Cluster>& out) {
    const double area = boxArea(cluster.tightBounds);
    const double fill = area > 0.0 ? cluster.areaSum / area : 1.0;
    if (depth >= 3U || cluster.cells.size() < 8U || fill >= fillThreshold) {
        out.push_back(cluster);
        return;
    }
    const double dx = cluster.tightBounds.max.x - cluster.tightBounds.min.x;
    const double dy = cluster.tightBounds.max.y - cluster.tightBounds.min.y;
    const bool splitX = dx >= dy;
    const double pivot = splitX
        ? 0.5 * (cluster.tightBounds.min.x + cluster.tightBounds.max.x)
        : 0.5 * (cluster.tightBounds.min.y + cluster.tightBounds.max.y);

    std::vector<std::size_t> left, right;
    left.reserve(cluster.cells.size());
    right.reserve(cluster.cells.size());
    for (const auto id : cluster.cells) {
        const auto c = fvMesh.cells[id].centre;
        const double coord = splitX ? c.x : c.y;
        if (coord < pivot) left.push_back(id);
        else right.push_back(id);
    }
    if (left.empty() || right.empty()) {
        out.push_back(cluster);
        return;
    }
    splitSparseCluster(makeCluster(left, cellBounds, fvMesh), cellBounds, fvMesh, fillThreshold, depth + 1U, out);
    splitSparseCluster(makeCluster(right, cellBounds, fvMesh), cellBounds, fvMesh, fillThreshold, depth + 1U, out);
}

std::size_t inferCellLevel(double domainSpan, double cellH) noexcept {
    if (!(domainSpan > 0.0) || !(cellH > 0.0)) return 0U;
    const double raw = std::log2(domainSpan / cellH);
    if (!std::isfinite(raw) || raw <= 0.0) return 0U;
    return static_cast<std::size_t>(std::clamp(std::llround(raw), 0LL, 28LL));
}

} // namespace

FlowAdaptiveRefinement2D suggestFlowAdaptiveRefinement2D(
    const TopologyMesh2D& topology,
    const FvMesh2D& fvMesh,
    const std::vector<double>& u,
    const std::vector<double>& v,
    const std::vector<double>& p,
    const FlowAdaptivePolicy2D& policy,
    const std::vector<double>& boundaryU,
    const std::vector<double>& boundaryV,
    const std::vector<double>& boundaryP,
    const std::vector<double>& scalar,
    const std::vector<double>& boundaryScalar) {
    FlowAdaptiveRefinement2D result;
    const std::size_t nCells = fvMesh.cells.size();
    const std::size_t nFaces = fvMesh.faces.size();

    if (!topology.valid() || topology.cells.empty() ||
        topology.cells.size() != nCells || topology.edges.size() != nFaces) {
        result.issues.push_back("flow adaptive sizing requires matching valid topology and FVM mesh");
        return result;
    }
    if (u.size() != nCells || v.size() != nCells || p.size() != nCells) {
        result.issues.push_back("velocity/pressure field lengths must match mesh cell count");
        return result;
    }
    if (!finiteValues(u) || !finiteValues(v) || !finiteValues(p)) {
        result.issues.push_back("velocity/pressure fields contain non-finite values");
        return result;
    }
    if (!scalar.empty() && (scalar.size() != nCells || !finiteValues(scalar))) {
        result.issues.push_back("optional scalar/temperature field must match cell count and be finite");
        return result;
    }
    if ((!boundaryU.empty() && (boundaryU.size() != nFaces || !finiteValues(boundaryU))) ||
        (!boundaryV.empty() && (boundaryV.size() != nFaces || !finiteValues(boundaryV))) ||
        (!boundaryP.empty() && (boundaryP.size() != nFaces || !finiteValues(boundaryP))) ||
        (!boundaryScalar.empty() && (boundaryScalar.size() != nFaces || !finiteValues(boundaryScalar)))) {
        result.issues.push_back("optional boundary field values must match face count and be finite");
        return result;
    }
    if (!std::isfinite(policy.refineFraction) || !(policy.refineFraction > 0.0) || policy.refineFraction > 0.5 ||
        !std::isfinite(policy.minRelativePeak) || policy.minRelativePeak < 0.0 || policy.minRelativePeak >= 1.0 ||
        policy.maxBoxes == 0U || policy.maxBoxes > 32U ||
        !std::isfinite(policy.padCells) || policy.padCells < 0.0 || policy.padCells > 8.0 ||
        !std::isfinite(policy.splitFillThreshold) || policy.splitFillThreshold < 0.0 || policy.splitFillThreshold > 1.0 ||
        !std::isfinite(policy.velocityWeight) || policy.velocityWeight < 0.0 ||
        !std::isfinite(policy.pressureWeight) || policy.pressureWeight < 0.0 ||
        !std::isfinite(policy.scalarWeight) || policy.scalarWeight < 0.0 ||
        (policy.velocityWeight == 0.0 && policy.pressureWeight == 0.0 &&
         (scalar.empty() || policy.scalarWeight == 0.0)) ||
        (policy.maxTargetLevel && *policy.maxTargetLevel > 28U) ||
        (policy.wallLevel && *policy.wallLevel > 28U)) {
        result.issues.push_back("flow adaptive sizing policy parameters are out of valid range");
        return result;
    }
    if (policy.bodyBounds) {
        const auto& b = *policy.bodyBounds;
        if (!std::isfinite(b.min.x) || !std::isfinite(b.min.y) ||
            !std::isfinite(b.max.x) || !std::isfinite(b.max.y) ||
            !(b.max.x > b.min.x) || !(b.max.y > b.min.y)) {
            result.issues.push_back("explicit body bounds must be finite with positive width and height");
            return result;
        }
    }

    result.cellCount = nCells;

    // 1. Domain bounds, cell AABBs, and body frame.
    std::vector<AABB2D> cellBounds(nCells);
    AABB2D domainBox = cellPolygonBounds(topology, 0U);
    double totalArea = 0.0;
    for (std::size_t i = 0; i < nCells; ++i) {
        cellBounds[i] = cellPolygonBounds(topology, i);
        domainBox = unionBox(domainBox, cellBounds[i]);
        totalArea += fvMesh.cells[i].area;
    }
    result.domainBounds = domainBox;
    const double domainSpan = std::max(domainBox.max.x - domainBox.min.x, domainBox.max.y - domainBox.min.y);

    AABB2D embeddedBox{
        {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()},
        {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}
    };
    bool hasEmbedded = false;
    std::vector<bool> wallCell(nCells, false);
    for (std::size_t f = 0; f < nFaces; ++f) {
        const auto& edge = topology.edges[f];
        if (!edge.neighbour && edge.patch == BoundaryPatch2D::EmbeddedBoundary) {
            hasEmbedded = true;
            wallCell[edge.owner] = true;
            const auto p0 = topology.vertices[edge.v0].point;
            const auto p1 = topology.vertices[edge.v1].point;
            embeddedBox.min.x = std::min({embeddedBox.min.x, p0.x, p1.x});
            embeddedBox.min.y = std::min({embeddedBox.min.y, p0.y, p1.y});
            embeddedBox.max.x = std::max({embeddedBox.max.x, p0.x, p1.x});
            embeddedBox.max.y = std::max({embeddedBox.max.y, p0.y, p1.y});
        }
    }

    AABB2D bodyBox = domainBox;
    if (policy.bodyBounds) {
        bodyBox = *policy.bodyBounds;
        result.bodyFromEmbeddedWall = false;
    } else if (hasEmbedded &&
               std::max(embeddedBox.max.x - embeddedBox.min.x, embeddedBox.max.y - embeddedBox.min.y) > 0.0) {
        bodyBox = embeddedBox;
        result.bodyFromEmbeddedWall = true;
    }
    result.bodyCentre = {0.5 * (bodyBox.min.x + bodyBox.max.x), 0.5 * (bodyBox.min.y + bodyBox.max.y)};
    result.bodySpan = std::max(bodyBox.max.x - bodyBox.min.x, bodyBox.max.y - bodyBox.min.y);
    if (!(result.bodySpan > 0.0)) result.bodySpan = std::max(domainSpan, 1.0);

    // 2. Infer local grid levels on the solved mesh using each cell's spatial
    // span (avoiding artificially inflated levels on small cut-cell slivers).
    std::vector<std::size_t> cellLevels(nCells, 0U);
    std::vector<std::size_t> wallLevels;
    std::vector<std::size_t> allLevels;
    allLevels.reserve(nCells);
    for (std::size_t i = 0; i < nCells; ++i) {
        const double spanX = cellBounds[i].max.x - cellBounds[i].min.x;
        const double spanY = cellBounds[i].max.y - cellBounds[i].min.y;
        const double cellSpan = std::max(spanX, spanY);
        const double h = cellSpan > 0.0 ? cellSpan : std::sqrt(std::max(fvMesh.cells[i].area, 0.0));
        cellLevels[i] = inferCellLevel(domainSpan, h);
        allLevels.push_back(cellLevels[i]);
        if (wallCell[i]) wallLevels.push_back(cellLevels[i]);
    }
    std::sort(allLevels.begin(), allLevels.end());
    if (!wallLevels.empty()) {
        std::sort(wallLevels.begin(), wallLevels.end());
        result.inferredWallLevel = wallLevels[wallLevels.size() / 2U];
    } else {
        result.inferredWallLevel = allLevels.back();
    }
    const std::size_t capLevel = policy.maxTargetLevel.value_or(result.inferredWallLevel);
    const std::size_t refWallLevel = policy.wallLevel.value_or(capLevel);
    result.maxTargetLevel = capLevel;

    // 3. Reconstruct gradients on actual polygon geometry using FvMesh2D.
    std::vector<double> bU = boundaryU;
    std::vector<double> bV = boundaryV;
    std::vector<double> bP = boundaryP;
    if (bU.empty()) {
        bU.resize(nFaces, 0.0);
        for (std::size_t f = 0; f < nFaces; ++f) {
            const auto& face = fvMesh.faces[f];
            if (!face.neighbour) {
                const bool noSlipWall = policy.embeddedWallsNoSlip && face.patch == BoundaryPatch2D::EmbeddedBoundary;
                bU[f] = noSlipWall ? 0.0 : u[face.owner];
            }
        }
    }
    if (bV.empty()) {
        bV.resize(nFaces, 0.0);
        for (std::size_t f = 0; f < nFaces; ++f) {
            const auto& face = fvMesh.faces[f];
            if (!face.neighbour) {
                const bool noSlipWall = policy.embeddedWallsNoSlip && face.patch == BoundaryPatch2D::EmbeddedBoundary;
                bV[f] = noSlipWall ? 0.0 : v[face.owner];
            }
        }
    }
    if (bP.empty()) {
        bP.resize(nFaces, 0.0);
        for (std::size_t f = 0; f < nFaces; ++f) {
            const auto& face = fvMesh.faces[f];
            if (!face.neighbour) bP[f] = p[face.owner];
        }
    }

    const auto gradU = reconstructGradient(fvMesh, u, bU);
    const auto gradV = reconstructGradient(fvMesh, v, bV);
    const auto gradP = reconstructGradient(fvMesh, p, bP);

    std::vector<Vector2D> gradS;
    double sSpan = 0.0;
    if (!scalar.empty()) {
        std::vector<double> bS = boundaryScalar;
        if (bS.empty()) {
            bS.resize(nFaces, 0.0);
            for (std::size_t f = 0; f < nFaces; ++f) {
                const auto& face = fvMesh.faces[f];
                if (!face.neighbour) bS[f] = scalar[face.owner];
            }
        }
        gradS = reconstructGradient(fvMesh, scalar, bS);
        double sMin = std::numeric_limits<double>::infinity();
        double sMax = -std::numeric_limits<double>::infinity();
        for (double val : scalar) {
            sMin = std::min(sMin, val);
            sMax = std::max(sMax, val);
        }
        sSpan = std::max(0.0, sMax - sMin);
    }

    double uRef = 0.0;
    double pMin = std::numeric_limits<double>::infinity();
    double pMax = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < nCells; ++i) {
        uRef = std::max(uRef, std::hypot(u[i], v[i]));
        pMin = std::min(pMin, p[i]);
        pMax = std::max(pMax, p[i]);
    }
    const double pSpan = std::max(0.0, pMax - pMin);
    const double pRef = std::max(pSpan, uRef * uRef);
    result.referenceVelocity = uRef;
    result.referencePressure = pRef;
    result.referenceScalarSpan = sSpan;

    result.cellIndicators.assign(nCells, 0.0);
    result.cellSelected.assign(nCells, false);
    for (std::size_t i = 0; i < nCells; ++i) {
        const double h = std::sqrt(std::max(fvMesh.cells[i].area, 0.0));
        double ind = 0.0;
        if (uRef > 1.0e-15 && policy.velocityWeight > 0.0) {
            const double velGradNorm = std::sqrt(
                gradU[i].x * gradU[i].x + gradU[i].y * gradU[i].y +
                gradV[i].x * gradV[i].x + gradV[i].y * gradV[i].y);
            ind += policy.velocityWeight * h * velGradNorm / uRef;
        }
        if (pRef > 1.0e-15 && policy.pressureWeight > 0.0) {
            const double pGradNorm = std::hypot(gradP[i].x, gradP[i].y);
            ind += policy.pressureWeight * h * pGradNorm / pRef;
        }
        if (!gradS.empty() && sSpan > 1.0e-15 && policy.scalarWeight > 0.0) {
            const double sGradNorm = std::hypot(gradS[i].x, gradS[i].y);
            ind += policy.scalarWeight * h * sGradNorm / sSpan;
        }
        if (std::isfinite(ind) && ind > 0.0) result.cellIndicators[i] = ind;
    }

    std::vector<double> sortedInd = result.cellIndicators;
    std::sort(sortedInd.begin(), sortedInd.end());
    result.indicatorP50 = percentileSorted(sortedInd, 0.50);
    result.indicatorP90 = percentileSorted(sortedInd, 0.90);
    result.indicatorMax = sortedInd.back();

    // Uniform / quiescent fields emit zero refinement boxes.
    constexpr double kAbsoluteIndicatorFloor = 1.0e-9;
    if (!(result.indicatorMax > kAbsoluteIndicatorFloor)) {
        return result;
    }

    // 4. Deterministic top-fraction cell selection with relative peak floor.
    const double cutoff = std::max(kAbsoluteIndicatorFloor, policy.minRelativePeak * result.indicatorMax);
    std::vector<std::size_t> ranked;
    ranked.reserve(nCells);
    for (std::size_t i = 0; i < nCells; ++i) {
        if (result.cellIndicators[i] >= cutoff) ranked.push_back(i);
    }
    if (ranked.empty()) return result;

    std::sort(ranked.begin(), ranked.end(), [&](std::size_t a, std::size_t b) {
        const double ia = result.cellIndicators[a];
        const double ib = result.cellIndicators[b];
        if (ia != ib) return ia > ib;
        return a < b;
    });
    const std::size_t targetSelect = std::max<std::size_t>(
        1U, static_cast<std::size_t>(std::ceil(policy.refineFraction * static_cast<double>(nCells))));
    if (ranked.size() > targetSelect) ranked.resize(targetSelect);

    std::vector<bool> selected(nCells, false);
    double selectedArea = 0.0;
    for (const auto id : ranked) {
        selected[id] = true;
        selectedArea += fvMesh.cells[id].area;
    }
    result.cellSelected = selected;
    result.selectedCellCount = ranked.size();
    result.selectedAreaFraction = totalArea > 0.0 ? selectedArea / totalArea : 0.0;

    // 4b. Emit compact flow-aligned cell-level AMR targets (one per selected
    // high-gradient cell) in deterministic cell-id order. Each cell upgrades by
    // at most +1 level (capped at capLevel) with a tight isotropic halo plus a
    // subtle downstream convective lead along the local velocity vector so
    // oblique wakes/shear layers stay crisp without dead bounding-box corners.
    result.cellTargets.reserve(ranked.size());
    for (std::size_t i = 0; i < nCells; ++i) {
        if (!selected[i]) continue;
        const std::size_t targetLevel = std::max<std::size_t>(
            1U, std::min(capLevel, cellLevels[i] + 1U));
        const double h = std::sqrt(std::max(fvMesh.cells[i].area, 0.0));
        const double isoPad = 0.5 * policy.padCells * h;
        AABB2D box = cellBounds[i];
        box.min.x -= isoPad;
        box.max.x += isoPad;
        box.min.y -= isoPad;
        box.max.y += isoPad;
        const double speed = std::hypot(u[i], v[i]);
        if (uRef > 1.0e-15 && speed > 1.0e-6 * uRef && policy.padCells > 0.0) {
            const double lead = 0.5 * policy.padCells * h;
            const double dxLead = (u[i] / speed) * lead;
            const double dyLead = (v[i] / speed) * lead;
            if (dxLead > 0.0) box.max.x += dxLead; else box.min.x += dxLead;
            if (dyLead > 0.0) box.max.y += dyLead; else box.min.y += dyLead;
        }
        box = clampBox(box, domainBox);
        if (box.max.x > box.min.x && box.max.y > box.min.y) {
            result.cellTargets.push_back({box, targetLevel});
        }
    }

    // 5. Face-adjacency connected components over selected cells, keeping
    // clusters within one grid level of their seed so fine near-wall cells do
    // not promote an entire coarse far-wake box to wall-level resolution.
    std::vector<bool> visited(nCells, false);
    std::vector<Cluster> clusters;
    for (std::size_t seed = 0; seed < nCells; ++seed) {
        if (!selected[seed] || visited[seed]) continue;
        const std::size_t seedLevel = cellLevels[seed];
        std::vector<std::size_t> comp;
        std::queue<std::size_t> q;
        visited[seed] = true;
        q.push(seed);
        while (!q.empty()) {
            const auto cur = q.front();
            q.pop();
            comp.push_back(cur);
            for (const auto fid : fvMesh.cells[cur].faces) {
                const auto& face = fvMesh.faces[fid];
                if (!face.neighbour) continue;
                const std::size_t other = (face.owner == cur) ? *face.neighbour : face.owner;
                if (!selected[other] || visited[other]) continue;
                const std::size_t diff = cellLevels[other] >= seedLevel
                    ? (cellLevels[other] - seedLevel) : (seedLevel - cellLevels[other]);
                if (diff > 1U) continue;
                visited[other] = true;
                q.push(other);
            }
        }
        std::sort(comp.begin(), comp.end());
        splitSparseCluster(makeCluster(comp, cellBounds, fvMesh), cellBounds, fvMesh,
                           policy.splitFillThreshold, 0U, clusters);
    }

    // 6. Convert clusters into padded, domain-clamped refinement boxes using
    // area-weighted median level so small fine cells cannot inflate a large box.
    std::vector<FlowAdaptiveBox2D> candidates;
    candidates.reserve(clusters.size());
    for (const auto& cluster : clusters) {
        std::vector<std::pair<std::size_t, double>> levelAreas;
        levelAreas.reserve(cluster.cells.size());
        double indSum = 0.0;
        double indMax = 0.0;
        double score = 0.0;
        for (const auto id : cluster.cells) {
            levelAreas.emplace_back(cellLevels[id], fvMesh.cells[id].area);
            const double ind = result.cellIndicators[id];
            indSum += ind;
            indMax = std::max(indMax, ind);
            score += ind * fvMesh.cells[id].area;
        }
        std::sort(levelAreas.begin(), levelAreas.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        double accumArea = 0.0;
        const double halfArea = 0.5 * cluster.areaSum;
        std::size_t areaMedianLevel = levelAreas.front().first;
        for (const auto& [lvl, a] : levelAreas) {
            accumArea += a;
            areaMedianLevel = lvl;
            if (accumArea >= halfArea) break;
        }
        const std::size_t targetLevel = std::min(capLevel, areaMedianLevel + 1U);

        const double localH = std::sqrt(std::max(cluster.areaSum / static_cast<double>(cluster.cells.size()), 0.0));
        const double pad = policy.padCells * localH;
        const AABB2D padded = clampBox(
            {{cluster.tightBounds.min.x - pad, cluster.tightBounds.min.y - pad},
             {cluster.tightBounds.max.x + pad, cluster.tightBounds.max.y + pad}},
            domainBox);
        if (!(padded.max.x > padded.min.x) || !(padded.max.y > padded.min.y)) continue;

        FlowAdaptiveBox2D box;
        box.region = {padded, targetLevel};
        box.cellCount = cluster.cells.size();
        box.meanIndicator = indSum / static_cast<double>(cluster.cells.size());
        box.maxIndicator = indMax;
        box.score = score;
        candidates.push_back(box);
    }

    // 7. Sort by integrated score and merge heavily overlapping same-level boxes.
    const auto boxOrder = [](const FlowAdaptiveBox2D& a, const FlowAdaptiveBox2D& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.region.targetLevel != b.region.targetLevel) return a.region.targetLevel > b.region.targetLevel;
        if (a.region.bounds.min.x != b.region.bounds.min.x) return a.region.bounds.min.x < b.region.bounds.min.x;
        if (a.region.bounds.min.y != b.region.bounds.min.y) return a.region.bounds.min.y < b.region.bounds.min.y;
        if (a.region.bounds.max.x != b.region.bounds.max.x) return a.region.bounds.max.x < b.region.bounds.max.x;
        return a.region.bounds.max.y < b.region.bounds.max.y;
    };
    std::sort(candidates.begin(), candidates.end(), boxOrder);

    std::vector<FlowAdaptiveBox2D> merged;
    for (const auto& cand : candidates) {
        bool absorbed = false;
        for (auto& existing : merged) {
            if (existing.region.targetLevel == cand.region.targetLevel &&
                overlapRatio(existing.region.bounds, cand.region.bounds) >= 0.5) {
                const double totalCells = static_cast<double>(existing.cellCount + cand.cellCount);
                existing.meanIndicator = (existing.meanIndicator * static_cast<double>(existing.cellCount) +
                                          cand.meanIndicator * static_cast<double>(cand.cellCount)) / totalCells;
                existing.maxIndicator = std::max(existing.maxIndicator, cand.maxIndicator);
                existing.score += cand.score;
                existing.cellCount += cand.cellCount;
                existing.region.bounds = unionBox(existing.region.bounds, cand.region.bounds);
                absorbed = true;
                break;
            }
        }
        if (!absorbed) merged.push_back(cand);
    }

    std::sort(merged.begin(), merged.end(), boxOrder);
    if (merged.size() > policy.maxBoxes) merged.resize(policy.maxBoxes);

    // 8. Finalize body-frame normalized bounds and levelsBelowWall.
    for (auto& box : merged) {
        box.levelsBelowWall = refWallLevel >= box.region.targetLevel
            ? (refWallLevel - box.region.targetLevel) : 0U;
        box.bodyRelativeBounds = {
            {(box.region.bounds.min.x - result.bodyCentre.x) / result.bodySpan,
             (box.region.bounds.min.y - result.bodyCentre.y) / result.bodySpan},
            {(box.region.bounds.max.x - result.bodyCentre.x) / result.bodySpan,
             (box.region.bounds.max.y - result.bodyCentre.y) / result.bodySpan}
        };
    }
    result.boxes = std::move(merged);
    return result;
}

std::string flowAdaptiveRefinementToJson2D(const FlowAdaptiveRefinement2D& suggestion, int indentSpaces) {
    const std::string i1(static_cast<std::size_t>(std::max(0, indentSpaces)), ' ');
    const std::string i2 = i1 + i1;
    std::ostringstream out;
    out << std::setprecision(17);
    out << "{\n";
    out << i1 << "\"format\": \"cartmesh2d-adaptive-refinement-v1\",\n";
    out << i1 << "\"valid\": " << (suggestion.valid() ? "true" : "false") << ",\n";
    out << i1 << "\"cell_count\": " << suggestion.cellCount << ",\n";
    out << i1 << "\"selected_cell_count\": " << suggestion.selectedCellCount << ",\n";
    out << i1 << "\"selected_area_fraction\": " << suggestion.selectedAreaFraction << ",\n";
    out << i1 << "\"inferred_wall_level\": " << suggestion.inferredWallLevel << ",\n";
    out << i1 << "\"max_target_level\": " << suggestion.maxTargetLevel << ",\n";
    out << i1 << "\"reference_velocity\": " << suggestion.referenceVelocity << ",\n";
    out << i1 << "\"reference_pressure\": " << suggestion.referencePressure << ",\n";
    out << i1 << "\"reference_scalar_span\": " << suggestion.referenceScalarSpan << ",\n";
    out << i1 << "\"indicator_p50\": " << suggestion.indicatorP50 << ",\n";
    out << i1 << "\"indicator_p90\": " << suggestion.indicatorP90 << ",\n";
    out << i1 << "\"indicator_max\": " << suggestion.indicatorMax << ",\n";
    out << i1 << "\"domain\": {\"xmin\": " << suggestion.domainBounds.min.x
        << ", \"ymin\": " << suggestion.domainBounds.min.y
        << ", \"xmax\": " << suggestion.domainBounds.max.x
        << ", \"ymax\": " << suggestion.domainBounds.max.y << "},\n";
    out << i1 << "\"body_frame\": {\"cx\": " << suggestion.bodyCentre.x
        << ", \"cy\": " << suggestion.bodyCentre.y
        << ", \"span\": " << suggestion.bodySpan
        << ", \"from_embedded_wall\": " << (suggestion.bodyFromEmbeddedWall ? "true" : "false") << "},\n";
    out << i1 << "\"boxes\": [\n";
    for (std::size_t k = 0; k < suggestion.boxes.size(); ++k) {
        const auto& b = suggestion.boxes[k];
        out << i2 << "{\"xmin\": " << b.region.bounds.min.x
            << ", \"ymin\": " << b.region.bounds.min.y
            << ", \"xmax\": " << b.region.bounds.max.x
            << ", \"ymax\": " << b.region.bounds.max.y
            << ", \"target_level\": " << b.region.targetLevel
            << ", \"levels_below_wall\": " << b.levelsBelowWall
            << ", \"body_xmin\": " << b.bodyRelativeBounds.min.x
            << ", \"body_ymin\": " << b.bodyRelativeBounds.min.y
            << ", \"body_xmax\": " << b.bodyRelativeBounds.max.x
            << ", \"body_ymax\": " << b.bodyRelativeBounds.max.y
            << ", \"cell_count\": " << b.cellCount
            << ", \"mean_indicator\": " << b.meanIndicator
            << ", \"max_indicator\": " << b.maxIndicator
            << ", \"score\": " << b.score << "}"
            << (k + 1U < suggestion.boxes.size() ? "," : "") << "\n";
    }
    out << i1 << "],\n";
    out << i1 << "\"cell_targets\": [\n";
    for (std::size_t k = 0; k < suggestion.cellTargets.size(); ++k) {
        const auto& t = suggestion.cellTargets[k];
        out << i2 << "{\"xmin\": " << t.bounds.min.x
            << ", \"ymin\": " << t.bounds.min.y
            << ", \"xmax\": " << t.bounds.max.x
            << ", \"ymax\": " << t.bounds.max.y
            << ", \"target_level\": " << t.targetLevel << "}"
            << (k + 1U < suggestion.cellTargets.size() ? "," : "") << "\n";
    }
    out << i1 << "],\n";
    out << i1 << "\"issues\": [\n";
    for (std::size_t k = 0; k < suggestion.issues.size(); ++k) {
        out << i2 << std::quoted(suggestion.issues[k])
            << (k + 1U < suggestion.issues.size() ? "," : "") << "\n";
    }
    out << i1 << "]\n";
    out << "}\n";
    return out.str();
}

std::string flowAdaptiveFieldsToJson2D(const FlowAdaptiveRefinement2D& suggestion) {
    std::ostringstream out;
    out << std::setprecision(10);
    out << "{\"format\":\"cartmesh2d-adaptive-fields-v1\",\"indicatorMax\":" << suggestion.indicatorMax
        << ",\"cells\":[";
    for (std::size_t i = 0; i < suggestion.cellIndicators.size(); ++i) {
        if (i > 0U) out << ',';
        const int sel = (i < suggestion.cellSelected.size() && suggestion.cellSelected[i]) ? 1 : 0;
        out << "{\"id\":" << i << ",\"indicator\":" << suggestion.cellIndicators[i]
            << ",\"selected\":" << sel << '}';
    }
    out << "]}\n";
    return out.str();
}

bool writeFlowAdaptiveVtk2D(
    const TopologyMesh2D& topology,
    const FlowAdaptiveRefinement2D& suggestion,
    const std::filesystem::path& vtkPath,
    std::string* error) {
    if (!suggestion.valid() || suggestion.cellIndicators.size() != topology.cells.size()) {
        if (error) *error = "adaptive VTK export requires valid suggestion matching topology cell count";
        return false;
    }
    const auto parent = vtkPath.parent_path();
    if (!parent.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(parent, ec);
    }
    if (!writeLegacyVtk2D(topology, vtkPath, error)) return false;
    std::ofstream vtk(vtkPath, std::ios::app);
    if (!vtk) {
        if (error) *error = "failed to reopen adaptive VTK file for scalar fields";
        return false;
    }
    vtk << std::setprecision(17)
        << "SCALARS adaptive_indicator double 1\nLOOKUP_TABLE default\n";
    for (double val : suggestion.cellIndicators) vtk << val << '\n';
    vtk << "SCALARS adaptive_selected int 1\nLOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < suggestion.cellIndicators.size(); ++i) {
        const int sel = (i < suggestion.cellSelected.size() && suggestion.cellSelected[i]) ? 1 : 0;
        vtk << sel << '\n';
    }
    if (!vtk.good()) {
        if (error) *error = "failed while writing adaptive VTK scalar fields";
        return false;
    }
    return true;
}

} // namespace cartmesh2d::fv
