#pragma once

#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <array>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

enum class FluidSide2D { Left, Right };
enum class WallFluidRegion2D { Exterior, Interior };
enum class WallChainOrientation2D { Open, Clockwise, CounterClockwise };
enum class WallVertexKind2D { Smooth, MildConvex, Concave, Sharp, Degenerate };

struct WallChain2D {
    std::size_t id = 0;
    std::string patchIdentity;
    std::vector<Point2D> vertices;
    std::vector<Segment2D> segments;
    bool closed = false;
    WallChainOrientation2D orientation = WallChainOrientation2D::Open;
    FluidSide2D fluidSide = FluidSide2D::Right;

    [[nodiscard]] std::size_t segmentCount() const noexcept {
        return segments.size();
    }
};

enum class WallChainFailureReason2D {
    None,
    InvalidBoundary,
    TooFewVertices,
    NonFiniteCoordinate,
    DuplicateVertex,
    DegenerateSegment,
    EmptyPatchIdentity
};

struct WallChainBuildResult2D {
    std::optional<WallChain2D> chain;
    WallChainFailureReason2D failureReason = WallChainFailureReason2D::None;
    std::size_t vertexId = 0;
    std::string message;

    [[nodiscard]] bool success() const noexcept { return chain.has_value(); }
};

// Closed exterior walls are canonicalized counter-clockwise so the solid is
// on the left and the fluid is on the right. Interior walls are canonicalized
// clockwise, again keeping the requested fluid region on the right.
[[nodiscard]] WallChainBuildResult2D makeClosedWallChain2D(
    const BoundaryLoop& boundary, std::size_t chainId,
    std::string patchIdentity,
    WallFluidRegion2D fluidRegion = WallFluidRegion2D::Exterior,
    const TolerancePolicy& tol = {});

// Open chains have no polygonal inside/outside. The caller must state the fluid
// side explicitly; input order is preserved after validation.
[[nodiscard]] WallChainBuildResult2D makeOpenWallChain2D(
    std::vector<Point2D> orderedVertices, std::size_t chainId,
    std::string patchIdentity, FluidSide2D fluidSide,
    const TolerancePolicy& tol = {});

enum class LayerThicknessMode2D { FirstLayerThickness, TotalThickness };

struct LayerParameters2D {
    std::size_t nLayers = 1;
    LayerThicknessMode2D thicknessMode =
        LayerThicknessMode2D::FirstLayerThickness;
    double thickness = 0.0;
    double growthRatio = 1.0;
};

struct ResolvedLayerParameters2D {
    std::size_t nLayers = 0;
    double firstLayerThickness = 0.0;
    double totalThickness = 0.0;
    double growthRatio = 0.0;
    std::vector<double> cumulativeNormalDistances;
};

struct LayerParameterResult2D {
    std::optional<ResolvedLayerParameters2D> parameters;
    std::string message;

    [[nodiscard]] bool success() const noexcept { return parameters.has_value(); }
};

[[nodiscard]] LayerParameterResult2D resolveLayerParameters2D(
    const LayerParameters2D& parameters);

// All H4-1 geometric thresholds live here. Length/area tolerances are derived
// from TolerancePolicy and the chain bounds. The dimensionless limits have a
// geometric meaning: maximum accepted solid-side turn, maximum tangential
// consumption of an adjacent segment, and retained collision clearance.
struct BoundaryLayerPolicy2D {
    TolerancePolicy tolerance{};
    double smoothTurnRadians = 0.08726646259971647;      // 5 degrees
    double maxConvexTurnRadians = 2.3561944901923448;    // 135 degrees
    double cornerLengthFraction = 0.45;
    double collisionClearanceFraction = 0.45;
};

enum class BoundaryLayerFailureReason2D {
    None,
    InvalidParameters,
    InvalidChain,
    DegenerateSegment,
    ConcaveCorner,
    SharpCorner,
    ConflictingHalfPlanes,
    ThicknessExceedsSafeLimit,
    EnvelopeSelfIntersection,
    EnvelopeWallIntersection,
    ChainCollision,
    NegativeAreaCell,
    SelfIntersectingCell,
    OverlappingCells,
    DuplicateGeometricVertex,
    TopologyConflict,
    IoFailure
};

struct BoundaryLayerFailure2D {
    BoundaryLayerFailureReason2D reason = BoundaryLayerFailureReason2D::None;
    std::string message;
    std::size_t chainId = 0;
    std::optional<std::size_t> vertexId;
    std::optional<std::size_t> edgeId;
    std::optional<std::size_t> cellId;
    double requestedThickness = 0.0;
    std::optional<double> safeThickness;
    std::size_t nLayers = 0;
    double growthRatio = 0.0;
};

struct BoundaryLayerVertex2D {
    std::size_t id = 0;
    std::size_t ring = 0;
    std::size_t chainVertex = 0;
    Point2D point;
};

struct BoundaryLayerCell2D {
    std::size_t id = 0;
    std::size_t layer = 0;
    std::size_t wallSegment = 0;
    std::array<std::size_t, 4> vertices{};
    double area = 0.0;
};

struct BoundaryLayerMetrics2D {
    std::size_t vertexCount = 0;
    std::size_t cellCount = 0;
    double minCellArea = 0.0;
    double maxCellArea = 0.0;
    double minLayerThickness = 0.0;
    double maxLayerThickness = 0.0;
    double requestedTotalThickness = 0.0;
    double usedTotalThickness = 0.0;
    double safeThicknessLimit = 0.0;
};

struct BoundaryLayerStrip2D {
    WallChain2D wallChain;
    ResolvedLayerParameters2D parameters;
    std::vector<WallVertexKind2D> wallVertexKinds;
    std::vector<Vector2D> marchingDirections;
    std::vector<std::vector<std::size_t>> ringVertexIds;
    std::vector<BoundaryLayerVertex2D> vertices;
    std::vector<BoundaryLayerCell2D> cells;
    BoundaryLayerMetrics2D metrics;

    [[nodiscard]] std::vector<Point2D> outerEnvelope() const;
};

enum class BoundaryLayerStatus2D { Success, Failed };

struct BoundaryLayerBuildResult2D {
    BoundaryLayerStatus2D status = BoundaryLayerStatus2D::Failed;
    std::vector<BoundaryLayerStrip2D> strips;
    BoundaryLayerFailure2D failure;

    [[nodiscard]] bool success() const noexcept {
        return status == BoundaryLayerStatus2D::Success;
    }
};

// Builds all strips as one transaction. Cross-chain clearance and envelope
// collisions are checked before any successful result is exposed.
[[nodiscard]] BoundaryLayerBuildResult2D buildBoundaryLayerStrips2D(
    const std::vector<WallChain2D>& wallChains,
    const LayerParameters2D& parameters,
    const BoundaryLayerPolicy2D& policy = {});

[[nodiscard]] BoundaryLayerBuildResult2D buildBoundaryLayerStrip2D(
    const WallChain2D& wallChain,
    const LayerParameters2D& parameters,
    const BoundaryLayerPolicy2D& policy = {});

[[nodiscard]] const char* boundaryLayerFailureReasonName(
    BoundaryLayerFailureReason2D reason) noexcept;
[[nodiscard]] const char* wallVertexKindName(WallVertexKind2D kind) noexcept;

// H4-1 debug evidence only. These writers do not produce hybrid solver meshes
// and are intentionally separate from the production OpenFOAM writer.
[[nodiscard]] bool writeBoundaryLayerLegacyVtk2D(
    const BoundaryLayerBuildResult2D& result,
    const std::filesystem::path& path, std::string* error = nullptr);
[[nodiscard]] bool writeBoundaryLayerReportJson2D(
    const BoundaryLayerBuildResult2D& result,
    const std::filesystem::path& path, std::string* error = nullptr);

} // namespace cartmesh2d
