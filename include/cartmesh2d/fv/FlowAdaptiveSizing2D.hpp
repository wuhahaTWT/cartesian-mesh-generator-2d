#pragma once

#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/quadtree/Quadtree2D.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d::fv {

// Flow-preview-driven adaptive mesh refinement (AMR) sizing controls.
// Computes a dimensionless, scale- and rotation-invariant feature indicator on
// the actual solved 2D polygon mesh, clusters high-gradient cells using face
// adjacency, splits sparse clusters into tight staircase boxes, and emits
// deterministic BoxRefinementRegion2D targets consumable by cartmesh2d_cli and
// the desktop app.
struct FlowAdaptivePolicy2D {
    // Fraction of fluid cells eligible for clustering, in (0, 0.5].
    double refineFraction = 0.15;
    // Minimum cell indicator relative to the peak indicator, in [0, 1).
    // Filters uniform freestream roundoff noise.
    double minRelativePeak = 0.05;
    // Maximum number of emitted refinement boxes, in [1, 32].
    std::size_t maxBoxes = 6;
    // Padding around each cluster in local cell widths, in [0, 8].
    double padCells = 1.0;
    // Area fill threshold below which clusters of >= 8 cells are recursively
    // bisected along their longest axis so curved wakes stay tight.
    double splitFillThreshold = 0.40;
    // Relative weights on velocity-gradient Frobenius norm, pressure gradient,
    // and optional passive scalar / temperature gradient.
    double velocityWeight = 1.0;
    double pressureWeight = 1.0;
    double scalarWeight = 1.0;
    // Optional ceiling on suggested targetLevel. Unset uses the inferred finest
    // level present in the solved mesh so re-meshing stays within the verified
    // wall level budget unless the caller explicitly raises it.
    std::optional<std::size_t> maxTargetLevel;
    // Optional reference wall level used for levelsBelowWall; defaults to the
    // effective maxTargetLevel.
    std::optional<std::size_t> wallLevel;
    // Optional body bounding box for body-frame normalized box coordinates.
    // Unset infers the box from EmbeddedBoundary vertices (or domain bounds if
    // the mesh has no embedded solid wall).
    std::optional<AABB2D> bodyBounds;
    // Treat EmbeddedBoundary faces as no-slip (u=0, v=0) when boundaryU/V are
    // not explicitly supplied by the caller.
    bool embeddedWallsNoSlip = true;
};

struct FlowAdaptiveBox2D {
    BoxRefinementRegion2D region;
    std::size_t levelsBelowWall = 0;
    // Body-frame normalized coordinates (about body centre, divided by body span),
    // matching the desktop app's refineBoxes convention.
    AABB2D bodyRelativeBounds;
    std::size_t cellCount = 0;
    double meanIndicator = 0.0;
    double maxIndicator = 0.0;
    double score = 0.0;
};

struct FlowAdaptiveRefinement2D {
    std::size_t cellCount = 0;
    std::size_t selectedCellCount = 0;
    double selectedAreaFraction = 0.0;
    std::size_t inferredWallLevel = 0;
    std::size_t maxTargetLevel = 0;
    double referenceVelocity = 0.0;
    double referencePressure = 0.0;
    double referenceScalarSpan = 0.0;
    double indicatorP50 = 0.0;
    double indicatorP90 = 0.0;
    double indicatorMax = 0.0;
    AABB2D domainBounds;
    Point2D bodyCentre;
    double bodySpan = 0.0;
    bool bodyFromEmbeddedWall = false;
    std::vector<double> cellIndicators;
    std::vector<bool> cellSelected;
    std::vector<FlowAdaptiveBox2D> boxes;
    // Per-cell compact flow-aligned refinement footprints (one per selected
    // high-gradient cell), enabling contour-hugging quadtree AMR without
    // axis-aligned bounding-box bloat on oblique wakes and curved shear layers.
    std::vector<BoxRefinementRegion2D> cellTargets;
    std::vector<std::string> issues;

    [[nodiscard]] bool valid() const noexcept { return issues.empty(); }
};

[[nodiscard]] FlowAdaptiveRefinement2D suggestFlowAdaptiveRefinement2D(
    const TopologyMesh2D& topology,
    const FvMesh2D& fvMesh,
    const std::vector<double>& u,
    const std::vector<double>& v,
    const std::vector<double>& p,
    const FlowAdaptivePolicy2D& policy = {},
    const std::vector<double>& boundaryU = {},
    const std::vector<double>& boundaryV = {},
    const std::vector<double>& boundaryP = {},
    const std::vector<double>& scalar = {},
    const std::vector<double>& boundaryScalar = {});

[[nodiscard]] std::string flowAdaptiveRefinementToJson2D(
    const FlowAdaptiveRefinement2D& suggestion,
    int indentSpaces = 2);

[[nodiscard]] std::string flowAdaptiveFieldsToJson2D(
    const FlowAdaptiveRefinement2D& suggestion);

[[nodiscard]] bool writeFlowAdaptiveVtk2D(
    const TopologyMesh2D& topology,
    const FlowAdaptiveRefinement2D& suggestion,
    const std::filesystem::path& vtkPath,
    std::string* error = nullptr);

} // namespace cartmesh2d::fv
