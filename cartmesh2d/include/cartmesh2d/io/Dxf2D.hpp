#pragma once

#include "cartmesh2d/geometry/Geometry2D.hpp"
#include "cartmesh2d/io/BoundaryMetadata2D.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

struct DxfImportOptions2D {
    // Maximum sagitta in output metres between an analytic curve and chord.
    double maximumChordError = 1.0e-3;
    // Absolute endpoint welding distance used only while assembling entities
    // into loops. Geometry predicates still use TolerancePolicy.
    double endpointWeldTolerance = 1.0e-10;
    // Missing/unitless DXF files fail closed unless an explicit code is supplied.
    std::optional<long long> sourceUnitsOverrideCode;
    bool rejectUnsupportedEntities = true;
};

enum class DxfIssueCode2D {
    CannotOpen,
    BinaryDxfUnsupported,
    MalformedGroupPair,
    MissingEntitiesSection,
    UnsupportedEntity,
    MissingRequiredGroup,
    InvalidNumericValue,
    NonPlanarEntity,
    UnsupportedExtrusion,
    NonZeroWidthOrThickness,
    UnknownOrUnitlessUnits,
    InvalidSplineDefinition,
    BoundaryMetadataConflict,
    OpenOrBranchedBoundary,
    DegenerateEntity,
    InvalidBoundaryRegion,
    OutputFailure
};

struct DxfIssue2D {
    DxfIssueCode2D code = DxfIssueCode2D::MalformedGroupPair;
    std::size_t line = 0;
    std::string entity;
    std::string layer;
    std::string message;
};

struct DxfImportReport2D {
    std::size_t lineEntityCount = 0;
    std::size_t lightweightPolylineCount = 0;
    std::size_t arcEntityCount = 0;
    std::size_t circleEntityCount = 0;
    std::size_t ellipseEntityCount = 0;
    std::size_t splineEntityCount = 0;
    std::size_t sourceEntityCount = 0;
    std::size_t outputLoopCount = 0;
    std::size_t outputVertexCount = 0;
    std::size_t sampledArcSegmentCount = 0;
    std::size_t sampledEllipseSegmentCount = 0;
    std::size_t sampledSplineSegmentCount = 0;
    double maximumChordError = 0.0;
    double endpointWeldTolerance = 0.0;
    double sourcePlaneZ = 0.0;
    bool sourcePlaneZDefined = false;
    long long insertionUnitsCode = 0;
    bool insertionUnitsCodeDefined = false;
    long long effectiveUnitsCode = 0;
    std::string effectiveUnitsName;
    double coordinateScaleToMetres = 0.0;
    bool unitsOverrideApplied = false;
    std::vector<std::string> layers;
};

struct DxfImportResult2D {
    std::optional<BoundaryRegion2D> boundary;
    std::vector<EmbeddedBoundaryPatch2D> embeddedPatches;
    DxfImportReport2D report;
    std::vector<DxfIssue2D> issues;

    [[nodiscard]] bool valid() const noexcept {
        return boundary.has_value() && issues.empty();
    }
};

[[nodiscard]] DxfImportResult2D readAsciiDxfBoundary2D(
    const std::filesystem::path& path,
    const DxfImportOptions2D& options = {},
    const TolerancePolicy& tol = {});

[[nodiscard]] bool writeBoundaryXy2D(
    const BoundaryRegion2D& boundary,
    const std::filesystem::path& path,
    std::string* error = nullptr,
    const std::vector<EmbeddedBoundaryPatch2D>& embeddedPatches = {});

[[nodiscard]] std::string dxfImportReportToJson(
    const DxfImportResult2D& result, int indentSpaces = 2);

} // namespace cartmesh2d
