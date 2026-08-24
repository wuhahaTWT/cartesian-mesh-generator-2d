#pragma once

#include "cartmesh2d/geometry/Geometry2D.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d {

struct DxfImportOptions2D {
    // Maximum sagitta between an analytic DXF arc and its emitted chord.
    double maximumChordError = 1.0e-3;
    // Absolute endpoint welding distance used only while assembling entities
    // into loops. Geometry predicates still use TolerancePolicy.
    double endpointWeldTolerance = 1.0e-10;
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
    std::size_t sourceEntityCount = 0;
    std::size_t outputLoopCount = 0;
    std::size_t outputVertexCount = 0;
    std::size_t sampledArcSegmentCount = 0;
    double maximumChordError = 0.0;
    double endpointWeldTolerance = 0.0;
    double sourcePlaneZ = 0.0;
    bool sourcePlaneZDefined = false;
    long long insertionUnitsCode = 0;
    bool insertionUnitsCodeDefined = false;
    std::vector<std::string> layers;
};

struct DxfImportResult2D {
    std::optional<BoundaryRegion2D> boundary;
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
    std::string* error = nullptr);

[[nodiscard]] std::string dxfImportReportToJson(
    const DxfImportResult2D& result, int indentSpaces = 2);

} // namespace cartmesh2d
