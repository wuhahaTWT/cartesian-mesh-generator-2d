#pragma once

#include "cartmesh2d/io/BoundaryMetadata2D.hpp"
#include "cartmesh2d/topology/Topology2D.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace cartmesh2d {

struct FluentZoneSummary2D {
    std::size_t zoneId = 0;
    std::string name;
    std::string bcType;
    BoundaryConditionRole2D role = BoundaryConditionRole2D::Wall;
    std::size_t faceCount = 0;
    std::size_t startFace = 0;
};

struct FluentMshWriteReport2D {
    std::size_t nodeCount = 0;
    std::size_t cellCount = 0;
    std::size_t faceCount = 0;
    std::size_t internalFaceCount = 0;
    std::size_t triangleCount = 0;
    std::size_t quadCount = 0;
    std::size_t polygonCount = 0;
    std::vector<FluentZoneSummary2D> zones;

    [[nodiscard]] bool valid() const noexcept {
        return nodeCount > 0 && cellCount > 0 && faceCount > 0;
    }
};

// Write a true native 2D ANSYS Fluent face-based ASCII mesh (.msh).
// Preserves exact 17-digit vertex coordinates (no wall movement), exact
// polygon connectivity (tri=1, quad=3, poly=7), and Fluent's 2D left-hand c0
// face orientation convention (c0 on the left of n0 -> n1; c1 = 0 on boundary).
[[nodiscard]] FluentMshWriteReport2D writeFluentMsh2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const std::filesystem::path& mshFile,
    std::string* error = nullptr,
    const TolerancePolicy& tol = {},
    const std::vector<EmbeddedBoundaryPatch2D>& embeddedPatches = {});

} // namespace cartmesh2d
