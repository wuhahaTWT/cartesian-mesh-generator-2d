#pragma once

#include "cartmesh2d/topology/Topology2D.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace cartmesh2d {

struct OpenFoamPatchSummary2D {
    std::string name;
    std::string type;
    std::size_t faceCount = 0;
    std::size_t startFace = 0;
};

struct OpenFoamWriteReport2D {
    std::size_t pointCount = 0;
    std::size_t faceCount = 0;
    std::size_t internalFaceCount = 0;
    std::size_t cellCount = 0;
    double thickness = 0.0;
    std::vector<OpenFoamPatchSummary2D> patches;

    [[nodiscard]] bool valid() const noexcept {
        return pointCount>0 && faceCount>0 && cellCount>0 && thickness>0.0;
    }
};

[[nodiscard]] OpenFoamWriteReport2D writeExtrudedOpenFoam2D(
    const TopologyMesh2D& topology,
    const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const std::filesystem::path& caseDirectory,
    double thickness,
    std::string* error = nullptr,
    const TolerancePolicy& tol = {});

} // namespace cartmesh2d
