#pragma once

#include "cartmesh2d/topology/Topology2D.hpp"

#include <filesystem>
#include <string>

namespace cartmesh2d {

struct MeshReadback2D {
    TopologyMesh2D topology;
    std::string error;

    [[nodiscard]] bool valid() const noexcept {
        return error.empty() && topology.valid();
    }
};

[[nodiscard]] bool writeLegacyVtk2D(
    const TopologyMesh2D& topology,
    const std::filesystem::path& path,
    std::string* error = nullptr);

[[nodiscard]] bool writeCm2dTopology(
    const TopologyMesh2D& topology,
    const std::filesystem::path& path,
    std::string* error = nullptr);

[[nodiscard]] MeshReadback2D readCm2dTopology(
    const std::filesystem::path& path);

} // namespace cartmesh2d
