#pragma once

#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"
#include "cartmesh/grid/LinearOctree.hpp"

#include <cstdint>
#include <string>

namespace cartmesh {

[[nodiscard]] std::uint64_t cut_cell_mesh_fingerprint_fnv1a64(
    const LinearOctree& tree, const ConvexCutCellMesh& mesh) noexcept;

[[nodiscard]] std::string fingerprint_hex(std::uint64_t fingerprint);

} // namespace cartmesh
