#pragma once

#include "cartmesh/grid/UniformCartesianGrid.hpp"
#include "cartmesh/scalable/CompactUniformCutCellMesh.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace cartmesh {

struct ScalableOpenFoamWriteStats {
    std::uint64_t point_count{};
    std::uint64_t face_count{};
    std::uint64_t internal_face_count{};
    std::uint64_t boundary_face_count{};
    std::uint64_t solver_cell_count{};
    std::uint64_t face_vertex_reference_count{};
    std::uint64_t background_interface_coverage_relaxation_count{};
    double maximum_background_interface_coverage_error{};
    std::uint64_t explicit_patch_coverage_relaxation_count{};
    double maximum_explicit_patch_coverage_error{};
    std::uint64_t topology_collapsed_face_count{};
    double topology_collapsed_face_area{};
    std::uint64_t topology_sealed_loop_count{};
    std::uint64_t topology_sealed_edge_count{};
    double topology_sealed_loop_area{};
    double maximum_topology_sealed_loop_area{};
    std::uint64_t written_bytes{};
    double preparation_seconds{};
    double writing_seconds{};
    double total_seconds{};
    std::uint64_t topology_hash_fnv1a64{};
};

// 将紧凑全域网格写成 OpenFOAM 2606 可直接读取的二进制
// constant/polyMesh。规则 Cartesian 区域按遍历顺序流式枚举；只有
// Cut-cell 邻域的公共细分面在内存中显式保存。
[[nodiscard]] ScalableOpenFoamWriteStats
write_scalable_openfoam_poly_mesh(
    const std::filesystem::path& case_directory,
    const UniformCartesianGrid& grid,
    const CompactUniformCutCellMesh& mesh,
    const std::vector<std::pair<std::uint64_t, std::string>>& boundary_names = {},
    double length_tolerance = 0.0);

} // 命名空间 cartmesh
