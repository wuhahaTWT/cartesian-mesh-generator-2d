#pragma once

#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"
#include "cartmesh/grid/LinearOctree.hpp"
#include "cartmesh/grid/UniformCartesianGrid.hpp"

#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace cartmesh {

struct OpenFoamCellSource {
    std::uint64_t background_cell_id{};
    std::uint64_t background_stable_id{};
    std::size_t component_id{};
    std::size_t local_piece_id{};
    double source_volume_fraction{1.0};
    bool full_cartesian{};
};

struct OpenFoamFace {
    std::vector<std::size_t> point_ids;
    std::size_t owner{};
    std::size_t neighbour{std::numeric_limits<std::size_t>::max()};
    std::uint64_t boundary_id{};
    bool farfield{};

    [[nodiscard]] bool internal() const noexcept {
        return neighbour != std::numeric_limits<std::size_t>::max();
    }
};

// reference writer 的最终内存拓扑。faces 采用 OpenFOAM 顺序：内部面在前，
// boundary faces 在后；点已经经过与实际 polyMesh 相同的确定性焊接。
struct OpenFoamMesh {
    std::string background_stable_id_kind;
    std::vector<Vec3> points;
    std::vector<OpenFoamFace> faces;
    std::vector<OpenFoamCellSource> cells;
    std::size_t internal_face_count{};
    double length_tolerance{};
};

[[nodiscard]] OpenFoamMesh build_openfoam_mesh(
    const UniformCartesianGrid& grid, const ConvexCutCellMesh& mesh,
    double length_tolerance = 0.0);

[[nodiscard]] OpenFoamMesh build_openfoam_mesh(
    const LinearOctree& tree, const ConvexCutCellMesh& mesh,
    double length_tolerance = 0.0);

void write_openfoam_poly_mesh(
    const std::filesystem::path& case_directory,
    const OpenFoamMesh& mesh,
    const std::vector<std::pair<std::uint64_t, std::string>>& boundary_names = {});

// 写出完整流体域 constant/polyMesh。普通流体背景单元与
// Cut-cell 的局部 fluid component 都是体单元；共享 Cartesian
// 面在两侧凸分解上取公共细分，不产生悬挂面。
void write_openfoam_poly_mesh(
    const std::filesystem::path& case_directory,
    const UniformCartesianGrid& grid,
    const ConvexCutCellMesh& mesh,
    const std::vector<std::pair<std::uint64_t, std::string>>& boundary_names = {},
    double length_tolerance = 0.0);

// 与均匀路径共用同一控制体组装和序列化实现。八叉树叶按稳定 Morton 顺序映射到
// background_cell_id；粗细交界以一个 coarse face 对一组 fine faces 整体守恒。
void write_openfoam_poly_mesh(
    const std::filesystem::path& case_directory,
    const LinearOctree& tree,
    const ConvexCutCellMesh& mesh,
    const std::vector<std::pair<std::uint64_t, std::string>>& boundary_names = {},
    double length_tolerance = 0.0);

} // 命名空间 cartmesh
