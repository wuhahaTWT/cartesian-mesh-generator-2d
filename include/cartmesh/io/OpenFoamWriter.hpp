#pragma once

#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"
#include "cartmesh/grid/LinearOctree.hpp"
#include "cartmesh/grid/UniformCartesianGrid.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace cartmesh {

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
