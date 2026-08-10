#pragma once

#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"

#include <filesystem>

namespace cartmesh {

// 写出所有 Cut-cell 的嵌入边界多边形。每个 VTK polygon 都带背景单元 ID、
// boundary ID、面积和流体侧外法向。
void write_embedded_boundary_vtp(const std::filesystem::path& path,
                                 const ConvexCutCellMesh& mesh);

// 写出体积 Cut-cell 的显式凸多面体分解，单元类型为 VTK_POLYHEDRON(42)。
void write_fluid_polyhedra_vtu(const std::filesystem::path& path,
                               const ConvexCutCellMesh& mesh);

// 将每个显式凸片以其内部质心和边界面三角扇分解为
// VTK_TETRA(10)，用于不依赖 VTK_POLYHEDRON 共面分区解释的外部检查。
void write_fluid_tetrahedra_vtu(const std::filesystem::path& path,
                                const ConvexCutCellMesh& mesh);

} // 命名空间 cartmesh
