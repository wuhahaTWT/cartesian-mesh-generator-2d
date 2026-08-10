#pragma once

#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"
#include "cartmesh/grid/UniformCartesianGrid.hpp"
#include "cartmesh/grid/LinearOctree.hpp"

#include <filesystem>

namespace cartmesh {

void write_cut_cell_geometry_json(const std::filesystem::path& path,
                                  const UniformCartesianGrid& grid,
                                  const ConvexCutCellMesh& mesh,
                                  bool solver_ready_cut_cell_mesh = true);

void write_cut_cell_geometry_json(const std::filesystem::path& path,
                                  const LinearOctree& tree,
                                  const ConvexCutCellMesh& mesh,
                                  bool solver_ready_cut_cell_mesh = true);

} // 命名空间 cartmesh
