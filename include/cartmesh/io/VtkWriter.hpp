#pragma once

#include "cartmesh/grid/UniformCartesianGrid.hpp"
#include "cartmesh/grid/LinearOctree.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace cartmesh {

struct VtkCellData {
    std::string name;
    std::vector<double> values;
};

void write_vtu(const std::filesystem::path& path, const UniformCartesianGrid& grid,
               const std::vector<VtkCellData>& cell_data = {});

void write_octree_vtu(const std::filesystem::path& path, const LinearOctree& tree,
                      const std::vector<VtkCellData>& cell_data = {});

} // 命名空间 cartmesh
