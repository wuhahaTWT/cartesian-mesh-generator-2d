#pragma once

#include "cartmesh/geometry/SurfaceDiagnostics.hpp"

#include <filesystem>

namespace cartmesh {

// 写出无效几何位置样本。issue_code：1退化面、2重复面、3边界边、4非流形边、
// 5方向冲突、6非流形顶点、7连通分量方向/嵌套不匹配、8共面面积重叠、
// 9非邻接自交、10非邻接接触、11极小连通分量。该文件是诊断标记，不是修复后的几何。
void write_surface_diagnostic_vtp(const std::filesystem::path& path, const SurfaceMesh& surface,
                                  const SurfaceDiagnostics& diagnostics);

} // 命名空间 cartmesh
