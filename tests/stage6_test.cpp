#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"
#include "cartmesh/io/ScalableOpenFoamWriter.hpp"
#include "cartmesh/scalable/CompactUniformCutCellMesh.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

[[nodiscard]] std::vector<cartmesh::Triangle> box_triangles(
    const cartmesh::Vec3& minimum, const cartmesh::Vec3& maximum,
    bool reverse = false) {
    const std::array<cartmesh::Vec3, 8> point = {
        cartmesh::Vec3{minimum.x, minimum.y, minimum.z},
        cartmesh::Vec3{maximum.x, minimum.y, minimum.z},
        cartmesh::Vec3{maximum.x, maximum.y, minimum.z},
        cartmesh::Vec3{minimum.x, maximum.y, minimum.z},
        cartmesh::Vec3{minimum.x, minimum.y, maximum.z},
        cartmesh::Vec3{maximum.x, minimum.y, maximum.z},
        cartmesh::Vec3{maximum.x, maximum.y, maximum.z},
        cartmesh::Vec3{minimum.x, maximum.y, maximum.z}};
    const std::array<std::array<std::size_t, 3>, 12> face = {{
        {0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
        {0, 1, 5}, {0, 5, 4}, {3, 7, 6}, {3, 6, 2},
        {0, 4, 7}, {0, 7, 3}, {1, 2, 6}, {1, 6, 5}}};
    std::vector<cartmesh::Triangle> result;
    result.reserve(face.size());
    for (const auto& indices : face) {
        result.emplace_back(
            point[indices[0]], point[reverse ? indices[2] : indices[1]],
            point[reverse ? indices[1] : indices[2]]);
    }
    return result;
}

void test_compact_cube_matches_reference() {
    const cartmesh::SurfaceMesh surface(
        box_triangles({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}),
        cartmesh::SurfaceFormat::ascii_stl, "stage6_cube");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 7);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-0.15, -0.15, -0.15}, {1.15, 1.15, 1.15}),
        10, 10, 10);
    const auto reference =
        cartmesh::build_triangulated_cut_cell_mesh(grid, cutter);
    const auto compact =
        cartmesh::build_compact_uniform_cut_cell_mesh(grid, cutter);
    expect(compact.invariants_pass(),
           "紧凑立方体网格必须通过完整不变量");
    expect(compact.background_cell_count == grid.cell_count() &&
               compact.full_fluid_cell_count == reference.fluid_cells.size() &&
               compact.full_solid_cell_count ==
                   reference.full_solid_cell_count &&
               compact.cut_cell_count == reference.cut_cell_count &&
               compact.global_fluid_region_count ==
                   reference.global_fluid_region_count,
           "紧凑立方体的单元和区域计数必须与基线一致");
    expect(std::abs(compact.total_fluid_volume -
                    reference.total_fluid_volume) < 1.0e-12 &&
               std::abs(compact.total_embedded_boundary_area -
                        reference.total_embedded_boundary_area) < 1.0e-12,
           "紧凑立方体的体积和边界面积必须与基线一致");
}

void test_compact_nested_shell_regions_and_determinism() {
    auto shell = box_triangles({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    auto cavity = box_triangles({0.2, 0.2, 0.2}, {0.8, 0.8, 0.8}, true);
    shell.insert(shell.end(), cavity.begin(), cavity.end());
    const cartmesh::SurfaceMesh surface(
        std::move(shell), cartmesh::SurfaceFormat::ascii_stl,
        "stage6_nested_shell");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 11);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-0.1, -0.1, -0.1}, {1.1, 1.1, 1.1}),
        12, 12, 12);
    const auto first =
        cartmesh::build_compact_uniform_cut_cell_mesh(grid, cutter);
    const auto second =
        cartmesh::build_compact_uniform_cut_cell_mesh(grid, cutter);
    expect(first.invariants_pass() && second.invariants_pass(),
           "紧凑嵌套薄壳必须通过完整不变量");
    expect(first.global_fluid_region_count == 2,
           "紧凑嵌套薄壳必须保留外部与内腔两个流体区");
    expect(first.result_hash_fnv1a64 == second.result_hash_fnv1a64 &&
               first.cell_states == second.cell_states,
           "紧凑网格状态、区域和几何 hash 必须确定");
}

void test_compact_storage_scales_with_surface_cells() {
    const cartmesh::SurfaceMesh surface(
        box_triangles({0.25, 0.25, 0.25}, {0.75, 0.75, 0.75}),
        cartmesh::SurfaceFormat::ascii_stl, "stage6_storage_cube");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 3);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}),
        48, 48, 48);
    const auto compact =
        cartmesh::build_compact_uniform_cut_cell_mesh(grid, cutter);
    expect(compact.invariants_pass(),
           "48^3 紧凑路径必须通过不变量");
    expect(compact.compact_storage_bytes <
               compact.background_cell_count * 128ULL,
           "紧凑稳态存储应显著低于旧的每单元显式面对象：bytes=" +
               std::to_string(compact.compact_storage_bytes) +
               " cells=" +
               std::to_string(compact.background_cell_count));
    expect(compact.explicit_surface_cell_count <
               compact.background_cell_count / 4ULL,
           "显式几何只能跟随表面单元，不能全域实体化");
}

void test_scalable_binary_openfoam_writer() {
    const cartmesh::SurfaceMesh surface(
        box_triangles({0.25, 0.25, 0.25}, {0.75, 0.75, 0.75}),
        cartmesh::SurfaceFormat::ascii_stl, "stage6_binary_cube");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 17);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}),
        9, 9, 9);
    const auto compact =
        cartmesh::build_compact_uniform_cut_cell_mesh(grid, cutter);
    expect(compact.invariants_pass(),
           "二进制导出测例的紧凑网格必须先通过不变量");
    const auto case_directory =
        std::filesystem::temp_directory_path() /
        "cartmesh-stage6-binary-writer-test";
    const auto stats = cartmesh::write_scalable_openfoam_poly_mesh(
        case_directory, grid, compact, {{17U, "cube_wall"}});
    expect(stats.solver_cell_count > 0U &&
               stats.solver_cell_count <= compact.solver_cell_count &&
               stats.point_count > 0U && stats.face_count > 0U &&
               stats.internal_face_count > 0U &&
               stats.boundary_face_count > 0U && stats.written_bytes > 0U,
           "二进制导出必须保留非空控制体，连通凸片并集数不得超过凸片数");
    std::ifstream faces(case_directory / "constant" / "polyMesh" / "faces",
                        std::ios::binary);
    std::string header(512, '\0');
    faces.read(header.data(), static_cast<std::streamsize>(header.size()));
    expect(header.find("format binary;") != std::string::npos &&
               header.find("class faceCompactList;") != std::string::npos,
           "阶段6 faces 必须使用 OpenFOAM faceCompactList 二进制格式");

}

} // 匿名命名空间

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"紧凑立方体与基线等价", test_compact_cube_matches_reference},
        {"嵌套薄壳区域与确定性",
         test_compact_nested_shell_regions_and_determinism},
        {"稳态存储跟随表面而非全域",
         test_compact_storage_scales_with_surface_cells},
        {"流式二进制 OpenFOAM 导出", test_scalable_binary_openfoam_writer},
    };
    std::size_t failures = 0;
    for (const auto& [name, test] : tests) {
        try {
            test();
            std::cout << "[通过] " << name << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "[失败] " << name << "：" << error.what() << '\n';
        }
    }
    std::cout << "阶段6紧凑核心测试数=" << tests.size()
              << " 失败数=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
