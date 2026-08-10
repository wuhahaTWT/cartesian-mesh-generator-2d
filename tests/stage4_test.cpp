#include "cartmesh/geometry/SurfaceDiagnostics.hpp"
#include "cartmesh/geometry/TriangleTriangleIntersection.hpp"
#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"
#include "cartmesh/grid/OctreeRefinement.hpp"

#include <functional>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <set>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct TestFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

void expect(bool condition, const std::string& message) {
    if (!condition) throw TestFailure(message);
}

[[nodiscard]] cartmesh::Triangle reversed(const cartmesh::Triangle& triangle) {
    const auto& vertex = triangle.vertices();
    return {vertex[0], vertex[2], vertex[1]};
}

[[nodiscard]] std::vector<cartmesh::Triangle> box_triangles(
    const cartmesh::Vec3& minimum, const cartmesh::Vec3& maximum,
    bool reverse_orientation = false) {
    const std::array<cartmesh::Vec3, 8> point = {
        cartmesh::Vec3{minimum.x, minimum.y, minimum.z},
        cartmesh::Vec3{maximum.x, minimum.y, minimum.z},
        cartmesh::Vec3{maximum.x, maximum.y, minimum.z},
        cartmesh::Vec3{minimum.x, maximum.y, minimum.z},
        cartmesh::Vec3{minimum.x, minimum.y, maximum.z},
        cartmesh::Vec3{maximum.x, minimum.y, maximum.z},
        cartmesh::Vec3{maximum.x, maximum.y, maximum.z},
        cartmesh::Vec3{minimum.x, maximum.y, maximum.z}};
    const std::array<std::array<std::size_t, 3>, 12> face = {
        std::array<std::size_t, 3>{0, 2, 1}, {0, 3, 2}, {4, 5, 6},
        {4, 6, 7}, {0, 1, 5}, {0, 5, 4}, {3, 7, 6}, {3, 6, 2},
        {0, 4, 7}, {0, 7, 3}, {1, 2, 6}, {1, 6, 5}};
    std::vector<cartmesh::Triangle> result;
    for (const auto& triangle : face) {
        result.emplace_back(
            point[triangle[0]],
            point[reverse_orientation ? triangle[2] : triangle[1]],
            point[reverse_orientation ? triangle[1] : triangle[2]]);
    }
    return result;
}

[[nodiscard]] std::vector<cartmesh::Triangle> tube_triangles(
    std::uint32_t segments, double outer_radius, double inner_radius,
    double half_height) {
    std::vector<cartmesh::Vec3> outer_bottom, outer_top, inner_bottom,
        inner_top;
    for (std::uint32_t segment = 0; segment < segments; ++segment) {
        const double angle = 2.0 * std::numbers::pi_v<double> *
                             static_cast<double>(segment) /
                             static_cast<double>(segments);
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        outer_bottom.push_back(
            {outer_radius * cosine, outer_radius * sine, -half_height});
        outer_top.push_back(
            {outer_radius * cosine, outer_radius * sine, half_height});
        inner_bottom.push_back(
            {inner_radius * cosine, inner_radius * sine, -half_height});
        inner_top.push_back(
            {inner_radius * cosine, inner_radius * sine, half_height});
    }
    std::vector<cartmesh::Triangle> triangles;
    triangles.reserve(static_cast<std::size_t>(segments) * 8U);
    for (std::uint32_t segment = 0; segment < segments; ++segment) {
        const auto next = (segment + 1U) % segments;
        triangles.emplace_back(outer_bottom[segment], outer_bottom[next],
                               outer_top[next]);
        triangles.emplace_back(outer_bottom[segment], outer_top[next],
                               outer_top[segment]);
        triangles.emplace_back(inner_bottom[segment], inner_top[next],
                               inner_bottom[next]);
        triangles.emplace_back(inner_bottom[segment], inner_top[segment],
                               inner_top[next]);
        triangles.emplace_back(outer_top[segment], outer_top[next],
                               inner_top[next]);
        triangles.emplace_back(outer_top[segment], inner_top[next],
                               inner_top[segment]);
        triangles.emplace_back(outer_bottom[segment], inner_bottom[next],
                               outer_bottom[next]);
        triangles.emplace_back(outer_bottom[segment], inner_bottom[segment],
                               inner_bottom[next]);
    }
    return triangles;
}

void test_triangle_triangle_relations() {
    const cartmesh::Triangle reference(
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
    const cartmesh::Triangle separated(
        {0.0, 0.0, 0.1}, {1.0, 0.0, 0.1}, {0.0, 1.0, 0.1});
    const cartmesh::Triangle overlap(
        {0.2, 0.2, 0.0}, {0.8, 0.2, 0.0}, {0.2, 0.8, 0.0});
    const cartmesh::Triangle adjacent(
        {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0});
    const cartmesh::Triangle crossing(
        {0.25, 0.1, -1.0}, {0.25, 0.1, 1.0}, {0.25, 0.8, 0.0});
    const cartmesh::Triangle point_touch(
        {0.25, 0.25, 0.0}, {0.25, 0.1, 1.0}, {0.25, 0.4, 1.0});

    expect(cartmesh::classify_triangle_triangle(reference, separated) ==
               cartmesh::TriangleTriangleRelation::disjoint,
           "平行分离三角形必须判为不相交");
    expect(cartmesh::classify_triangle_triangle(reference, overlap) ==
               cartmesh::TriangleTriangleRelation::coplanar_area_overlap,
           "部分共面覆盖必须判为面积重叠");
    expect(cartmesh::classify_triangle_triangle(reference, adjacent) ==
               cartmesh::TriangleTriangleRelation::boundary_contact,
           "合法共边三角形的几何关系必须是边界接触");
    expect(cartmesh::classify_triangle_triangle(reference, crossing) ==
               cartmesh::TriangleTriangleRelation::proper_intersection,
           "穿过三角形内部的非共面片必须判为真自交");
    expect(cartmesh::classify_triangle_triangle(reference, point_touch) ==
               cartmesh::TriangleTriangleRelation::boundary_contact,
           "单点接触必须与穿透相交区分");
    expect(cartmesh::classify_triangle_triangle(reversed(reference),
                                                reversed(crossing)) ==
               cartmesh::TriangleTriangleRelation::proper_intersection,
           "翻转三角形绕序不得改变自交关系");
    const cartmesh::Vec3 shift{1.0e9, 1.0e9, 1.0e9};
    const auto translate = [&](const cartmesh::Triangle& triangle) {
        const auto& vertex = triangle.vertices();
        return cartmesh::Triangle(vertex[0] + shift, vertex[1] + shift,
                                  vertex[2] + shift);
    };
    expect(cartmesh::classify_triangle_triangle(
               translate(reference), translate(adjacent)) ==
               cartmesh::TriangleTriangleRelation::boundary_contact,
           "大坐标平移不得把合法共边误判为面积重叠");
}

void test_surface_pair_diagnostics() {
    const cartmesh::Triangle reference(
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
    const cartmesh::Triangle overlap(
        {0.2, 0.2, 0.0}, {0.8, 0.2, 0.0}, {0.2, 0.8, 0.0});
    const cartmesh::Triangle crossing(
        {0.25, 0.1, -1.0}, {0.25, 0.1, 1.0}, {0.25, 0.8, 0.0});
    const cartmesh::Triangle adjacent(
        {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0});

    const cartmesh::SurfaceMesh overlapping(
        {reference, overlap}, cartmesh::SurfaceFormat::ascii_stl,
        "partial_coplanar_overlap");
    const auto overlap_diagnostics = cartmesh::diagnose_surface(overlapping);
    expect(overlap_diagnostics.overlapping_triangle_pair_count == 1 &&
               overlap_diagnostics.self_intersection_pair_count == 0 &&
               !overlap_diagnostics.overlapping_triangle_examples.empty(),
           "表面诊断必须报告部分共面覆盖的三角形对和位置");

    const cartmesh::SurfaceMesh self_intersecting(
        {reference, crossing}, cartmesh::SurfaceFormat::ascii_stl,
        "proper_self_intersection");
    const auto intersection_diagnostics =
        cartmesh::diagnose_surface(self_intersecting);
    expect(intersection_diagnostics.self_intersection_pair_count == 1 &&
               !intersection_diagnostics.self_intersection_examples.empty(),
           "表面诊断必须报告真自交三角形对和位置");

    const cartmesh::SurfaceMesh legal_adjacent(
        {reference, adjacent}, cartmesh::SurfaceFormat::ascii_stl,
        "legal_shared_edge");
    const auto adjacent_diagnostics = cartmesh::diagnose_surface(legal_adjacent);
    expect(adjacent_diagnostics.overlapping_triangle_pair_count == 0 &&
               adjacent_diagnostics.self_intersection_pair_count == 0 &&
               adjacent_diagnostics.non_adjacent_contact_pair_count == 0,
           "合法共享边不得被工业几何诊断误报为自交或重叠");
}

void test_multiple_components_and_nested_cavity_cut_cells() {
    auto two_boxes = box_triangles({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    auto second_box = box_triangles({2.0, 0.0, 0.0}, {3.0, 1.0, 1.0});
    two_boxes.insert(two_boxes.end(), second_box.begin(), second_box.end());
    const cartmesh::SurfaceMesh disjoint_surface(
        std::move(two_boxes), cartmesh::SurfaceFormat::ascii_stl,
        "two_disjoint_solids");
    const auto disjoint_diagnostics =
        cartmesh::diagnose_surface(disjoint_surface);
    expect(disjoint_diagnostics.connected_component_count == 2 &&
               disjoint_diagnostics.component_orientation_mismatch_count == 0,
           "两个分离外壳必须识别为两个正向材料分量");
    const cartmesh::TriangulatedSurfaceCutter disjoint_cutter(
        disjoint_surface, 31);
    const auto disjoint_cut = disjoint_cutter.cut_box(
        cartmesh::AABB({-0.5, -0.5, -0.5}, {3.5, 1.5, 1.5}));
    expect(std::abs(disjoint_cut.solid_volume - 2.0) < 2.0e-12 &&
               std::abs(disjoint_cut.fluid_volume - 14.0) < 2.0e-12 &&
               std::abs(disjoint_cut.embedded_boundary_area - 12.0) < 2.0e-12,
           "多部件 Cut-cell 必须守恒累计两个实体的体积和边界面积");
    const cartmesh::UniformCartesianGrid disjoint_grid(
        cartmesh::AABB({-0.5, -0.5, -0.5}, {3.5, 1.5, 1.5}), 8, 4, 4);
    const auto disjoint_mesh =
        cartmesh::build_triangulated_cut_cell_mesh(disjoint_grid,
                                                   disjoint_cutter);
    expect(disjoint_mesh.global_fluid_region_count == 1 &&
               disjoint_mesh.global_fluid_region_volumes.size() == 1 &&
               std::abs(disjoint_mesh.global_fluid_region_volumes.front() -
                        disjoint_mesh.total_fluid_volume) < 2.0e-11,
           "两个分离实体之外的流体必须连成同一个全局 region");

    auto nested = box_triangles({0.0, 0.0, 0.0}, {2.0, 2.0, 2.0});
    auto cavity = box_triangles({0.75, 0.75, 0.75}, {1.25, 1.25, 1.25}, true);
    nested.insert(nested.end(), cavity.begin(), cavity.end());
    std::vector<std::uint64_t> boundary_ids(nested.size(), 10);
    std::fill(boundary_ids.begin() + 12, boundary_ids.end(), 20);
    const cartmesh::SurfaceMesh cavity_surface(
        std::move(nested), cartmesh::SurfaceFormat::ascii_stl,
        "solid_with_internal_cavity");
    const auto cavity_diagnostics = cartmesh::diagnose_surface(cavity_surface);
    expect(cavity_diagnostics.connected_component_count == 2 &&
               cavity_diagnostics.components[1].nesting_depth == 1 &&
               cavity_diagnostics.components[1].expected_orientation_sign == -1 &&
               cavity_diagnostics.component_orientation_mismatch_count == 0,
           "内腔必须由嵌套深度和反向壳层识别");
    const cartmesh::TriangulatedSurfaceCutter cavity_cutter(
        cavity_surface, boundary_ids);
    const cartmesh::AABB cavity_domain(
        {-0.5, -0.5, -0.5}, {2.5, 2.5, 2.5});
    const auto cavity_cut = cavity_cutter.cut_box(cavity_domain);
    expect(std::abs(cavity_cut.solid_volume - 7.875) < 3.0e-12 &&
               std::abs(cavity_cut.fluid_volume - 19.125) < 3.0e-12 &&
               std::abs(cavity_cut.embedded_boundary_area - 25.5) < 3.0e-12,
           "嵌套内腔必须从材料体积中扣除并保留内外边界面积");
    bool outer_boundary_seen = false;
    bool cavity_boundary_seen = false;
    for (const auto& face : cavity_cut.embedded_boundary_faces) {
        outer_boundary_seen = outer_boundary_seen || face.boundary_id == 10;
        cavity_boundary_seen = cavity_boundary_seen || face.boundary_id == 20;
    }
    expect(outer_boundary_seen && cavity_boundary_seen,
           "多壳层 Cut-cell 必须保留外壁和内腔不同 boundary ID");

    const cartmesh::UniformCartesianGrid grid(cavity_domain, 8, 8, 8);
    const auto mesh =
        cartmesh::build_triangulated_cut_cell_mesh(grid, cavity_cutter);
    expect(std::abs(mesh.total_fluid_volume - 19.125) < 2.0e-11 &&
               std::abs(mesh.total_embedded_boundary_area - 25.5) < 2.0e-11 &&
               mesh.nonclosed_cell_count == 0 &&
               mesh.shared_face_mismatch_count == 0 &&
               mesh.component_analysis_pending_cell_count == 0,
           "嵌套内腔均匀网格必须保持体积、边界和邻接闭合");
    expect(mesh.global_fluid_region_count == 2 &&
               mesh.global_fluid_region_volumes.size() == 2,
           "密闭内腔和外部流体必须生成两个全局 region");
    auto region_volumes = mesh.global_fluid_region_volumes;
    std::sort(region_volumes.begin(), region_volumes.end());
    expect(std::abs(region_volumes[0] - 0.125) < 2.0e-11 &&
               std::abs(region_volumes[1] - 19.0) < 2.0e-11 &&
               std::abs(region_volumes[0] + region_volumes[1] -
                        mesh.total_fluid_volume) < 2.0e-11,
           "全局 region 体积必须分别等于内腔与外部流体解析真值");
    for (const auto& connection : mesh.component_internal_faces) {
        const auto& first = mesh.fluid_cells[connection.first_fluid_cell_index];
        const auto& second = mesh.fluid_cells[connection.second_fluid_cell_index];
        expect(connection.global_region_id ==
                   first.fluid_component_region_ids[connection.first_component_id] &&
                   connection.global_region_id ==
                   second.fluid_component_region_ids[connection.second_component_id],
               "跨单元分量连接必须只连接同一全局 region");
    }
}

void test_small_disconnected_component_diagnostics() {
    auto triangles = box_triangles({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    auto fragment = box_triangles({2.0, 0.0, 0.0},
                                  {2.0000001, 0.0000001, 0.0000001});
    triangles.insert(triangles.end(), fragment.begin(), fragment.end());
    const cartmesh::SurfaceMesh surface(
        std::move(triangles), cartmesh::SurfaceFormat::ascii_stl,
        "unit_cube_with_tiny_fragment");
    const auto diagnostics = cartmesh::diagnose_surface(surface);
    expect(diagnostics.connected_component_count == 2 &&
               diagnostics.small_component_count == 1 &&
               diagnostics.minimum_component_bounding_diagonal <
                   diagnostics.small_component_diagonal_threshold &&
               diagnostics.minimum_component_absolute_volume < 2.0e-21,
           "极小独立碎片必须按包围盒尺度和绝对体积显式统计");
    expect(diagnostics.components[1].surface_area > 0.0 &&
               diagnostics.components[1].bounds.center().x > 1.5,
           "极小碎片诊断必须保留面积、包围盒和位置");
}

void test_subcell_thin_wall_and_small_hole_pipeline() {
    // 壁厚 0.05，低于背景单元边长 0.1；仍必须保留内外两个流体区。
    auto thin_shell = box_triangles({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    auto thin_cavity = box_triangles({0.05, 0.05, 0.05},
                                     {0.95, 0.95, 0.95}, true);
    thin_shell.insert(thin_shell.end(), thin_cavity.begin(), thin_cavity.end());
    const cartmesh::SurfaceMesh shell_surface(
        std::move(thin_shell), cartmesh::SurfaceFormat::ascii_stl,
        "subcell_thin_shell");
    const cartmesh::TriangulatedSurfaceCutter shell_cutter(shell_surface, 41);
    const cartmesh::AABB shell_domain({-0.1, -0.1, -0.1},
                                     {1.1, 1.1, 1.1});
    const cartmesh::UniformCartesianGrid shell_grid(shell_domain, 12, 12, 12);
    const auto shell_mesh =
        cartmesh::build_triangulated_cut_cell_mesh(shell_grid, shell_cutter);
    const double cavity_volume = 0.9 * 0.9 * 0.9;
    const double solid_volume = 1.0 - cavity_volume;
    expect(shell_mesh.global_fluid_region_count == 2 &&
               std::abs(shell_domain.volume() - shell_mesh.total_fluid_volume -
                        solid_volume) < 3.0e-11 &&
               std::abs(shell_mesh.total_embedded_boundary_area -
                        (6.0 + 6.0 * 0.9 * 0.9)) < 3.0e-11,
           "低于背景单元尺寸的薄壁必须保留内腔、体积和两侧边界");

    constexpr double inner_radius = 0.06;
    const cartmesh::SurfaceMesh tube_surface(
        tube_triangles(24, 0.5, inner_radius, 0.5),
        cartmesh::SurfaceFormat::ascii_stl, "stage4_small_hole_tube");
    const auto tube_diagnostics = cartmesh::diagnose_surface(tube_surface);
    expect(tube_diagnostics.valid_for_stage1_classification(),
           "小孔端到端案例必须是封闭定向圆管");
    const cartmesh::AABB tube_domain({-0.6, -0.6, -0.6},
                                    {0.6, 0.6, 0.6});
    cartmesh::LinearOctree tree(tube_domain, 2, 6);
    cartmesh::OctreeRefinementConfiguration refinement;
    refinement.surface_target_level = 4;
    refinement.gap = cartmesh::GapRefinementRule{0.14, -0.9, 0.8, 4};
    const cartmesh::TriangleBvh tube_bvh(tube_surface, 8);
    const auto adaptation =
        cartmesh::OctreeRefinementEngine(refinement, &tube_bvh).apply(tree);
    expect(adaptation.gap_resolution_failure_count == 0 &&
               adaptation.gap_rule_hits > 0 &&
               adaptation.maximum_required_gap_level == 6,
           "内径0.12的小孔必须自适应到至少四单元跨孔");
    std::set<cartmesh::OctreeNodeCode> hole_leaves;
    const std::uint32_t maximum_coordinate = 1U << tree.maximum_level();
    for (std::uint32_t sample = 0; sample < 1000; ++sample) {
        const double x = -0.059 + 0.118 * static_cast<double>(sample) / 999.0;
        const auto coordinate = [&](double value, double minimum,
                                    double extent) {
            return static_cast<std::uint32_t>(std::clamp(
                std::floor((value - minimum) / extent * maximum_coordinate),
                0.0, static_cast<double>(maximum_coordinate - 1U)));
        };
        const auto leaf = tree.find_leaf_covering_maximum_level_cell(
            coordinate(x, -0.6, 1.2), coordinate(0.0, -0.6, 1.2),
            coordinate(0.0, -0.6, 1.2));
        expect(leaf.has_value(), "小孔直径采样必须落入唯一叶单元");
        hole_leaves.insert(tree.leaf_code(*leaf));
    }
    expect(hole_leaves.size() >= 4,
           "小孔直径必须实际穿过至少四个自适应叶单元");
    const cartmesh::TriangulatedSurfaceCutter tube_cutter(tube_surface, 50);
    const auto tube_mesh =
        cartmesh::build_triangulated_cut_cell_mesh(tree, tube_cutter, 1.0e-11);
    bool center_hole_fluid = false;
    for (const auto& cell : tube_mesh.fluid_cells) {
        const auto leaf_box = tree.cell_bounds(tree.leaf_code(cell.background_cell_id));
        if (leaf_box.contains({0.0, 0.0, 0.0})) {
            center_hole_fluid = cell.volume > 0.0;
            break;
        }
    }
    expect(center_hole_fluid && tube_mesh.global_fluid_region_count == 1 &&
               tube_mesh.nonclosed_cell_count == 0 &&
               tube_mesh.shared_face_mismatch_count == 0 &&
               tube_mesh.component_analysis_pending_cell_count == 0,
           "自适应 Cut-cell 必须保留小孔中心流体并与外部连通：center=" +
               std::to_string(center_hole_fluid) + " regions=" +
               std::to_string(tube_mesh.global_fluid_region_count) +
               " nonclosed=" + std::to_string(tube_mesh.nonclosed_cell_count) +
               " mismatch=" +
               std::to_string(tube_mesh.shared_face_mismatch_count) +
               " maxArea=" +
               std::to_string(tube_mesh.maximum_shared_face_area_mismatch) +
               " maxCentroid=" +
               std::to_string(tube_mesh.maximum_shared_face_centroid_mismatch) +
               " pending=" +
               std::to_string(tube_mesh.component_analysis_pending_cell_count));
}

void test_large_tolerance_reports_unresolved_components() {
    auto thin_shell = box_triangles({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    auto thin_cavity = box_triangles({0.05, 0.05, 0.05},
                                     {0.95, 0.95, 0.95}, true);
    thin_shell.insert(thin_shell.end(), thin_cavity.begin(), thin_cavity.end());
    const cartmesh::SurfaceMesh surface(
        std::move(thin_shell), cartmesh::SurfaceFormat::ascii_stl,
        "thin_shell_large_tolerance");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 41);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-0.1, -0.1, -0.1}, {1.1, 1.1, 1.1}), 12, 12, 12);
    const auto mesh =
        cartmesh::build_triangulated_cut_cell_mesh(grid, cutter, 0.04);
    expect(mesh.component_analysis_pending_cell_count > 0,
           "过大的几何容差丢弃流体分片时必须报告未解析分量，而不是崩溃");
    for (const auto& connection : mesh.component_internal_faces) {
        const auto& first = mesh.fluid_cells[connection.first_fluid_cell_index];
        const auto& second = mesh.fluid_cells[connection.second_fluid_cell_index];
        expect(connection.first_component_id <
                   first.fluid_component_region_ids.size() &&
                   connection.second_component_id <
                   second.fluid_component_region_ids.size(),
               "分量连接不得引用被容差丢弃的空分量");
    }
}

void test_thin_shell_translation_and_scale_invariance() {
    const auto build_shell = [](const cartmesh::Vec3& origin, double scale) {
        auto shell = box_triangles(origin,
                                   origin + cartmesh::Vec3{scale, scale, scale});
        const double inset = 0.05 * scale;
        auto cavity = box_triangles(
            origin + cartmesh::Vec3{inset, inset, inset},
            origin + cartmesh::Vec3{scale - inset, scale - inset,
                                    scale - inset},
            true);
        shell.insert(shell.end(), cavity.begin(), cavity.end());
        const cartmesh::SurfaceMesh surface(
            std::move(shell), cartmesh::SurfaceFormat::ascii_stl,
            "translated_scaled_thin_shell");
        const cartmesh::Vec3 padding{0.1 * scale, 0.1 * scale, 0.1 * scale};
        const cartmesh::UniformCartesianGrid grid(
            cartmesh::AABB(origin - padding,
                           origin + cartmesh::Vec3{scale, scale, scale} +
                               padding),
            12, 12, 12);
        return cartmesh::build_triangulated_cut_cell_mesh(
            grid, cartmesh::TriangulatedSurfaceCutter(surface, 41));
    };

    const auto translated =
        build_shell({1.0e6, -2.0e6, 3.0e6}, 1.0);
    expect(translated.cut_cell_count == 488 &&
               translated.global_fluid_region_count == 2 &&
               translated.nonclosed_cell_count == 0 &&
               translated.shared_face_mismatch_count == 0 &&
               translated.component_analysis_pending_cell_count == 0 &&
               std::abs(translated.total_fluid_volume - 1.457) < 2.0e-9 &&
               std::abs(translated.total_embedded_boundary_area - 10.86) <
                   2.0e-8,
           "大平移薄壁必须保持体积、面积、两流体区和闭合拓扑：cut=" +
               std::to_string(translated.cut_cell_count) +
               " regions=" +
               std::to_string(translated.global_fluid_region_count) +
               " nonclosed=" +
               std::to_string(translated.nonclosed_cell_count) +
               " shared=" +
               std::to_string(translated.shared_face_mismatch_count) +
               " pending=" +
               std::to_string(translated.component_analysis_pending_cell_count) +
               " volume=" + std::to_string(translated.total_fluid_volume) +
               " area=" +
               std::to_string(translated.total_embedded_boundary_area));

    constexpr double scale = 1.0e-7;
    const auto tiny = build_shell({0.0, 0.0, 0.0}, scale);
    expect(tiny.cut_cell_count == 488 && tiny.global_fluid_region_count == 2 &&
               tiny.nonclosed_cell_count == 0 &&
               tiny.shared_face_mismatch_count == 0 &&
               tiny.component_analysis_pending_cell_count == 0 &&
               std::abs(tiny.total_fluid_volume /
                            (scale * scale * scale) -
                        1.457) < 2.0e-11 &&
               std::abs(tiny.total_embedded_boundary_area /
                            (scale * scale) -
                        10.86) < 2.0e-10,
           "小尺度薄壁必须保持归一化体积、面积、两流体区和闭合拓扑");
}

} // 匿名命名空间

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"三角形对窄相关系", test_triangle_triangle_relations},
        {"表面重叠与自交诊断", test_surface_pair_diagnostics},
        {"多部件与嵌套内腔 Cut-cell", test_multiple_components_and_nested_cavity_cut_cells},
        {"极小独立碎片诊断", test_small_disconnected_component_diagnostics},
        {"亚单元薄壁与小孔全流程", test_subcell_thin_wall_and_small_hole_pipeline},
        {"大容差未解析分量", test_large_tolerance_reports_unresolved_components},
        {"薄壁平移与缩放不变量", test_thin_shell_translation_and_scale_invariance},
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
    std::cout << "阶段4工业几何测试数=" << tests.size()
              << " 失败数=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
