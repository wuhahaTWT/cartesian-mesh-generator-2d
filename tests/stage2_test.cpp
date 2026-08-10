#include "cartmesh/grid/LinearOctree.hpp"
#include "cartmesh/grid/OctreeRefinement.hpp"
#include "cartmesh/geometry/SurfaceMesh.hpp"
#include "cartmesh/classify/SurfaceClassifier.hpp"
#include "cartmesh/spatial/TriangleBvh.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct TestFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

void expect(bool condition, const std::string& message) {
    if (!condition) {
        throw TestFailure(message);
    }
}

template <typename Function>
void expect_throw(Function&& function, const std::string& message) {
    try {
        function();
    } catch (const std::exception&) {
        return;
    }
    throw TestFailure(message);
}

cartmesh::SurfaceMesh make_unit_cube_surface() {
    using cartmesh::Triangle;
    using cartmesh::Vec3;
    const Vec3 p000{0.0, 0.0, 0.0};
    const Vec3 p100{1.0, 0.0, 0.0};
    const Vec3 p010{0.0, 1.0, 0.0};
    const Vec3 p110{1.0, 1.0, 0.0};
    const Vec3 p001{0.0, 0.0, 1.0};
    const Vec3 p101{1.0, 0.0, 1.0};
    const Vec3 p011{0.0, 1.0, 1.0};
    const Vec3 p111{1.0, 1.0, 1.0};
    return cartmesh::SurfaceMesh(
        {Triangle(p000, p010, p110), Triangle(p000, p110, p100),
         Triangle(p001, p101, p111), Triangle(p001, p111, p011),
         Triangle(p000, p100, p101), Triangle(p000, p101, p001),
         Triangle(p010, p011, p111), Triangle(p010, p111, p110),
         Triangle(p000, p001, p011), Triangle(p000, p011, p010),
         Triangle(p100, p110, p111), Triangle(p100, p111, p101)},
        cartmesh::SurfaceFormat::ascii_stl, "unit_cube");
}

std::uint8_t leaf_level_at_point(const cartmesh::LinearOctree& tree,
                                 const cartmesh::Vec3& point) {
    for (const auto code : tree.leaf_codes()) {
        const auto bounds = tree.cell_bounds(code);
        const auto minimum = bounds.minimum();
        const auto maximum = bounds.maximum();
        if (point.x >= minimum.x && point.x < maximum.x && point.y >= minimum.y &&
            point.y < maximum.y && point.z >= minimum.z && point.z < maximum.z) {
            return cartmesh::decode_octree_node(code).level;
        }
    }
    throw TestFailure("找不到覆盖内部查询点的八叉树叶子");
}

void test_level_aware_morton_code() {
    const std::array<std::uint8_t, 10> tested_levels{0, 1, 2, 3, 4, 5, 6, 7, 8, 21};
    for (const std::uint8_t level : tested_levels) {
        const std::uint32_t maximum = (1U << level) - 1U;
        const auto code = cartmesh::encode_octree_node(level, maximum, maximum / 2U, 0);
        const auto decoded = cartmesh::decode_octree_node(code);
        expect(decoded.level == level && decoded.x == maximum &&
                   decoded.y == maximum / 2U && decoded.z == 0,
               "层级 Morton 节点码必须无损往返");
    }
    const auto root = cartmesh::encode_octree_node(0, 0, 0, 0);
    for (std::uint8_t child = 0; child < 8U; ++child) {
        const auto code = cartmesh::octree_child(root, child);
        expect(cartmesh::octree_parent(code) == root, "子节点的父节点必须回到根");
        expect(code != root, "层级哨兵必须避免父子 Morton 冲突");
    }
    expect_throw([] { static_cast<void>(cartmesh::decode_octree_node(0)); },
                 "零节点码必须失败");
    expect_throw([] { static_cast<void>(cartmesh::encode_octree_node(2, 4, 0, 0)); },
                 "越界层级坐标必须失败");
}

void test_initial_partition_and_bounds() {
    const cartmesh::LinearOctree tree(cartmesh::AABB({-1.0, -2.0, -3.0}, {1.0, 2.0, 3.0}),
                                      2, 5);
    expect(tree.leaf_count() == 64, "基础层 2 必须包含 4^3 个叶子");
    expect(tree.validate_partition(), "初始线性八叉树必须完整无重叠分割根域");
    expect(tree.compact_storage_bytes() >= tree.leaf_count() * 8U,
           "紧凑叶数组必须按单个 64 位节点码计量");
    const auto first = cartmesh::decode_octree_node(tree.leaf_code(0));
    expect(first.level == 2 && first.x == 0 && first.y == 0 && first.z == 0,
           "第一个叶子必须是 Morton 原点");
    const auto bounds = tree.cell_bounds(tree.leaf_code(0));
    expect(bounds.minimum().x == -1.0 && bounds.minimum().y == -2.0 &&
               bounds.minimum().z == -3.0,
           "首叶最小角必须等于根域最小角");
    expect(bounds.maximum().x == -0.5 && bounds.maximum().y == -1.0 &&
               bounds.maximum().z == -1.5,
           "层 2 叶边长必须是根域各轴的四分之一");
}

void test_refine_and_coarsen() {
    cartmesh::LinearOctree tree(cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 0, 4);
    const auto root = tree.leaf_code(0);
    expect(tree.refine_leaf(root), "根叶必须可细分");
    expect(tree.leaf_count() == 8 && tree.validate_partition(),
           "一次细分必须用八个子叶替代父叶");
    const auto first_child = cartmesh::octree_child(root, 0);
    expect(tree.refine_leaf(first_child), "一级叶必须继续细分");
    expect(tree.leaf_count() == 15 && tree.validate_partition(),
           "第二次细分净增加七个叶子");
    expect(tree.coarsen_parent(first_child), "完整八兄弟必须可粗化");
    expect(tree.leaf_count() == 8 && tree.validate_partition(),
           "粗化必须恢复一级八叶分区");
    expect(!tree.coarsen_parent(first_child), "子孙不全时不得粗化");
}

void test_face_neighbors_across_levels() {
    cartmesh::LinearOctree tree(cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 1, 4);
    const auto coarse = cartmesh::encode_octree_node(1, 0, 0, 0);
    const auto refined = cartmesh::encode_octree_node(1, 1, 0, 0);
    expect(tree.refine_leaf(refined), "相邻一级叶必须可局部细分");
    const auto neighbors = tree.face_neighbors(coarse, cartmesh::FaceDirection::positive_x);
    expect(neighbors.size() == 4, "粗叶一个面必须枚举四个细一级邻居");
    for (const auto neighbor : neighbors) {
        expect(cartmesh::decode_octree_node(neighbor).level == 2,
               "粗细界面的四个邻居必须都位于下一层");
        const auto reverse =
            tree.face_neighbors(neighbor, cartmesh::FaceDirection::negative_x);
        expect(reverse.size() == 1 && reverse.front() == coarse,
               "细叶反向查询必须返回同一个粗邻居");
    }
    expect(tree.face_neighbors(coarse, cartmesh::FaceDirection::negative_x).empty(),
           "根域边界之外不得伪造邻居");
}

void test_deterministic_refinement_and_balance() {
    const auto build = [] {
        cartmesh::LinearOctree tree(cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 1,
                                    5);
        const auto refinement = tree.refine_to_desired_levels(
            [](cartmesh::OctreeNodeCode, const cartmesh::AABB& bounds) {
                const auto center = bounds.center();
                return center.x < 0.26 && center.y < 0.26 && center.z < 0.26
                           ? std::uint8_t{5}
                           : std::uint8_t{1};
            });
        expect(refinement.split_count > 0 && tree.validate_partition(),
               "局部目标层级必须产生有效自适应分区");
        const auto before = tree.check_face_balance();
        expect(!before.balanced && before.maximum_level_difference > 1,
               "深层局部细化在平衡前必须暴露层级跳变");
        const auto balance = tree.balance_faces_2_to_1();
        expect(balance.split_count > 0 && balance.iteration_count > 0,
               "2:1 平衡必须实际扩散细化");
        const auto after = tree.check_face_balance();
        expect(after.balanced && after.maximum_level_difference <= 1,
               "平衡后所有共享面层级差必须不超过一");
        expect(tree.validate_partition(), "平衡不得破坏根域分区");
        return tree;
    };
    const auto first = build();
    const auto second = build();
    expect(first.leaf_codes().size() == second.leaf_codes().size(),
           "相同细化规则必须得到相同叶数");
    expect(std::equal(first.leaf_codes().begin(), first.leaf_codes().end(),
                      second.leaf_codes().begin()),
           "相同细化规则必须得到逐字节相同的 Morton 叶顺序");
    expect(first.result_hash_fnv1a64() == second.result_hash_fnv1a64(),
           "相同八叉树结果哈希必须稳定");
}

void test_bvh_nearest_surface_distance() {
    const auto surface = make_unit_cube_surface();
    const cartmesh::TriangleBvh bvh(surface, 2);
    expect(std::abs(bvh.distance_to_surface({0.5, 0.5, 0.5}) - 0.5) < 1.0e-14,
           "立方体中心到表面的 BVH 最短距离必须为 0.5");
    expect(std::abs(bvh.distance_to_surface({1.25, 0.5, 0.5}) - 0.25) < 1.0e-14,
           "立方体外点的 BVH 最短距离必须正确");
    expect(bvh.distance_to_surface({0.3, 0.7, 0.0}) < 1.0e-14,
           "表面点的最短距离必须为零");
    expect_throw([&] {
        static_cast<void>(bvh.distance_to_surface(
            {0.0, 0.0, std::numeric_limits<double>::quiet_NaN()}));
    },
                 "非有限距离查询必须失败");
}

void test_user_region_refinement() {
    cartmesh::LinearOctree tree(cartmesh::AABB({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}),
                                1, 5);
    cartmesh::OctreeRefinementConfiguration configuration;
    configuration.boxes.push_back(
        {cartmesh::AABB({-0.9, -0.9, -0.9}, {-0.7, -0.7, -0.7}), 4});
    configuration.spheres.push_back({{0.5, 0.5, 0.5}, 0.14, 5});
    configuration.cylinders.push_back(
        {{-0.5, 0.5, -0.5}, {-0.5, 0.5, 0.5}, 0.12, 4});
    const cartmesh::OctreeRefinementEngine engine(configuration);
    const auto statistics = engine.apply(tree);
    expect(statistics.user_region_rule_hits > 0,
           "用户 box/sphere/cylinder 区域必须实际命中叶子");
    expect(leaf_level_at_point(tree, {-0.8, -0.8, -0.8}) == 4,
           "box 内部点必须达到指定层级");
    expect(leaf_level_at_point(tree, {0.5, 0.5, 0.5}) == 5,
           "sphere 圆心必须达到指定层级");
    expect(leaf_level_at_point(tree, {-0.5, 0.5, 0.0}) == 4,
           "cylinder 轴线内部点必须达到指定层级");
    expect(tree.validate_partition() && tree.check_face_balance().balanced,
           "用户区域细化后必须保持完整分区和面 2:1 平衡");
}

void test_surface_distance_and_curvature_refinement() {
    const auto surface = make_unit_cube_surface();
    const cartmesh::TriangleBvh bvh(surface, 2);
    cartmesh::LinearOctree tree(cartmesh::AABB({-0.5, -0.5, -0.5}, {1.5, 1.5, 1.5}),
                                1, 5);
    cartmesh::OctreeRefinementConfiguration configuration;
    configuration.surface_target_level = 5;
    configuration.distance_bands.push_back({0.2, 4});
    configuration.curvature = cartmesh::CurvatureRefinementRule{35.0, 1.25, 5};
    const cartmesh::OctreeRefinementEngine engine(configuration, &bvh);
    const auto statistics = engine.apply(tree);
    expect(statistics.surface_rule_hits > 0 && statistics.distance_rule_hits > 0 &&
               statistics.curvature_rule_hits > 0,
           "表面、距离和曲率规则必须分别命中叶子");
    for (const auto code : tree.leaf_codes()) {
        if (bvh.intersects_surface(tree.cell_bounds(code))) {
            expect(cartmesh::decode_octree_node(code).level == 5,
                   "所有精确表面相交叶必须达到表面目标层级");
        }
    }
    expect(leaf_level_at_point(tree, {0.01, 0.01, 0.01}) == 5,
           "立方体锐角邻域必须被曲率规则细化");
    expect(tree.validate_partition() && tree.check_face_balance().balanced,
           "几何规则细化后必须保持完整分区和 2:1 平衡");

    cartmesh::LinearOctree distance_tree(
        cartmesh::AABB({-0.5, -0.5, -0.5}, {1.5, 1.5, 1.5}), 1, 5);
    cartmesh::OctreeRefinementConfiguration distance_only;
    distance_only.distance_bands.push_back({0.2, 4});
    static_cast<void>(cartmesh::OctreeRefinementEngine(distance_only, &bvh).apply(distance_tree));
    expect(leaf_level_at_point(distance_tree, {0.5, 0.5, 1.1}) >= 4 &&
               leaf_level_at_point(distance_tree, {-0.4, -0.4, -0.4}) < 4,
           "仅启用距离规则时必须区分近表面点和远离表面的点");

    cartmesh::LinearOctree curvature_tree(
        cartmesh::AABB({-0.5, -0.5, -0.5}, {1.5, 1.5, 1.5}), 1, 5);
    cartmesh::OctreeRefinementConfiguration curvature_only;
    curvature_only.curvature = cartmesh::CurvatureRefinementRule{35.0, 1.25, 5};
    static_cast<void>(
        cartmesh::OctreeRefinementEngine(curvature_only, &bvh).apply(curvature_tree));
    expect(leaf_level_at_point(curvature_tree, {0.01, 0.01, 0.01}) == 5 &&
               leaf_level_at_point(curvature_tree, {0.5, 0.5, -0.1}) < 5,
           "仅启用法向变化规则时锐角必须比平面中部更细");
}

void test_gap_protection_refinement() {
    using cartmesh::Triangle;
    const cartmesh::SurfaceMesh opposing_planes(
        {Triangle({0.4, 0.0, 0.0}, {0.4, 1.0, 0.0}, {0.4, 0.0, 1.0}),
         Triangle({0.6, 0.0, 0.0}, {0.6, 0.0, 1.0}, {0.6, 1.0, 0.0})},
        cartmesh::SurfaceFormat::ascii_stl, "opposing_gap");
    const cartmesh::TriangleBvh bvh(opposing_planes, 1);
    cartmesh::LinearOctree tree(cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 1,
                                5);
    cartmesh::OctreeRefinementConfiguration configuration;
    configuration.gap = cartmesh::GapRefinementRule{0.3, -0.9, 0.9, 4};
    const cartmesh::OctreeRefinementEngine engine(configuration, &bvh);
    const auto statistics = engine.apply(tree);
    expect(statistics.gap_rule_hits > 0, "相对表面必须触发狭缝保护");
    expect(statistics.gap_resolution_failure_count == 0 &&
               statistics.maximum_required_gap_level == 5,
           "可满足的 0.2 间隙必须报告所需层级 5 且无解析失败");
    expect(leaf_level_at_point(tree, {0.5, 0.3, 0.3}) == 5,
           "0.2 间隙按至少 4 个单元的要求必须达到最大层");
    expect(tree.validate_partition() && tree.check_face_balance().balanced,
           "狭缝保护后必须保持完整分区和 2:1 平衡");
    cartmesh::OctreeRefinementConfiguration invalid;
    invalid.gap = cartmesh::GapRefinementRule{0.0, -0.8, 0.5, 4};
    const cartmesh::OctreeRefinementEngine invalid_engine(invalid, &bvh);
    cartmesh::LinearOctree invalid_tree(cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}),
                                        1, 5);
    expect_throw([&] { static_cast<void>(invalid_engine.apply(invalid_tree)); },
                 "非正狭缝搜索距离必须被拒绝");

    cartmesh::OctreeRefinementConfiguration unresolved_configuration;
    unresolved_configuration.gap = cartmesh::GapRefinementRule{0.3, -0.9, 0.9, 16};
    const cartmesh::OctreeRefinementEngine unresolved_engine(unresolved_configuration, &bvh);
    cartmesh::LinearOctree unresolved_tree(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 1, 4);
    const auto unresolved = unresolved_engine.apply(unresolved_tree);
    expect(unresolved.gap_resolution_failure_count > 0 &&
               unresolved.maximum_required_gap_level == 7,
           "最大层级无法在间隙内放入指定单元数时必须显式报告");

    const cartmesh::Vec3 oblique_normal =
        cartmesh::Vec3{1.0, 1.0, 1.0} / std::sqrt(3.0);
    const cartmesh::Vec3 tangent_u =
        cartmesh::Vec3{1.0, -1.0, 0.0} / std::sqrt(2.0);
    const cartmesh::Vec3 tangent_v = cartmesh::cross(oblique_normal, tangent_u);
    const cartmesh::Vec3 middle{0.49, 0.49, 0.49};
    const double oblique_width = 0.13;
    const auto plane_triangle = [&](const cartmesh::Vec3& center, bool reverse) {
        const auto first = center - tangent_u * 2.0 - tangent_v * 2.0;
        const auto second = center + tangent_u * 2.0 - tangent_v * 2.0;
        const auto third = center + tangent_v * 2.0;
        return reverse ? cartmesh::Triangle(first, third, second)
                       : cartmesh::Triangle(first, second, third);
    };
    const cartmesh::SurfaceMesh oblique_planes(
        {plane_triangle(middle - oblique_normal * (0.5 * oblique_width), false),
         plane_triangle(middle + oblique_normal * (0.5 * oblique_width), true)},
        cartmesh::SurfaceFormat::ascii_stl, "oblique_gap");
    const cartmesh::TriangleBvh oblique_bvh(oblique_planes, 1);
    cartmesh::LinearOctree oblique_tree(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 1, 6);
    cartmesh::OctreeRefinementConfiguration oblique_configuration;
    oblique_configuration.gap = cartmesh::GapRefinementRule{0.2, -0.9, 0.9, 4};
    const auto oblique_statistics =
        cartmesh::OctreeRefinementEngine(oblique_configuration, &oblique_bvh)
            .apply(oblique_tree);
    expect(oblique_statistics.maximum_required_gap_level == 6 &&
               oblique_statistics.gap_resolution_failure_count == 0,
           "斜向间隙必须按单元在法向上的投影厚度计算所需层级");
    std::set<cartmesh::OctreeNodeCode> crossed_leaves;
    const std::uint32_t maximum_coordinate = 1U << oblique_tree.maximum_level();
    for (std::uint32_t sample = 1; sample < 1000; ++sample) {
        const double fraction = static_cast<double>(sample) / 1000.0;
        const auto point =
            middle + oblique_normal * ((fraction - 0.5) * oblique_width);
        const auto coordinate = [&](double value) {
            return static_cast<std::uint32_t>(std::clamp(
                std::floor(value * maximum_coordinate), 0.0,
                static_cast<double>(maximum_coordinate - 1U)));
        };
        const auto leaf = oblique_tree.find_leaf_covering_maximum_level_cell(
            coordinate(point.x), coordinate(point.y), coordinate(point.z));
        expect(leaf.has_value(), "斜向采样线上的点必须落入唯一八叉树叶");
        crossed_leaves.insert(oblique_tree.leaf_code(*leaf));
    }
    expect(crossed_leaves.size() >= 4,
           "声称已解析的斜向间隙必须实际穿过至少请求数量的不同叶单元");

    std::vector<cartmesh::Triangle> small_hole_wall;
    constexpr std::size_t segments = 24;
    constexpr double radius = 0.08;
    for (std::size_t segment = 0; segment < segments; ++segment) {
        const auto angle = 2.0 * std::numbers::pi * static_cast<double>(segment) /
                           static_cast<double>(segments);
        const auto next_angle =
            2.0 * std::numbers::pi * static_cast<double>((segment + 1U) % segments) /
            static_cast<double>(segments);
        const cartmesh::Vec3 bottom{radius * std::cos(angle), radius * std::sin(angle), -0.5};
        const cartmesh::Vec3 top{radius * std::cos(angle), radius * std::sin(angle), 0.5};
        const cartmesh::Vec3 next_bottom{radius * std::cos(next_angle),
                                         radius * std::sin(next_angle), -0.5};
        const cartmesh::Vec3 next_top{radius * std::cos(next_angle),
                                      radius * std::sin(next_angle), 0.5};
        small_hole_wall.emplace_back(bottom, next_top, next_bottom);
        small_hole_wall.emplace_back(bottom, top, next_top);
    }
    const cartmesh::SurfaceMesh small_hole_surface(
        std::move(small_hole_wall), cartmesh::SurfaceFormat::ascii_stl,
        "small_hole_inner_wall");
    const cartmesh::TriangleBvh small_hole_bvh(small_hole_surface, 2);
    cartmesh::LinearOctree small_hole_tree(
        cartmesh::AABB({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}), 2, 7);
    cartmesh::OctreeRefinementConfiguration small_hole_configuration;
    small_hole_configuration.gap = cartmesh::GapRefinementRule{0.18, -0.9, 0.8, 4};
    const auto small_hole_statistics =
        cartmesh::OctreeRefinementEngine(small_hole_configuration, &small_hole_bvh)
            .apply(small_hole_tree);
    expect(small_hole_statistics.gap_rule_hits > 0 &&
               small_hole_statistics.gap_resolution_failure_count == 0 &&
               leaf_level_at_point(small_hole_tree, {0.0, 0.0, 0.0}) >= 6,
           "直径远小于根单元的解析小孔必须触发至少四单元跨孔细化");
}

void test_adaptive_leaf_classification() {
    const auto surface = make_unit_cube_surface();
    const cartmesh::TriangleBvh bvh(surface, 2);
    const cartmesh::SurfaceClassifier classifier(bvh);
    const auto build_and_classify = [&] {
        cartmesh::LinearOctree tree(
            cartmesh::AABB({-0.25, -0.25, -0.25}, {1.25, 1.25, 1.25}), 2, 5);
        cartmesh::OctreeRefinementConfiguration configuration;
        configuration.surface_target_level = 5;
        configuration.distance_bands.push_back({0.1, 4});
        const cartmesh::OctreeRefinementEngine engine(configuration, &bvh);
        static_cast<void>(engine.apply(tree));
        return std::pair{tree, cartmesh::classify_octree_leaves(tree, classifier)};
    };
    const auto [first_tree, first] = build_and_classify();
    const auto [second_tree, second] = build_and_classify();
    expect(first.outside_count + first.inside_count + first.intersected_count +
                   first.conflict_count ==
               first_tree.leaf_count(),
           "自适应叶的四类计数必须完整分割所有叶子");
    expect(first.conflict_count == 0 && first.intersected_count > 0 &&
               first.inside_count > 0 && first.outside_count > 0,
           "封闭立方体自适应分类必须有内外相交叶且无冲突");
    expect(first.inside_volume <= 1.0 && first.inside_plus_intersected_volume >= 1.0,
           "自适应叶体积上下界必须包含单位立方体真体积");
    expect(first_tree.result_hash_fnv1a64() == second_tree.result_hash_fnv1a64() &&
               first.cell_classification == second.cell_classification &&
               first.center_point_classification == second.center_point_classification,
           "相同自适应输入必须产生完全相同的叶顺序与分类字节");
}

} // 匿名命名空间

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"层级 Morton 节点码", test_level_aware_morton_code},
        {"初始分区与包围盒", test_initial_partition_and_bounds},
        {"细化与粗化", test_refine_and_coarsen},
        {"跨层面邻居", test_face_neighbors_across_levels},
        {"确定性细化与 2:1 平衡", test_deterministic_refinement_and_balance},
        {"BVH 最近表面距离", test_bvh_nearest_surface_distance},
        {"用户区域细化", test_user_region_refinement},
        {"表面、距离与曲率细化", test_surface_distance_and_curvature_refinement},
        {"狭缝保护细化", test_gap_protection_refinement},
        {"自适应叶几何分类", test_adaptive_leaf_classification},
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
    std::cout << "阶段2核心测试数=" << tests.size() << " 失败数=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
