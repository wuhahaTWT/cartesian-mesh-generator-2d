#include "cartmesh/quality/SolverMeshStabilizer.hpp"
#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

struct MeshBuilder {
    cartmesh::OpenFoamMesh mesh;
    std::map<std::tuple<double, double, double>, std::size_t> point_ids;
    std::vector<cartmesh::OpenFoamFace> internal;
    std::vector<cartmesh::OpenFoamFace> boundary;

    std::size_t point(double x, double y, double z) {
        const auto key = std::tuple(x, y, z);
        const auto found = point_ids.find(key);
        if (found != point_ids.end()) return found->second;
        const std::size_t id = mesh.points.size();
        mesh.points.push_back({x, y, z});
        point_ids.emplace(key, id);
        return id;
    }

    void face(std::initializer_list<std::tuple<double, double, double>> vertices,
              std::size_t owner,
              std::size_t neighbour = std::numeric_limits<std::size_t>::max(),
              std::uint64_t boundary_id = 0U, bool farfield = true) {
        cartmesh::OpenFoamFace face;
        for (const auto& [x, y, z] : vertices) {
            face.point_ids.push_back(point(x, y, z));
        }
        face.owner = owner;
        face.neighbour = neighbour;
        face.boundary_id = boundary_id;
        face.farfield = farfield;
        (face.internal() ? internal : boundary).push_back(std::move(face));
    }

    cartmesh::OpenFoamMesh finish() {
        mesh.background_stable_id_kind = "stage63_test";
        mesh.internal_face_count = internal.size();
        mesh.faces = std::move(internal);
        mesh.faces.insert(mesh.faces.end(),
                          std::make_move_iterator(boundary.begin()),
                          std::make_move_iterator(boundary.end()));
        return std::move(mesh);
    }
};

cartmesh::OpenFoamCellSource source(std::uint64_t id, double fraction,
                                    double background_volume,
                                    std::uint64_t region = 1U) {
    cartmesh::OpenFoamCellSource result{
        id, id, 0U, 0U, fraction, false, region, {}};
    result.members.push_back({id, id, 0U, 0U, region, background_volume});
    return result;
}

[[nodiscard]] cartmesh::OpenFoamMesh two_boxes(double split = 1.0e-3,
                                               double tiny_fraction = 1.0e-3) {
    MeshBuilder builder;
    builder.mesh.cells = {source(101U, tiny_fraction, 1.0),
                          source(102U, 1.0 - split, 1.0)};
    builder.face({{split,0,0},{split,1,0},{split,1,1},{split,0,1}}, 0, 1);
    builder.face({{0,0,0},{0,0,1},{0,1,1},{0,1,0}}, 0);
    builder.face({{0,0,0},{split,0,0},{split,0,1},{0,0,1}}, 0);
    builder.face({{0,1,0},{0,1,1},{split,1,1},{split,1,0}}, 0);
    builder.face({{0,0,0},{0,1,0},{split,1,0},{split,0,0}}, 0);
    builder.face({{0,0,1},{split,0,1},{split,1,1},{0,1,1}}, 0);
    builder.face({{1,0,0},{1,1,0},{1,1,1},{1,0,1}}, 1);
    builder.face({{split,0,0},{1,0,0},{1,0,1},{split,0,1}}, 1);
    builder.face({{split,1,0},{split,1,1},{1,1,1},{1,1,0}}, 1);
    builder.face({{split,0,0},{split,1,0},{1,1,0},{1,0,0}}, 1);
    builder.face({{split,0,1},{1,0,1},{1,1,1},{split,1,1}}, 1);
    return builder.finish();
}

[[nodiscard]] cartmesh::OpenFoamMesh t_junction_mesh() {
    MeshBuilder builder;
    builder.mesh.cells = {source(201U, 1.0, 1.0),
                          source(202U, 1.0e-4, 1.0),
                          source(203U, 1.0, 0.5)};
    builder.face({{1,0,0},{1,.5,0},{1,.5,1},{1,0,1}}, 0, 1);
    builder.face({{1,.5,0},{1,1,0},{1,1,1},{1,.5,1}}, 0, 2);
    builder.face({{1,.5,0},{1,.5,1},{2,.5,1},{2,.5,0}}, 1, 2);
    builder.face({{0,0,0},{0,0,1},{0,1,1},{0,1,0}}, 0);
    builder.face({{0,0,0},{1,0,0},{1,0,1},{0,0,1}}, 0);
    builder.face({{0,1,0},{0,1,1},{1,1,1},{1,1,0}}, 0);
    builder.face({{0,0,0},{0,1,0},{1,1,0},{1,.5,0},{1,0,0}}, 0);
    builder.face({{0,0,1},{1,0,1},{1,.5,1},{1,1,1},{0,1,1}}, 0);
    builder.face({{2,0,0},{2,.5,0},{2,.5,1},{2,0,1}}, 1);
    builder.face({{1,0,0},{2,0,0},{2,0,1},{1,0,1}}, 1);
    builder.face({{1,0,0},{1,.5,0},{2,.5,0},{2,0,0}}, 1);
    builder.face({{1,0,1},{2,0,1},{2,.5,1},{1,.5,1}}, 1);
    builder.face({{2,.5,0},{2,1,0},{2,1,1},{2,.5,1}}, 2);
    builder.face({{1,1,0},{1,1,1},{2,1,1},{2,1,0}}, 2);
    builder.face({{1,.5,0},{1,1,0},{2,1,0},{2,.5,0}}, 2);
    builder.face({{1,.5,1},{2,.5,1},{2,1,1},{1,1,1}}, 2);
    return builder.finish();
}

[[nodiscard]] std::vector<cartmesh::Triangle> box_triangles(
    const cartmesh::Vec3& minimum, const cartmesh::Vec3& maximum) {
    const std::array<cartmesh::Vec3, 8> point = {{
        {minimum.x, minimum.y, minimum.z},
        {maximum.x, minimum.y, minimum.z},
        {maximum.x, maximum.y, minimum.z},
        {minimum.x, maximum.y, minimum.z},
        {minimum.x, minimum.y, maximum.z},
        {maximum.x, minimum.y, maximum.z},
        {maximum.x, maximum.y, maximum.z},
        {minimum.x, maximum.y, maximum.z}}};
    const std::array<std::array<std::size_t, 3>, 12> face = {{
        {0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},
        {3,7,6},{3,6,2},{0,4,7},{0,7,3},{1,2,6},{1,6,5}}};
    std::vector<cartmesh::Triangle> result;
    for (const auto& indices : face) {
        result.emplace_back(point[indices[0]], point[indices[1]],
                            point[indices[2]]);
    }
    return result;
}

void test_tiny_sliver_agglomerates_conservatively() {
    const auto before = two_boxes();
    const auto before_quality = cartmesh::evaluate_solver_mesh_quality(before);
    const auto result = cartmesh::stabilize_solver_mesh(before);
    const auto after_quality = cartmesh::evaluate_solver_mesh_quality(result.mesh);
    expect(result.report.pass() && result.report.agglomeration_count == 1U &&
               result.mesh.cells.size() == 1U,
           "tiny sliver 必须通过一次保守聚并消除");
    expect(std::abs(result.report.initial_volume - result.report.final_volume) <
               1.0e-14 &&
               cartmesh::norm(result.report.initial_first_moment -
                              result.report.final_first_moment) < 1.0e-14,
           "聚并前后体积和一阶矩必须保守");
    expect(before_quality.topology_pass() && after_quality.quality_pass(),
           "tiny sliver 聚并不得破坏拓扑或几何质量");
    expect(result.mesh.cells.front().members.size() == 2U,
           "聚并控制体必须保留两个来源记录");
}

void test_t_junction_survives_agglomeration() {
    const auto result = cartmesh::stabilize_solver_mesh(t_junction_mesh());
    const auto quality = cartmesh::evaluate_solver_mesh_quality(result.mesh);
    expect(result.report.pass() && result.report.agglomeration_count == 1U &&
               result.mesh.cells.size() == 2U,
           "T-junction 上的小单元必须以真实邻接聚并");
    expect(quality.topology_pass() &&
               quality.summary.minimum_face_pyramid_volume > 0.0,
           "T-junction 聚并后必须保持 owner/neighbour、闭合和正棱锥");
}

void test_extreme_volume_fraction_is_not_deleted() {
    const auto result = cartmesh::stabilize_solver_mesh(two_boxes(1.0e-6, 1.0e-12));
    expect(result.report.pass() && result.mesh.cells.size() == 1U &&
               result.mesh.cells.front().members.front().background_stable_id == 101U,
           "极端体积分数必须保留来源并聚并，不得删除单元");
}

void test_coarse_fine_cut_interface_stays_conformal() {
    const cartmesh::SurfaceMesh surface(
        box_triangles({0.2,0.2,0.2}, {0.65,0.65,0.65}),
        cartmesh::SurfaceFormat::ascii_stl, "stage63_coarse_fine_cut");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 9U);
    cartmesh::LinearOctree tree(
        cartmesh::AABB({0,0,0}, {1,1,1}), 1, 3);
    expect(tree.refine_leaf(cartmesh::encode_octree_node(1,0,0,0)),
           "coarse-fine+cut 测例必须建立一个局部细化叶");
    static_cast<void>(tree.balance_faces_2_to_1());
    const auto cut = cartmesh::build_triangulated_cut_cell_mesh(tree, cutter);
    const auto raw = cartmesh::build_openfoam_mesh(tree, cut);
    expect(cartmesh::evaluate_solver_mesh_quality(raw).quality_pass(),
           "coarse-fine+cut 基线必须先通过原生质量评估");
    bool exercised = false;
    for (std::size_t cell = 0; cell < raw.cells.size() && !exercised; ++cell) {
        if (raw.cells[cell].full_cartesian) continue;
        const bool has_neighbour = std::any_of(
            raw.faces.begin(), raw.faces.begin() +
                static_cast<std::ptrdiff_t>(raw.internal_face_count),
            [&](const auto& face) {
                return face.owner == cell || face.neighbour == cell;
            });
        if (!has_neighbour) continue;
        auto triggered = raw;
        triggered.cells[cell].source_volume_fraction = 1.0e-6;
        const auto stabilized = cartmesh::stabilize_solver_mesh(triggered);
        if (!stabilized.report.pass() ||
            stabilized.report.agglomeration_count == 0U) {
            continue;
        }
        const auto quality =
            cartmesh::evaluate_solver_mesh_quality(stabilized.mesh);
        expect(quality.quality_pass() &&
                   stabilized.report.initial_cell_count ==
                       stabilized.report.final_cell_count + 1U,
               "coarse-fine+cut 聚并必须保持共形拓扑和正质量");
        exercised = true;
    }
    expect(exercised,
           "coarse-fine+cut 测例必须实际覆盖一次安全聚并");
}

void test_non_star_candidate_is_rejected() {
    auto mesh = two_boxes();
    std::reverse(mesh.faces[1].point_ids.begin(), mesh.faces[1].point_ids.end());
    const auto result = cartmesh::stabilize_solver_mesh(mesh);
    const bool rejected = std::any_of(
        result.report.actions.begin(), result.report.actions.end(),
        [](const auto& action) {
            return action.kind ==
                   cartmesh::StabilizationActionKind::rejected_non_star;
        });
    expect(rejected && !result.report.pass() && result.mesh.cells.size() == 2U &&
               result.report.refinement_requested_stable_ids ==
                   std::vector<std::uint64_t>{101U},
           "合并后非星形的候选必须拒绝并请求来源细化");
}

void test_adaptive_refinement_and_max_level_stop() {
    cartmesh::LinearOctree tree(
        cartmesh::AABB({0,0,0}, {1,1,1}), 1, 3);
    const auto target = cartmesh::encode_octree_node(1, 0, 0, 0);
    const auto refined = cartmesh::refine_stabilization_sources(tree, {target});
    expect(refined.pass() && refined.refined_leaf_count == 1U &&
               !tree.find_leaf(target).has_value() &&
               tree.validate_partition() && tree.check_face_balance().balanced,
           "局部回退必须细分原 Morton 叶并恢复 2:1 平衡");
    cartmesh::LinearOctree maxed(
        cartmesh::AABB({0,0,0}, {1,1,1}), 1, 1);
    const auto stopped =
        cartmesh::refine_stabilization_sources(maxed, {target});
    expect(!stopped.pass() && stopped.refined_leaf_count == 0U &&
               stopped.unresolved_stable_ids ==
                   std::vector<std::uint64_t>{target},
           "达到最大层级时必须显式停止，不得假装稳定");
    cartmesh::LinearOctree mixed(
        cartmesh::AABB({0,0,0}, {1,1,1}), 1, 2);
    const auto refinable = cartmesh::encode_octree_node(1, 1, 0, 0);
    expect(mixed.refine_leaf(target),
           "混合原子性测例必须先构造 max-level 叶");
    const auto max_level_child = cartmesh::octree_child(target, 0U);
    const auto before_count = mixed.leaf_count();
    const auto mixed_stop = cartmesh::refine_stabilization_sources(
        mixed, {refinable, max_level_child});
    expect(!mixed_stop.pass() && mixed.leaf_count() == before_count &&
               mixed.find_leaf(refinable).has_value(),
           "局部细化请求必须全部可执行才原子修改 tree");
}

void test_result_and_report_are_deterministic() {
    const auto first = cartmesh::stabilize_solver_mesh(t_junction_mesh());
    const auto second = cartmesh::stabilize_solver_mesh(t_junction_mesh());
    const auto first_path = std::filesystem::temp_directory_path() /
                            "cartmesh-stage63-first.json";
    const auto second_path = std::filesystem::temp_directory_path() /
                             "cartmesh-stage63-second.json";
    cartmesh::write_solver_mesh_stabilization_json(first_path, first.report);
    cartmesh::write_solver_mesh_stabilization_json(second_path, second.report);
    const auto read = [](const std::filesystem::path& path) {
        std::ifstream input(path);
        return std::string(std::istreambuf_iterator<char>(input),
                           std::istreambuf_iterator<char>());
    };
    const auto first_text = read(first_path);
    const auto second_text = read(second_path);
    std::filesystem::remove(first_path);
    std::filesystem::remove(second_path);
    const bool same_faces = first.mesh.faces.size() == second.mesh.faces.size() &&
        std::equal(first.mesh.faces.begin(), first.mesh.faces.end(),
                   second.mesh.faces.begin(), [](const auto& lhs, const auto& rhs) {
            return lhs.point_ids == rhs.point_ids && lhs.owner == rhs.owner &&
                   lhs.neighbour == rhs.neighbour &&
                   lhs.boundary_id == rhs.boundary_id &&
                   lhs.farfield == rhs.farfield;
        });
    expect(first_text == second_text && same_faces,
           "聚并选择、拓扑和 JSON 报告必须逐字节确定");
}

} // namespace

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"tiny sliver 保守聚并", test_tiny_sliver_agglomerates_conservatively},
        {"T-junction 拓扑", test_t_junction_survives_agglomeration},
        {"极端体积分数", test_extreme_volume_fraction_is_not_deleted},
        {"coarse-fine+cut 共形性",
         test_coarse_fine_cut_interface_stays_conformal},
        {"非星形候选拒绝", test_non_star_candidate_is_rejected},
        {"局部细化与最大层级停止",
         test_adaptive_refinement_and_max_level_stop},
        {"确定性", test_result_and_report_are_deterministic},
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
    std::cout << "Stage 6.3 稳定化测试数=" << tests.size()
              << " 失败数=" << failures << '\n';
    return failures == 0U ? 0 : 1;
}
