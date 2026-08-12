#include "cartmesh/quality/SolverMeshQuality.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

[[nodiscard]] bool has_issue(const cartmesh::MeshQualityReport& report,
                             cartmesh::MeshQualityIssueKind kind) {
    return std::any_of(report.issues.begin(), report.issues.end(),
                       [&](const auto& issue) { return issue.kind == kind; });
}

[[nodiscard]] cartmesh::OpenFoamMesh unit_cube_mesh() {
    cartmesh::OpenFoamMesh mesh;
    mesh.background_stable_id_kind = "test_stable_id";
    mesh.points = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}};
    mesh.cells.push_back({7U, 42U, 0U, 0U, 1.0, true, 0U, {}});
    mesh.faces = {
        {{0, 4, 7, 3}, 0}, {{1, 2, 6, 5}, 0},
        {{0, 1, 5, 4}, 0}, {{3, 7, 6, 2}, 0},
        {{0, 3, 2, 1}, 0}, {{4, 5, 6, 7}, 0}};
    return mesh;
}

void test_valid_cube_metrics() {
    const auto report = cartmesh::evaluate_solver_mesh_quality(unit_cube_mesh());
    expect(report.topology_pass() && report.quality_pass(),
           "单位立方体不应产生质量问题");
    expect(std::abs(report.summary.minimum_cell_volume - 1.0) < 1.0e-14,
           "单位立方体体积必须为 1");
    expect(report.summary.maximum_cell_closure_ratio < 1.0e-14,
           "单位立方体面积向量必须闭合");
    expect(std::abs(report.summary.minimum_face_pyramid_volume - 1.0 / 6.0) <
               1.0e-14,
           "单位立方体 face pyramid 体积必须为 1/6");
}

void test_nonclosed_and_wrong_pyramid_are_located() {
    auto mesh = unit_cube_mesh();
    std::reverse(mesh.faces.back().point_ids.begin(),
                 mesh.faces.back().point_ids.end());
    const auto report = cartmesh::evaluate_solver_mesh_quality(mesh);
    expect(has_issue(report, cartmesh::MeshQualityIssueKind::cell_not_closed),
           "翻转一张面必须定位非闭合单元");
    expect(has_issue(report, cartmesh::MeshQualityIssueKind::wrong_face_pyramid),
           "翻转一张面必须定位负 face pyramid");
    expect(has_issue(report, cartmesh::MeshQualityIssueKind::non_star_shaped_cell),
           "负 face pyramid 必须汇总为非星形单元");
    const auto issue = std::find_if(
        report.issues.begin(), report.issues.end(), [](const auto& candidate) {
            return candidate.kind ==
                       cartmesh::MeshQualityIssueKind::wrong_face_pyramid &&
                   candidate.face_id == 5U;
        });
    expect(issue != report.issues.end() && issue->cell_id == 0U &&
               issue->face_id == 5U && issue->background_stable_id == 42U &&
               issue->source_type == "full_cartesian",
           "失败项必须包含 solver cell/face、稳定背景 ID 和来源类型");
}

void test_tiny_face_edge_duplicate_and_fraction() {
    auto mesh = unit_cube_mesh();
    mesh.cells.front().source_volume_fraction = 1.0e-4;
    mesh.faces.push_back({{0, 0, 0}, 0});
    mesh.faces.push_back(mesh.faces.front());
    const auto report = cartmesh::evaluate_solver_mesh_quality(mesh);
    expect(has_issue(report, cartmesh::MeshQualityIssueKind::zero_or_tiny_face),
           "零面积面必须可定位");
    expect(has_issue(report, cartmesh::MeshQualityIssueKind::zero_or_tiny_edge),
           "零长度边必须可定位");
    expect(has_issue(report,
                     cartmesh::MeshQualityIssueKind::baffle_like_duplicate),
           "重复 boundary face 必须识别为 baffle-like duplicate");
    expect(has_issue(report, cartmesh::MeshQualityIssueKind::tiny_volume_fraction),
           "极小来源体积分数必须可定位");
}

void test_concave_face_is_detected() {
    auto mesh = unit_cube_mesh();
    mesh.points.push_back({0.5, 0.25, 1.0});
    mesh.faces.back().point_ids = {4, 5, 6, 8, 7};
    const auto report = cartmesh::evaluate_solver_mesh_quality(mesh);
    expect(has_issue(report, cartmesh::MeshQualityIssueKind::concave_face),
           "含内凹顶点的面必须识别为 concave face");
}

void test_non_orthogonality_and_skewness_are_located() {
    auto mesh = unit_cube_mesh();
    mesh.points.insert(mesh.points.end(), {
        {2.0, 0.5, 0.0}, {2.0, 1.5, 0.0},
        {2.0, 0.5, 1.0}, {2.0, 1.5, 1.0}});
    mesh.cells.push_back({8U, 43U, 0U, 0U, 1.0, true, 0U, {}});
    mesh.faces = {
        {{1, 2, 6, 5}, 0, 1},
        {{0, 4, 7, 3}, 0}, {{0, 1, 5, 4}, 0},
        {{3, 7, 6, 2}, 0}, {{0, 3, 2, 1}, 0},
        {{4, 5, 6, 7}, 0},
        {{8, 9, 11, 10}, 1}, {{1, 8, 10, 5}, 1},
        {{2, 6, 11, 9}, 1}, {{1, 2, 9, 8}, 1},
        {{5, 10, 11, 6}, 1}};
    mesh.internal_face_count = 1U;
    cartmesh::MeshQualityThresholds thresholds;
    thresholds.maximum_non_orthogonality_degrees = 5.0;
    thresholds.maximum_internal_skewness = 0.05;
    const auto report =
        cartmesh::evaluate_solver_mesh_quality(mesh, thresholds);
    expect(has_issue(
               report,
               cartmesh::MeshQualityIssueKind::excessive_non_orthogonality),
           "斜置相邻单元必须定位内部面非正交");
    expect(has_issue(report,
                     cartmesh::MeshQualityIssueKind::excessive_skewness),
           "斜置相邻单元必须定位内部面偏斜");
    expect(report.summary.maximum_non_orthogonality_degrees > 5.0 &&
               report.summary.maximum_skewness > 0.05,
           "非正交和偏斜 worst value 必须进入汇总");
}

void test_report_is_byte_deterministic() {
    auto mesh = unit_cube_mesh();
    std::reverse(mesh.faces.back().point_ids.begin(),
                 mesh.faces.back().point_ids.end());
    const auto report = cartmesh::evaluate_solver_mesh_quality(mesh);
    const auto first = std::filesystem::temp_directory_path() /
                       "cartmesh-stage62-quality-first.json";
    const auto second = std::filesystem::temp_directory_path() /
                        "cartmesh-stage62-quality-second.json";
    cartmesh::write_solver_mesh_quality_json(first, report);
    cartmesh::write_solver_mesh_quality_json(second, report);
    const auto read = [](const std::filesystem::path& path) {
        std::ifstream input(path);
        return std::string(std::istreambuf_iterator<char>(input),
                           std::istreambuf_iterator<char>());
    };
    const auto first_text = read(first);
    const auto second_text = read(second);
    std::filesystem::remove(first);
    std::filesystem::remove(second);
    expect(first_text == second_text &&
               first_text.find("\"schema\": \"cartmesh-solver-mesh-quality-v1\"") !=
                   std::string::npos &&
               first_text.find("\"backgroundStableId\":\"42\"") !=
                   std::string::npos,
           "质量报告必须逐字节确定且保存可追溯稳定 ID");
}

void test_nonfinite_diagnostic_json_stays_valid() {
    auto mesh = unit_cube_mesh();
    mesh.points.front().x = std::numeric_limits<double>::quiet_NaN();
    const auto report = cartmesh::evaluate_solver_mesh_quality(mesh);
    expect(has_issue(report,
                     cartmesh::MeshQualityIssueKind::non_finite_geometry),
           "非有限坐标必须显式报告");
    const auto path = std::filesystem::temp_directory_path() /
                      "cartmesh-stage62-nonfinite-quality.json";
    cartmesh::write_solver_mesh_quality_json(path, report);
    std::ifstream input(path);
    const std::string text{std::istreambuf_iterator<char>(input),
                           std::istreambuf_iterator<char>()};
    std::filesystem::remove(path);
    expect(text.find(":nan") == std::string::npos &&
               text.find(":inf") == std::string::npos &&
               text.find("null") != std::string::npos,
           "非有限诊断值必须序列化为 JSON null，不能输出 nan/inf");
}

} // 匿名命名空间

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"有效单位立方体", test_valid_cube_metrics},
        {"非闭合与负 face pyramid 定位",
         test_nonclosed_and_wrong_pyramid_are_located},
        {"微小面边重复面与体积分数",
         test_tiny_face_edge_duplicate_and_fraction},
        {"凹面检测", test_concave_face_is_detected},
        {"非正交与偏斜定位",
         test_non_orthogonality_and_skewness_are_located},
        {"报告确定性", test_report_is_byte_deterministic},
        {"非有限诊断 JSON", test_nonfinite_diagnostic_json_stays_valid},
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
    std::cout << "Stage 6.2 原生质量测试数=" << tests.size()
              << " 失败数=" << failures << '\n';
    return failures == 0U ? 0 : 1;
}
