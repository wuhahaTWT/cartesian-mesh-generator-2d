#include "cartmesh/classify/SurfaceClassifier.hpp"
#include "cartmesh/geometry/SurfaceDiagnostics.hpp"
#include "cartmesh/geometry/TriangleBoxIntersection.hpp"
#include "cartmesh/io/StlReader.hpp"
#include "cartmesh/io/DiagnosticVtkWriter.hpp"
#include "cartmesh/spatial/TriangleBvh.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <numbers>
#include <sstream>
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

void expect_near(double actual, double expected, double tolerance, const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        std::ostringstream detail;
        detail << message << "：实际值=" << actual << " 期望值=" << expected
               << " 容差=" << tolerance;
        throw TestFailure(detail.str());
    }
}

template <typename Function> void expect_throw(Function&& function, const std::string& message) {
    try {
        function();
    } catch (const std::exception&) {
        return;
    }
    throw TestFailure(message);
}

[[nodiscard]] std::uint64_t classification_hash(
    const cartmesh::UniformCartesianGrid& grid,
    const cartmesh::UniformClassification& classification) noexcept {
    std::uint64_t hash = 14695981039346656037ULL;
    const auto hash_byte = [&](std::uint8_t value) {
        hash ^= value;
        hash *= 1099511628211ULL;
    };
    const auto hash_u64 = [&](std::uint64_t value) {
        for (std::uint32_t shift = 0; shift < 64; shift += 8) {
            hash_byte(static_cast<std::uint8_t>((value >> shift) & 0xffULL));
        }
    };
    hash_u64(grid.nx());
    hash_u64(grid.ny());
    hash_u64(grid.nz());
    for (std::size_t index = 0; index < classification.cell_classification.size(); ++index) {
        hash_byte(classification.cell_classification[index]);
        hash_byte(classification.center_point_classification[index]);
    }
    return hash;
}

[[nodiscard]] std::vector<cartmesh::Triangle> cube_triangles() {
    const std::array<cartmesh::Vec3, 8> point = {
        cartmesh::Vec3{0.0, 0.0, 0.0}, cartmesh::Vec3{1.0, 0.0, 0.0},
        cartmesh::Vec3{1.0, 1.0, 0.0}, cartmesh::Vec3{0.0, 1.0, 0.0},
        cartmesh::Vec3{0.0, 0.0, 1.0}, cartmesh::Vec3{1.0, 0.0, 1.0},
        cartmesh::Vec3{1.0, 1.0, 1.0}, cartmesh::Vec3{0.0, 1.0, 1.0}};
    const std::array<std::array<std::size_t, 3>, 12> faces = {
        std::array<std::size_t, 3>{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
        {0, 1, 5}, {0, 5, 4}, {3, 7, 6}, {3, 6, 2}, {0, 4, 7}, {0, 7, 3},
        {1, 2, 6}, {1, 6, 5}};
    std::vector<cartmesh::Triangle> triangles;
    triangles.reserve(faces.size());
    for (const auto& face : faces) {
        triangles.emplace_back(point[face[0]], point[face[1]], point[face[2]]);
    }
    return triangles;
}

[[nodiscard]] cartmesh::SurfaceMesh cube_surface() {
    return cartmesh::SurfaceMesh(cube_triangles(), cartmesh::SurfaceFormat::ascii_stl,
                                 "unit_cube");
}

[[nodiscard]] std::vector<cartmesh::Triangle> transformed_cube_triangles(
    const cartmesh::Vec3& shift, double scale, bool reverse = false) {
    std::vector<cartmesh::Triangle> result;
    result.reserve(12);
    for (const auto& triangle : cube_triangles()) {
        const auto& vertex = triangle.vertices();
        const auto transform = [&](const cartmesh::Vec3& point) {
            return shift + point * scale;
        };
        result.emplace_back(transform(vertex[0]), transform(reverse ? vertex[2] : vertex[1]),
                            transform(reverse ? vertex[1] : vertex[2]));
    }
    return result;
}

[[nodiscard]] std::vector<cartmesh::Triangle> tetrahedron_triangles(
    const std::array<cartmesh::Vec3, 4>& point) {
    const std::array<std::array<std::size_t, 3>, 4> faces = {
        std::array<std::size_t, 3>{0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2}};
    const auto center = (point[0] + point[1] + point[2] + point[3]) / 4.0;
    std::vector<cartmesh::Triangle> triangles;
    for (const auto& face : faces) {
        cartmesh::Triangle triangle(point[face[0]], point[face[1]], point[face[2]]);
        if (cartmesh::dot(triangle.area_vector(), triangle.centroid() - center) < 0.0) {
            triangle = cartmesh::Triangle(point[face[0]], point[face[2]], point[face[1]]);
        }
        triangles.push_back(triangle);
    }
    return triangles;
}

void write_ascii_stl(const std::filesystem::path& path,
                     const std::vector<cartmesh::Triangle>& triangles) {
    std::ofstream output(path, std::ios::trunc);
    output.precision(17);
    output << "solid unit_cube\n";
    for (const auto& triangle : triangles) {
        const auto normal = triangle.area_vector();
        output << "  facet normal " << normal.x << ' ' << normal.y << ' ' << normal.z << "\n"
               << "    outer loop\n";
        for (const auto& vertex : triangle.vertices()) {
            output << "      vertex " << vertex.x << ' ' << vertex.y << ' ' << vertex.z << "\n";
        }
        output << "    endloop\n  endfacet\n";
    }
    output << "endsolid unit_cube\n";
}

[[nodiscard]] std::vector<cartmesh::Triangle> sphere_triangles(std::uint32_t latitude_count,
                                                              std::uint32_t longitude_count) {
    const cartmesh::Vec3 top{0.0, 0.0, 1.0};
    const cartmesh::Vec3 bottom{0.0, 0.0, -1.0};
    std::vector<std::vector<cartmesh::Vec3>> rings;
    for (std::uint32_t latitude = 1; latitude < latitude_count; ++latitude) {
        const double phi = std::numbers::pi_v<double> * static_cast<double>(latitude) /
                           static_cast<double>(latitude_count);
        std::vector<cartmesh::Vec3> ring;
        ring.reserve(longitude_count);
        for (std::uint32_t longitude = 0; longitude < longitude_count; ++longitude) {
            const double theta = 2.0 * std::numbers::pi_v<double> *
                                 static_cast<double>(longitude) /
                                 static_cast<double>(longitude_count);
            ring.push_back(
                {std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi)});
        }
        rings.push_back(std::move(ring));
    }
    std::vector<cartmesh::Triangle> triangles;
    triangles.reserve(2U * longitude_count * (latitude_count - 1U));
    for (std::uint32_t longitude = 0; longitude < longitude_count; ++longitude) {
        const auto next = (longitude + 1U) % longitude_count;
        triangles.emplace_back(top, rings.front()[longitude], rings.front()[next]);
        for (std::size_t ring = 0; ring + 1 < rings.size(); ++ring) {
            triangles.emplace_back(rings[ring][longitude], rings[ring + 1][longitude],
                                   rings[ring + 1][next]);
            triangles.emplace_back(rings[ring][longitude], rings[ring + 1][next],
                                   rings[ring][next]);
        }
        triangles.emplace_back(bottom, rings.back()[next], rings.back()[longitude]);
    }
    return triangles;
}

[[nodiscard]] std::vector<cartmesh::Triangle> cylinder_triangles(std::uint32_t segments) {
    const cartmesh::Vec3 bottom_center{0.0, 0.0, -1.0};
    const cartmesh::Vec3 top_center{0.0, 0.0, 1.0};
    std::vector<cartmesh::Vec3> bottom;
    std::vector<cartmesh::Vec3> top;
    for (std::uint32_t segment = 0; segment < segments; ++segment) {
        const double theta = 2.0 * std::numbers::pi_v<double> * static_cast<double>(segment) /
                             static_cast<double>(segments);
        bottom.push_back({std::cos(theta), std::sin(theta), -1.0});
        top.push_back({std::cos(theta), std::sin(theta), 1.0});
    }
    std::vector<cartmesh::Triangle> triangles;
    triangles.reserve(static_cast<std::size_t>(segments) * 4U);
    for (std::uint32_t segment = 0; segment < segments; ++segment) {
        const auto next = (segment + 1U) % segments;
        triangles.emplace_back(bottom_center, bottom[next], bottom[segment]);
        triangles.emplace_back(top_center, top[segment], top[next]);
        triangles.emplace_back(bottom[segment], bottom[next], top[next]);
        triangles.emplace_back(bottom[segment], top[next], top[segment]);
    }
    return triangles;
}

[[nodiscard]] std::vector<cartmesh::Triangle> tube_triangles(std::uint32_t segments,
                                                            double inner_radius) {
    std::vector<cartmesh::Vec3> outer_bottom;
    std::vector<cartmesh::Vec3> outer_top;
    std::vector<cartmesh::Vec3> inner_bottom;
    std::vector<cartmesh::Vec3> inner_top;
    for (std::uint32_t segment = 0; segment < segments; ++segment) {
        const double theta = 2.0 * std::numbers::pi_v<double> * static_cast<double>(segment) /
                             static_cast<double>(segments);
        const double cosine = std::cos(theta);
        const double sine = std::sin(theta);
        outer_bottom.push_back({cosine, sine, -1.0});
        outer_top.push_back({cosine, sine, 1.0});
        inner_bottom.push_back({inner_radius * cosine, inner_radius * sine, -1.0});
        inner_top.push_back({inner_radius * cosine, inner_radius * sine, 1.0});
    }
    std::vector<cartmesh::Triangle> triangles;
    triangles.reserve(static_cast<std::size_t>(segments) * 8U);
    for (std::uint32_t segment = 0; segment < segments; ++segment) {
        const auto next = (segment + 1U) % segments;
        triangles.emplace_back(outer_bottom[segment], outer_bottom[next], outer_top[next]);
        triangles.emplace_back(outer_bottom[segment], outer_top[next], outer_top[segment]);
        triangles.emplace_back(inner_bottom[segment], inner_top[next], inner_bottom[next]);
        triangles.emplace_back(inner_bottom[segment], inner_top[segment], inner_top[next]);
        triangles.emplace_back(outer_top[segment], outer_top[next], inner_top[next]);
        triangles.emplace_back(outer_top[segment], inner_top[next], inner_top[segment]);
        triangles.emplace_back(outer_bottom[segment], inner_bottom[next], outer_bottom[next]);
        triangles.emplace_back(outer_bottom[segment], inner_bottom[segment], inner_bottom[next]);
    }
    return triangles;
}

struct ClassificationMetrics {
    double center_sample_volume{};
    double definitely_inside_volume{};
    double inside_plus_intersected_volume{};
    std::uint64_t intersected_count{};
};

[[nodiscard]] ClassificationMetrics classified_metrics(const cartmesh::SurfaceMesh& surface,
                                                        std::uint32_t resolution) {
    const auto diagnostics = cartmesh::diagnose_surface(surface);
    expect(diagnostics.valid_for_stage1_classification(), "收敛案例表面必须通过诊断");
    const cartmesh::TriangleBvh bvh(surface);
    const cartmesh::SurfaceClassifier classifier(bvh);
    const auto bounds = surface.bounds();
    const auto extent = bounds.extent();
    const double padding = std::max({extent.x, extent.y, extent.z}) * 0.05;
    const cartmesh::Vec3 delta{padding, padding, padding};
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB(bounds.minimum() - delta, bounds.maximum() + delta), resolution,
        resolution, resolution);
    const auto classification = cartmesh::classify_uniform_cells(grid, classifier);
    expect(classification.conflict_count == 0, "收敛案例不应存在分类冲突");
    expect(classification.center_conflict_count == 0, "收敛案例中心射线不应存在冲突");
    expect(classification.center_on_surface_count == 0,
           "收敛案例网格中心不应恰好落在表面");
    return {
        static_cast<double>(classification.center_inside_count) * grid.cell_volume(),
        static_cast<double>(classification.inside_count) * grid.cell_volume(),
        static_cast<double>(classification.inside_count + classification.intersected_count) *
            grid.cell_volume(),
        classification.intersected_count,
    };
}

void write_little_u32(std::ofstream& output, std::uint32_t value) {
    const std::array<char, 4> bytes = {static_cast<char>(value & 0xffU),
                                       static_cast<char>((value >> 8U) & 0xffU),
                                       static_cast<char>((value >> 16U) & 0xffU),
                                       static_cast<char>((value >> 24U) & 0xffU)};
    output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
}

void write_little_f32(std::ofstream& output, float value) {
    write_little_u32(output, std::bit_cast<std::uint32_t>(value));
}

void write_binary_stl(const std::filesystem::path& path,
                      const std::vector<cartmesh::Triangle>& triangles) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    std::array<char, 80> header{};
    const std::string text = "solid binary_header_must_not_force_ascii";
    std::copy(text.begin(), text.end(), header.begin());
    output.write(header.data(), static_cast<std::streamsize>(header.size()));
    write_little_u32(output, static_cast<std::uint32_t>(triangles.size()));
    for (const auto& triangle : triangles) {
        const auto normal = triangle.area_vector();
        write_little_f32(output, static_cast<float>(normal.x));
        write_little_f32(output, static_cast<float>(normal.y));
        write_little_f32(output, static_cast<float>(normal.z));
        for (const auto& vertex : triangle.vertices()) {
            write_little_f32(output, static_cast<float>(vertex.x));
            write_little_f32(output, static_cast<float>(vertex.y));
            write_little_f32(output, static_cast<float>(vertex.z));
        }
        const std::array<char, 2> attribute{};
        output.write(attribute.data(), static_cast<std::streamsize>(attribute.size()));
    }
}

void test_closed_oriented_manifold_diagnostics() {
    const auto diagnostics = cartmesh::diagnose_surface(cube_surface());
    expect(diagnostics.triangle_count == 12, "立方体三角形数量");
    expect(diagnostics.unique_vertex_count == 8, "立方体唯一顶点数量");
    expect(diagnostics.unique_edge_count == 18, "立方体含面对角线的唯一边数量");
    expect(diagnostics.boundary_edge_count == 0, "封闭立方体没有边界边");
    expect(diagnostics.non_manifold_edge_count == 0, "封闭立方体没有非流形边");
    expect(diagnostics.orientation_conflict_edge_count == 0, "立方体共享边方向一致");
    expect(diagnostics.connected_component_count == 1, "立方体只有一个连通分量");
    expect(diagnostics.valid_for_stage1_classification(), "封闭定向流形立方体可以分类");
    expect_near(diagnostics.signed_volume, 1.0, 1.0e-14, "立方体有向体积");
}

void test_invalid_surface_diagnostics_remain_visible() {
    auto open = cube_triangles();
    open.pop_back();
    const auto open_diagnostics = cartmesh::diagnose_surface(
        cartmesh::SurfaceMesh(open, cartmesh::SurfaceFormat::ascii_stl, "open_cube"));
    expect(open_diagnostics.boundary_edge_count == 3, "缺一个三角形暴露三条边界边");
    expect(!open_diagnostics.closed, "缺面立方体不得标记为封闭");
    expect(!open_diagnostics.boundary_edge_examples.empty(), "边界边必须保留位置证据");

    auto flipped = cube_triangles();
    const auto vertices = flipped[0].vertices();
    flipped[0] = cartmesh::Triangle(vertices[0], vertices[2], vertices[1]);
    const auto flipped_diagnostics = cartmesh::diagnose_surface(
        cartmesh::SurfaceMesh(flipped, cartmesh::SurfaceFormat::ascii_stl, "flipped_face"));
    expect(flipped_diagnostics.orientation_conflict_edge_count == 3,
           "翻转一个三角形必须暴露三条方向冲突边");
    expect(!flipped_diagnostics.consistently_oriented, "方向冲突不得隐藏");

    auto duplicate = cube_triangles();
    duplicate.push_back(duplicate.front());
    const auto duplicate_diagnostics = cartmesh::diagnose_surface(
        cartmesh::SurfaceMesh(duplicate, cartmesh::SurfaceFormat::ascii_stl, "duplicate_face"));
    expect(duplicate_diagnostics.duplicate_triangle_count == 1, "重复三角形必须被检测");
    expect(duplicate_diagnostics.non_manifold_edge_count == 3,
           "重复三角形使其三条边成为非流形边");
    expect(!duplicate_diagnostics.valid_for_stage1_classification(),
           "重复面输入不得进入阶段 1 分类");
}

void test_non_manifold_vertex_diagnostic() {
    auto triangles = tetrahedron_triangles(
        {cartmesh::Vec3{0.0, 0.0, 0.0}, cartmesh::Vec3{1.0, 0.0, 0.0},
         cartmesh::Vec3{0.0, 1.0, 0.0}, cartmesh::Vec3{0.0, 0.0, 1.0}});
    auto second = tetrahedron_triangles(
        {cartmesh::Vec3{0.0, 0.0, 0.0}, cartmesh::Vec3{-1.0, 0.0, 0.0},
         cartmesh::Vec3{0.0, -1.0, 0.0}, cartmesh::Vec3{0.0, 0.0, -1.0}});
    triangles.insert(triangles.end(), second.begin(), second.end());
    const auto diagnostics = cartmesh::diagnose_surface(cartmesh::SurfaceMesh(
        triangles, cartmesh::SurfaceFormat::ascii_stl, "pinched_vertex"));
    expect(diagnostics.boundary_edge_count == 0,
           "只共享顶点的两个闭合四面体不产生边界边");
    expect(diagnostics.non_manifold_edge_count == 0,
           "只共享顶点的两个闭合四面体不产生非流形边");
    expect(diagnostics.non_manifold_vertex_count == 1,
           "顶点 link 由两个不相连圆环组成时必须标记非流形顶点");
    expect(!diagnostics.manifold && !diagnostics.closed,
           "非流形顶点不得被封闭边检查掩盖");
    expect(!diagnostics.non_manifold_vertex_examples.empty(),
           "非流形顶点必须保留位置证据");
}

void test_ascii_and_binary_stl_import() {
    const auto directory = std::filesystem::temp_directory_path();
    const auto ascii_path = directory / "cartmesh_stage1_cube_ascii.stl";
    const auto binary_path = directory / "cartmesh_stage1_cube_binary_solid_header.stl";
    write_ascii_stl(ascii_path, cube_triangles());
    write_binary_stl(binary_path, cube_triangles());

    const auto ascii = cartmesh::read_stl(ascii_path);
    const auto binary = cartmesh::read_stl(binary_path);
    expect(ascii.format() == cartmesh::SurfaceFormat::ascii_stl, "ASCII STL 格式识别");
    expect(binary.format() == cartmesh::SurfaceFormat::binary_stl,
           "以 solid 开头的二进制 STL 仍按长度识别");
    expect(ascii.triangles().size() == 12 && binary.triangles().size() == 12,
           "ASCII/二进制 STL 三角形数量一致");
    expect_near(cartmesh::diagnose_surface(ascii).signed_volume,
                cartmesh::diagnose_surface(binary).signed_volume, 1.0e-14,
                "ASCII/二进制导入的有向体积一致");

    std::filesystem::remove(ascii_path);
    std::filesystem::remove(binary_path);
}

void test_malformed_stl_rejected() {
    const auto path = std::filesystem::temp_directory_path() / "cartmesh_stage1_bad_ascii.stl";
    {
        std::ofstream output(path, std::ios::trunc);
        output << "solid bad\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nendloop\n"
                  "endfacet\nendsolid bad\n";
    }
    expect_throw([&] { static_cast<void>(cartmesh::read_stl(path)); },
                 "缺少两个顶点的 ASCII STL 必须失败");
    std::filesystem::remove(path);
}

void test_invalid_surface_marker_output() {
    auto open = cube_triangles();
    open.pop_back();
    const cartmesh::SurfaceMesh surface(open, cartmesh::SurfaceFormat::ascii_stl, "open_cube");
    const auto diagnostics = cartmesh::diagnose_surface(surface);
    const auto path =
        std::filesystem::temp_directory_path() / "cartmesh_stage1_geometry_issues.vtp";
    cartmesh::write_surface_diagnostic_vtp(path, surface, diagnostics);
    std::ifstream input(path);
    const std::string contents((std::istreambuf_iterator<char>(input)),
                               std::istreambuf_iterator<char>());
    expect(contents.find("Name=\"issue_code\"") != std::string::npos,
           "无效几何 VTP 必须包含问题类型字段");
    expect(contents.find("NumberOfVerts=\"3\"") != std::string::npos,
           "缺一个三角形的立方体必须输出三个边界边位置标记");
    std::filesystem::remove(path);
}

void test_bvh_query_and_classification() {
    const auto surface = cube_surface();
    const cartmesh::TriangleBvh bvh(surface, 2);
    expect(bvh.statistics().triangle_count == 12, "BVH 三角形统计");
    expect(bvh.statistics().leaf_count > 1, "小叶容量应产生多个 BVH 叶节点");

    const cartmesh::AABB query_box({-0.01, 0.2, 0.2}, {0.01, 0.8, 0.8});
    auto expected = std::vector<std::uint64_t>{};
    for (std::uint64_t index = 0; index < surface.triangles().size(); ++index) {
        if (surface.triangles()[static_cast<std::size_t>(index)].bounds().intersects(query_box)) {
            expected.push_back(index);
        }
    }
    expect(bvh.query(query_box) == expected, "BVH AABB 查询必须等于全三角形真值");

    const cartmesh::SurfaceClassifier classifier(bvh);
    expect(classifier.classify({0.5, 0.5, 0.5}).classification ==
               cartmesh::PointClassification::inside,
           "立方体中心必须分类为内部");
    expect(classifier.classify({1.5, 0.5, 0.5}).classification ==
               cartmesh::PointClassification::outside,
           "立方体外部点必须分类为外部");
    expect(classifier.classify({0.0, 0.5, 0.5}).classification ==
               cartmesh::PointClassification::on_surface,
           "位于三角曲面的点必须单独标记");
    const cartmesh::AABB crossing_cell({-0.1, 0.2, 0.2}, {0.05, 0.8, 0.8});
    expect(classifier.classify(crossing_cell.center()).classification ==
               cartmesh::PointClassification::outside,
           "回归案例的单元中心必须位于立方体外部");
    expect(classifier.classify_cell(crossing_cell) ==
               cartmesh::CellClassification::intersected,
           "即使中心在外，表面穿过的单元也必须标记为相交");

    const cartmesh::AABB surface_centered_cell({-0.1, 0.4, 0.4}, {0.1, 0.6, 0.6});
    const auto surface_center = classifier.classify(surface_centered_cell.center());
    expect(surface_center.classification == cartmesh::PointClassification::on_surface,
           "单元中心恰在表面时，点分类必须保留 on_surface");
    expect(classifier.classify_cell(surface_centered_cell, surface_center.classification) ==
               cartmesh::CellClassification::intersected,
           "单元中心恰在表面时，几何单元分类必须优先标记 intersected");
    const cartmesh::UniformCartesianGrid surface_centered_grid(surface_centered_cell, 1, 1, 1);
    const auto surface_centered =
        cartmesh::classify_uniform_cells(surface_centered_grid, classifier);
    expect(surface_centered.center_on_surface_count == 1 &&
               surface_centered.intersected_count == 1 && surface_centered.conflict_count == 0,
           "均匀分类路径必须同时保留表面中心证据和相交单元结论");

    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-0.5, -0.5, -0.5}, {1.5, 1.5, 1.5}), 20, 20, 20);
    const auto classified = cartmesh::classify_uniform_cells(grid, classifier);
    expect(classified.intersected_count == 1216,
           "闭集接触策略应把单位立方体表面接触的两侧单元都标记为相交");
    expect(classified.inside_count == 512, "严格位于单位立方体内部的单元应为 8^3");
    expect(classified.outside_count == 6272, "相交壳层之外的单元应分类为外部");
    expect(classified.conflict_count == 0, "封闭立方体的多射线分类冲突必须为零");
    expect(classified.outside_count + classified.inside_count +
                   classified.intersected_count + classified.conflict_count ==
               grid.cell_count(),
           "四类单元计数必须完整分割背景网格");
    expect(classified.center_inside_count == 1000,
           "20^3 背景网格中的单位立方体应有 10^3 个内部中心");
    expect(classified.center_outside_count == 7000, "其余中心应在单位立方体外");
    expect(classified.center_on_surface_count == 0, "此错位网格没有中心落在表面");
}

void test_exact_triangle_box_intersection() {
    const cartmesh::AABB box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    expect(cartmesh::triangle_intersects_aabb(
               cartmesh::Triangle({0.2, 0.2, 0.5}, {0.8, 0.2, 0.5}, {0.2, 0.8, 0.5}), box),
           "完全位于盒内的三角形必须相交");
    expect(cartmesh::triangle_intersects_aabb(
               cartmesh::Triangle({0.0, 0.2, 0.2}, {0.0, 0.8, 0.2}, {0.0, 0.2, 0.8}), box),
           "接触盒面的三角形按闭集语义必须相交");
    expect(cartmesh::triangle_intersects_aabb(
               cartmesh::Triangle({1.0, 1.0, 1.0}, {1.2, 1.0, 1.0}, {1.0, 1.2, 1.0}), box),
           "接触盒顶点的三角形按闭集语义必须相交");
    expect(!cartmesh::triangle_intersects_aabb(
               cartmesh::Triangle({1.1, 0.2, 0.2}, {1.1, 0.8, 0.2}, {1.1, 0.2, 0.8}), box),
           "与盒分离的三角形不得误报相交");
    expect(!cartmesh::triangle_intersects_aabb(
               cartmesh::Triangle({1.0 + 1.0e-8, 0.2, 0.2},
                                  {1.0 + 1.0e-8, 0.8, 0.2},
                                  {1.0 + 1.0e-8, 0.2, 0.8}),
               box),
           "具有明确正间隔的近接触三角形不得误报相交");

    const cartmesh::Triangle bounds_false_positive(
        {2.0, 0.9, 0.9}, {0.9, 2.0, 0.9}, {0.9, 0.9, 2.0});
    expect(bounds_false_positive.bounds().intersects(box),
           "最小回归案例的三角形包围盒必须与盒重叠");
    expect(!cartmesh::triangle_intersects_aabb(bounds_false_positive, box),
           "三角形包围盒重叠不能冒充精确三角形-盒相交");

    const cartmesh::Triangle point_triangle({0.5, 0.5, 0.5}, {0.5, 0.5, 0.5},
                                            {0.5, 0.5, 0.5});
    expect(cartmesh::triangle_intersects_aabb(point_triangle, box),
           "退化为盒内点的三角形仍应按几何集合报告相交");

    const std::array<cartmesh::Vec3, 3> base = {
        cartmesh::Vec3{-0.2, 0.5, 0.5}, cartmesh::Vec3{1.2, 0.5, 0.5},
        cartmesh::Vec3{0.5, 1.2, 0.5}};
    expect(std::none_of(base.begin(), base.end(), [&](const auto& vertex) {
               return box.contains(vertex);
           }),
           "穿透回归三角形的三个顶点都必须位于盒外");
    expect(cartmesh::triangle_intersects_aabb(
               cartmesh::Triangle(base[0], base[1], base[2]), box),
           "即使没有任何三角形顶点在盒内，穿过盒子的三角形仍必须相交");
    const std::array<std::array<std::size_t, 3>, 6> permutations = {
        std::array<std::size_t, 3>{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
        {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (const auto& permutation : permutations) {
        expect(cartmesh::triangle_intersects_aabb(
                   cartmesh::Triangle(base[permutation[0]], base[permutation[1]],
                                      base[permutation[2]]),
                   box),
               "三角形顶点排列和朝向不得改变相交关系");
    }

    const cartmesh::Vec3 shift{1.0e6, -2.0e6, 3.0e6};
    const double scale = 1.0e-3;
    const auto transform = [&](const cartmesh::Vec3& value) { return shift + value * scale; };
    expect(cartmesh::triangle_intersects_aabb(
               cartmesh::Triangle(transform(base[0]), transform(base[1]), transform(base[2])),
               cartmesh::AABB(transform({0.0, 0.0, 0.0}), transform({1.0, 1.0, 1.0}))),
           "平移并缩放后不得改变三角形-盒相交关系");
    const double translated_gap = 2.0e-6;
    const cartmesh::Triangle translated_near_miss(
        transform({1.0 + translated_gap, 0.2, 0.2}),
        transform({1.0 + translated_gap, 0.8, 0.2}),
        transform({1.0 + translated_gap, 0.2, 0.8}));
    expect(!cartmesh::triangle_intersects_aabb(
               translated_near_miss,
               cartmesh::AABB(transform({0.0, 0.0, 0.0}), transform({1.0, 1.0, 1.0}))),
           "大平移后仍可由多个坐标 ULP 分辨的正间隙不得被容差吞掉");

    const auto surface = cube_surface();
    const cartmesh::TriangleBvh bvh(surface, 2);
    const std::array<cartmesh::AABB, 4> boxes = {
        cartmesh::AABB({-0.1, 0.2, 0.2}, {0.1, 0.8, 0.8}),
        cartmesh::AABB({0.2, 0.2, 0.2}, {0.8, 0.8, 0.8}),
        cartmesh::AABB({1.1, 1.1, 1.1}, {1.2, 1.2, 1.2}),
        cartmesh::AABB({0.9, 0.9, 0.9}, {1.1, 1.1, 1.1})};
    for (const auto& candidate_box : boxes) {
        const bool brute_force = std::any_of(
            surface.triangles().begin(), surface.triangles().end(),
            [&](const auto& triangle) {
                return cartmesh::triangle_intersects_aabb(triangle, candidate_box);
            });
        expect(bvh.intersects_surface(candidate_box) == brute_force,
               "BVH 精确表面相交必须等于逐三角 SAT 真值");
    }
}

void test_surface_order_and_cyclic_vertex_invariance() {
    const auto baseline_triangles = cube_triangles();
    const cartmesh::SurfaceMesh baseline_surface(
        baseline_triangles, cartmesh::SurfaceFormat::ascii_stl, "baseline_cube");

    std::vector<cartmesh::Triangle> cyclic_triangles;
    cyclic_triangles.reserve(baseline_triangles.size());
    for (const auto& triangle : baseline_triangles) {
        const auto& vertex = triangle.vertices();
        cyclic_triangles.emplace_back(vertex[1], vertex[2], vertex[0]);
    }
    const cartmesh::SurfaceMesh cyclic_surface(
        cyclic_triangles, cartmesh::SurfaceFormat::ascii_stl, "cyclic_cube");

    const std::array<std::size_t, 12> shuffled_order = {7, 2, 10, 0, 5, 11,
                                                        3, 8, 1, 9, 4, 6};
    std::vector<cartmesh::Triangle> shuffled_cyclic_triangles;
    shuffled_cyclic_triangles.reserve(cyclic_triangles.size());
    for (const auto index : shuffled_order) {
        shuffled_cyclic_triangles.push_back(cyclic_triangles[index]);
    }
    const cartmesh::SurfaceMesh shuffled_surface(
        std::move(shuffled_cyclic_triangles), cartmesh::SurfaceFormat::ascii_stl,
        "shuffled_cyclic_cube");

    const auto baseline_diagnostics = cartmesh::diagnose_surface(baseline_surface);
    const auto cyclic_diagnostics = cartmesh::diagnose_surface(cyclic_surface);
    const auto shuffled_diagnostics = cartmesh::diagnose_surface(shuffled_surface);
    expect(baseline_diagnostics.valid_for_stage1_classification() &&
               cyclic_diagnostics.valid_for_stage1_classification() &&
               shuffled_diagnostics.valid_for_stage1_classification(),
           "三角形顺序和循环顶点重排不得改变封闭、定向、流形有效性");
    expect_near(cyclic_diagnostics.signed_volume, baseline_diagnostics.signed_volume,
                1.0e-14, "循环顶点重排不得改变有向体积");
    expect_near(shuffled_diagnostics.signed_volume, baseline_diagnostics.signed_volume,
                1.0e-14, "三角形顺序和循环顶点重排不得改变有向体积");

    const cartmesh::TriangleBvh baseline_bvh(baseline_surface, 2);
    const cartmesh::TriangleBvh cyclic_bvh(cyclic_surface, 2);
    const cartmesh::TriangleBvh shuffled_bvh(shuffled_surface, 2);
    const cartmesh::SurfaceClassifier baseline_classifier(baseline_bvh);
    const cartmesh::SurfaceClassifier cyclic_classifier(cyclic_bvh);
    const cartmesh::SurfaceClassifier shuffled_classifier(shuffled_bvh);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-0.5, -0.5, -0.5}, {1.5, 1.5, 1.5}), 8, 8, 8);
    const auto baseline = cartmesh::classify_uniform_cells(grid, baseline_classifier);
    const auto cyclic = cartmesh::classify_uniform_cells(grid, cyclic_classifier);
    const auto shuffled = cartmesh::classify_uniform_cells(grid, shuffled_classifier);

    expect(cyclic.cell_classification == baseline.cell_classification &&
               cyclic.center_point_classification == baseline.center_point_classification,
           "单独循环顶点重排不得改变单元或中心点分类字节");
    expect(shuffled.cell_classification == baseline.cell_classification,
           "三角形顺序和循环顶点重排不得改变物理单元分类字节");
    expect(shuffled.center_point_classification == baseline.center_point_classification,
           "三角形顺序和循环顶点重排不得改变中心点分类字节");
    expect(shuffled.outside_count == baseline.outside_count &&
               shuffled.inside_count == baseline.inside_count &&
               shuffled.intersected_count == baseline.intersected_count &&
               shuffled.conflict_count == baseline.conflict_count,
           "三角形顺序和循环顶点重排不得改变物理分类计数");
    expect(classification_hash(grid, cyclic) == classification_hash(grid, baseline) &&
               classification_hash(grid, shuffled) == classification_hash(grid, baseline),
           "三角形顺序和循环顶点重排后分类结果哈希必须稳定");
}

void test_full_pipeline_scale_translation_and_component_orientation() {
    const auto classify_cube = [](const cartmesh::Vec3& shift, double scale) {
        const cartmesh::SurfaceMesh surface(
            transformed_cube_triangles(shift, scale), cartmesh::SurfaceFormat::ascii_stl,
            "transformed_cube");
        const auto diagnostics = cartmesh::diagnose_surface(surface);
        expect(diagnostics.valid_for_stage1_classification(),
               "正向平移/缩放立方体必须通过逐分量方向诊断");
        expect_near(diagnostics.signed_volume, scale * scale * scale,
                    std::abs(scale * scale * scale) * 1.0e-12,
                    "平移/缩放后的稳定有向体积");
        const cartmesh::Vec3 padding{0.05 * scale, 0.05 * scale, 0.05 * scale};
        const cartmesh::UniformCartesianGrid grid(
            cartmesh::AABB(shift - padding,
                           shift + cartmesh::Vec3{scale, scale, scale} + padding),
            20, 20, 20);
        const cartmesh::TriangleBvh bvh(surface, 2);
        return cartmesh::classify_uniform_cells(grid, cartmesh::SurfaceClassifier(bvh));
    };

    const auto unit = classify_cube({0.0, 0.0, 0.0}, 1.0);
    const auto tiny = classify_cube({0.0, 0.0, 0.0}, 1.0e-7);
    const auto translated = classify_cube({1.0e9, 1.0e9, 1.0e9}, 1.0);
    for (const auto* classification : {&unit, &tiny, &translated}) {
        expect(classification->inside_count == 5832 && classification->outside_count == 0 &&
                   classification->intersected_count == 2168 &&
                   classification->conflict_count == 0,
               "完整分类流水线必须在常规、小尺度和大平移下保持解析立方体计数");
    }

    const cartmesh::SurfaceMesh reversed(
        transformed_cube_triangles({0.0, 0.0, 0.0}, 1.0, true),
        cartmesh::SurfaceFormat::ascii_stl, "reversed_cube");
    const auto reversed_diagnostics = cartmesh::diagnose_surface(reversed);
    expect(reversed_diagnostics.component_orientation_mismatch_count == 1 &&
               reversed_diagnostics.valid_for_stage1_classification(),
           "整体反向外壳应报告方向警告，但不阻断奇偶分类");
    expect_near(reversed_diagnostics.material_volume, 1.0, 1.0e-13,
                "反向外壳不应使材料体积变为负数");

    auto mixed_triangles = transformed_cube_triangles({0.0, 0.0, 0.0}, 1.0);
    const auto reversed_second = transformed_cube_triangles({2.0, 0.0, 0.0}, 1.0, true);
    mixed_triangles.insert(mixed_triangles.end(), reversed_second.begin(), reversed_second.end());
    const cartmesh::SurfaceMesh mixed(std::move(mixed_triangles),
                                      cartmesh::SurfaceFormat::ascii_stl,
                                      "two_cubes_one_reversed");
    const auto mixed_diagnostics = cartmesh::diagnose_surface(mixed);
    expect(mixed_diagnostics.connected_component_count == 2 &&
               mixed_diagnostics.component_orientation_mismatch_count == 1 &&
               mixed_diagnostics.valid_for_stage1_classification(),
           "多分量方向问题应显式报告但不阻断分类");
    expect_near(mixed_diagnostics.material_volume, 2.0, 1.0e-13,
                "多分量的材料体积不得被相反法向相消");

    auto nested_triangles = transformed_cube_triangles({0.0, 0.0, 0.0}, 2.0);
    const auto cavity = transformed_cube_triangles({0.75, 0.75, 0.75}, 0.5, true);
    nested_triangles.insert(nested_triangles.end(), cavity.begin(), cavity.end());
    const cartmesh::SurfaceMesh nested(std::move(nested_triangles),
                                       cartmesh::SurfaceFormat::ascii_stl,
                                       "nested_cavity");
    const auto nested_diagnostics = cartmesh::diagnose_surface(nested);
    expect(nested_diagnostics.valid_for_stage1_classification() &&
               nested_diagnostics.components.size() == 2 &&
               nested_diagnostics.components[1].nesting_depth == 1 &&
               nested_diagnostics.components[1].expected_orientation_sign == -1 &&
               nested_diagnostics.components[1].orientation_matches_nesting,
           "嵌套内腔分量必须要求与嵌套深度一致的负向方向");
    expect_near(nested_diagnostics.signed_volume, 7.875, 1.0e-13,
                "外壳减内腔的稳定有向体积");
    expect_near(nested_diagnostics.material_volume, 7.875, 1.0e-13,
                "嵌套内腔的材料体积");

    const std::array<cartmesh::Vec3, 4> sliver_points = {
        cartmesh::Vec3{0.0, 0.0, 0.0}, cartmesh::Vec3{1.0, 0.0, 0.0},
        cartmesh::Vec3{0.0, 1.0e-20, 0.0}, cartmesh::Vec3{0.0, 0.0, 1.0}};
    const cartmesh::SurfaceMesh sliver(tetrahedron_triangles(sliver_points),
                                       cartmesh::SurfaceFormat::ascii_stl,
                                       "sliver_tetrahedron");
    const auto sliver_diagnostics = cartmesh::diagnose_surface(sliver);
    expect(sliver_diagnostics.degenerate_triangle_count > 0 &&
               !sliver_diagnostics.valid_for_stage1_classification(),
           "长边正常但高度近零的极瘦三角形必须在几何诊断阶段被拒绝");
}

void test_stl_sphere_and_cylinder_convergence() {
    const auto directory = std::filesystem::temp_directory_path();
    const auto sphere_path = directory / "cartmesh_stage1_faceted_sphere.stl";
    const auto cylinder_path = directory / "cartmesh_stage1_faceted_cylinder.stl";
    const auto tube_path = directory / "cartmesh_stage1_faceted_tube.stl";
    write_ascii_stl(sphere_path, sphere_triangles(12, 24));
    write_ascii_stl(cylinder_path, cylinder_triangles(48));
    write_ascii_stl(tube_path, tube_triangles(48, 0.5));
    const auto sphere = cartmesh::read_stl(sphere_path);
    const auto cylinder = cartmesh::read_stl(cylinder_path);
    const auto tube = cartmesh::read_stl(tube_path);

    const double sphere_exact = std::abs(cartmesh::diagnose_surface(sphere).signed_volume);
    const auto sphere_coarse = classified_metrics(sphere, 12);
    const auto sphere_fine = classified_metrics(sphere, 32);
    const double sphere_coarse_error =
        std::abs(sphere_coarse.center_sample_volume - sphere_exact);
    const double sphere_fine_error = std::abs(sphere_fine.center_sample_volume - sphere_exact);
    expect(sphere_fine_error < sphere_coarse_error,
           "STL 球体中心分类体积误差从 12^3 到 32^3 应下降");
    expect(sphere_coarse.definitely_inside_volume <= sphere_exact &&
               sphere_exact <= sphere_coarse.inside_plus_intersected_volume,
           "粗球体的真实多面体体积必须落在 inside/intersected 分类体积区间");
    expect(sphere_fine.definitely_inside_volume <= sphere_exact &&
               sphere_exact <= sphere_fine.inside_plus_intersected_volume,
           "细球体的真实多面体体积必须落在 inside/intersected 分类体积区间");
    expect(sphere_fine.inside_plus_intersected_volume -
                   sphere_fine.definitely_inside_volume <
               sphere_coarse.inside_plus_intersected_volume -
                   sphere_coarse.definitely_inside_volume,
           "球体分类体积区间应随网格加密收缩");

    const double cylinder_exact = std::abs(cartmesh::diagnose_surface(cylinder).signed_volume);
    const auto cylinder_coarse = classified_metrics(cylinder, 12);
    const auto cylinder_fine = classified_metrics(cylinder, 36);
    const double cylinder_coarse_error =
        std::abs(cylinder_coarse.center_sample_volume - cylinder_exact);
    const double cylinder_fine_error =
        std::abs(cylinder_fine.center_sample_volume - cylinder_exact);
    expect(cylinder_fine_error < cylinder_coarse_error,
           "STL 圆柱中心分类体积误差从 12^3 到 36^3 应下降");
    expect(cylinder_coarse.definitely_inside_volume <= cylinder_exact &&
               cylinder_exact <= cylinder_coarse.inside_plus_intersected_volume,
           "粗圆柱的真实多面体体积必须落在 inside/intersected 分类体积区间");
    expect(cylinder_fine.definitely_inside_volume <= cylinder_exact &&
               cylinder_exact <= cylinder_fine.inside_plus_intersected_volume,
           "细圆柱的真实多面体体积必须落在 inside/intersected 分类体积区间");
    expect(cylinder_fine.inside_plus_intersected_volume -
                   cylinder_fine.definitely_inside_volume <
               cylinder_coarse.inside_plus_intersected_volume -
                   cylinder_coarse.definitely_inside_volume,
           "圆柱分类体积区间应随网格加密收缩");

    const double tube_exact = std::abs(cartmesh::diagnose_surface(tube).signed_volume);
    const auto tube_coarse = classified_metrics(tube, 12);
    const auto tube_fine = classified_metrics(tube, 36);
    expect(tube_coarse.definitely_inside_volume <= tube_exact &&
               tube_exact <= tube_coarse.inside_plus_intersected_volume,
           "粗圆管的真实多面体体积必须落在 inside/intersected 分类体积区间");
    expect(tube_fine.definitely_inside_volume <= tube_exact &&
               tube_exact <= tube_fine.inside_plus_intersected_volume,
           "细圆管的真实多面体体积必须落在 inside/intersected 分类体积区间");
    expect(tube_fine.inside_plus_intersected_volume - tube_fine.definitely_inside_volume <
               tube_coarse.inside_plus_intersected_volume -
                   tube_coarse.definitely_inside_volume,
           "圆管分类体积区间应随网格加密收缩");

    std::filesystem::remove(sphere_path);
    std::filesystem::remove(cylinder_path);
    std::filesystem::remove(tube_path);
}

} // 匿名命名空间

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"封闭定向流形诊断", test_closed_oriented_manifold_diagnostics},
        {"无效表面诊断可见", test_invalid_surface_diagnostics_remain_visible},
        {"非流形顶点诊断", test_non_manifold_vertex_diagnostic},
        {"ASCII 与二进制 STL 导入", test_ascii_and_binary_stl_import},
        {"畸形 STL 拒绝", test_malformed_stl_rejected},
        {"无效表面位置标记", test_invalid_surface_marker_output},
        {"BVH 查询与均匀分类", test_bvh_query_and_classification},
        {"精确三角形与盒相交", test_exact_triangle_box_intersection},
        {"表面顺序与循环顶点不变性", test_surface_order_and_cyclic_vertex_invariance},
        {"完整尺度平移与逐分量方向", test_full_pipeline_scale_translation_and_component_orientation},
        {"STL 球体与圆柱分类收敛", test_stl_sphere_and_cylinder_convergence},
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
    std::cout << "阶段1测试数=" << tests.size() << " 失败数=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
