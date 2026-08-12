#include "cartmesh/cutcell/ConvexPolyhedron.hpp"
#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"
#include "cartmesh/cutcell/ConvexSurfaceCutter.hpp"
#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"
#include "cartmesh/io/CutCellVtkWriter.hpp"
#include "cartmesh/io/CutCellJsonWriter.hpp"
#include "cartmesh/io/OpenFoamWriter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <limits>
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

void expect_near(double actual, double expected, double tolerance,
                 const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        std::ostringstream detail;
        detail << std::setprecision(17) << message << "：实际值=" << actual << " 期望值=" << expected
               << " 容差=" << tolerance;
        throw TestFailure(detail.str());
    }
}

void expect_vec_near(const cartmesh::Vec3& actual, const cartmesh::Vec3& expected,
                     double tolerance, const std::string& message) {
    expect_near(actual.x, expected.x, tolerance, message + " x");
    expect_near(actual.y, expected.y, tolerance, message + " y");
    expect_near(actual.z, expected.z, tolerance, message + " z");
}

[[nodiscard]] std::string read_all(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw TestFailure("无法读取测试产物：" + path.string());
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

[[nodiscard]] std::vector<std::uint64_t> read_foam_label_list(
    const std::filesystem::path& path) {
    const std::string text = read_all(path);
    const std::size_t begin_marker = text.find("\n(\n");
    const std::size_t end_marker = text.rfind("\n)\n");
    if (begin_marker == std::string::npos || end_marker == std::string::npos ||
        end_marker <= begin_marker) {
        throw TestFailure("OpenFOAM labelList 结构无效：" + path.string());
    }
    std::istringstream values(
        text.substr(begin_marker + 3U, end_marker - begin_marker - 3U));
    std::vector<std::uint64_t> result;
    std::uint64_t value = 0;
    while (values >> value) result.push_back(value);
    return result;
}

void expect_openfoam_cell_labels(const std::filesystem::path& case_directory,
                                 std::uint64_t expected_cell_count) {
    const auto poly_mesh = case_directory / "constant/polyMesh";
    const auto owner = read_foam_label_list(poly_mesh / "owner");
    const auto neighbour = read_foam_label_list(poly_mesh / "neighbour");
    expect(!owner.empty(), "OpenFOAM owner 不得为空");
    std::uint64_t maximum_cell = 0;
    for (const auto id : owner) maximum_cell = std::max(maximum_cell, id);
    for (const auto id : neighbour) maximum_cell = std::max(maximum_cell, id);
    expect(maximum_cell + 1U == expected_cell_count,
           "OpenFOAM owner/neighbour 必须覆盖全部 solver cell");
    expect(neighbour.size() <= owner.size(),
           "OpenFOAM neighbour 数不得超过 owner 数");
    for (std::size_t face = 0; face < neighbour.size(); ++face) {
        expect(owner[face] < neighbour[face],
               "OpenFOAM 内部面必须满足 owner < neighbour");
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

[[nodiscard]] std::size_t embedded_face_index(const cartmesh::ConvexPolyhedron& polyhedron) {
    std::size_t result = polyhedron.faces.size();
    for (std::size_t index = 0; index < polyhedron.faces.size(); ++index) {
        if (polyhedron.faces[index].kind == cartmesh::PolyhedronFaceKind::embedded_boundary) {
            expect(result == polyhedron.faces.size(), "单平面裁剪只能产生一个几何边界面");
            result = index;
        }
    }
    expect(result != polyhedron.faces.size(), "裁剪结果必须包含几何边界面");
    return result;
}

[[nodiscard]] std::vector<cartmesh::Triangle> cube_triangles(bool reverse = false) {
    const std::vector<cartmesh::Vec3> point = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}};
    const std::vector<std::array<std::size_t, 3>> face = {
        {0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
        {3, 7, 6}, {3, 6, 2}, {0, 4, 7}, {0, 7, 3}, {1, 2, 6}, {1, 6, 5}};
    std::vector<cartmesh::Triangle> result;
    for (const auto& indices : face) {
        result.emplace_back(point[indices[0]], point[reverse ? indices[2] : indices[1]],
                            point[reverse ? indices[1] : indices[2]]);
    }
    return result;
}

[[nodiscard]] std::vector<cartmesh::Triangle> axis_aligned_box_triangles(
    const cartmesh::Vec3& minimum, const cartmesh::Vec3& maximum) {
    auto result = cube_triangles();
    const cartmesh::Vec3 extent = maximum - minimum;
    for (auto& triangle : result) {
        const auto vertices = triangle.vertices();
        triangle = cartmesh::Triangle(
            minimum + cartmesh::Vec3{vertices[0].x * extent.x,
                                     vertices[0].y * extent.y,
                                     vertices[0].z * extent.z},
            minimum + cartmesh::Vec3{vertices[1].x * extent.x,
                                     vertices[1].y * extent.y,
                                     vertices[1].z * extent.z},
            minimum + cartmesh::Vec3{vertices[2].x * extent.x,
                                     vertices[2].y * extent.y,
                                     vertices[2].z * extent.z});
    }
    return result;
}

[[nodiscard]] std::vector<cartmesh::Triangle> octahedron_triangles() {
    const cartmesh::Vec3 positive_x{1.0, 0.0, 0.0};
    const cartmesh::Vec3 negative_x{-1.0, 0.0, 0.0};
    const cartmesh::Vec3 positive_y{0.0, 1.0, 0.0};
    const cartmesh::Vec3 negative_y{0.0, -1.0, 0.0};
    const cartmesh::Vec3 positive_z{0.0, 0.0, 1.0};
    const cartmesh::Vec3 negative_z{0.0, 0.0, -1.0};
    return {{positive_z, positive_x, positive_y},
            {positive_z, positive_y, negative_x},
            {positive_z, negative_x, negative_y},
            {positive_z, negative_y, positive_x},
            {negative_z, positive_y, positive_x},
            {negative_z, negative_x, positive_y},
            {negative_z, negative_y, negative_x},
            {negative_z, positive_x, negative_y}};
}

[[nodiscard]] std::vector<cartmesh::Triangle> subdivided_unit_sphere(
    std::uint32_t levels) {
    auto triangles = octahedron_triangles();
    for (std::uint32_t level = 0; level < levels; ++level) {
        std::vector<cartmesh::Triangle> refined;
        refined.reserve(triangles.size() * 4U);
        for (const auto& triangle : triangles) {
            const auto& vertex = triangle.vertices();
            const cartmesh::Vec3 ab =
                (vertex[0] + vertex[1]) / cartmesh::norm(vertex[0] + vertex[1]);
            const cartmesh::Vec3 bc =
                (vertex[1] + vertex[2]) / cartmesh::norm(vertex[1] + vertex[2]);
            const cartmesh::Vec3 ca =
                (vertex[2] + vertex[0]) / cartmesh::norm(vertex[2] + vertex[0]);
            refined.emplace_back(vertex[0], ab, ca);
            refined.emplace_back(ab, vertex[1], bc);
            refined.emplace_back(ca, bc, vertex[2]);
            refined.emplace_back(ab, bc, ca);
        }
        triangles = std::move(refined);
    }
    return triangles;
}

[[nodiscard]] std::vector<cartmesh::Triangle> l_prism_triangles() {
    std::vector<cartmesh::Triangle> result;
    const auto point = [](double x, double y, double z) {
        return cartmesh::Vec3{x, y, z};
    };
    const std::array<std::array<double, 4>, 3> squares = {
        std::array<double, 4>{0.0, 0.0, 1.0, 1.0},
        std::array<double, 4>{0.0, 1.0, 1.0, 2.0},
        std::array<double, 4>{1.0, 0.0, 2.0, 1.0}};
    for (const auto& square : squares) {
        const auto a0 = point(square[0], square[1], 0.0);
        const auto b0 = point(square[2], square[1], 0.0);
        const auto c0 = point(square[2], square[3], 0.0);
        const auto d0 = point(square[0], square[3], 0.0);
        const auto a1 = point(square[0], square[1], 1.0);
        const auto b1 = point(square[2], square[1], 1.0);
        const auto c1 = point(square[2], square[3], 1.0);
        const auto d1 = point(square[0], square[3], 1.0);
        result.emplace_back(a1, b1, c1);
        result.emplace_back(a1, c1, d1);
        result.emplace_back(a0, c0, b0);
        result.emplace_back(a0, d0, c0);
    }
    const std::array<std::array<double, 4>, 8> boundary_edges = {
        std::array<double, 4>{0.0, 0.0, 1.0, 0.0},
        std::array<double, 4>{1.0, 0.0, 2.0, 0.0},
        std::array<double, 4>{2.0, 0.0, 2.0, 1.0},
        std::array<double, 4>{2.0, 1.0, 1.0, 1.0},
        std::array<double, 4>{1.0, 1.0, 1.0, 2.0},
        std::array<double, 4>{1.0, 2.0, 0.0, 2.0},
        std::array<double, 4>{0.0, 2.0, 0.0, 1.0},
        std::array<double, 4>{0.0, 1.0, 0.0, 0.0}};
    for (const auto& edge : boundary_edges) {
        const auto a0 = point(edge[0], edge[1], 0.0);
        const auto b0 = point(edge[2], edge[3], 0.0);
        const auto a1 = point(edge[0], edge[1], 1.0);
        const auto b1 = point(edge[2], edge[3], 1.0);
        result.emplace_back(a0, b0, b1);
        result.emplace_back(a0, b1, a1);
    }
    return result;
}

void expect_closed_positive(const cartmesh::PolyhedronGeometry& geometry,
                            const std::string& message) {
    expect(geometry.closed, message + " 必须闭合");
    expect(geometry.positive_volume, message + " 必须具有正体积");
    expect(cartmesh::norm(geometry.oriented_area_vector_sum) < 1.0e-12,
           message + " 的各面有向面积和必须为零");
}

void test_box_geometry() {
    const auto polyhedron = cartmesh::make_box_polyhedron(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
    const auto geometry = cartmesh::measure_polyhedron(polyhedron);
    expect(polyhedron.vertices.size() == 8 && polyhedron.faces.size() == 6,
           "背景盒必须由 8 顶点和 6 个四边形面组成");
    expect_closed_positive(geometry, "单位背景盒");
    expect_near(geometry.volume, 1.0, 1.0e-14, "单位背景盒体积");
    expect_vec_near(geometry.centroid, {0.5, 0.5, 0.5}, 1.0e-14,
                    "单位背景盒质心");
    for (const auto& face : geometry.faces) {
        expect_near(face.area, 1.0, 1.0e-14, "单位背景盒面积");
    }
}

void test_tetrahedron_constructor_orientation() {
    const auto positive = cartmesh::make_tetrahedron_polyhedron(
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0});
    const auto negative_input = cartmesh::make_tetrahedron_polyhedron(
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
        {0.0, 1.0, 0.0});
    const auto positive_geometry = cartmesh::measure_polyhedron(positive);
    const auto negative_geometry = cartmesh::measure_polyhedron(negative_input);
    expect_closed_positive(positive_geometry, "正序四面体");
    expect_closed_positive(negative_geometry, "任意输入序四面体");
    expect_near(positive_geometry.volume, 1.0 / 6.0, 1.0e-14,
                "标准四面体体积");
    expect_near(negative_geometry.volume, positive_geometry.volume, 1.0e-14,
                "构造器必须统一四面体面绕序");
}

void test_axis_aligned_cut_and_complement() {
    const auto box = cartmesh::make_box_polyhedron(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
    const auto lower = cartmesh::clip_convex_polyhedron(
        box, cartmesh::OrientedHalfSpace({0.25, 0.0, 0.0}, {1.0, 0.0, 0.0}, 42));
    const auto upper = cartmesh::clip_convex_polyhedron(
        box, cartmesh::OrientedHalfSpace({0.25, 0.0, 0.0}, {-1.0, 0.0, 0.0}, 42));
    const auto lower_geometry = cartmesh::measure_polyhedron(lower);
    const auto upper_geometry = cartmesh::measure_polyhedron(upper);
    expect_closed_positive(lower_geometry, "x<=0.25 Cut-cell");
    expect_closed_positive(upper_geometry, "x>=0.25 互补 Cut-cell");
    expect_near(lower_geometry.volume, 0.25, 1.0e-14, "x<=0.25 体积");
    expect_near(upper_geometry.volume, 0.75, 1.0e-14, "x>=0.25 体积");
    expect_near(lower_geometry.volume + upper_geometry.volume, 1.0, 1.0e-14,
                "流体与互补子体积必须闭合");
    expect_vec_near(lower_geometry.centroid, {0.125, 0.5, 0.5}, 1.0e-14,
                    "x<=0.25 质心");
    expect_vec_near(upper_geometry.centroid, {0.625, 0.5, 0.5}, 1.0e-14,
                    "x>=0.25 质心");

    const auto lower_cut = embedded_face_index(lower);
    const auto upper_cut = embedded_face_index(upper);
    expect(lower.faces[lower_cut].source_id == 42 && upper.faces[upper_cut].source_id == 42,
           "几何边界 ID 必须穿过裁剪保留");
    expect_near(lower_geometry.faces[lower_cut].area, 1.0, 1.0e-14,
                "x=0.25 切面面积");
    expect_vec_near(lower_geometry.faces[lower_cut].centroid, {0.25, 0.5, 0.5},
                    1.0e-14, "x=0.25 切面质心");
    expect_vec_near(lower_geometry.faces[lower_cut].outward_normal, {1.0, 0.0, 0.0},
                    1.0e-14, "x<=0.25 切面外法向");
    expect_vec_near(upper_geometry.faces[upper_cut].outward_normal, {-1.0, 0.0, 0.0},
                    1.0e-14, "互补切面外法向");
}

void test_oblique_tetrahedron_cut() {
    const auto box = cartmesh::make_box_polyhedron(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
    const auto cut = cartmesh::clip_convex_polyhedron(
        box, cartmesh::OrientedHalfSpace({1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 7));
    const auto geometry = cartmesh::measure_polyhedron(cut);
    expect(cut.vertices.size() == 4 && cut.faces.size() == 4,
           "x+y+z<=1 必须产生四面体");
    expect_closed_positive(geometry, "倾斜平面四面体 Cut-cell");
    expect_near(geometry.volume, 1.0 / 6.0, 1.0e-14,
                "x+y+z<=1 四面体体积");
    expect_vec_near(geometry.centroid, {0.25, 0.25, 0.25}, 1.0e-14,
                    "x+y+z<=1 四面体质心");
    const auto cut_face = embedded_face_index(cut);
    expect_near(geometry.faces[cut_face].area, std::sqrt(3.0) * 0.5, 1.0e-14,
                "x+y+z=1 三角切面面积");
    expect_vec_near(geometry.faces[cut_face].centroid,
                    {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, 1.0e-14,
                    "x+y+z=1 三角切面质心");
    expect_vec_near(geometry.faces[cut_face].outward_normal,
                    {1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0),
                     1.0 / std::sqrt(3.0)},
                    1.0e-14, "x+y+z<=1 切面外法向");
}

void test_small_and_translated_cuts() {
    constexpr double epsilon = 1.0e-6;
    const auto box = cartmesh::make_box_polyhedron(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
    const auto tiny = cartmesh::clip_convex_polyhedron(
        box, cartmesh::OrientedHalfSpace({epsilon, 0.0, 0.0}, {1.0, 1.0, 1.0}));
    const auto tiny_geometry = cartmesh::measure_polyhedron(tiny);
    expect_closed_positive(tiny_geometry, "近顶点微小 Cut-cell");
    expect_near(tiny_geometry.volume, epsilon * epsilon * epsilon / 6.0, 2.0e-29,
                "近顶点四面体体积");
    expect_vec_near(tiny_geometry.centroid,
                    {epsilon / 4.0, epsilon / 4.0, epsilon / 4.0}, 2.0e-17,
                    "近顶点四面体质心");

    constexpr double shift = 1.0e9;
    const auto translated_box = cartmesh::make_box_polyhedron(cartmesh::AABB(
        {shift, shift, shift}, {shift + 1.0, shift + 1.0, shift + 1.0}));
    const auto translated = cartmesh::clip_convex_polyhedron(
        translated_box,
        cartmesh::OrientedHalfSpace({shift + 0.25, shift, shift}, {1.0, 0.0, 0.0}));
    const auto translated_geometry = cartmesh::measure_polyhedron(translated);
    expect_closed_positive(translated_geometry, "大平移 Cut-cell");
    expect_near(translated_geometry.volume, 0.25, 1.0e-13,
                "大平移后体积必须不变");
    expect_vec_near(translated_geometry.centroid,
                    {shift + 0.125, shift + 0.5, shift + 0.5}, 2.0e-7,
                    "大平移后质心");
}

void test_empty_and_unchanged_clips() {
    const auto box = cartmesh::make_box_polyhedron(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
    const auto unchanged = cartmesh::clip_convex_polyhedron(
        box, cartmesh::OrientedHalfSpace({2.0, 0.0, 0.0}, {1.0, 0.0, 0.0}));
    const auto empty = cartmesh::clip_convex_polyhedron(
        box, cartmesh::OrientedHalfSpace({-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}));
    const auto face_contact = cartmesh::clip_convex_polyhedron(
        box, cartmesh::OrientedHalfSpace({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}));
    const auto vertex_contact = cartmesh::clip_convex_polyhedron(
        box, cartmesh::OrientedHalfSpace({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
    expect_near(cartmesh::measure_polyhedron(unchanged).volume, 1.0, 1.0e-14,
                "完全在保留侧的盒不得被修改");
    expect(empty.empty() && cartmesh::measure_polyhedron(empty).volume == 0.0,
           "完全在丢弃侧的盒必须返回空集");
    expect(face_contact.empty() && vertex_contact.empty(),
           "仅与盒面或顶点接触的零体积保留集必须返回空集");
    expect_throw([] {
        static_cast<void>(cartmesh::OrientedHalfSpace({0.0, 0.0, 0.0},
                                                       {0.0, 0.0, 0.0}));
    }, "零法向半空间必须拒绝");
}

void test_convex_stl_cube_cut_cells() {
    const cartmesh::SurfaceMesh surface(cube_triangles(),
                                        cartmesh::SurfaceFormat::ascii_stl,
                                        "unit_cube");
    const cartmesh::SurfaceMesh reversed_surface(
        cube_triangles(true), cartmesh::SurfaceFormat::ascii_stl,
        "reversed_unit_cube");
    const cartmesh::ConvexSurfaceCutter cutter(surface, 9);
    const cartmesh::ConvexSurfaceCutter reversed_cutter(reversed_surface, 9);
    expect(cutter.plane_count() == 6 && reversed_cutter.plane_count() == 6,
           "立方体 12 三角面必须确定性合并为 6 个支撑平面");
    expect(!cutter.input_orientation_reversed() &&
               reversed_cutter.input_orientation_reversed(),
           "凸 STL 工作表面必须显式记录整壳翻向");

    const cartmesh::AABB crossing_cell({0.75, 0.25, 0.25}, {1.25, 0.75, 0.75});
    const auto cut = cutter.cut_box(crossing_cell);
    const auto reversed_cut = reversed_cutter.cut_box(crossing_cell);
    expect(cut.cut && cut.solid_geometry.closed && cut.solid_geometry.positive_volume,
           "穿过 x=1 外壳的单元必须产生闭合正体积 Cut-cell");
    expect_near(cut.solid_geometry.volume, 0.0625, 1.0e-14,
                "立方体相交单元的固体体积");
    expect_near(cut.solid_volume_fraction, 0.5, 1.0e-14,
                "立方体相交单元的固体体积分数");
    expect_near(cut.fluid_volume_fraction, 0.5, 1.0e-14,
                "立方体相交单元的流体体积分数");
    expect_vec_near(cut.solid_geometry.centroid, {0.875, 0.5, 0.5}, 1.0e-14,
                    "立方体相交单元的固体质心");
    expect_vec_near(cut.fluid_centroid, {1.125, 0.5, 0.5}, 1.0e-14,
                    "立方体相交单元的流体质心");
    expect_near(cut.embedded_boundary_area, 0.25, 1.0e-14,
                "立方体相交单元的几何边界面积");
    expect_near(cut.volume_conservation_residual, 0.0, 1.0e-15,
                "单元固体加流体体积守恒");
    expect_near(reversed_cut.solid_geometry.volume, cut.solid_geometry.volume,
                1.0e-14, "整壳反向不得改变 Cut-cell 体积");

    const auto inside = cutter.cut_box(
        cartmesh::AABB({0.25, 0.25, 0.25}, {0.75, 0.75, 0.75}));
    const auto outside = cutter.cut_box(
        cartmesh::AABB({1.1, 0.25, 0.25}, {1.6, 0.75, 0.75}));
    expect(!inside.cut && inside.solid_volume_fraction == 1.0 &&
               inside.fluid_volume == 0.0,
           "完全固体内单元不是 Cut-cell");
    expect(!outside.cut && outside.solid_volume_fraction == 0.0 &&
               outside.fluid_volume_fraction == 1.0,
           "完全固体外单元应保留全部流体体积");
}

void test_convex_stl_multiple_surface_patches() {
    const cartmesh::SurfaceMesh octahedron(
        octahedron_triangles(), cartmesh::SurfaceFormat::ascii_stl,
        "unit_octahedron");
    const cartmesh::ConvexSurfaceCutter cutter(octahedron, 17);
    const auto cut = cutter.cut_box(
        cartmesh::AABB({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}));
    expect(cutter.plane_count() == 8,
           "八面体的 8 个不共面三角片必须全部参与重建");
    expect(cut.cut && cut.solid_geometry.closed && cut.solid_geometry.positive_volume,
           "单元内多三角表面片必须重建为闭合多面体");
    expect_near(cut.solid_geometry.volume, 4.0 / 3.0, 1.0e-13,
                "单位八面体体积");
    expect_vec_near(cut.solid_geometry.centroid, {0.0, 0.0, 0.0}, 1.0e-14,
                    "单位八面体质心");
    expect_near(cut.embedded_boundary_area, 4.0 * std::sqrt(3.0), 1.0e-13,
                "单位八面体的 8 个几何边界面面积和");
    expect_near(cut.solid_volume_fraction, 1.0 / 6.0, 1.0e-14,
                "八面体在 2x2x2 背景盒中的体积分数");
    expect_near(cut.fluid_volume_fraction, 5.0 / 6.0, 1.0e-14,
                "八面体外部流体体积分数");
    expect_near(cut.volume_conservation_residual, 0.0, 1.0e-14,
                "八面体固体/流体体积守恒");
}

void test_uniform_grid_cut_cell_topology() {
    const cartmesh::SurfaceMesh surface(
        axis_aligned_box_triangles({0.25, 0.25, 0.25}, {0.75, 0.75, 0.75}),
        cartmesh::SurfaceFormat::ascii_stl, "centered_half_cube");
    const cartmesh::ConvexSurfaceCutter cutter(surface, 23);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 2, 1, 1);
    const auto mesh = cartmesh::build_convex_cut_cell_mesh(grid, cutter);

    expect(mesh.fluid_cells.size() == 2 && mesh.cut_cell_count == 2,
           "两个背景单元都必须重建为 Cut-cell");
    expect(mesh.full_fluid_cell_count == 0 && mesh.full_solid_cell_count == 0,
           "该解析案例不含完整流体或完整固体单元");
    expect_near(mesh.total_fluid_volume, 0.875, 1.0e-14,
                "网格总流体体积");
    expect_near(mesh.total_embedded_boundary_area, 1.5, 1.0e-14,
                "跨单元 STL 嵌入边界总面积");
    expect(mesh.nonclosed_cell_count == 0,
           "所有 Cut-cell 的有向面积向量必须闭合");
    expect(mesh.shared_face_mismatch_count == 0,
           "相邻 Cut-cell 共享面的面积和质心必须一致");
    expect(mesh.maximum_cell_area_closure_residual < 1.0e-14,
           "Cut-cell 面积向量闭合残差");

    for (const auto& cell : mesh.fluid_cells) {
        expect_near(cell.volume, 0.4375, 1.0e-14,
                    "每个 Cut-cell 的流体体积");
        expect_near(cell.volume_fraction, 0.875, 1.0e-14,
                    "每个 Cut-cell 的流体体积分数");
        expect(cell.embedded_boundary_faces.size() == 5,
               "每个半盒 Cut-cell 必须保留五个 STL 边界片");
        for (const auto& face : cell.embedded_boundary_faces) {
            expect(face.boundary_id == 23, "嵌入边界 patch ID 必须保留");
            expect(face.vertices.size() == 4,
                   "轴对齐立方体的嵌入边界片必须保留四个显式顶点");
        }
        expect(cell.cartesian_faces[0].oriented_boundary_loops.size() >= 1 &&
                   cell.cartesian_faces[1].oriented_boundary_loops.size() >= 1,
               "每个 Cartesian 开口面必须保留显式有向边界环");
    }

    expect(mesh.internal_faces.size() == 1,
           "两个流体单元之间必须只有一个共享连接面");
    const auto& connection = mesh.internal_faces.front();
    expect(connection.first_background_cell_id == 0 &&
               connection.second_background_cell_id == 1,
           "共享面连接必须按背景单元 ID 确定性排序");
    expect(connection.first_local_face == 1 && connection.second_local_face == 0,
           "共享面局部编号必须为 x+ / x-");
    expect_near(connection.area, 0.75, 1.0e-14,
                "共享面的流体开口面积");
    expect_vec_near(connection.centroid, {0.5, 0.5, 0.5}, 1.0e-14,
                    "共享面流体开口质心");
    expect_vec_near(connection.normal, {1.0, 0.0, 0.0}, 1.0e-14,
                    "共享面法向必须从低 ID 指向高 ID");

    const auto small = cartmesh::analyze_small_cut_cells(mesh, 0.9);
    expect_near(small.minimum_cut_cell_volume_fraction, 0.875, 1.0e-14,
                "小 Cut-cell 统计必须给出最小流体体积分数");
    expect(small.cells.size() == 2 &&
               small.cells[0].background_cell_id == 0 &&
               small.cells[1].background_cell_id == 1,
           "小 Cut-cell 必须按背景单元 ID 确定性列出");
    for (const auto& entry : small.cells) {
        expect(entry.boundary_ids == std::vector<std::uint64_t>{23},
               "小 Cut-cell 必须记录所属边界 patch ID");
    }
}

void test_uniform_grid_endpoint_fractions() {
    const cartmesh::SurfaceMesh surface(cube_triangles(),
                                        cartmesh::SurfaceFormat::ascii_stl,
                                        "unit_cube");
    const cartmesh::ConvexSurfaceCutter cutter(surface);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-0.1, -0.1, -0.1}, {1.1, 1.1, 1.1}), 8, 8, 8);
    const auto mesh = cartmesh::build_convex_cut_cell_mesh(grid, cutter);
    expect(mesh.full_solid_cell_count == 216,
           "完全位于单位立方体内的 6^3 单元必须归一化为完整固体");
    expect(mesh.cut_cell_count == 296 && mesh.full_fluid_cell_count == 0,
           "外层 8^3-6^3 单元必须全部是 Cut-cell，不得出现机器精度伪流体单元");
    expect_near(mesh.total_fluid_volume, 0.728, 5.0e-15,
                "端点归一化后域内流体体积");
}

void test_disconnected_fluid_pieces_are_detected() {
    const cartmesh::SurfaceMesh slab(
        axis_aligned_box_triangles({0.4, -1.0, -1.0}, {0.6, 2.0, 2.0}),
        cartmesh::SurfaceFormat::ascii_stl, "cell_spanning_slab");
    const cartmesh::ConvexSurfaceCutter cutter(slab);
    const cartmesh::AABB cell({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const auto cut = cutter.cut_box(cell);
    expect(cut.cut && cut.fluid_decomposition_polyhedra.size() >= 2,
           "跨越单元的固体薄板必须产生显式流体多面体分解");
    expect(cut.fluid_component_count == 2,
           "薄板左右两侧必须识别为两个不连通流体分量");
    expect_near(cut.fluid_volume, 0.8, 1.0e-14,
                "两个流体分量的总体积");
    const cartmesh::UniformCartesianGrid grid(cell, 1, 1, 1);
    const auto mesh = cartmesh::build_convex_cut_cell_mesh(grid, cutter);
    expect(mesh.disconnected_fluid_cell_count == 1,
           "网格统计必须显式报告包含多个流体分量的 Cut-cell");
}

void test_grid_aligned_surface_becomes_coplanar_wall() {
    const cartmesh::SurfaceMesh surface(cube_triangles(),
                                        cartmesh::SurfaceFormat::ascii_stl,
                                        "aligned_unit_cube");
    const cartmesh::ConvexSurfaceCutter cutter(surface, 31);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-0.5, -0.5, -0.5}, {1.5, 1.5, 1.5}), 4, 4, 4);
    const auto mesh = cartmesh::build_convex_cut_cell_mesh(grid, cutter);
    expect(mesh.full_solid_cell_count == 8 && mesh.cut_cell_count == 0 &&
               mesh.full_fluid_cell_count == 56,
           "网格对齐单位立方体必须包含 8 个完整固体和 56 个完整体积流体单元");
    expect(mesh.shared_face_mismatch_count == 0 && mesh.nonclosed_cell_count == 0,
           "共面壁面必须关闭固体/流体共享面并保持流体控制体闭合");
    expect_near(mesh.total_embedded_boundary_area, 6.0, 1.0e-14,
                "共面嵌入壁面的总面积");
    expect(mesh.boundary_cell_count == 24,
           "共面边界必须显式落在 24 个面相邻流体单元上");
    expect_near(mesh.total_fluid_volume, 7.0, 1.0e-14,
                "对齐案例的流体总体积");
}

void test_general_triangle_chain_matches_convex_cube() {
    const cartmesh::SurfaceMesh surface(cube_triangles(),
                                        cartmesh::SurfaceFormat::ascii_stl,
                                        "general_cube");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 41);
    const auto cut = cutter.cut_box(
        cartmesh::AABB({0.75, 0.25, 0.25}, {1.25, 0.75, 0.75}));
    expect(cut.cut, "通用有向四面体链必须识别立方体穿透单元");
    expect_near(cut.solid_volume, 0.0625, 2.0e-14,
                "通用立方体相交固体体积");
    expect_near(cut.fluid_volume, 0.0625, 2.0e-14,
                "通用立方体相交流体体积");
    expect_vec_near(cut.solid_centroid, {0.875, 0.5, 0.5}, 2.0e-14,
                    "通用立方体固体质心");
    expect_vec_near(cut.fluid_centroid, {1.125, 0.5, 0.5}, 2.0e-14,
                    "通用立方体流体质心");
    expect_near(cut.embedded_boundary_area, 0.25, 2.0e-14,
                "通用立方体嵌入边界面积");
}

void test_general_multiple_boundary_patches_in_one_cell() {
    const auto triangles = cube_triangles();
    const cartmesh::SurfaceMesh surface(triangles,
                                        cartmesh::SurfaceFormat::ascii_stl,
                                        "multi_patch_cube");
    std::vector<std::uint64_t> patch_ids(triangles.size(), 7);
    patch_ids[10] = 110;
    patch_ids[11] = 111;
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, patch_ids);
    const auto cut = cutter.cut_box(
        cartmesh::AABB({0.75, 0.25, 0.25}, {1.25, 0.75, 0.75}));
    std::vector<std::uint64_t> found;
    for (const auto& face : cut.embedded_boundary_faces) {
        found.push_back(face.boundary_id);
    }
    std::sort(found.begin(), found.end());
    expect(found == std::vector<std::uint64_t>({110, 111}),
           "同一 Cut-cell 内两个三角 patch ID 必须分别保留");
    expect_near(cut.embedded_boundary_area, 0.25, 2.0e-14,
                "多 patch 划分不得改变几何边界总面积");
}

void test_general_nonconvex_l_prism() {
    const cartmesh::SurfaceMesh surface(l_prism_triangles(),
                                        cartmesh::SurfaceFormat::ascii_stl,
                                        "nonconvex_l_prism");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 43);
    const auto cut = cutter.cut_box(
        cartmesh::AABB({0.0, 0.0, 0.0}, {2.0, 2.0, 1.0}));
    expect(cut.cut, "非凸 L 形棱柱必须形成真实流体/固体分割");
    expect_near(cut.solid_volume, 3.0, 5.0e-13,
                "非凸 L 形棱柱固体体积");
    expect_near(cut.fluid_volume, 1.0, 5.0e-13,
                "非凸 L 形棱柱缺角流体体积");
    expect_vec_near(cut.solid_centroid, {5.0 / 6.0, 5.0 / 6.0, 0.5},
                    5.0e-13, "非凸 L 形棱柱固体质心");
    expect_vec_near(cut.fluid_centroid, {1.5, 1.5, 0.5}, 5.0e-13,
                    "非凸 L 形棱柱流体质心");
    expect_near(cut.embedded_boundary_area, 2.0, 5.0e-13,
                "只计入与单元内流体相邻的两个凹角壁面");
    expect_near(cut.solid_cartesian_faces[0].area, 2.0, 5.0e-13,
                "L 棱柱 x- 固体占据面积");
    expect_near(cut.solid_cartesian_faces[1].area, 1.0, 5.0e-13,
                "L 棱柱 x+ 固体占据面积");
    expect_near(cut.solid_cartesian_faces[4].area, 3.0, 5.0e-13,
                "L 棱柱 z- 固体占据面积");
    expect_near(cut.volume_conservation_residual, 0.0, 1.0e-14,
                "非凸 L 棱柱单元体积闭合");
}

void test_general_nonconvex_uniform_mesh() {
    const cartmesh::SurfaceMesh surface(l_prism_triangles(),
                                        cartmesh::SurfaceFormat::ascii_stl,
                                        "nonconvex_l_prism_grid");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 47);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({0.0, 0.0, 0.0}, {2.0, 2.0, 1.0}), 1, 1, 1);
    const auto mesh = cartmesh::build_triangulated_cut_cell_mesh(grid, cutter);
    expect(mesh.fluid_cells.size() == 1 && mesh.cut_cell_count == 1,
           "非凸 L 棱柱必须在均匀网格中产生一个真实 Cut-cell");
    expect_near(mesh.total_fluid_volume, 1.0, 5.0e-13,
                "非凸网格总流体体积");
    expect_near(mesh.total_embedded_boundary_area, 2.0, 5.0e-13,
                "非凸网格嵌入壁面积");
    expect(mesh.nonclosed_cell_count == 0 &&
               mesh.shared_face_mismatch_count == 0,
           "非凸 Cut-cell 的边界表示必须闭合且拓扑一致");
    expect(mesh.component_analysis_pending_cell_count == 0 &&
               mesh.disconnected_fluid_cell_count == 0 &&
               mesh.fluid_cells.front().fluid_component_count == 1,
           "局部三角平面排列必须证明 L 形缺角流体只有一个连通分量");
    const auto& cell = mesh.fluid_cells.front();
    expect_near(cell.cartesian_faces[1].area, 1.0, 5.0e-13,
                "非凸 L 缺角的 x+ 开口面积");
    expect_near(cell.cartesian_faces[3].area, 1.0, 5.0e-13,
                "非凸 L 缺角的 y+ 开口面积");
    expect_near(cell.cartesian_faces[4].area, 1.0, 5.0e-13,
                "非凸 L 缺角的 z- 开口面积");
}

void test_general_surface_triangle_reordering_is_byte_stable() {
    auto reordered = l_prism_triangles();
    std::reverse(reordered.begin(), reordered.end());
    for (std::size_t index = 0; index < reordered.size(); ++index) {
        if ((index & 1U) == 0U) {
            const auto vertices = reordered[index].vertices();
            reordered[index] =
                cartmesh::Triangle(vertices[1], vertices[2], vertices[0]);
        }
    }
    const cartmesh::SurfaceMesh first_surface(
        l_prism_triangles(), cartmesh::SurfaceFormat::ascii_stl,
        "l_prism_original_order");
    const cartmesh::SurfaceMesh second_surface(
        std::move(reordered), cartmesh::SurfaceFormat::ascii_stl,
        "l_prism_reordered");
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-0.2, -0.2, -0.2}, {2.2, 2.2, 1.2}), 4, 4, 4);
    const auto first = cartmesh::build_triangulated_cut_cell_mesh(
        grid, cartmesh::TriangulatedSurfaceCutter(first_surface, 47));
    const auto second = cartmesh::build_triangulated_cut_cell_mesh(
        grid, cartmesh::TriangulatedSurfaceCutter(second_surface, 47));
    const auto temporary = std::filesystem::temp_directory_path();
    const auto first_path = temporary / "cartmesh_stage3_order_first.json";
    const auto second_path = temporary / "cartmesh_stage3_order_second.json";
    cartmesh::write_cut_cell_geometry_json(first_path, grid, first, true);
    cartmesh::write_cut_cell_geometry_json(second_path, grid, second, true);
    const auto read_all = [](const std::filesystem::path& path) {
        std::ifstream input(path, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(input),
                           std::istreambuf_iterator<char>());
    };
    const auto first_bytes = read_all(first_path);
    const auto second_bytes = read_all(second_path);
    std::filesystem::remove(first_path);
    std::filesystem::remove(second_path);
    expect(first_bytes == second_bytes,
           "三角片顺序和循环顶点编号不得改变通用 Cut-cell 几何字节");
}

void test_general_arrangement_detects_disconnected_slab() {
    const cartmesh::SurfaceMesh slab(
        axis_aligned_box_triangles({0.4, -1.0, -1.0}, {0.6, 2.0, 2.0}),
        cartmesh::SurfaceFormat::ascii_stl, "general_cell_spanning_slab");
    const cartmesh::TriangulatedSurfaceCutter cutter(slab);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 1, 1, 1);
    const auto mesh = cartmesh::build_triangulated_cut_cell_mesh(grid, cutter);
    expect(mesh.component_analysis_pending_cell_count == 0,
           "薄板案例的局部平面排列不得留下未决分量");
    expect(mesh.disconnected_fluid_cell_count == 1 &&
               mesh.fluid_cells.front().fluid_component_count == 2,
           "通用局部平面排列必须检测薄板左右两个流体分量");
    expect_near(mesh.total_fluid_volume, 0.8, 2.0e-13,
                "通用薄板两个流体分量总体积");
    std::array<bool, 2> component_seen{false, false};
    double piece_volume = 0.0;
    for (const auto& piece : mesh.fluid_cells.front().fluid_polyhedron_pieces) {
        expect(piece.component_id < component_seen.size(),
               "薄板流体分解的分量 ID 必须连续");
        component_seen[piece.component_id] = true;
        expect(piece.geometry.closed && piece.geometry.positive_volume,
               "每个不连通流体片必须输出闭合正体积多面体");
        piece_volume += piece.geometry.volume;
    }
    expect(component_seen[0] && component_seen[1],
           "不连通薄板案例必须显式输出两个流体分量");
    expect_near(piece_volume, mesh.total_fluid_volume, 2.0e-13,
                "不连通流体多面体片体积和");
}

void test_general_adaptive_octree_topology() {
    const cartmesh::SurfaceMesh surface(
        axis_aligned_box_triangles({0.6, 0.6, 0.6}, {1.4, 1.4, 1.4}),
        cartmesh::SurfaceFormat::ascii_stl, "adaptive_embedded_cube");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 59);
    cartmesh::LinearOctree tree(
        cartmesh::AABB({0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}), 1, 3);
    const auto refinement = tree.refine_to_desired_levels(
        [&](cartmesh::OctreeNodeCode, const cartmesh::AABB& box) {
            return static_cast<std::uint8_t>(
                cutter.bvh().intersects_surface(box) ? 3U : 1U);
        });
    const auto balancing = tree.balance_faces_2_to_1();
    expect(refinement.split_count > 0 && balancing.final_leaf_count == tree.leaf_count(),
           "自适应案例必须实际执行表面细化与 2:1 平衡");
    const auto levels = tree.level_statistics();
    expect(levels.minimum_leaf_level < levels.maximum_leaf_level &&
               tree.check_face_balance().balanced && tree.validate_partition(),
           "自适应案例必须保留多个层级、2:1 平衡和完整空间分区");

    const auto mesh = cartmesh::build_triangulated_cut_cell_mesh(tree, cutter);
    expect_near(mesh.total_fluid_volume, 8.0 - 0.8 * 0.8 * 0.8, 2.0e-12,
                "自适应网格总流体体积守恒");
    expect_near(mesh.total_embedded_boundary_area, 6.0 * 0.8 * 0.8, 2.0e-12,
                "自适应网格嵌入边界面积守恒");
    expect(mesh.nonclosed_cell_count == 0 &&
               mesh.negative_volume_cell_count == 0 &&
               mesh.shared_face_mismatch_count == 0 &&
               mesh.classification_conflict_count == 0 &&
               mesh.component_analysis_pending_cell_count == 0,
           "自适应 Cut-cell 必须保持闭合、正体积、共享面一致和确定分类");

    std::size_t coarse_fine_connection_count = 0;
    for (const auto& connection : mesh.internal_faces) {
        const auto first = cartmesh::decode_octree_node(
            tree.leaf_code(connection.first_background_cell_id));
        const auto second = cartmesh::decode_octree_node(
            tree.leaf_code(connection.second_background_cell_id));
        if (first.level != second.level) ++coarse_fine_connection_count;
    }
    expect(coarse_fine_connection_count > 0,
           "自适应 cell-face-neighbor 拓扑必须显式包含粗细层级连接");

    const auto first_case = std::filesystem::temp_directory_path() /
                            "cartmesh_stage61_adaptive_cut_first";
    const auto second_case = std::filesystem::temp_directory_path() /
                             "cartmesh_stage61_adaptive_cut_second";
    std::filesystem::remove_all(first_case);
    std::filesystem::remove_all(second_case);
    cartmesh::write_openfoam_poly_mesh(first_case, tree, mesh);
    cartmesh::write_openfoam_poly_mesh(second_case, tree, mesh);
    std::uint64_t solver_cell_count = 0;
    for (const auto& cell : mesh.fluid_cells) {
        solver_cell_count += cell.cut ? cell.fluid_polyhedron_pieces.size() : 1U;
    }
    expect_openfoam_cell_labels(first_case, solver_cell_count);
    for (const std::string name : {"points", "faces", "owner", "neighbour",
                                   "boundary", "cartmeshCellMapping.json"}) {
        expect(read_all(first_case / "constant/polyMesh" / name) ==
                   read_all(second_case / "constant/polyMesh" / name),
               "自适应 Cut-cell OpenFOAM 核心五文件及映射必须确定性一致");
    }
    std::filesystem::remove_all(first_case);
    std::filesystem::remove_all(second_case);
}

void test_openfoam_adaptive_coarse_fine_without_cut() {
    const cartmesh::SurfaceMesh outside_surface(
        axis_aligned_box_triangles({2.0, 2.0, 2.0}, {3.0, 3.0, 3.0}),
        cartmesh::SurfaceFormat::ascii_stl, "outside_domain_box");
    const cartmesh::TriangulatedSurfaceCutter cutter(outside_surface, 0);
    cartmesh::LinearOctree tree(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 1, 3);
    expect(tree.refine_leaf(cartmesh::encode_octree_node(1, 0, 0, 0)),
           "粗细交界测试必须细分一个基础叶");
    expect(tree.check_face_balance().balanced && tree.validate_partition(),
           "粗细交界测试树必须保持 2:1 与完整分区");
    const auto mesh = cartmesh::build_triangulated_cut_cell_mesh(tree, cutter);
    expect(mesh.cut_cell_count == 0 &&
               mesh.full_fluid_cell_count == tree.leaf_count(),
           "无切割粗细案例必须保留所有 fluid leaves");
    std::size_t coarse_fine_connections = 0;
    for (const auto& connection : mesh.internal_faces) {
        if (cartmesh::decode_octree_node(
                tree.leaf_code(connection.first_background_cell_id)).level !=
            cartmesh::decode_octree_node(
                tree.leaf_code(connection.second_background_cell_id)).level) {
            ++coarse_fine_connections;
        }
    }
    expect(coarse_fine_connections > 0,
           "无切割案例必须真实包含 coarse-fine connections");
    const auto case_directory = std::filesystem::temp_directory_path() /
                                "cartmesh_stage61_coarse_fine_no_cut";
    std::filesystem::remove_all(case_directory);
    cartmesh::write_openfoam_poly_mesh(case_directory, tree, mesh);
    expect_openfoam_cell_labels(case_directory, tree.leaf_count());
    const auto boundary = read_all(
        case_directory / "constant/polyMesh/boundary");
    expect(boundary.find("farfield") != std::string::npos,
           "无切割自适应 OpenFOAM 必须写出远场 patch");
    const auto mapping = read_all(
        case_directory / "constant/polyMesh/cartmeshCellMapping.json");
    expect(mapping.find("\"backgroundStableIdKind\": \"octree_node_code\"") !=
               std::string::npos &&
               mapping.find("\"solverCellCount\": 15") != std::string::npos,
           "自适应 OpenFOAM 必须保存 Morton leaf 到 solver cell 的映射");
    std::filesystem::remove_all(case_directory);
}

void test_nearly_blocked_shared_face_uses_conserved_moments() {
    constexpr double gap = 1.0e-7;
    const cartmesh::SurfaceMesh surface(
        axis_aligned_box_triangles({0.25, gap, gap},
                                   {0.75, 1.0 - gap, 1.0 - gap}),
        cartmesh::SurfaceFormat::ascii_stl,
        "nearly_blocked_shared_face");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 2, 1, 1);
    const auto mesh =
        cartmesh::build_triangulated_cut_cell_mesh(grid, cutter);

    expect(mesh.internal_faces.size() == 1,
           "几乎封死的共享面仍必须保留唯一流体连接");
    expect(mesh.shared_face_mismatch_count == 0,
           "微小开口应依据面积和一阶矩守恒，不能被病态派生质心误拒");
    expect(mesh.maximum_shared_face_area_mismatch < 1.0e-14,
           "微小开口两侧面积必须一致");
    expect(mesh.maximum_shared_face_first_moment_mismatch < 1.0e-14,
           "微小开口两侧一阶矩必须一致");
}

void test_polyhedron_vtu_compacts_unused_vertices() {
    auto box = cartmesh::make_box_polyhedron(
        cartmesh::AABB({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}));
    const auto box_geometry = cartmesh::measure_polyhedron(box);
    auto polyhedron = cartmesh::make_tetrahedron_polyhedron(
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0});
    const auto geometry = cartmesh::measure_polyhedron(polyhedron);
    polyhedron.vertices.push_back({9.0, 9.0, 9.0});
    cartmesh::FluidCellGeometry box_cell;
    box_cell.background_cell_id = 6;
    box_cell.fluid_polyhedron_pieces.push_back(
        {std::move(box), box_geometry, 0});
    cartmesh::FluidCellGeometry tetra_cell;
    tetra_cell.background_cell_id = 7;
    tetra_cell.fluid_polyhedron_pieces.push_back(
        {std::move(polyhedron), geometry, 0});
    cartmesh::ConvexCutCellMesh mesh;
    mesh.fluid_cells.push_back(std::move(box_cell));
    mesh.fluid_cells.push_back(std::move(tetra_cell));
    const auto path = std::filesystem::temp_directory_path() /
                      "cartmesh_stage3_unused_vertex_polyhedron.vtu";
    cartmesh::write_fluid_polyhedra_vtu(path, mesh);
    std::ifstream input(path);
    const std::string xml((std::istreambuf_iterator<char>(input)),
                          std::istreambuf_iterator<char>());
    std::filesystem::remove(path);
    expect(xml.find("NumberOfPoints=\"12\"") != std::string::npos,
           "显式 polyhedron VTU 不得把未被面引用的旧顶点写入 connectivity");
    expect(xml.find("Name=\"offsets\" format=\"ascii\">\n4\n12\n") !=
               std::string::npos,
           "不同顶点数的 polyhedron 必须按外部 meshio 分块顺序稳定写出");
}

void test_sphere_volume_area_convergence() {
    const cartmesh::AABB domain({-1.1, -1.1, -1.1}, {1.1, 1.1, 1.1});
    const double exact_volume = 4.0 * std::numbers::pi / 3.0;
    const double exact_area = 4.0 * std::numbers::pi;
    double previous_volume_error = std::numeric_limits<double>::infinity();
    double previous_area_error = std::numeric_limits<double>::infinity();
    for (std::uint32_t level = 0; level <= 2; ++level) {
        const cartmesh::SurfaceMesh surface(
            subdivided_unit_sphere(level), cartmesh::SurfaceFormat::ascii_stl,
            "subdivided_unit_sphere");
        const cartmesh::ConvexSurfaceCutter cutter(surface);
        const auto cut = cutter.cut_box(domain);
        const double volume_error = std::abs(cut.solid_geometry.volume - exact_volume);
        const double area_error = std::abs(cut.embedded_boundary_area - exact_area);
        expect(volume_error < previous_volume_error,
               "球面三角片细化必须单调减小体积误差");
        expect(area_error < previous_area_error,
               "球面三角片细化必须单调减小表面积误差");
        previous_volume_error = volume_error;
        previous_area_error = area_error;
    }

    const cartmesh::SurfaceMesh surface(
        subdivided_unit_sphere(1), cartmesh::SurfaceFormat::ascii_stl,
        "grid_invariant_unit_sphere");
    const cartmesh::ConvexSurfaceCutter cutter(surface);
    const auto reference = cutter.cut_box(domain);
    for (const std::uint32_t resolution : {4U, 6U}) {
        const cartmesh::UniformCartesianGrid grid(
            domain, resolution, resolution, resolution);
        const auto mesh = cartmesh::build_convex_cut_cell_mesh(grid, cutter);
        const double reconstructed_solid = domain.volume() - mesh.total_fluid_volume;
        expect_near(reconstructed_solid, reference.solid_geometry.volume, 2.0e-12,
                    "固定 STL 的总 Cut-cell 体积不得依赖背景网格分辨率");
        expect_near(mesh.total_embedded_boundary_area,
                    reference.embedded_boundary_area, 2.0e-12,
                    "固定 STL 的嵌入边界总面积不得依赖背景网格分辨率");
        expect(mesh.nonclosed_cell_count == 0 &&
                   mesh.shared_face_mismatch_count == 0,
               "球体多面体网格必须保持逐单元闭合和共享面一致");
    }
}

void test_openfoam_merges_coplanar_wall_triangles() {
    const cartmesh::SurfaceMesh surface(
        cube_triangles(), cartmesh::SurfaceFormat::ascii_stl,
        "openfoam_coplanar_cube");
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 0);
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-0.1, -0.1, -0.1}, {1.1, 1.1, 1.1}), 8, 8, 8);
    const auto mesh = cartmesh::build_triangulated_cut_cell_mesh(grid, cutter);
    const auto case_directory = std::filesystem::temp_directory_path() /
                                "cartmesh_stage3_openfoam_coplanar_merge";
    std::filesystem::remove_all(case_directory);
    cartmesh::write_openfoam_poly_mesh(case_directory, grid, mesh);
    const auto read_all = [](const std::filesystem::path& path) {
        std::ifstream input(path);
        return std::string(std::istreambuf_iterator<char>(input),
                           std::istreambuf_iterator<char>());
    };
    const auto faces = read_all(case_directory / "constant/polyMesh/faces");
    const auto boundary =
        read_all(case_directory / "constant/polyMesh/boundary");
    std::filesystem::remove_all(case_directory);
    expect(faces.find("\n1956\n(\n") != std::string::npos,
           "共面 STL 三角片必须合并为单个 OpenFOAM 壁面多边形");
    expect(boundary.find("boundary_0\n{\n    type wall;\n    nFaces 384;") !=
               std::string::npos,
           "立方体壁面不得保留导致伪凹单元的三角形内部分割");
}

} // 匿名命名空间

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"背景盒多面体几何", test_box_geometry},
        {"四面体构造绕序", test_tetrahedron_constructor_orientation},
        {"轴对齐切割与互补体积", test_axis_aligned_cut_and_complement},
        {"倾斜平面四面体切割", test_oblique_tetrahedron_cut},
        {"微小与大平移 Cut-cell", test_small_and_translated_cuts},
        {"空集与不变裁剪", test_empty_and_unchanged_clips},
        {"凸 STL 立方体 Cut-cell", test_convex_stl_cube_cut_cells},
        {"凸 STL 多表面片重建", test_convex_stl_multiple_surface_patches},
        {"均匀网格 Cut-cell 拓扑", test_uniform_grid_cut_cell_topology},
        {"均匀网格体积分数端点", test_uniform_grid_endpoint_fractions},
        {"不连通流体片检测", test_disconnected_fluid_pieces_are_detected},
        {"网格共面壁面", test_grid_aligned_surface_becomes_coplanar_wall},
        {"通用三角链凸立方体", test_general_triangle_chain_matches_convex_cube},
        {"通用单元多 patch ID", test_general_multiple_boundary_patches_in_one_cell},
        {"通用非凸 L 形棱柱", test_general_nonconvex_l_prism},
        {"通用非凸均匀网格", test_general_nonconvex_uniform_mesh},
        {"通用三角片重排确定性", test_general_surface_triangle_reordering_is_byte_stable},
        {"通用平面排列不连通检测", test_general_arrangement_detects_disconnected_slab},
        {"通用自适应八叉树拓扑", test_general_adaptive_octree_topology},
        {"自适应粗细交界 OpenFOAM", test_openfoam_adaptive_coarse_fine_without_cut},
        {"微小开口共享面守恒量", test_nearly_blocked_shared_face_uses_conserved_moments},
        {"多面体 VTU 紧凑顶点", test_polyhedron_vtu_compacts_unused_vertices},
        {"球体体积与面积收敛", test_sphere_volume_area_convergence},
        {"OpenFOAM 共面壁面合并", test_openfoam_merges_coplanar_wall_triangles},
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
    std::cout << "阶段3几何内核测试数=" << tests.size() << " 失败数=" << failures
              << '\n';
    return failures == 0 ? 0 : 1;
}
