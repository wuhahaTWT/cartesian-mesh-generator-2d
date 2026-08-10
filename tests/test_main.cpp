#include "cartmesh/geometry/AnalyticShapes.hpp"
#include "cartmesh/geometry/Triangle.hpp"
#include "cartmesh/grid/MortonCode.hpp"
#include "cartmesh/grid/UniformCartesianGrid.hpp"
#include "cartmesh/io/VtkWriter.hpp"

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

double center_sampled_sphere_volume(std::uint32_t resolution) {
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-1.25, -1.25, -1.25}, {1.25, 1.25, 1.25}), resolution, resolution,
        resolution);
    const cartmesh::AnalyticSphere sphere({0.0, 0.0, 0.0}, 1.0);
    std::uint64_t inside = 0;
    for (std::uint64_t id = 0; id < grid.cell_count(); ++id) {
        inside += static_cast<std::uint64_t>(sphere.contains(grid.cell_center(grid.cell_key(id))));
    }
    return static_cast<double>(inside) * grid.cell_volume();
}

void test_vec3_and_aabb() {
    const cartmesh::Vec3 a{1.0, 2.0, 3.0};
    const cartmesh::Vec3 b{-2.0, 1.0, 4.0};
    expect_near(cartmesh::dot(a, b), 12.0, 1.0e-14, "点积");
    const auto crossed = cartmesh::cross(a, b);
    expect_near(crossed.x, 5.0, 1.0e-14, "叉积 x 分量");
    expect_near(crossed.y, -10.0, 1.0e-14, "叉积 y 分量");
    expect_near(crossed.z, 5.0, 1.0e-14, "叉积 z 分量");

    const cartmesh::AABB box({-1.0, -2.0, -3.0}, {3.0, 2.0, 1.0});
    expect_near(box.volume(), 64.0, 1.0e-14, "AABB 体积");
    expect(box.contains({0.0, 0.0, 0.0}), "AABB 包含原点");
    expect(box.intersects(cartmesh::AABB({2.0, 1.0, 0.0}, {4.0, 3.0, 2.0})),
           "AABB 相交");
    expect_throw([] { cartmesh::AABB({1.0, 0.0, 0.0}, {0.0, 1.0, 1.0}); },
                 "坐标颠倒的 AABB 必须失败");
}

void test_triangle() {
    const cartmesh::Triangle triangle({0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 3.0, 0.0});
    expect_near(triangle.area(), 3.0, 1.0e-14, "三角形面积");
    const auto centroid = triangle.centroid();
    expect_near(centroid.x, 2.0 / 3.0, 1.0e-14, "三角形质心 x 分量");
    expect_near(centroid.y, 1.0, 1.0e-14, "三角形质心 y 分量");
    expect_near(triangle.bounds().volume(), 0.0, 1.0e-14, "平面三角形包围盒体积");
}

void test_analytic_shapes() {
    const cartmesh::AnalyticCube cube(cartmesh::AABB({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}));
    expect_near(cube.volume(), 8.0, 1.0e-14, "立方体体积");
    expect_near(cube.surface_area(), 24.0, 1.0e-14, "立方体表面积");
    expect(cube.contains({1.0, 0.0, 0.0}), "立方体包含其闭合边界");

    const cartmesh::AnalyticSphere sphere({1.0, 2.0, 3.0}, 2.0);
    expect_near(sphere.volume(), 32.0 * std::numbers::pi_v<double> / 3.0, 1.0e-13,
                "球体体积");
    expect_near(sphere.surface_area(), 16.0 * std::numbers::pi_v<double>, 1.0e-13,
                "球体表面积");
    expect_near(sphere.signed_distance({3.0, 2.0, 3.0}), 0.0, 1.0e-14,
                "球体有符号距离");
    expect_throw([] { cartmesh::AnalyticSphere({0.0, 0.0, 0.0}, 0.0); },
                 "零半径球体必须失败");
}

void test_morton_round_trip() {
    const std::vector<cartmesh::MortonCoordinates> cases = {
        {0, 0, 0}, {1, 2, 3}, {17, 31, 9}, {1024, 65535, 77},
        {cartmesh::max_morton_coordinate, cartmesh::max_morton_coordinate,
         cartmesh::max_morton_coordinate}};
    for (const auto coordinates : cases) {
        const auto decoded = cartmesh::decode_morton_3d(
            cartmesh::encode_morton_3d(coordinates.x, coordinates.y, coordinates.z));
        expect(decoded.x == coordinates.x && decoded.y == coordinates.y &&
                   decoded.z == coordinates.z,
               "Morton 编解码往返");
    }
    expect(cartmesh::encode_morton_3d(1, 0, 0) == 1, "Morton x 位");
    expect(cartmesh::encode_morton_3d(0, 1, 0) == 2, "Morton y 位");
    expect(cartmesh::encode_morton_3d(0, 0, 1) == 4, "Morton z 位");
    expect_throw(
        [] {
            static_cast<void>(
                cartmesh::encode_morton_3d(cartmesh::max_morton_coordinate + 1U, 0, 0));
        },
        "Morton 坐标溢出必须失败");
}

void test_uniform_grid() {
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({-1.0, -3.0, 2.0}, {3.0, 3.0, 10.0}), 2, 3, 4);
    expect(grid.cell_count() == 24, "网格单元数");
    expect(grid.point_count() == 60, "网格点数");
    expect_near(grid.cell_volume(), 8.0, 1.0e-14, "网格单元体积");
    expect_near(grid.cell_volume() * static_cast<double>(grid.cell_count()),
                grid.domain().volume(), 1.0e-13, "网格体积守恒");
    for (std::uint64_t id = 0; id < grid.cell_count(); ++id) {
        expect(grid.linear_id(grid.cell_key(id)) == id, "线性 ID 往返转换");
    }
    const auto first_center = grid.cell_center({0, 0, 0, 0});
    expect_near(first_center.x, 0.0, 1.0e-14, "首个单元中心 x 分量");
    expect_near(first_center.y, -2.0, 1.0e-14, "首个单元中心 y 分量");
    expect_near(first_center.z, 3.0, 1.0e-14, "首个单元中心 z 分量");
    expect_throw(
        [] {
            cartmesh::UniformCartesianGrid invalid(
                cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 0, 1, 1);
        },
        "网格任一方向为零时必须失败");
    expect_throw([&] { static_cast<void>(grid.cell_key(grid.cell_count())); },
                 "超出范围的线性 ID 必须失败");
}

void test_sphere_convergence() {
    const auto exact = 4.0 * std::numbers::pi_v<double> / 3.0;
    const auto coarse_error = std::abs(center_sampled_sphere_volume(12) - exact);
    const auto fine_error = std::abs(center_sampled_sphere_volume(48) - exact);
    expect(fine_error < coarse_error, "球体中心采样误差从 12³ 到 48³ 应当减小");
    expect(fine_error / exact < 0.015, "48³ 球体相对体积误差应低于 1.5%");
}

void test_vtu_writer() {
    const cartmesh::UniformCartesianGrid grid(
        cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 2, 1, 1);
    const auto path = std::filesystem::temp_directory_path() / "cartmesh_stage0_writer_test.vtu";
    cartmesh::write_vtu(path, grid, {{"region", {1.0, 0.0}}});
    std::ifstream input(path);
    const std::string contents((std::istreambuf_iterator<char>(input)),
                               std::istreambuf_iterator<char>());
    expect(contents.find("NumberOfPoints=\"12\"") != std::string::npos,
           "VTU 点数元数据");
    expect(contents.find("NumberOfCells=\"2\"") != std::string::npos,
           "VTU 单元数元数据");
    expect(contents.find("Name=\"connectivity\"") != std::string::npos,
           "VTU 连接数组");
    expect(contents.find("Name=\"region\"") != std::string::npos, "VTU 区域数据");
    std::filesystem::remove(path);
    expect_throw([&] { cartmesh::write_vtu(path, grid, {{"bad", {1.0}}}); },
                 "单元数据长度不一致时必须失败");
}

} // 匿名命名空间

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"三维向量与包围盒", test_vec3_and_aabb},
        {"三角形", test_triangle},
        {"解析几何", test_analytic_shapes},
        {"Morton 编解码", test_morton_round_trip},
        {"均匀网格", test_uniform_grid},
        {"球体收敛", test_sphere_convergence},
        {"VTU 写出器", test_vtu_writer},
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
    std::cout << "测试数=" << tests.size() << " 失败数=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
