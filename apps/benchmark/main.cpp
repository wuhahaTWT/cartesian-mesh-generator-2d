#include "cartmesh/geometry/AnalyticShapes.hpp"
#include "cartmesh/grid/UniformCartesianGrid.hpp"
#include "cartmesh/util/Performance.hpp"

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

[[nodiscard]] std::uint32_t parse_resolution(int argc, char** argv) {
    if (argc == 1) {
        return 128;
    }
    if (argc != 2) {
        throw std::invalid_argument("用法：cartmesh_benchmark [分辨率]");
    }
    std::size_t parsed = 0;
    const auto value = std::stoull(argv[1], &parsed);
    if (parsed != std::string(argv[1]).size() || value == 0 ||
        value > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument("分辨率必须是 32 位正整数");
    }
    return static_cast<std::uint32_t>(value);
}

} // 匿名命名空间

int main(int argc, char** argv) {
    try {
        const auto resolution = parse_resolution(argc, argv);
        const auto start = std::chrono::steady_clock::now();
        const cartmesh::UniformCartesianGrid grid(
            cartmesh::AABB({-1.25, -1.25, -1.25}, {1.25, 1.25, 1.25}), resolution,
            resolution, resolution);
        const cartmesh::AnalyticSphere sphere({0.0, 0.0, 0.0}, 1.0);
        std::uint64_t inside = 0;
        std::uint64_t checksum = 14695981039346656037ULL;
        for (std::uint64_t id = 0; id < grid.cell_count(); ++id) {
            const bool selected = sphere.contains(grid.cell_center(grid.cell_key(id)));
            inside += static_cast<std::uint64_t>(selected);
            checksum ^= static_cast<std::uint8_t>(selected);
            checksum *= 1099511628211ULL;
        }
        const auto end = std::chrono::steady_clock::now();
        const auto elapsed = std::chrono::duration<double>(end - start).count();
        const auto estimated_volume = static_cast<double>(inside) * grid.cell_volume();
        std::cout << std::setprecision(17) << "{\"stage\":0,\"case\":\"sphere_center_sample\""
                  << ",\"resolution\":" << resolution << ",\"cellCount\":"
                  << grid.cell_count() << ",\"insideCellCount\":" << inside
                  << ",\"uniformGridObjectBytes\":" << sizeof(grid)
                  << ",\"cellKeyBytes\":" << sizeof(cartmesh::CellKey)
                  << ",\"mortonCodeBytes\":" << sizeof(cartmesh::MortonCode)
                  << ",\"perCellCoreStorageBytes\":0"
                  << ",\"estimatedVolume\":" << estimated_volume << ",\"exactVolume\":"
                  << sphere.volume() << ",\"wallSeconds\":" << elapsed
                  << ",\"peakRssBytes\":" << cartmesh::peak_rss_bytes()
                  << ",\"threads\":1,\"checksum\":" << checksum << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "cartmesh_benchmark 错误：" << error.what() << '\n';
        return 1;
    }
}
