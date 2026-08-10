#include "cartmesh/io/VtkWriter.hpp"

#include <fstream>
#include <iomanip>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cartmesh {
namespace {

[[nodiscard]] std::uint64_t point_id(std::uint32_t i, std::uint32_t j, std::uint32_t k,
                                     std::uint32_t nx, std::uint32_t ny) noexcept {
    return static_cast<std::uint64_t>(i) +
           (static_cast<std::uint64_t>(nx) + 1) *
               (static_cast<std::uint64_t>(j) +
                (static_cast<std::uint64_t>(ny) + 1) * static_cast<std::uint64_t>(k));
}

[[nodiscard]] std::string xml_escape(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (const char character : value) {
        switch (character) {
        case '&':
            escaped += "&amp;";
            break;
        case '<':
            escaped += "&lt;";
            break;
        case '>':
            escaped += "&gt;";
            break;
        case '\"':
            escaped += "&quot;";
            break;
        case '\'':
            escaped += "&apos;";
            break;
        default:
            escaped += character;
        }
    }
    return escaped;
}

struct LatticePoint {
    std::uint32_t x{};
    std::uint32_t y{};
    std::uint32_t z{};

    [[nodiscard]] bool operator==(const LatticePoint&) const noexcept = default;
};

struct LatticePointHash {
    [[nodiscard]] std::size_t operator()(const LatticePoint& point) const noexcept {
        std::uint64_t hash = 14695981039346656037ULL;
        const auto add = [&](std::uint32_t value) {
            for (std::uint32_t shift = 0; shift < 32U; shift += 8U) {
                hash ^= static_cast<std::uint8_t>((value >> shift) & 0xffU);
                hash *= 1099511628211ULL;
            }
        };
        add(point.x);
        add(point.y);
        add(point.z);
        return static_cast<std::size_t>(hash);
    }
};

} // 匿名命名空间

void write_vtu(const std::filesystem::path& path, const UniformCartesianGrid& grid,
               const std::vector<VtkCellData>& cell_data) {
    for (const auto& field : cell_data) {
        if (field.name.empty()) {
            throw std::invalid_argument("VTK 单元数据名称不得为空");
        }
        if (field.values.size() != grid.cell_count()) {
            throw std::invalid_argument("VTK 单元数据长度与网格单元数不一致");
        }
    }

    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("无法打开 VTK 输出文件：" + path.string());
    }
    output << std::setprecision(17);
    output << "<?xml version=\"1.0\"?>\n";
    output << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
              "byte_order=\"LittleEndian\">\n";
    output << "  <UnstructuredGrid>\n";
    output << "    <Piece NumberOfPoints=\"" << grid.point_count() << "\" NumberOfCells=\""
           << grid.cell_count() << "\">\n";

    output << "      <Points>\n";
    output << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n";
    const auto origin = grid.domain().minimum();
    const auto step = grid.spacing();
    for (std::uint32_t k = 0; k <= grid.nz(); ++k) {
        for (std::uint32_t j = 0; j <= grid.ny(); ++j) {
            for (std::uint32_t i = 0; i <= grid.nx(); ++i) {
                output << origin.x + step.x * static_cast<double>(i) << ' '
                       << origin.y + step.y * static_cast<double>(j) << ' '
                       << origin.z + step.z * static_cast<double>(k) << '\n';
            }
        }
    }
    output << "        </DataArray>\n";
    output << "      </Points>\n";

    output << "      <Cells>\n";
    output << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::uint32_t k = 0; k < grid.nz(); ++k) {
        for (std::uint32_t j = 0; j < grid.ny(); ++j) {
            for (std::uint32_t i = 0; i < grid.nx(); ++i) {
                output << point_id(i, j, k, grid.nx(), grid.ny()) << ' '
                       << point_id(i + 1, j, k, grid.nx(), grid.ny()) << ' '
                       << point_id(i + 1, j + 1, k, grid.nx(), grid.ny()) << ' '
                       << point_id(i, j + 1, k, grid.nx(), grid.ny()) << ' '
                       << point_id(i, j, k + 1, grid.nx(), grid.ny()) << ' '
                       << point_id(i + 1, j, k + 1, grid.nx(), grid.ny()) << ' '
                       << point_id(i + 1, j + 1, k + 1, grid.nx(), grid.ny()) << ' '
                       << point_id(i, j + 1, k + 1, grid.nx(), grid.ny()) << '\n';
            }
        }
    }
    output << "        </DataArray>\n";
    output << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::uint64_t cell = 1; cell <= grid.cell_count(); ++cell) {
        output << cell * 8 << '\n';
    }
    output << "        </DataArray>\n";
    output << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::uint64_t cell = 0; cell < grid.cell_count(); ++cell) {
        output << "12\n";
    }
    output << "        </DataArray>\n";
    output << "      </Cells>\n";

    if (!cell_data.empty()) {
        output << "      <CellData>\n";
        for (const auto& field : cell_data) {
            output << "        <DataArray type=\"Float64\" Name=\"" << xml_escape(field.name)
                   << "\" format=\"ascii\">\n";
            for (const double value : field.values) {
                output << value << '\n';
            }
            output << "        </DataArray>\n";
        }
        output << "      </CellData>\n";
    }

    output << "    </Piece>\n";
    output << "  </UnstructuredGrid>\n";
    output << "</VTKFile>\n";
    output.flush();
    if (!output) {
        throw std::runtime_error("写入 VTK 输出文件失败：" + path.string());
    }
}

void write_octree_vtu(const std::filesystem::path& path, const LinearOctree& tree,
                      const std::vector<VtkCellData>& cell_data) {
    for (const auto& field : cell_data) {
        if (field.name.empty()) {
            throw std::invalid_argument("VTK 单元数据名称不得为空");
        }
        if (field.values.size() != tree.leaf_count()) {
            throw std::invalid_argument("VTK 单元数据长度与八叉树叶数不一致");
        }
    }
    if (tree.leaf_count() >
        static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max() / 8U)) {
        throw std::length_error("八叉树 VTK 连接数组超出容器范围");
    }
    std::vector<LatticePoint> points;
    std::vector<std::array<std::uint64_t, 8>> connectivity;
    const auto maximum_point_capacity = static_cast<std::size_t>(tree.leaf_count()) * 8U;
    points.reserve(maximum_point_capacity);
    connectivity.reserve(static_cast<std::size_t>(tree.leaf_count()));
    std::unordered_map<LatticePoint, std::uint64_t, LatticePointHash> point_ids;
    point_ids.reserve(maximum_point_capacity);
    constexpr std::array<std::array<std::uint8_t, 3>, 8> corners = {
        std::array<std::uint8_t, 3>{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    for (const auto code : tree.leaf_codes()) {
        const auto node = decode_octree_node(code);
        const auto shift = static_cast<std::uint8_t>(tree.maximum_level() - node.level);
        const std::uint32_t x0 = node.x << shift;
        const std::uint32_t y0 = node.y << shift;
        const std::uint32_t z0 = node.z << shift;
        const std::uint32_t span = 1U << shift;
        std::array<std::uint64_t, 8> cell{};
        for (std::size_t corner = 0; corner < corners.size(); ++corner) {
            const LatticePoint point{x0 + static_cast<std::uint32_t>(corners[corner][0]) * span,
                                     y0 + static_cast<std::uint32_t>(corners[corner][1]) * span,
                                     z0 + static_cast<std::uint32_t>(corners[corner][2]) * span};
            const auto [iterator, inserted] =
                point_ids.try_emplace(point, static_cast<std::uint64_t>(points.size()));
            if (inserted) {
                points.push_back(point);
            }
            cell[corner] = iterator->second;
        }
        connectivity.push_back(cell);
    }

    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("无法打开八叉树 VTK 输出文件：" + path.string());
    }
    output << std::setprecision(17);
    output << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
              "byte_order=\"LittleEndian\">\n"
           << "  <UnstructuredGrid>\n"
           << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\""
           << tree.leaf_count() << "\">\n"
           << "      <Points>\n"
           << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n";
    const auto root_minimum = tree.domain().minimum();
    const auto root_extent = tree.domain().extent();
    const double denominator = std::ldexp(1.0, tree.maximum_level());
    for (const auto point : points) {
        output << root_minimum.x + root_extent.x * static_cast<double>(point.x) / denominator
               << ' '
               << root_minimum.y + root_extent.y * static_cast<double>(point.y) / denominator
               << ' '
               << root_minimum.z + root_extent.z * static_cast<double>(point.z) / denominator
               << '\n';
    }
    output << "        </DataArray>\n"
           << "      </Points>\n"
           << "      <Cells>\n"
           << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& cell : connectivity) {
        for (std::size_t corner = 0; corner < cell.size(); ++corner) {
            output << cell[corner] << (corner + 1U == cell.size() ? '\n' : ' ');
        }
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::uint64_t cell = 1; cell <= tree.leaf_count(); ++cell) {
        output << cell * 8U << '\n';
    }
    output << "        </DataArray>\n"
           << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::uint64_t cell = 0; cell < tree.leaf_count(); ++cell) {
        output << "12\n";
    }
    output << "        </DataArray>\n"
           << "      </Cells>\n";
    if (!cell_data.empty()) {
        output << "      <CellData>\n";
        for (const auto& field : cell_data) {
            output << "        <DataArray type=\"Float64\" Name=\""
                   << xml_escape(field.name) << "\" format=\"ascii\">\n";
            for (const double value : field.values) {
                output << value << '\n';
            }
            output << "        </DataArray>\n";
        }
        output << "      </CellData>\n";
    }
    output << "    </Piece>\n"
           << "  </UnstructuredGrid>\n"
           << "</VTKFile>\n";
    output.flush();
    if (!output) {
        throw std::runtime_error("写入八叉树 VTK 输出文件失败：" + path.string());
    }
}

} // 命名空间 cartmesh
