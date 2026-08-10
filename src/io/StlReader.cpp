#include "cartmesh/io/StlReader.hpp"

#include <array>
#include <bit>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace cartmesh {
namespace {

[[nodiscard]] std::vector<std::uint8_t> read_bytes(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) {
        throw std::runtime_error("无法打开 STL 文件：" + path.string());
    }
    const auto end = input.tellg();
    if (end < 0) {
        throw std::runtime_error("无法读取 STL 文件大小：" + path.string());
    }
    const auto size = static_cast<std::uint64_t>(end);
    if (size > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("STL 文件大于当前进程可寻址内存");
    }
    std::vector<std::uint8_t> bytes(static_cast<std::size_t>(size));
    input.seekg(0);
    if (!bytes.empty()) {
        input.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    }
    if (!input) {
        throw std::runtime_error("读取 STL 文件失败：" + path.string());
    }
    return bytes;
}

[[nodiscard]] std::uint32_t little_u32(const std::uint8_t* bytes) noexcept {
    return static_cast<std::uint32_t>(bytes[0]) |
           (static_cast<std::uint32_t>(bytes[1]) << 8U) |
           (static_cast<std::uint32_t>(bytes[2]) << 16U) |
           (static_cast<std::uint32_t>(bytes[3]) << 24U);
}

[[nodiscard]] float little_f32(const std::uint8_t* bytes) noexcept {
    return std::bit_cast<float>(little_u32(bytes));
}

[[nodiscard]] bool has_binary_layout(const std::vector<std::uint8_t>& bytes) noexcept {
    if (bytes.size() < 84) {
        return false;
    }
    const auto count = static_cast<std::uint64_t>(little_u32(bytes.data() + 80));
    return count <= (std::numeric_limits<std::uint64_t>::max() - 84ULL) / 50ULL &&
           84ULL + count * 50ULL == bytes.size();
}

[[nodiscard]] std::string trimmed(std::string_view line) {
    std::size_t begin = 0;
    while (begin < line.size() && std::isspace(static_cast<unsigned char>(line[begin])) != 0) {
        ++begin;
    }
    std::size_t end = line.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(line[end - 1])) != 0) {
        --end;
    }
    return std::string(line.substr(begin, end - begin));
}

[[nodiscard]] bool starts_with_keyword(const std::string& line, std::string_view keyword) {
    return line == keyword ||
           (line.size() > keyword.size() && line.compare(0, keyword.size(), keyword) == 0 &&
            std::isspace(static_cast<unsigned char>(line[keyword.size()])) != 0);
}

[[nodiscard]] Vec3 parse_vec3_line(const std::string& line, std::string_view keyword,
                                   std::size_t line_number) {
    std::istringstream parser(line);
    std::string parsed_keyword;
    Vec3 value;
    std::string extra;
    if (!(parser >> parsed_keyword >> value.x >> value.y >> value.z) || parsed_keyword != keyword ||
        (parser >> extra)) {
        throw std::runtime_error("ASCII STL 第 " + std::to_string(line_number) +
                                 " 行格式无效，期望 " + std::string(keyword) + " x y z");
    }
    if (!is_finite(value)) {
        throw std::runtime_error("ASCII STL 第 " + std::to_string(line_number) +
                                 " 行包含非有限坐标");
    }
    return value;
}

[[nodiscard]] SurfaceMesh parse_ascii(const std::vector<std::uint8_t>& bytes) {
    if (bytes.empty()) {
        throw std::runtime_error("STL 文件为空");
    }
    if (std::memchr(bytes.data(), '\0', bytes.size()) != nullptr) {
        throw std::runtime_error("STL 文件既不符合二进制长度，也不是有效 ASCII 文本");
    }
    const std::string text(bytes.begin(), bytes.end());
    std::istringstream input(text);
    std::vector<Triangle> triangles;
    std::string name;
    std::string line;
    std::size_t line_number = 0;
    bool saw_solid = false;
    bool saw_endsolid = false;
    while (std::getline(input, line)) {
        ++line_number;
        line = trimmed(line);
        if (line.empty()) {
            continue;
        }
        if (!saw_solid) {
            if (!starts_with_keyword(line, "solid")) {
                throw std::runtime_error("ASCII STL 必须以 solid 行开始");
            }
            saw_solid = true;
            name = trimmed(std::string_view(line).substr(5));
            continue;
        }
        if (starts_with_keyword(line, "endsolid")) {
            saw_endsolid = true;
            while (std::getline(input, line)) {
                ++line_number;
                if (!trimmed(line).empty()) {
                    throw std::runtime_error("ASCII STL 的 endsolid 之后存在非空内容");
                }
            }
            break;
        }
        if (!starts_with_keyword(line, "facet")) {
            throw std::runtime_error("ASCII STL 第 " + std::to_string(line_number) +
                                     " 行期望 facet normal");
        }
        {
            std::istringstream parser(line);
            std::string facet;
            std::string normal;
            Vec3 ignored_normal;
            std::string extra;
            if (!(parser >> facet >> normal >> ignored_normal.x >> ignored_normal.y >>
                  ignored_normal.z) ||
                facet != "facet" || normal != "normal" || (parser >> extra) ||
                !is_finite(ignored_normal)) {
                throw std::runtime_error("ASCII STL 第 " + std::to_string(line_number) +
                                         " 行 facet normal 格式无效");
            }
        }
        auto next_nonempty = [&]() -> std::string {
            while (std::getline(input, line)) {
                ++line_number;
                auto value = trimmed(line);
                if (!value.empty()) {
                    return value;
                }
            }
            throw std::runtime_error("ASCII STL 在三角形记录中意外结束");
        };
        if (next_nonempty() != "outer loop") {
            throw std::runtime_error("ASCII STL 第 " + std::to_string(line_number) +
                                     " 行期望 outer loop");
        }
        const auto a = parse_vec3_line(next_nonempty(), "vertex", line_number);
        const auto b = parse_vec3_line(next_nonempty(), "vertex", line_number);
        const auto c = parse_vec3_line(next_nonempty(), "vertex", line_number);
        if (next_nonempty() != "endloop") {
            throw std::runtime_error("ASCII STL 第 " + std::to_string(line_number) +
                                     " 行期望 endloop");
        }
        if (next_nonempty() != "endfacet") {
            throw std::runtime_error("ASCII STL 第 " + std::to_string(line_number) +
                                     " 行期望 endfacet");
        }
        triangles.emplace_back(a, b, c);
    }
    if (!saw_solid || !saw_endsolid) {
        throw std::runtime_error("ASCII STL 缺少 solid 或 endsolid");
    }
    return SurfaceMesh(std::move(triangles), SurfaceFormat::ascii_stl, std::move(name));
}

[[nodiscard]] SurfaceMesh parse_binary(const std::vector<std::uint8_t>& bytes) {
    const auto count = little_u32(bytes.data() + 80);
    std::vector<Triangle> triangles;
    triangles.reserve(count);
    for (std::uint32_t triangle_index = 0; triangle_index < count; ++triangle_index) {
        const auto* record = bytes.data() + 84 + static_cast<std::size_t>(triangle_index) * 50;
        std::array<Vec3, 3> vertices;
        for (std::size_t vertex = 0; vertex < vertices.size(); ++vertex) {
            const auto* coordinates = record + 12 + vertex * 12;
            vertices[vertex] = {static_cast<double>(little_f32(coordinates)),
                                static_cast<double>(little_f32(coordinates + 4)),
                                static_cast<double>(little_f32(coordinates + 8))};
            if (!is_finite(vertices[vertex])) {
                throw std::runtime_error("二进制 STL 三角形 " +
                                         std::to_string(triangle_index) + " 包含非有限坐标");
            }
        }
        triangles.emplace_back(vertices[0], vertices[1], vertices[2]);
    }
    std::string header(reinterpret_cast<const char*>(bytes.data()), 80);
    while (!header.empty() && (header.back() == '\0' || header.back() == ' ')) {
        header.pop_back();
    }
    return SurfaceMesh(std::move(triangles), SurfaceFormat::binary_stl, std::move(header));
}

} // 匿名命名空间

SurfaceMesh read_stl(const std::filesystem::path& path) {
    const auto bytes = read_bytes(path);
    return has_binary_layout(bytes) ? parse_binary(bytes) : parse_ascii(bytes);
}

} // 命名空间 cartmesh
