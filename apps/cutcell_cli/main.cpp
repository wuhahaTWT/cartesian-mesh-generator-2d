#include "cartmesh/cutcell/ConvexCutCellMesh.hpp"
#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"
#include "cartmesh/geometry/SurfaceDiagnostics.hpp"
#include "cartmesh/grid/OctreeRefinement.hpp"
#include "cartmesh/io/DiagnosticVtkWriter.hpp"
#include "cartmesh/io/CutCellVtkWriter.hpp"
#include "cartmesh/io/CutCellJsonWriter.hpp"
#include "cartmesh/io/OpenFoamWriter.hpp"
#include "cartmesh/io/StlReader.hpp"
#include "cartmesh/io/VtkWriter.hpp"
#include "cartmesh/quality/SolverMeshQuality.hpp"
#include "cartmesh/quality/SolverMeshStabilizer.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace {

struct BoundaryRange {
    std::uint64_t first_triangle{};
    std::uint64_t last_triangle{};
    std::uint64_t boundary_id{};
    std::string name;
};

struct RegionName {
    std::uint64_t region_id{};
    std::string name;
};

struct Options {
    std::filesystem::path stl;
    std::filesystem::path output{"artifacts/stage3_cut_cells.vtu"};
    std::filesystem::path boundary_output{"artifacts/stage3_embedded_boundary.vtp"};
    std::filesystem::path geometry_output{"artifacts/stage3_cut_cell_geometry.json"};
    std::filesystem::path polyhedra_output{"artifacts/stage3_fluid_polyhedra.vtu"};
    std::filesystem::path tetrahedra_output{"artifacts/stage4_fluid_tetrahedra.vtu"};
    std::filesystem::path report{"artifacts/stage3_cut_cells.json"};
    std::optional<std::filesystem::path> openfoam_case;
    std::optional<std::filesystem::path> quality_output;
    std::optional<std::filesystem::path> stabilization_output;
    std::uint32_t resolution{24};
    double padding_fraction{0.1};
    double small_cell_threshold{0.01};
    double geometric_tolerance{};
    bool adaptive{};
    std::uint8_t base_level{2};
    std::uint8_t maximum_level{6};
    std::optional<std::uint8_t> surface_level;
    cartmesh::OctreeRefinementConfiguration refinement;
    bool strict_gap_resolution{};
    bool stabilize{};
    std::uint32_t maximum_stabilization_rounds{2};
    std::vector<BoundaryRange> boundary_ranges;
    std::vector<RegionName> region_names;
    bool write_vtk{true};
};

[[nodiscard]] std::string json_escape(std::string_view value) {
    std::string result;
    for (const char character : value) {
        switch (character) {
        case '\\': result += "\\\\"; break;
        case '"': result += "\\\""; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default: result += character; break;
        }
    }
    return result;
}

[[nodiscard]] std::uint32_t parse_u32(std::string_view text) {
    std::size_t used = 0;
    const auto value = std::stoull(std::string(text), &used);
    if (used != text.size() || value == 0 || value > 0xffffffffULL) {
        throw std::invalid_argument("resolution 必须是正的 32 位整数");
    }
    return static_cast<std::uint32_t>(value);
}

[[nodiscard]] double parse_nonnegative(std::string_view text) {
    std::size_t used = 0;
    const double value = std::stod(std::string(text), &used);
    if (used != text.size() || !std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument("padding-fraction 必须是非负有限数");
    }
    return value;
}

[[nodiscard]] std::uint8_t parse_level(std::string_view text,
                                       std::string_view option) {
    std::size_t used = 0;
    const auto value = std::stoul(std::string(text), &used);
    if (used != text.size() || value > cartmesh::maximum_octree_level) {
        throw std::invalid_argument(std::string(option) +
                                    " 需要 0..21 的整数层级");
    }
    return static_cast<std::uint8_t>(value);
}

[[nodiscard]] std::uint32_t parse_positive_u32(std::string_view text,
                                               std::string_view option) {
    std::size_t used = 0;
    const auto value = std::stoull(std::string(text), &used);
    if (used != text.size() || value == 0 ||
        value > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument(std::string(option) + " 需要 32 位正整数");
    }
    return static_cast<std::uint32_t>(value);
}

[[nodiscard]] std::vector<double> parse_csv(std::string_view text,
                                            std::size_t count,
                                            std::string_view option) {
    std::vector<double> result;
    std::size_t begin = 0;
    while (begin <= text.size()) {
        const auto comma = text.find(',', begin);
        const auto token = text.substr(
            begin, comma == std::string_view::npos ? text.size() - begin
                                                   : comma - begin);
        std::size_t used = 0;
        const double value = std::stod(std::string(token), &used);
        if (used != token.size() || !std::isfinite(value)) {
            throw std::invalid_argument(std::string(option) + " 包含无效数值");
        }
        result.push_back(value);
        if (comma == std::string_view::npos) break;
        begin = comma + 1U;
    }
    if (result.size() != count) {
        throw std::invalid_argument(std::string(option) + " 数值个数不正确");
    }
    return result;
}

[[nodiscard]] std::pair<std::string_view, std::uint8_t> split_level_suffix(
    std::string_view text, std::string_view option) {
    const auto colon = text.rfind(':');
    if (colon == std::string_view::npos || colon == 0 ||
        colon + 1U == text.size()) {
        throw std::invalid_argument(std::string(option) +
                                    " 需要 values:level 格式");
    }
    return {text.substr(0, colon),
            parse_level(text.substr(colon + 1U), option)};
}

[[nodiscard]] std::uint64_t parse_u64(std::string_view text,
                                      std::string_view field) {
    std::size_t used = 0;
    const auto value = std::stoull(std::string(text), &used);
    if (used != text.size()) {
        throw std::invalid_argument(std::string(field) + " 必须是非负整数");
    }
    return value;
}

[[nodiscard]] BoundaryRange parse_boundary_range(std::string_view text) {
    std::array<std::string_view, 4> field{};
    std::size_t begin = 0;
    for (std::size_t index = 0; index < 3; ++index) {
        const std::size_t separator = text.find(':', begin);
        if (separator == std::string_view::npos) {
            throw std::invalid_argument(
                "boundary-range 需要 FIRST:LAST:ID:NAME 格式");
        }
        field[index] = text.substr(begin, separator - begin);
        begin = separator + 1U;
    }
    field[3] = text.substr(begin);
    BoundaryRange result{parse_u64(field[0], "boundary FIRST"),
                         parse_u64(field[1], "boundary LAST"),
                         parse_u64(field[2], "boundary ID"),
                         std::string(field[3])};
    if (result.first_triangle >= result.last_triangle || result.name.empty()) {
        throw std::invalid_argument(
            "boundary-range 必须是非空半开三角片范围且 NAME 非空");
    }
    return result;
}

[[nodiscard]] RegionName parse_region_name(std::string_view text) {
    const auto separator = text.find(':');
    if (separator == std::string_view::npos || separator + 1U == text.size()) {
        throw std::invalid_argument("region-name 需要 ID:NAME 格式");
    }
    return {parse_u64(text.substr(0, separator), "region ID"),
            std::string(text.substr(separator + 1U))};
}

[[nodiscard]] Options parse_options(int argc, char** argv) {
    Options options;
    for (int argument = 1; argument < argc; ++argument) {
        const std::string value = argv[argument];
        const auto next = [&]() -> std::string {
            if (++argument >= argc) {
                throw std::invalid_argument(value + " 缺少参数");
            }
            return argv[argument];
        };
        if (value == "--stl") {
            options.stl = next();
        } else if (value == "--resolution") {
            options.resolution = parse_u32(next());
        } else if (value == "--padding-fraction") {
            options.padding_fraction = parse_nonnegative(next());
        } else if (value == "--small-cell-threshold") {
            options.small_cell_threshold = parse_nonnegative(next());
        } else if (value == "--geometric-tolerance") {
            options.geometric_tolerance = parse_nonnegative(next());
        } else if (value == "--adaptive") {
            options.adaptive = true;
        } else if (value == "--base-level") {
            options.adaptive = true;
            options.base_level = parse_level(next(), value);
        } else if (value == "--max-level") {
            options.adaptive = true;
            options.maximum_level = parse_level(next(), value);
        } else if (value == "--surface-level") {
            options.adaptive = true;
            options.surface_level = parse_level(next(), value);
        } else if (value == "--distance") {
            options.adaptive = true;
            const auto specification = next();
            const auto [distance, level] =
                split_level_suffix(specification, value);
            options.refinement.distance_bands.push_back(
                {parse_nonnegative(distance), level});
        } else if (value == "--curvature") {
            options.adaptive = true;
            const auto specification = next();
            const auto [degrees, level] =
                split_level_suffix(specification, value);
            options.refinement.curvature = cartmesh::CurvatureRefinementRule{
                parse_nonnegative(degrees), 1.5, level};
        } else if (value == "--gap") {
            options.adaptive = true;
            const auto specification = next();
            const auto text = std::string_view(specification);
            const auto colon = text.rfind(':');
            if (colon == std::string_view::npos) {
                throw std::invalid_argument("--gap 需要 D:N 格式");
            }
            options.refinement.gap = cartmesh::GapRefinementRule{
                parse_nonnegative(text.substr(0, colon)), -0.75, 0.5,
                parse_positive_u32(text.substr(colon + 1U), value)};
        } else if (value == "--box") {
            options.adaptive = true;
            const auto specification = next();
            const auto [coordinates, level] =
                split_level_suffix(specification, value);
            const auto v = parse_csv(coordinates, 6, value);
            options.refinement.boxes.push_back(
                {cartmesh::AABB({v[0], v[1], v[2]}, {v[3], v[4], v[5]}),
                 level});
        } else if (value == "--sphere") {
            options.adaptive = true;
            const auto specification = next();
            const auto [coordinates, level] =
                split_level_suffix(specification, value);
            const auto v = parse_csv(coordinates, 4, value);
            options.refinement.spheres.push_back(
                {{v[0], v[1], v[2]}, v[3], level});
        } else if (value == "--cylinder") {
            options.adaptive = true;
            const auto specification = next();
            const auto [coordinates, level] =
                split_level_suffix(specification, value);
            const auto v = parse_csv(coordinates, 7, value);
            options.refinement.cylinders.push_back(
                {{v[0], v[1], v[2]}, {v[3], v[4], v[5]}, v[6], level});
        } else if (value == "--no-balance") {
            options.adaptive = true;
            options.refinement.enforce_face_2_to_1_balance = false;
        } else if (value == "--strict-gaps") {
            options.adaptive = true;
            options.strict_gap_resolution = true;
        } else if (value == "--boundary-range") {
            options.boundary_ranges.push_back(parse_boundary_range(next()));
        } else if (value == "--region-name") {
            options.region_names.push_back(parse_region_name(next()));
        } else if (value == "--output") {
            options.output = next();
        } else if (value == "--boundary-output") {
            options.boundary_output = next();
        } else if (value == "--geometry-output") {
            options.geometry_output = next();
        } else if (value == "--polyhedra-output") {
            options.polyhedra_output = next();
        } else if (value == "--tetrahedra-output") {
            options.tetrahedra_output = next();
        } else if (value == "--report") {
            options.report = next();
        } else if (value == "--openfoam-case") {
            options.openfoam_case = next();
        } else if (value == "--quality-output") {
            options.quality_output = next();
        } else if (value == "--stabilize") {
            options.stabilize = true;
        } else if (value == "--stabilization-output") {
            options.stabilize = true;
            options.stabilization_output = next();
        } else if (value == "--max-stabilization-rounds") {
            options.stabilize = true;
            options.maximum_stabilization_rounds =
                parse_positive_u32(next(), value);
        } else if (value == "--no-vtk") {
            options.write_vtk = false;
        } else if (value == "--help") {
            std::cout
                << "用法：cartmesh_cutcell_cli --stl FILE [选项]\n"
                << "  --resolution N          三个方向的均匀单元数（默认 24）\n"
                << "  --padding-fraction F    相对最大 STL 尺寸的域外扩比例（默认 0.1）\n"
                << "  --small-cell-threshold F 小流体单元体积分数阈值（默认 0.01）\n"
                << "  --geometric-tolerance F 显式几何长度容差（默认自动）\n"
                << "  --adaptive             使用阶段2线性八叉树背景网格\n"
                << "  --base-level/--max-level/--surface-level L\n"
                << "  --distance D:L --curvature DEG:L --gap D:N\n"
                << "  --box ...:L --sphere ...:L --cylinder ...:L\n"
                << "  --no-balance --strict-gaps\n"
                << "  --boundary-range FIRST:LAST:ID:NAME  三角片半开范围的边界命名\n"
                << "  --region-name ID:NAME  确定性全局流体 region 命名\n"
                << "  --output FILE.vtu       背景网格与流体体积分数\n"
                << "  --boundary-output FILE.vtp  嵌入边界多边形\n"
                << "  --geometry-output FILE.json 完整单元面与邻接几何\n"
                << "  --polyhedra-output FILE.vtu 显式流体凸多面体片\n"
                << "  --tetrahedra-output FILE.vtu 凸片的外部检查四面体分解\n"
              << "  --report FILE.json      阶段四 Cut-cell 机器报告\n"
              << "  --openfoam-case DIR     写出完整流体域 constant/polyMesh\n"
              << "  --quality-output FILE   写出同一 solver mesh 的原生质量诊断\n"
              << "  --stabilize             执行 Stage 6.3 保守聚并/局部细化闭环\n"
              << "  --stabilization-output FILE 写出 Stage 6.3 稳定化报告\n"
              << "  --max-stabilization-rounds N 最多局部细化轮数（默认 2）\n"
                << "  --no-vtk                只计算并写 JSON\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("未知参数：" + value);
        }
    }
    if (options.stl.empty()) {
        throw std::invalid_argument("必须提供 --stl FILE");
    }
    if (options.base_level > options.maximum_level) {
        throw std::invalid_argument("基础层级不得超过最大层级");
    }
    if (options.adaptive) {
        options.refinement.surface_target_level =
            options.surface_level.value_or(options.maximum_level);
    }
    return options;
}

[[nodiscard]] std::uint64_t peak_rss_bytes() noexcept {
#if defined(__unix__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0;
    }
#if defined(__APPLE__)
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
#else
    return 0;
#endif
}

void write_vec3(std::ostream& output, const cartmesh::Vec3& value) {
    output << '[' << value.x << ", " << value.y << ", " << value.z << ']';
}

void write_refinement_rules(std::ostream& output, const Options& options) {
    output << "  \"refinementRules\": {\"surfaceTargetLevel\":"
           << static_cast<unsigned>(*options.refinement.surface_target_level)
           << ",\"distanceBands\":[";
    for (std::size_t index = 0;
         index < options.refinement.distance_bands.size(); ++index) {
        if (index != 0) output << ',';
        const auto& band = options.refinement.distance_bands[index];
        output << "{\"maximumDistance\":" << band.maximum_distance
               << ",\"targetLevel\":"
               << static_cast<unsigned>(band.target_level) << '}';
    }
    output << "],\"curvature\":";
    if (options.refinement.curvature) {
        const auto& rule = *options.refinement.curvature;
        output << "{\"minimumNormalAngleDegrees\":"
               << rule.minimum_normal_angle_degrees
               << ",\"neighborhoodCellDiagonals\":"
               << rule.neighborhood_cell_diagonals << ",\"targetLevel\":"
               << static_cast<unsigned>(rule.target_level) << '}';
    } else {
        output << "null";
    }
    output << ",\"gap\":";
    if (options.refinement.gap) {
        const auto& rule = *options.refinement.gap;
        output << "{\"maximumSearchDistance\":"
               << rule.maximum_search_distance
               << ",\"maximumOpposingNormalDot\":"
               << rule.maximum_opposing_normal_dot
               << ",\"minimumFacingDot\":" << rule.minimum_facing_dot
               << ",\"minimumCellsAcrossGap\":"
               << rule.minimum_cells_across_gap << '}';
    } else {
        output << "null";
    }
    output << ",\"boxes\":[";
    for (std::size_t index = 0; index < options.refinement.boxes.size(); ++index) {
        if (index != 0) output << ',';
        const auto& region = options.refinement.boxes[index];
        output << "{\"minimum\":";
        write_vec3(output, region.bounds.minimum());
        output << ",\"maximum\":";
        write_vec3(output, region.bounds.maximum());
        output << ",\"targetLevel\":"
               << static_cast<unsigned>(region.target_level) << '}';
    }
    output << "],\"spheres\":[";
    for (std::size_t index = 0; index < options.refinement.spheres.size();
         ++index) {
        if (index != 0) output << ',';
        const auto& region = options.refinement.spheres[index];
        output << "{\"center\":";
        write_vec3(output, region.center);
        output << ",\"radius\":" << region.radius << ",\"targetLevel\":"
               << static_cast<unsigned>(region.target_level) << '}';
    }
    output << "],\"cylinders\":[";
    for (std::size_t index = 0; index < options.refinement.cylinders.size();
         ++index) {
        if (index != 0) output << ',';
        const auto& region = options.refinement.cylinders[index];
        output << "{\"firstAxisPoint\":";
        write_vec3(output, region.first_axis_point);
        output << ",\"secondAxisPoint\":";
        write_vec3(output, region.second_axis_point);
        output << ",\"radius\":" << region.radius << ",\"targetLevel\":"
               << static_cast<unsigned>(region.target_level) << '}';
    }
    output << "]},\n";
}

void ensure_parent(const std::filesystem::path& path) {
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
}

[[nodiscard]] std::filesystem::path diagnostic_marker_path(
    const std::filesystem::path& report_path) {
    auto path = report_path;
    path.replace_extension();
    path += "_geometry_issues.vtp";
    return path;
}

void write_rejected_surface_report(
    const Options& options, const cartmesh::SurfaceMesh& surface,
    const cartmesh::SurfaceDiagnostics& diagnostics,
    const std::string& input_hash,
    const std::filesystem::path& marker_path) {
    ensure_parent(options.report);
    std::ofstream report(options.report, std::ios::trunc);
    if (!report) {
        throw std::runtime_error("无法打开阶段四无效几何报告：" +
                                 options.report.string());
    }
    report << std::setprecision(17)
           << "{\n"
           << "  \"schema\": \"cartmesh-stage4-cutcell-v1\",\n"
           << "  \"status\": \"rejected_invalid_surface\",\n"
           << "  \"stage3Complete\": false,\n"
           << "  \"stage4Complete\": false,\n"
           << "  \"solverReadyCutCellMesh\": false,\n"
           << "  \"input\": \"" << json_escape(options.stl.string())
           << "\",\n"
           << "  \"inputHashFnv1a64\": \"" << input_hash << "\",\n"
           << "  \"surfaceFormat\": \""
           << cartmesh::surface_format_name(surface.format()) << "\",\n"
           << "  \"diagnosticMarker\": \""
           << json_escape(marker_path.string()) << "\",\n"
           << "  \"surfaceDiagnostics\": {\n"
           << "    \"triangleCount\": " << diagnostics.triangle_count << ",\n"
           << "    \"degenerateTriangleCount\": "
           << diagnostics.degenerate_triangle_count << ",\n"
           << "    \"duplicateTriangleCount\": "
           << diagnostics.duplicate_triangle_count << ",\n"
           << "    \"boundaryEdgeCount\": "
           << diagnostics.boundary_edge_count << ",\n"
           << "    \"nonManifoldEdgeCount\": "
           << diagnostics.non_manifold_edge_count << ",\n"
           << "    \"nonManifoldVertexCount\": "
           << diagnostics.non_manifold_vertex_count << ",\n"
           << "    \"orientationConflictEdgeCount\": "
           << diagnostics.orientation_conflict_edge_count << ",\n"
           << "    \"overlappingTrianglePairCount\": "
           << diagnostics.overlapping_triangle_pair_count << ",\n"
           << "    \"selfIntersectionPairCount\": "
           << diagnostics.self_intersection_pair_count << ",\n"
           << "    \"nonAdjacentContactPairCount\": "
           << diagnostics.non_adjacent_contact_pair_count << ",\n"
           << "    \"closed\": " << (diagnostics.closed ? "true" : "false")
           << ",\n"
           << "    \"manifold\": "
           << (diagnostics.manifold ? "true" : "false") << ",\n"
           << "    \"consistentlyOriented\": "
           << (diagnostics.consistently_oriented ? "true" : "false") << "\n"
           << "  }\n"
           << "}\n";
    if (!report) {
        throw std::runtime_error("写入阶段四无效几何报告失败：" +
                                 options.report.string());
    }
}

void hash_byte(std::uint64_t& hash, std::uint8_t value) noexcept {
    hash ^= value;
    hash *= 1099511628211ULL;
}

void hash_u64(std::uint64_t& hash, std::uint64_t value) noexcept {
    for (std::uint32_t shift = 0; shift < 64U; shift += 8U) {
        hash_byte(hash, static_cast<std::uint8_t>((value >> shift) & 0xffU));
    }
}

void hash_double(std::uint64_t& hash, double value) noexcept {
    hash_u64(hash, std::bit_cast<std::uint64_t>(value));
}

void hash_vec3(std::uint64_t& hash, const cartmesh::Vec3& value) noexcept {
    hash_double(hash, value.x);
    hash_double(hash, value.y);
    hash_double(hash, value.z);
}

[[nodiscard]] std::string hex_hash(std::uint64_t hash) {
    std::ostringstream output;
    output << std::hex << std::setfill('0') << std::setw(16) << hash;
    return output.str();
}

[[nodiscard]] std::string file_hash(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("无法计算输入哈希：" + path.string());
    std::uint64_t hash = 14695981039346656037ULL;
    std::array<char, 65536> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        for (std::streamsize index = 0; index < input.gcount(); ++index) {
            hash_byte(hash, static_cast<std::uint8_t>(buffer[static_cast<std::size_t>(index)]));
        }
    }
    return hex_hash(hash);
}

[[nodiscard]] constexpr std::string_view build_type() noexcept {
#if defined(NDEBUG)
    return "Release";
#else
    return "Debug";
#endif
}

[[nodiscard]] constexpr std::string_view compiler() noexcept {
#if defined(__clang__)
    return "Clang " __clang_version__;
#elif defined(__GNUC__)
    return "GCC " __VERSION__;
#else
    return "unknown";
#endif
}

void hash_mesh_payload(std::uint64_t& hash,
                       const cartmesh::ConvexCutCellMesh& mesh) {
    for (const auto& cell : mesh.fluid_cells) {
        hash_u64(hash, cell.background_cell_id);
        hash_double(hash, cell.volume);
        hash_double(hash, cell.volume_fraction);
        hash_vec3(hash, cell.centroid);
        hash_u64(hash, cell.fluid_component_count);
        for (const auto& face : cell.cartesian_faces) {
            hash_double(hash, face.area);
            hash_vec3(hash, face.centroid);
            hash_vec3(hash, face.outward_normal);
            for (const auto& loop : face.oriented_boundary_loops) {
                hash_u64(hash, loop.size());
                for (const auto& vertex : loop) hash_vec3(hash, vertex);
            }
        }
        for (const auto& face : cell.embedded_boundary_faces) {
            hash_u64(hash, face.boundary_id);
            hash_double(hash, face.area);
            hash_vec3(hash, face.centroid);
            hash_vec3(hash, face.outward_normal);
            for (const auto& vertex : face.vertices) hash_vec3(hash, vertex);
        }
        for (const auto& piece : cell.fluid_polyhedron_pieces) {
            hash_u64(hash, piece.component_id);
            hash_u64(hash, piece.global_region_id);
            hash_double(hash, piece.geometry.volume);
            for (const auto& vertex : piece.polyhedron.vertices) hash_vec3(hash, vertex);
            for (const auto& face : piece.polyhedron.faces) {
                hash_byte(hash, static_cast<std::uint8_t>(face.kind));
                hash_u64(hash, face.source_id);
                for (const auto vertex : face.vertex_indices) hash_u64(hash, vertex);
            }
        }
        for (const auto region_id : cell.fluid_component_region_ids) {
            hash_u64(hash, region_id);
        }
    }
    for (const auto& face : mesh.internal_faces) {
        hash_u64(hash, face.first_background_cell_id);
        hash_u64(hash, face.second_background_cell_id);
        hash_double(hash, face.area);
        hash_vec3(hash, face.centroid);
        hash_vec3(hash, face.normal);
    }
    for (const auto& face : mesh.component_internal_faces) {
        hash_u64(hash, face.first_background_cell_id);
        hash_u64(hash, face.second_background_cell_id);
        hash_u64(hash, face.first_component_id);
        hash_u64(hash, face.second_component_id);
        hash_u64(hash, face.global_region_id);
        hash_double(hash, face.area);
        hash_vec3(hash, face.centroid);
    }
    for (const auto volume : mesh.global_fluid_region_volumes) {
        hash_double(hash, volume);
    }
}

[[nodiscard]] std::string mesh_hash(const cartmesh::UniformCartesianGrid& grid,
                                    const cartmesh::ConvexCutCellMesh& mesh) {
    std::uint64_t hash = 14695981039346656037ULL;
    hash_vec3(hash, grid.domain().minimum());
    hash_vec3(hash, grid.domain().maximum());
    hash_u64(hash, grid.nx());
    hash_u64(hash, grid.ny());
    hash_u64(hash, grid.nz());
    hash_mesh_payload(hash, mesh);
    return hex_hash(hash);
}

[[nodiscard]] std::string mesh_hash(const cartmesh::LinearOctree& tree,
                                    const cartmesh::ConvexCutCellMesh& mesh) {
    std::uint64_t hash = 14695981039346656037ULL;
    hash_vec3(hash, tree.domain().minimum());
    hash_vec3(hash, tree.domain().maximum());
    hash_u64(hash, tree.base_level());
    hash_u64(hash, tree.maximum_level());
    for (const auto code : tree.leaf_codes()) hash_u64(hash, code);
    hash_mesh_payload(hash, mesh);
    return hex_hash(hash);
}

} // 匿名命名空间

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        const auto start = std::chrono::steady_clock::now();
        const cartmesh::SurfaceMesh surface = cartmesh::read_stl(options.stl);
        const std::string input_hash = file_hash(options.stl);
        const auto diagnostics = cartmesh::diagnose_surface(surface);
        if (!diagnostics.valid_for_stage1_classification()) {
            const auto marker_path = diagnostic_marker_path(options.report);
            ensure_parent(marker_path);
            cartmesh::write_surface_diagnostic_vtp(marker_path, surface,
                                                   diagnostics);
            write_rejected_surface_report(options, surface, diagnostics,
                                          input_hash, marker_path);
            std::cerr << "错误：输入 STL 未通过封闭、定向、流形和相交诊断；详情见 "
                      << options.report << '\n';
            return 2;
        }
        const cartmesh::Vec3 surface_extent = surface.bounds().extent();
        const double maximum_extent =
            std::max({surface_extent.x, surface_extent.y, surface_extent.z});
        const double padding = options.padding_fraction * maximum_extent;
        const cartmesh::Vec3 padding_vector{padding, padding, padding};
        const cartmesh::AABB domain(surface.bounds().minimum() - padding_vector,
                                    surface.bounds().maximum() + padding_vector);
        std::optional<cartmesh::UniformCartesianGrid> grid;
        std::optional<cartmesh::LinearOctree> tree;
        std::unique_ptr<cartmesh::TriangulatedSurfaceCutter>
            triangulated_cutter;
        cartmesh::OctreeAdaptationStatistics adaptation;
        cartmesh::ConvexCutCellMesh mesh;
        std::string kernel;
        std::string supported_geometry;
        std::size_t supporting_plane_count = 0;
        std::optional<std::vector<std::uint64_t>> named_boundary_ids;
        if (!options.boundary_ranges.empty()) {
            std::vector<std::uint64_t> triangle_boundary_ids(
                surface.triangles().size(), 0);
            std::vector<bool> boundary_assigned(surface.triangles().size(), false);
            for (const auto& range : options.boundary_ranges) {
                if (range.last_triangle > surface.triangles().size()) {
                    throw std::invalid_argument(
                        "boundary-range 超出 STL 三角片数量");
                }
                for (std::uint64_t triangle = range.first_triangle;
                     triangle < range.last_triangle; ++triangle) {
                    const auto index = static_cast<std::size_t>(triangle);
                    if (boundary_assigned[index]) {
                        throw std::invalid_argument(
                            "boundary-range 之间不得重叠");
                    }
                    boundary_assigned[index] = true;
                    triangle_boundary_ids[index] = range.boundary_id;
                }
            }
            if (std::find(boundary_assigned.begin(), boundary_assigned.end(),
                          false) != boundary_assigned.end()) {
                throw std::invalid_argument(
                    "提供 boundary-range 时必须恰好覆盖全部 STL 三角片");
            }
            named_boundary_ids = std::move(triangle_boundary_ids);
        }
        if (options.adaptive) {
            tree.emplace(domain, options.base_level, options.maximum_level);
            const auto build_adaptive = [&](const auto& cutter) {
                adaptation = cartmesh::OctreeRefinementEngine(
                                 options.refinement, &cutter.bvh())
                                 .apply(*tree);
                mesh = cartmesh::build_triangulated_cut_cell_mesh(
                    *tree, cutter, options.geometric_tolerance);
            };
            if (named_boundary_ids) {
                triangulated_cutter =
                    std::make_unique<cartmesh::TriangulatedSurfaceCutter>(
                        surface, std::move(*named_boundary_ids));
                build_adaptive(*triangulated_cutter);
                supported_geometry = "named_closed_oriented_shells";
            } else {
                triangulated_cutter =
                    std::make_unique<cartmesh::TriangulatedSurfaceCutter>(
                        surface, 0);
                build_adaptive(*triangulated_cutter);
                supported_geometry =
                    "closed_oriented_nested_or_disjoint_shells";
            }
            kernel = "adaptive_oriented_tetrahedral_chain";
        } else {
            grid.emplace(domain, options.resolution, options.resolution,
                         options.resolution);
            if (named_boundary_ids) {
                triangulated_cutter =
                    std::make_unique<cartmesh::TriangulatedSurfaceCutter>(
                        surface, std::move(*named_boundary_ids));
                mesh = cartmesh::build_triangulated_cut_cell_mesh(
                    *grid, *triangulated_cutter,
                    options.geometric_tolerance);
                kernel = "oriented_tetrahedral_chain";
                supported_geometry = "named_closed_oriented_shells";
            } else if (options.openfoam_case || options.quality_output ||
                       options.stabilize) {
                triangulated_cutter =
                    std::make_unique<cartmesh::TriangulatedSurfaceCutter>(
                        surface, 0);
                mesh = cartmesh::build_triangulated_cut_cell_mesh(
                    *grid, *triangulated_cutter,
                    options.geometric_tolerance);
                kernel = "oriented_tetrahedral_chain";
                supported_geometry =
                    "closed_oriented_nested_or_disjoint_shells";
            } else try {
                const cartmesh::ConvexSurfaceCutter cutter(surface, 0);
                supporting_plane_count = cutter.plane_count();
                mesh = cartmesh::build_convex_cut_cell_mesh(
                    *grid, cutter, options.geometric_tolerance);
                kernel = "convex_halfspace";
                supported_geometry = "single_component_closed_convex_stl";
            } catch (const std::invalid_argument&) {
                triangulated_cutter =
                    std::make_unique<cartmesh::TriangulatedSurfaceCutter>(
                        surface, 0);
                mesh = cartmesh::build_triangulated_cut_cell_mesh(
                    *grid, *triangulated_cutter,
                    options.geometric_tolerance);
                kernel = "oriented_tetrahedral_chain";
                supported_geometry =
                    "closed_oriented_nested_or_disjoint_shells";
            }
        }
        std::optional<std::filesystem::path> quality_output =
            options.quality_output;
        if (!quality_output && options.openfoam_case) {
            quality_output = *options.openfoam_case / "cartmeshQuality.json";
        }
        std::optional<std::filesystem::path> stabilization_output =
            options.stabilization_output;
        if (options.stabilize && !stabilization_output) {
            if (options.openfoam_case) {
                stabilization_output =
                    *options.openfoam_case / "cartmeshStabilization.json";
            } else {
                stabilization_output = options.report;
                stabilization_output->replace_extension();
                *stabilization_output += "_stabilization.json";
            }
        }
        std::optional<cartmesh::OpenFoamMesh> solver_mesh;
        std::optional<cartmesh::MeshQualityReport> quality_report;
        std::optional<cartmesh::SolverMeshStabilizationReport>
            stabilization_report;
        std::uint64_t stabilization_refined_leaf_count = 0U;
        std::uint64_t stabilization_balance_split_count = 0U;
        if (options.stabilize) {
            if (!triangulated_cutter) {
                throw std::runtime_error(
                    "Stage 6.3 稳定化需要可重建的三角化曲面 cutter");
            }
            cartmesh::MeshQualityThresholds thresholds;
            thresholds.minimum_volume_fraction = options.small_cell_threshold;
            cartmesh::SolverMeshStabilizationOptions stabilization_options;
            stabilization_options.minimum_volume_fraction =
                options.small_cell_threshold;
            const double writer_tolerance = std::max(
                options.geometric_tolerance,
                diagnostics.suggested_length_tolerance);
            cartmesh::SolverMeshStabilizationReport aggregate;
            bool first_round = true;
            for (std::uint32_t round = 0;
                 round <= options.maximum_stabilization_rounds; ++round) {
                auto raw_solver_mesh = tree
                                           ? cartmesh::build_openfoam_mesh(
                                                 *tree, mesh, writer_tolerance)
                                           : cartmesh::build_openfoam_mesh(
                                                 *grid, mesh, writer_tolerance);
                auto stabilized = cartmesh::stabilize_solver_mesh(
                    raw_solver_mesh, thresholds, stabilization_options);
                if (first_round) {
                    aggregate.initial_cell_count =
                        stabilized.report.initial_cell_count;
                    aggregate.initial_volume = stabilized.report.initial_volume;
                    aggregate.initial_first_moment =
                        stabilized.report.initial_first_moment;
                    first_round = false;
                }
                aggregate.final_cell_count = stabilized.report.final_cell_count;
                aggregate.final_volume = stabilized.report.final_volume;
                aggregate.final_first_moment =
                    stabilized.report.final_first_moment;
                aggregate.agglomeration_count +=
                    stabilized.report.agglomeration_count;
                aggregate.rejected_candidate_count +=
                    stabilized.report.rejected_candidate_count;
                aggregate.actions.insert(
                    aggregate.actions.end(), stabilized.report.actions.begin(),
                    stabilized.report.actions.end());
                aggregate.refinement_requested_stable_ids.insert(
                    aggregate.refinement_requested_stable_ids.end(),
                    stabilized.report.refinement_requested_stable_ids.begin(),
                    stabilized.report.refinement_requested_stable_ids.end());
                solver_mesh = std::move(stabilized.mesh);
                if (stabilized.report.pass()) {
                    aggregate.unresolved_stable_ids.clear();
                    break;
                }
                aggregate.unresolved_stable_ids =
                    stabilized.report.unresolved_stable_ids;
                if (!tree || round == options.maximum_stabilization_rounds) {
                    break;
                }
                const auto refinement =
                    cartmesh::refine_stabilization_sources(
                        *tree, stabilized.report.refinement_requested_stable_ids);
                stabilization_refined_leaf_count +=
                    refinement.refined_leaf_count;
                stabilization_balance_split_count +=
                    refinement.balance_split_count;
                for (const auto stable_id : refinement.refined_stable_ids) {
                    aggregate.actions.push_back(
                        {cartmesh::StabilizationActionKind::conformal_refined,
                         stable_id, 0U, 0.0, 0.0,
                         "refined source Morton leaf and restored 2:1 balance"});
                }
                if (!refinement.unresolved_stable_ids.empty()) {
                    aggregate.unresolved_stable_ids =
                        refinement.unresolved_stable_ids;
                    for (const auto stable_id :
                         refinement.unresolved_stable_ids) {
                        aggregate.actions.push_back(
                            {cartmesh::StabilizationActionKind::unresolved_max_level,
                             stable_id, 0U, 0.0, 0.0,
                             "source is not a current refinable leaf or reached maximum level"});
                    }
                    break;
                }
                mesh = cartmesh::build_triangulated_cut_cell_mesh(
                    *tree, *triangulated_cutter,
                    options.geometric_tolerance);
            }
            std::sort(aggregate.refinement_requested_stable_ids.begin(),
                      aggregate.refinement_requested_stable_ids.end());
            aggregate.refinement_requested_stable_ids.erase(
                std::unique(aggregate.refinement_requested_stable_ids.begin(),
                            aggregate.refinement_requested_stable_ids.end()),
                aggregate.refinement_requested_stable_ids.end());
            std::sort(aggregate.unresolved_stable_ids.begin(),
                      aggregate.unresolved_stable_ids.end());
            aggregate.unresolved_stable_ids.erase(
                std::unique(aggregate.unresolved_stable_ids.begin(),
                            aggregate.unresolved_stable_ids.end()),
                aggregate.unresolved_stable_ids.end());
            const double conservation_scale = std::max(
                {1.0, std::abs(aggregate.initial_volume),
                 std::abs(aggregate.final_volume),
                 cartmesh::norm(aggregate.initial_first_moment),
                 cartmesh::norm(aggregate.final_first_moment)});
            aggregate.conservation_pass =
                std::abs(aggregate.initial_volume - aggregate.final_volume) <=
                    1.0e-10 * conservation_scale &&
                cartmesh::norm(aggregate.initial_first_moment -
                               aggregate.final_first_moment) <=
                    1.0e-10 * conservation_scale;
            stabilization_report = std::move(aggregate);
            cartmesh::write_solver_mesh_stabilization_json(
                *stabilization_output, *stabilization_report);
        }
        const bool stabilization_pass =
            !stabilization_report || stabilization_report->pass();
        const double solid_volume = domain.volume() - mesh.total_fluid_volume;
        for (const auto& region : options.region_names) {
            if (region.region_id >= mesh.global_fluid_region_count) {
                throw std::invalid_argument("region-name ID 超出生成的全局流体区数");
            }
        }
        const auto region_name = [&](std::size_t id) -> std::string {
            for (auto iterator = options.region_names.rbegin();
                 iterator != options.region_names.rend(); ++iterator) {
                if (iterator->region_id == id) return iterator->name;
            }
            return {};
        };
        double triangle_area = 0.0;
        for (const auto& triangle : surface.triangles()) {
            triangle_area += triangle.area();
        }
        const double conservation_tolerance =
            1.0e-10 * std::max(1.0, triangle_area);
        const double region_volume_sum = std::accumulate(
            mesh.global_fluid_region_volumes.begin(),
            mesh.global_fluid_region_volumes.end(), 0.0);
        const double volume_tolerance =
            1.0e-10 * std::max(1.0, domain.volume());
        const bool gap_constraint_pass =
            !options.adaptive ||
            adaptation.gap_resolution_failure_count == 0 ||
            !options.strict_gap_resolution;
        const bool invariants_pass =
            mesh.nonclosed_cell_count == 0 &&
            mesh.shared_face_mismatch_count == 0 &&
            mesh.component_analysis_pending_cell_count == 0 &&
            mesh.classification_conflict_count == 0 &&
            mesh.negative_volume_cell_count == 0 && solid_volume >= 0.0 &&
            std::abs(region_volume_sum - mesh.total_fluid_volume) <=
                volume_tolerance && gap_constraint_pass &&
            std::abs(mesh.total_embedded_boundary_area - triangle_area) <=
                conservation_tolerance;
        const bool stage4_complete =
            invariants_pass &&
            (!options.adaptive ||
             adaptation.gap_resolution_failure_count == 0);
        if (options.openfoam_case || quality_output) {
            std::vector<std::pair<std::uint64_t, std::string>> boundary_names;
            boundary_names.reserve(options.boundary_ranges.size());
            for (const auto& range : options.boundary_ranges) {
                boundary_names.emplace_back(range.boundary_id, range.name);
            }
            const double writer_tolerance = std::max(
                options.geometric_tolerance,
                diagnostics.suggested_length_tolerance);
            if (!solver_mesh) {
                solver_mesh = tree
                                  ? cartmesh::build_openfoam_mesh(
                                        *tree, mesh, writer_tolerance)
                                  : cartmesh::build_openfoam_mesh(
                                        *grid, mesh, writer_tolerance);
            }
            cartmesh::MeshQualityThresholds quality_thresholds;
            quality_thresholds.minimum_volume_fraction =
                options.small_cell_threshold;
            quality_report = cartmesh::evaluate_solver_mesh_quality(
                *solver_mesh, quality_thresholds);
            if (quality_output) {
                cartmesh::write_solver_mesh_quality_json(
                    *quality_output, *quality_report);
            }
            if (options.openfoam_case && stabilization_pass) {
                cartmesh::write_openfoam_poly_mesh(
                    *options.openfoam_case, *solver_mesh, boundary_names);
            }
        }
        const std::uint64_t background_cell_count =
            tree ? tree->leaf_count() : grid->cell_count();
        if (tree) {
            cartmesh::write_cut_cell_geometry_json(
                options.geometry_output, *tree, mesh, false);
        } else {
            cartmesh::write_cut_cell_geometry_json(
                options.geometry_output, *grid, mesh, false);
        }
        const std::string result_hash =
            tree ? mesh_hash(*tree, mesh) : mesh_hash(*grid, mesh);

        const auto small_cells = cartmesh::analyze_small_cut_cells(
            mesh, options.small_cell_threshold);
        std::uint64_t explicit_piece_count = 0;
        double explicit_piece_volume = 0.0;
        for (const auto& cell : mesh.fluid_cells) {
            explicit_piece_count += cell.fluid_polyhedron_pieces.size();
            for (const auto& piece : cell.fluid_polyhedron_pieces) {
                explicit_piece_volume += piece.geometry.volume;
            }
        }

        if (options.write_vtk) {
            std::vector<double> fluid_fraction(
                static_cast<std::size_t>(background_cell_count), 0.0);
            std::vector<double> cut_cell(
                static_cast<std::size_t>(background_cell_count), 0.0);
            std::vector<double> small_cut_cell(
                static_cast<std::size_t>(background_cell_count), 0.0);
            std::vector<double> fluid_region_id(
                static_cast<std::size_t>(background_cell_count), -1.0);
            std::vector<double> multiple_fluid_components(
                static_cast<std::size_t>(background_cell_count), 0.0);
            for (const auto& cell : mesh.fluid_cells) {
                fluid_fraction[static_cast<std::size_t>(cell.background_cell_id)] =
                    cell.volume_fraction;
                cut_cell[static_cast<std::size_t>(cell.background_cell_id)] =
                    cell.cut ? 1.0 : 0.0;
                small_cut_cell[static_cast<std::size_t>(cell.background_cell_id)] =
                    cell.cut && cell.volume_fraction < options.small_cell_threshold
                        ? 1.0
                        : 0.0;
                if (!cell.fluid_component_region_ids.empty()) {
                    fluid_region_id[static_cast<std::size_t>(
                        cell.background_cell_id)] =
                        static_cast<double>(cell.fluid_component_region_ids.front());
                }
                multiple_fluid_components[static_cast<std::size_t>(
                    cell.background_cell_id)] =
                    cell.fluid_component_count > 1 ? 1.0 : 0.0;
            }
            std::vector<cartmesh::VtkCellData> fields{
                {"fluid_volume_fraction", std::move(fluid_fraction)},
                {"cut_cell", std::move(cut_cell)},
                {"small_cut_cell", std::move(small_cut_cell)},
                {"fluid_region_id", std::move(fluid_region_id)},
                {"multiple_fluid_components",
                 std::move(multiple_fluid_components)}};
            if (tree) {
                std::vector<double> level(tree->leaf_count());
                std::vector<double> code_low(tree->leaf_count());
                std::vector<double> code_high(tree->leaf_count());
                for (std::uint64_t leaf = 0; leaf < tree->leaf_count(); ++leaf) {
                    const auto code = tree->leaf_code(leaf);
                    const auto index = static_cast<std::size_t>(leaf);
                    level[index] = static_cast<double>(
                        cartmesh::decode_octree_node(code).level);
                    code_low[index] =
                        static_cast<double>(code & 0xffffffffULL);
                    code_high[index] = static_cast<double>(code >> 32U);
                }
                fields.push_back({"octree_level", std::move(level)});
                fields.push_back(
                    {"octree_node_code_low32", std::move(code_low)});
                fields.push_back(
                    {"octree_node_code_high32", std::move(code_high)});
                cartmesh::write_octree_vtu(options.output, *tree, fields);
            } else {
                cartmesh::write_vtu(options.output, *grid, fields);
            }
            cartmesh::write_embedded_boundary_vtp(options.boundary_output, mesh);
            cartmesh::write_fluid_polyhedra_vtu(options.polyhedra_output, mesh);
            cartmesh::write_fluid_tetrahedra_vtu(options.tetrahedra_output, mesh);
        }

        const auto end = std::chrono::steady_clock::now();
        const double wall_seconds = std::chrono::duration<double>(end - start).count();
        ensure_parent(options.report);
        std::ofstream report(options.report, std::ios::trunc);
        if (!report) {
            throw std::runtime_error("无法打开阶段四报告：" + options.report.string());
        }
        report << std::setprecision(17)
               << "{\n"
               << "  \"schema\": \"cartmesh-stage4-cutcell-v1\",\n"
               << "  \"status\": \""
               << (!stabilization_pass
                       ? "failed_solver_mesh_stabilization"
                       : invariants_pass
                       ? (options.adaptive &&
                                  adaptation.gap_resolution_failure_count > 0
                              ? "pass_with_gap_resolution_warning"
                              : "geometry_pass_external_cfd_pending")
                       : (options.adaptive && !gap_constraint_pass
                              ? "failed_gap_resolution"
                              : "fail"))
               << "\",\n"
               << "  \"stage3Complete\": "
               << "false,\n"
               << "  \"stage4Complete\": "
               << "false,\n"
               << "  \"stage3GeometryTopologyComplete\": "
               << (invariants_pass ? "true" : "false") << ",\n"
               << "  \"stage4GeometryRobustnessComplete\": "
               << (stage4_complete ? "true" : "false") << ",\n"
               << "  \"externalCfdCheckerAccepted\": false,\n"
               << "  \"solverReadyCutCellMesh\": "
               << "false,\n"
               << "  \"acceptanceBlockers\": "
               << (!stabilization_pass
                       ? "[\"solver_mesh_stabilization\",\"external_cfd_checker\"]"
                       : options.openfoam_case
                             ? "[\"external_cfd_checker\"]"
                             : "[\"complete_volume_mesh_export_and_external_cfd_checker\"]")
               << ",\n"
               << "  \"stage4Capabilities\": {\"multipleComponents\":true,"
                  "\"nestedCavities\":true,\"globalFluidRegions\":true,"
                  "\"namedBoundaries\":true,\"thinWallCutCells\":true,"
                  "\"adaptiveOctreeCutCells\":true},\n"
               << "  \"geometryRepairApplied\": false,\n"
               << "  \"geometricTolerance\": "
               << options.geometric_tolerance << ",\n"
               << "  \"surfaceLengthTolerance\": "
               << diagnostics.suggested_length_tolerance << ",\n"
               << "  \"toleranceActions\": "
               << (options.geometric_tolerance > 0.0
                       ? "[\"explicit_geometric_length_tolerance\"]"
                       : "[]")
               << ",\n"
               << "  \"kernel\": \"" << kernel << "\",\n"
               << "  \"meshKind\": \""
               << (tree ? "adaptive_linear_octree_cut_cells"
                        : "uniform_cartesian_cut_cells")
               << "\",\n"
               << "  \"supportedGeometry\": \"" << supported_geometry << "\",\n"
               << "  \"boundaries\": [";
        for (std::size_t index = 0; index < options.boundary_ranges.size(); ++index) {
            if (index != 0) report << ',';
            const auto& range = options.boundary_ranges[index];
            report << "{\"firstTriangle\":" << range.first_triangle
                   << ",\"lastTriangleExclusive\":" << range.last_triangle
                   << ",\"boundaryId\":" << range.boundary_id
                   << ",\"name\":\"" << json_escape(range.name) << "\"}";
        }
        report << "],\n"
               << "  \"input\": \"" << json_escape(options.stl.string()) << "\",\n"
               << "  \"inputHashFnv1a64\": \"" << input_hash << "\",\n"
               << "  \"resultHashFnv1a64\": \"" << result_hash << "\",\n"
               << "  \"surfaceFormat\": \""
               << cartmesh::surface_format_name(surface.format()) << "\",\n"
               << "  \"triangleCount\": " << surface.triangles().size() << ",\n"
               << "  \"surfaceDiagnostics\": {"
               << "\"closed\":" << (diagnostics.closed ? "true" : "false")
               << ",\"manifold\":"
               << (diagnostics.manifold ? "true" : "false")
               << ",\"consistentlyOriented\":"
               << (diagnostics.consistently_oriented ? "true" : "false")
               << ",\"degenerateTriangleCount\":"
               << diagnostics.degenerate_triangle_count
               << ",\"duplicateTriangleCount\":"
               << diagnostics.duplicate_triangle_count
               << ",\"overlappingTrianglePairCount\":"
               << diagnostics.overlapping_triangle_pair_count
               << ",\"selfIntersectionPairCount\":"
               << diagnostics.self_intersection_pair_count
               << ",\"nonAdjacentContactPairCount\":"
               << diagnostics.non_adjacent_contact_pair_count
               << ",\"connectedComponentCount\":"
               << diagnostics.connected_component_count
               << ",\"smallComponentCount\":"
               << diagnostics.small_component_count << "},\n"
               << "  \"supportingPlaneCount\": " << supporting_plane_count << ",\n"
               << "  \"inputOrientationReversed\": "
               << (diagnostics.signed_volume < 0.0 ? "true" : "false") << ",\n"
               << "  \"domain\": {\"minimum\": ";
        write_vec3(report, domain.minimum());
        report << ", \"maximum\": ";
        write_vec3(report, domain.maximum());
        report << "},\n";
        if (tree) {
            const auto balance = tree->check_face_balance();
            report << "  \"baseLevel\": "
                   << static_cast<unsigned>(tree->base_level()) << ",\n"
                   << "  \"maximumLevel\": "
                   << static_cast<unsigned>(tree->maximum_level()) << ",\n"
                   << "  \"partitionValid\": "
                   << (tree->validate_partition() ? "true" : "false") << ",\n"
                   << "  \"faceBalanceValid\": "
                   << (balance.balanced ? "true" : "false") << ",\n"
                   << "  \"refinement\": {\"ruleSplitCount\":"
                   << adaptation.rule_refinement.split_count
                   << ",\"balanceSplitCount\":"
                   << adaptation.balance.split_count
                   << ",\"surfaceHits\":" << adaptation.surface_rule_hits
                   << ",\"distanceHits\":" << adaptation.distance_rule_hits
                   << ",\"curvatureHits\":"
                   << adaptation.curvature_rule_hits << ",\"gapHits\":"
                   << adaptation.gap_rule_hits
                   << ",\"gapResolutionFailureCount\":"
                   << adaptation.gap_resolution_failure_count
                   << ",\"maximumRequiredGapLevel\":"
                   << adaptation.maximum_required_gap_level
                   << ",\"userRegionHits\":"
                   << adaptation.user_region_rule_hits << "},\n"
                   << "  \"strictGapResolution\": "
                   << (options.strict_gap_resolution ? "true" : "false")
                   << ",\n";
            write_refinement_rules(report, options);
        } else {
            report << "  \"dimensions\": [" << grid->nx() << ", "
                   << grid->ny() << ", " << grid->nz() << "],\n"
                   << "  \"spacing\": ";
            write_vec3(report, grid->spacing());
            report << ",\n";
        }
        report << "  \"backgroundCellCount\": " << background_cell_count << ",\n"
               << "  \"fluidCellCount\": " << mesh.fluid_cells.size() << ",\n"
               << "  \"fullFluidCellCount\": " << mesh.full_fluid_cell_count << ",\n"
               << "  \"fullSolidCellCount\": " << mesh.full_solid_cell_count << ",\n"
               << "  \"cutCellCount\": " << mesh.cut_cell_count << ",\n"
               << "  \"volumeCutCellCount\": " << mesh.cut_cell_count << ",\n"
               << "  \"explicitFluidPieceCount\": "
               << explicit_piece_count << ",\n"
               << "  \"explicitFluidPieceVolume\": "
               << explicit_piece_volume << ",\n"
               << "  \"smallCellThreshold\": " << options.small_cell_threshold
               << ",\n"
               << "  \"minimumCutCellVolumeFraction\": "
               << small_cells.minimum_cut_cell_volume_fraction << ",\n"
               << "  \"smallCutCellCount\": " << small_cells.cells.size() << ",\n"
               << "  \"smallCutCells\": [";
        for (std::size_t index = 0; index < small_cells.cells.size(); ++index) {
            if (index != 0) report << ',';
            const auto& cell = small_cells.cells[index];
            report << "{\"backgroundCellId\":" << cell.background_cell_id
                   << ",\"volumeFraction\":" << cell.volume_fraction
                   << ",\"centroid\":";
            write_vec3(report, cell.centroid);
            report << ",\"boundaryIds\":[";
            for (std::size_t boundary = 0; boundary < cell.boundary_ids.size(); ++boundary) {
                if (boundary != 0) report << ',';
                report << cell.boundary_ids[boundary];
            }
            report << "]}";
        }
        report << "],\n"
               << "  \"boundaryCellCount\": " << mesh.boundary_cell_count << ",\n"
               << "  \"internalFluidFaceCount\": " << mesh.internal_faces.size() << ",\n"
               << "  \"componentInternalFaceCount\": "
               << mesh.component_internal_faces.size() << ",\n"
               << "  \"globalFluidRegionCount\": "
               << mesh.global_fluid_region_count << ",\n"
               << "  \"globalFluidRegions\": [";
        for (std::size_t region = 0;
             region < mesh.global_fluid_region_volumes.size(); ++region) {
            if (region != 0) report << ',';
            report << "{\"regionId\":" << region << ",\"volume\":"
                   << mesh.global_fluid_region_volumes[region]
                   << ",\"name\":\""
                   << json_escape(region_name(region)) << "\"}";
        }
        report << "],\n"
               << "  \"totalFluidVolume\": " << mesh.total_fluid_volume << ",\n"
               << "  \"totalSolidVolume\": " << solid_volume << ",\n"
               << "  \"domainVolume\": " << domain.volume() << ",\n"
               << "  \"volumeConservationResidual\": "
               << (mesh.total_fluid_volume + solid_volume - domain.volume()) << ",\n"
               << "  \"globalRegionVolumeSum\": " << region_volume_sum << ",\n"
               << "  \"totalEmbeddedBoundaryArea\": "
               << mesh.total_embedded_boundary_area << ",\n"
               << "  \"inputTriangleArea\": " << triangle_area << ",\n"
               << "  \"nonclosedCellCount\": " << mesh.nonclosed_cell_count << ",\n"
               << "  \"negativeVolumeCellCount\": "
               << mesh.negative_volume_cell_count << ",\n"
               << "  \"disconnectedFluidCellCount\": "
               << mesh.disconnected_fluid_cell_count << ",\n"
               << "  \"componentAnalysisPendingCellCount\": "
               << mesh.component_analysis_pending_cell_count << ",\n"
               << "  \"classificationConflictCount\": "
               << mesh.classification_conflict_count << ",\n"
               << "  \"sharedFaceMismatchCount\": "
               << mesh.shared_face_mismatch_count << ",\n"
               << "  \"discardedNumericalPieceCount\": "
               << mesh.discarded_numerical_piece_count << ",\n"
               << "  \"discardedNumericalPieceVolume\": "
               << mesh.discarded_numerical_piece_volume << ",\n"
               << "  \"maximumCellAreaClosureResidual\": "
               << mesh.maximum_cell_area_closure_residual << ",\n"
               << "  \"maximumSharedFaceAreaMismatch\": "
               << mesh.maximum_shared_face_area_mismatch << ",\n"
               << "  \"maximumSharedFaceCentroidMismatch\": "
               << mesh.maximum_shared_face_centroid_mismatch << ",\n"
               << "  \"maximumSharedFaceFirstMomentMismatch\": "
               << mesh.maximum_shared_face_first_moment_mismatch << ",\n"
               << "  \"wallSeconds\": " << wall_seconds << ",\n"
               << "  \"peakRssBytes\": " << peak_rss_bytes() << ",\n"
               << "  \"runtimeThreads\": 1,\n"
               << "  \"buildType\": \"" << build_type() << "\",\n"
               << "  \"compiler\": \"" << json_escape(compiler()) << "\",\n"
               << "  \"vtkWritten\": " << (options.write_vtk ? "true" : "false")
               << ",\n  \"tetrahedraOutput\": \""
               << json_escape(options.tetrahedra_output.string()) << "\",\n"
               << "  \"tetrahedraOutputScope\": "
                  "\"cut_cell_fluid_polyhedron_pieces_only\",\n"
               << "  \"polyhedraOutputScope\": "
                  "\"debug_cut_cell_fluid_polyhedron_pieces_only\",\n"
               << "  \"completeSolverVolumeMeshWritten\": "
               << (options.openfoam_case && stabilization_pass ? "true"
                                                               : "false")
               << ",\n"
               << "  \"solverMeshStabilizationEnabled\": "
               << (options.stabilize ? "true" : "false") << ",\n"
               << "  \"solverMeshStabilizationPass\": "
               << (stabilization_pass ? "true" : "false") << ",\n"
               << "  \"solverMeshAgglomerationCount\": "
               << (stabilization_report
                       ? stabilization_report->agglomeration_count
                       : 0U)
               << ",\n"
               << "  \"solverMeshRejectedCandidateCount\": "
               << (stabilization_report
                       ? stabilization_report->rejected_candidate_count
                       : 0U)
               << ",\n"
               << "  \"solverMeshRefinedLeafCount\": "
               << stabilization_refined_leaf_count << ",\n"
               << "  \"solverMeshBalanceSplitCount\": "
               << stabilization_balance_split_count << ",\n"
               << "  \"solverMeshUnresolvedSourceCount\": "
               << (stabilization_report
                       ? stabilization_report->unresolved_stable_ids.size()
                       : 0U)
               << ",\n"
               << "  \"solverMeshStabilizationReport\": ";
        if (stabilization_output) {
            report << "\"" << json_escape(stabilization_output->string())
                   << "\",\n";
        } else {
            report << "null,\n";
        }
        report
               << "  \"nativeSolverQualityEvaluated\": "
               << (quality_report ? "true" : "false") << ",\n"
               << "  \"nativeSolverTopologyPass\": "
               << (quality_report && quality_report->topology_pass()
                       ? "true" : "false") << ",\n"
               << "  \"nativeSolverQualityPass\": "
               << (quality_report && quality_report->quality_pass()
                       ? "true" : "false") << ",\n"
               << "  \"nativeSolverQualityIssueCount\": "
               << (quality_report ? quality_report->summary.issue_count : 0U)
               << ",\n"
               << "  \"nativeSolverMaximumNonOrthogonalityDegrees\": "
               << (quality_report
                       ? quality_report->summary.maximum_non_orthogonality_degrees
                       : 0.0) << ",\n"
               << "  \"nativeSolverMaximumSkewness\": "
               << (quality_report ? quality_report->summary.maximum_skewness
                                  : 0.0) << ",\n"
               << "  \"nativeSolverQualityReport\": ";
        if (quality_output) {
            report << "\"" << json_escape(quality_output->string()) << "\",\n";
        } else {
            report << "null,\n";
        }
        report << "  \"completeSolverVolumeMeshOutput\": ";
        if (options.openfoam_case && stabilization_pass) {
            report << "\"" << json_escape(options.openfoam_case->string())
                   << "/constant/polyMesh\"\n";
        } else {
            report << "null\n";
        }
        report << "}\n";
        report.flush();
        if (!report) {
            throw std::runtime_error("写入阶段四报告失败：" + options.report.string());
        }
        std::cout << "阶段4 STL Cut-cell：kernel=" << kernel
                  << " 背景单元=" << background_cell_count
                  << " cut=" << mesh.cut_cell_count
                  << " solid=" << mesh.full_solid_cell_count
                  << " fluidVolume=" << std::setprecision(17)
                  << mesh.total_fluid_volume << " closureFailures="
                  << mesh.nonclosed_cell_count << " sharedFaceFailures="
                  << mesh.shared_face_mismatch_count << " hash=" << result_hash << '\n';
        return invariants_pass && stabilization_pass ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "错误：" << error.what() << '\n';
        return 1;
    }
}
