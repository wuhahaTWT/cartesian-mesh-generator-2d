#include "cartmesh/classify/SurfaceClassifier.hpp"
#include "cartmesh/geometry/AnalyticShapes.hpp"
#include "cartmesh/geometry/SurfaceDiagnostics.hpp"
#include "cartmesh/grid/UniformCartesianGrid.hpp"
#include "cartmesh/io/StlReader.hpp"
#include "cartmesh/io/DiagnosticVtkWriter.hpp"
#include "cartmesh/io/VtkWriter.hpp"
#include "cartmesh/spatial/TriangleBvh.hpp"
#include "cartmesh/util/Performance.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
    std::string case_name{"cube"};
    std::optional<std::filesystem::path> stl;
    std::uint32_t nx{32};
    std::uint32_t ny{32};
    std::uint32_t nz{32};
    std::optional<std::filesystem::path> output{"artifacts/cube_32.vtu"};
    std::filesystem::path report{"artifacts/cube_32.json"};
    double padding_fraction{0.05};
};

[[nodiscard]] double seconds_between(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double>(end - start).count();
}

[[nodiscard]] std::uint32_t parse_dimension(std::string_view value, std::string_view option) {
    std::size_t parsed = 0;
    const auto number = std::stoull(std::string(value), &parsed);
    if (parsed != value.size() || number == 0 ||
        number > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument(std::string(option) + " 需要 [1, 2^32-1] 范围内的整数");
    }
    return static_cast<std::uint32_t>(number);
}

[[nodiscard]] double parse_nonnegative_double(std::string_view value, std::string_view option) {
    std::size_t parsed = 0;
    const double number = std::stod(std::string(value), &parsed);
    if (parsed != value.size() || !std::isfinite(number) || number < 0.0) {
        throw std::invalid_argument(std::string(option) + " 需要非负有限数");
    }
    return number;
}

[[nodiscard]] std::string json_escape(std::string_view value) {
    std::ostringstream output;
    for (const char character : value) {
        switch (character) {
        case '\\':
            output << "\\\\";
            break;
        case '\"':
            output << "\\\"";
            break;
        case '\n':
            output << "\\n";
            break;
        case '\r':
            output << "\\r";
            break;
        case '\t':
            output << "\\t";
            break;
        default:
            output << character;
        }
    }
    return output.str();
}

void print_help() {
    std::cout
        << "用法：cartmesh_cli [选项]\n"
        << "  --case cube|sphere       阶段 0 解析案例（默认：cube）\n"
        << "  --stl FILE.stl           阶段 1 封闭、定向、流形 STL 输入\n"
        << "  --padding-fraction X     STL 包围盒相对留白（默认：0.05）\n"
        << "  --resolution N           设置 nx=ny=nz=N\n"
        << "  --nx N --ny N --nz N     分别设置三个方向的网格数\n"
        << "  --output FILE.vtu         VTK XML 输出路径\n"
        << "  --report FILE.json        机器可读的测量报告\n"
        << "  --no-vtk                  运行但不写出可视化文件\n"
        << "  --help                    显示本帮助信息\n";
}

[[nodiscard]] Options parse_options(int argc, char** argv) {
    Options options;
    bool output_was_set = false;
    bool report_was_set = false;
    bool case_was_set = false;
    for (int argument = 1; argument < argc; ++argument) {
        const std::string_view current(argv[argument]);
        const auto require_value = [&](std::string_view option) -> std::string_view {
            if (argument + 1 >= argc) {
                throw std::invalid_argument(std::string(option) + " 需要提供一个值");
            }
            ++argument;
            return argv[argument];
        };
        if (current == "--help") {
            print_help();
            std::exit(0);
        }
        if (current == "--case") {
            options.case_name = require_value(current);
            case_was_set = true;
        } else if (current == "--stl") {
            options.stl = std::filesystem::path(require_value(current));
        } else if (current == "--padding-fraction") {
            options.padding_fraction = parse_nonnegative_double(require_value(current), current);
        } else if (current == "--resolution") {
            const auto value = parse_dimension(require_value(current), current);
            options.nx = value;
            options.ny = value;
            options.nz = value;
        } else if (current == "--nx") {
            options.nx = parse_dimension(require_value(current), current);
        } else if (current == "--ny") {
            options.ny = parse_dimension(require_value(current), current);
        } else if (current == "--nz") {
            options.nz = parse_dimension(require_value(current), current);
        } else if (current == "--output") {
            options.output = std::filesystem::path(require_value(current));
            output_was_set = true;
        } else if (current == "--report") {
            options.report = std::filesystem::path(require_value(current));
            report_was_set = true;
        } else if (current == "--no-vtk") {
            options.output.reset();
            output_was_set = true;
        } else {
            throw std::invalid_argument("未知选项：" + std::string(current));
        }
    }
    if (options.case_name != "cube" && options.case_name != "sphere") {
        throw std::invalid_argument("--case 必须是 cube 或 sphere");
    }
    if (options.stl && case_was_set) {
        throw std::invalid_argument("--stl 与 --case 不能同时使用");
    }
    const auto suffix = std::to_string(options.nx) + "x" + std::to_string(options.ny) + "x" +
                        std::to_string(options.nz);
    const auto output_stem = options.stl ? options.stl->stem().string() : options.case_name;
    if (!output_was_set) {
        options.output = std::filesystem::path("artifacts") /
                         (output_stem + "_" + suffix + ".vtu");
    }
    if (!report_was_set) {
        options.report = std::filesystem::path("artifacts") /
                         (output_stem + "_" + suffix + ".json");
    }
    return options;
}

void hash_byte(std::uint64_t& hash, std::uint8_t value) noexcept {
    hash ^= value;
    hash *= 1099511628211ULL;
}

void hash_u64(std::uint64_t& hash, std::uint64_t value) noexcept {
    for (unsigned shift = 0; shift < 64; shift += 8) {
        hash_byte(hash, static_cast<std::uint8_t>((value >> shift) & 0xffULL));
    }
}

void hash_double(std::uint64_t& hash, double value) noexcept {
    hash_u64(hash, std::bit_cast<std::uint64_t>(value));
}

[[nodiscard]] std::string hash_hex(std::uint64_t hash) {
    std::ostringstream output;
    output << std::hex << std::setfill('0') << std::setw(16) << hash;
    return output.str();
}

[[nodiscard]] std::string file_hash_fnv1a64(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("无法为 STL 计算输入哈希：" + path.string());
    }
    std::uint64_t hash = 14695981039346656037ULL;
    std::array<char, 64 * 1024> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const auto count = input.gcount();
        for (std::streamsize index = 0; index < count; ++index) {
            hash_byte(hash, static_cast<std::uint8_t>(buffer[static_cast<std::size_t>(index)]));
        }
    }
    if (!input.eof()) {
        throw std::runtime_error("读取 STL 输入哈希失败：" + path.string());
    }
    return hash_hex(hash);
}

[[nodiscard]] constexpr std::string_view build_type_name() noexcept {
#if defined(NDEBUG)
    return "Release";
#else
    return "Debug";
#endif
}

[[nodiscard]] constexpr std::string_view compiler_name() noexcept {
#if defined(__clang__)
    return "Clang " __clang_version__;
#elif defined(__GNUC__)
    return "GCC " __VERSION__;
#elif defined(_MSC_VER)
    return "MSVC";
#else
    return "unknown";
#endif
}

void write_vec3_json(std::ostream& output, const cartmesh::Vec3& value) {
    output << '[' << value.x << ',' << value.y << ',' << value.z << ']';
}

void write_id_examples(std::ostream& output, const std::vector<std::uint64_t>& examples) {
    output << '[';
    for (std::size_t index = 0; index < examples.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        output << examples[index];
    }
    output << ']';
}

void write_edge_examples(
    std::ostream& output,
    const std::vector<cartmesh::SurfaceDiagnosticLocation>& examples) {
    output << '[';
    for (std::size_t index = 0; index < examples.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        output << "{\"firstTriangle\":" << examples[index].first_triangle
               << ",\"secondTriangle\":" << examples[index].second_triangle
               << ",\"start\":";
        write_vec3_json(output, examples[index].edge_start);
        output << ",\"end\":";
        write_vec3_json(output, examples[index].edge_end);
        output << '}';
    }
    output << ']';
}

void write_vertex_examples(
    std::ostream& output,
    const std::vector<cartmesh::SurfaceDiagnosticVertex>& examples) {
    output << '[';
    for (std::size_t index = 0; index < examples.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        output << "{\"position\":";
        write_vec3_json(output, examples[index].position);
        output << ",\"incidentTriangleCount\":" << examples[index].incident_triangle_count
               << '}';
    }
    output << ']';
}

void write_triangle_pair_examples(
    std::ostream& output,
    const std::vector<cartmesh::SurfaceDiagnosticTrianglePair>& examples) {
    output << '[';
    for (std::size_t index = 0; index < examples.size(); ++index) {
        if (index != 0) output << ',';
        const auto& example = examples[index];
        output << "{\"firstTriangle\":" << example.first_triangle
               << ",\"secondTriangle\":" << example.second_triangle
               << ",\"firstPosition\":";
        write_vec3_json(output, example.first_position);
        output << ",\"secondPosition\":";
        write_vec3_json(output, example.second_position);
        output << '}';
    }
    output << ']';
}

void write_component_diagnostics(
    std::ostream& output,
    const std::vector<cartmesh::SurfaceDiagnosticComponent>& components) {
    output << '[';
    for (std::size_t index = 0; index < components.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        const auto& component = components[index];
        output << "{\"componentId\":" << component.component_id
               << ",\"triangleCount\":" << component.triangle_count
               << ",\"surfaceArea\":" << component.surface_area
               << ",\"bounds\":{\"minimum\":";
        write_vec3_json(output, component.bounds.minimum());
        output << ",\"maximum\":";
        write_vec3_json(output, component.bounds.maximum());
        output << "},\"boundingDiagonal\":" << component.bounding_diagonal
               << ",\"signedVolume\":" << component.signed_volume
               << ",\"nestingDepth\":" << component.nesting_depth
               << ",\"expectedOrientationSign\":"
               << component.expected_orientation_sign
               << ",\"orientationChecked\":" << std::boolalpha
               << component.orientation_checked
               << ",\"orientationMatchesNesting\":"
               << component.orientation_matches_nesting << ",\"samplePosition\":";
        write_vec3_json(output, component.sample_position);
        output << '}';
    }
    output << ']';
}

void write_surface_diagnostics(std::ostream& output,
                               const cartmesh::SurfaceDiagnostics& diagnostics,
                               std::string_view indentation) {
    output << indentation << "{\n"
           << indentation << "  \"triangleCount\": " << diagnostics.triangle_count << ",\n"
           << indentation << "  \"uniqueVertexCount\": " << diagnostics.unique_vertex_count
           << ",\n"
           << indentation << "  \"uniqueEdgeCount\": " << diagnostics.unique_edge_count
           << ",\n"
           << indentation << "  \"degenerateTriangleCount\": "
           << diagnostics.degenerate_triangle_count << ",\n"
           << indentation << "  \"duplicateTriangleCount\": "
           << diagnostics.duplicate_triangle_count << ",\n"
           << indentation << "  \"boundaryEdgeCount\": " << diagnostics.boundary_edge_count
           << ",\n"
           << indentation << "  \"nonManifoldEdgeCount\": "
           << diagnostics.non_manifold_edge_count << ",\n"
           << indentation << "  \"nonManifoldVertexCount\": "
           << diagnostics.non_manifold_vertex_count << ",\n"
           << indentation << "  \"orientationConflictEdgeCount\": "
           << diagnostics.orientation_conflict_edge_count << ",\n"
           << indentation << "  \"connectedComponentCount\": "
           << diagnostics.connected_component_count << ",\n"
           << indentation << "  \"componentOrientationMismatchCount\": "
           << diagnostics.component_orientation_mismatch_count << ",\n"
           << indentation << "  \"overlappingTrianglePairCount\": "
           << diagnostics.overlapping_triangle_pair_count << ",\n"
           << indentation << "  \"selfIntersectionPairCount\": "
           << diagnostics.self_intersection_pair_count << ",\n"
           << indentation << "  \"nonAdjacentContactPairCount\": "
           << diagnostics.non_adjacent_contact_pair_count << ",\n"
           << indentation << "  \"smallComponentCount\": "
           << diagnostics.small_component_count << ",\n"
           << indentation << "  \"suggestedLengthTolerance\": "
           << diagnostics.suggested_length_tolerance << ",\n"
           << indentation << "  \"degenerateAreaTolerance\": "
           << diagnostics.degenerate_area_tolerance << ",\n"
           << indentation << "  \"smallComponentDiagonalThreshold\": "
           << diagnostics.small_component_diagonal_threshold << ",\n"
           << indentation << "  \"minimumComponentBoundingDiagonal\": "
           << diagnostics.minimum_component_bounding_diagonal << ",\n"
           << indentation << "  \"minimumComponentAbsoluteVolume\": "
           << diagnostics.minimum_component_absolute_volume << ",\n"
           << indentation << "  \"signedVolume\": " << diagnostics.signed_volume << ",\n"
           << indentation << "  \"materialVolume\": " << diagnostics.material_volume << ",\n"
           << indentation << "  \"closed\": " << std::boolalpha << diagnostics.closed << ",\n"
           << indentation << "  \"manifold\": " << diagnostics.manifold << ",\n"
           << indentation << "  \"consistentlyOriented\": "
           << diagnostics.consistently_oriented << ",\n"
           << indentation << "  \"degenerateTriangleExamples\": ";
    write_id_examples(output, diagnostics.degenerate_triangle_examples);
    output << ",\n" << indentation << "  \"duplicateTriangleExamples\": ";
    write_id_examples(output, diagnostics.duplicate_triangle_examples);
    output << ",\n" << indentation << "  \"boundaryEdgeExamples\": ";
    write_edge_examples(output, diagnostics.boundary_edge_examples);
    output << ",\n" << indentation << "  \"nonManifoldEdgeExamples\": ";
    write_edge_examples(output, diagnostics.non_manifold_edge_examples);
    output << ",\n" << indentation << "  \"orientationConflictExamples\": ";
    write_edge_examples(output, diagnostics.orientation_conflict_examples);
    output << ",\n" << indentation << "  \"nonManifoldVertexExamples\": ";
    write_vertex_examples(output, diagnostics.non_manifold_vertex_examples);
    output << ",\n" << indentation << "  \"overlappingTriangleExamples\": ";
    write_triangle_pair_examples(output, diagnostics.overlapping_triangle_examples);
    output << ",\n" << indentation << "  \"selfIntersectionExamples\": ";
    write_triangle_pair_examples(output, diagnostics.self_intersection_examples);
    output << ",\n" << indentation << "  \"nonAdjacentContactExamples\": ";
    write_triangle_pair_examples(output, diagnostics.non_adjacent_contact_examples);
    output << ",\n" << indentation << "  \"components\": ";
    write_component_diagnostics(output, diagnostics.components);
    output << ",\n"
           << indentation << "  \"validForStage1Classification\": "
           << diagnostics.valid_for_stage1_classification() << "\n"
           << indentation << "}";
}

[[nodiscard]] cartmesh::AABB padded_domain(const cartmesh::AABB& surface_bounds,
                                           double padding_fraction) {
    const auto extent = surface_bounds.extent();
    const double reference = std::max({extent.x, extent.y, extent.z});
    if (!(reference > 0.0) || !std::isfinite(reference)) {
        throw std::runtime_error("STL 包围盒必须至少在一个方向具有正尺度");
    }
    const double padding = reference * padding_fraction;
    const double minimum_padding =
        std::max(reference * 1.0e-9, 64.0 * std::numeric_limits<double>::epsilon() * reference);
    const double applied = std::max(padding, minimum_padding);
    const cartmesh::Vec3 delta{applied, applied, applied};
    return cartmesh::AABB(surface_bounds.minimum() - delta, surface_bounds.maximum() + delta);
}

void ensure_report_parent(const std::filesystem::path& path) {
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
}

[[nodiscard]] std::filesystem::path diagnostic_marker_path(const Options& options) {
    auto path = options.report;
    path.replace_extension();
    path += "_geometry_issues.vtp";
    return path;
}

void write_rejected_stl_report(const Options& options, const cartmesh::SurfaceMesh& surface,
                               const cartmesh::SurfaceDiagnostics& diagnostics,
                               const std::filesystem::path& marker_path,
                               std::string_view input_hash, double import_seconds,
                               double diagnostic_seconds, double total_seconds) {
    ensure_report_parent(options.report);
    std::ofstream report(options.report, std::ios::trunc);
    if (!report) {
        throw std::runtime_error("无法打开报告输出文件：" + options.report.string());
    }
    report << std::setprecision(17);
    report << "{\n"
           << "  \"schemaVersion\": 3,\n"
           << "  \"projectStage\": 1,\n"
           << "  \"status\": \"rejected_invalid_surface\",\n"
           << "  \"input\": \"" << json_escape(options.stl->string()) << "\",\n"
           << "  \"inputFormat\": \"" << cartmesh::surface_format_name(surface.format())
           << "\",\n"
           << "  \"inputHashFnv1a64\": \"" << input_hash << "\",\n"
           << "  \"solverReadyCutCellMesh\": false,\n"
           << "  \"diagnosticMarkerVtp\": \"" << json_escape(marker_path.string())
           << "\",\n"
           << "  \"surfaceDiagnostics\": ";
    write_surface_diagnostics(report, diagnostics, "  ");
    report << ",\n"
           << "  \"threads\": 1,\n"
           << "  \"buildType\": \"" << build_type_name() << "\",\n"
           << "  \"compiler\": \"" << json_escape(compiler_name()) << "\",\n"
           << "  \"peakRssBytes\": " << cartmesh::peak_rss_bytes() << ",\n"
           << "  \"timingsSeconds\": {\"stlImport\": " << import_seconds
           << ", \"surfaceDiagnostics\": " << diagnostic_seconds << ", \"total\": "
           << total_seconds << "}\n"
           << "}\n";
}

void write_report(const Options& options, const cartmesh::UniformCartesianGrid& grid,
                  std::uint64_t selected_cells, double exact_volume, double estimated_volume,
                  std::string_view classification_method, std::string_view result_hash,
                  double grid_seconds, double classification_seconds, double vtk_seconds,
                  double total_seconds) {
    if (options.report.has_parent_path()) {
        std::filesystem::create_directories(options.report.parent_path());
    }
    std::ofstream report(options.report, std::ios::trunc);
    if (!report) {
        throw std::runtime_error("无法打开报告输出文件：" + options.report.string());
    }
    const double error = std::abs(estimated_volume - exact_volume);
    report << std::setprecision(17);
    report << "{\n"
           << "  \"schemaVersion\": 1,\n"
           << "  \"projectStage\": 0,\n"
           << "  \"case\": \"" << json_escape(options.case_name) << "\",\n"
           << "  \"classificationMethod\": \"" << classification_method << "\",\n"
           << "  \"solverReadyCutCellMesh\": false,\n"
           << "  \"dimensions\": {\"nx\": " << grid.nx() << ", \"ny\": " << grid.ny()
           << ", \"nz\": " << grid.nz() << "},\n"
           << "  \"domain\": {\"minimum\": ";
    write_vec3_json(report, grid.domain().minimum());
    report << ", \"maximum\": ";
    write_vec3_json(report, grid.domain().maximum());
    report << "},\n"
           << "  \"cellSpacing\": ";
    write_vec3_json(report, grid.spacing());
    report << ",\n"
           << "  \"pointCount\": " << grid.point_count() << ",\n"
           << "  \"cellCount\": " << grid.cell_count() << ",\n"
           << "  \"selectedCellCount\": " << selected_cells << ",\n"
           << "  \"exactVolume\": " << exact_volume << ",\n"
           << "  \"estimatedVolume\": " << estimated_volume << ",\n"
           << "  \"absoluteVolumeError\": " << error << ",\n"
           << "  \"relativeVolumeError\": " << error / exact_volume << ",\n"
           << "  \"resultHashFnv1a64\": \"" << result_hash << "\",\n"
           << "  \"threads\": 1,\n"
           << "  \"peakRssBytes\": " << cartmesh::peak_rss_bytes() << ",\n"
           << "  \"timingsSeconds\": {\n"
           << "    \"gridConstruction\": " << grid_seconds << ",\n"
           << "    \"classification\": " << classification_seconds << ",\n"
           << "    \"vtkWrite\": " << vtk_seconds << ",\n"
           << "    \"total\": " << total_seconds << "\n"
           << "  },\n"
           << "  \"vtkOutput\": ";
    if (options.output) {
        report << "\"" << json_escape(options.output->string()) << "\"\n";
    } else {
        report << "null\n";
    }
    report << "}\n";
    if (!report) {
        throw std::runtime_error("写入报告失败：" + options.report.string());
    }
}

void write_stage1_report(const Options& options, const cartmesh::SurfaceMesh& surface,
                         const cartmesh::SurfaceDiagnostics& diagnostics,
                         const cartmesh::TriangleBvh& bvh,
                         const cartmesh::UniformCartesianGrid& grid,
                         const cartmesh::UniformClassification& classification,
                         std::string_view input_hash, std::string_view result_hash,
                         double estimated_volume, double surface_tolerance, double import_seconds,
                         double diagnostic_seconds, double bvh_seconds, double grid_seconds,
                         double classification_seconds, double vtk_seconds, double total_seconds) {
    ensure_report_parent(options.report);
    std::ofstream report(options.report, std::ios::trunc);
    if (!report) {
        throw std::runtime_error("无法打开报告输出文件：" + options.report.string());
    }
    const double exact_polyhedral_volume = diagnostics.material_volume;
    const double error = std::abs(estimated_volume - exact_polyhedral_volume);
    const auto unresolved = classification.conflict_count;
    const double definitely_inside_volume =
        static_cast<double>(classification.inside_count) * grid.cell_volume();
    const double inside_plus_intersected_volume =
        static_cast<double>(classification.inside_count + classification.intersected_count) *
        grid.cell_volume();
    const bool volume_bracket_pass =
        definitely_inside_volume <= exact_polyhedral_volume &&
        exact_polyhedral_volume <= inside_plus_intersected_volume;
    const auto& bvh_statistics = bvh.statistics();
    report << std::setprecision(17);
    report << "{\n"
           << "  \"schemaVersion\": 3,\n"
           << "  \"projectStage\": 1,\n"
           << "  \"status\": \""
           << (unresolved != 0
                   ? "classification_unresolved"
                   : (volume_bracket_pass ? "pass" : "failed_invariant"))
           << "\",\n"
           << "  \"input\": \"" << json_escape(options.stl->string()) << "\",\n"
           << "  \"inputFormat\": \"" << cartmesh::surface_format_name(surface.format())
           << "\",\n"
           << "  \"inputHashFnv1a64\": \"" << input_hash << "\",\n"
           << "  \"classificationMethod\": "
              "\"bvh_exact_triangle_aabb_then_three_direction_parity\",\n"
           << "  \"intersectedUsesExactTriangleAabbSat\": true,\n"
           << "  \"insideOutsideUsesCellCenterSample\": true,\n"
           << "  \"solverReadyCutCellMesh\": false,\n"
           << "  \"warnings\": ["
           << (diagnostics.component_orientation_mismatch_count == 0
                   ? ""
                   : "\"componentOrientationDoesNotMatchNesting\"")
           << "],\n"
           << "  \"surfaceDiagnostics\": ";
    write_surface_diagnostics(report, diagnostics, "  ");
    report << ",\n"
           << "  \"bvh\": {\"triangleCount\": " << bvh_statistics.triangle_count
           << ", \"nodeCount\": " << bvh_statistics.node_count << ", \"leafCount\": "
           << bvh_statistics.leaf_count << ", \"maximumDepth\": "
           << bvh_statistics.maximum_depth << ", \"maximumLeafTriangles\": "
           << bvh_statistics.maximum_leaf_triangles << "},\n"
           << "  \"dimensions\": {\"nx\": " << grid.nx() << ", \"ny\": " << grid.ny()
           << ", \"nz\": " << grid.nz() << "},\n"
           << "  \"domain\": {\"minimum\": ";
    write_vec3_json(report, grid.domain().minimum());
    report << ", \"maximum\": ";
    write_vec3_json(report, grid.domain().maximum());
    report << "},\n"
           << "  \"cellSpacing\": ";
    write_vec3_json(report, grid.spacing());
    report << ",\n"
           << "  \"paddingFraction\": " << options.padding_fraction << ",\n"
           << "  \"surfaceTolerance\": " << surface_tolerance << ",\n"
           << "  \"pointCount\": " << grid.point_count() << ",\n"
           << "  \"cellCount\": " << grid.cell_count() << ",\n"
           << "  \"classificationCounts\": {\"outside\": " << classification.outside_count
           << ", \"inside\": " << classification.inside_count << ", \"intersected\": "
           << classification.intersected_count << ", \"conflict\": "
           << classification.conflict_count << "},\n"
           << "  \"classificationValueLegend\": {\"outside\": 0, \"inside\": 1, "
              "\"intersected\": 2, \"conflict\": 3},\n"
           << "  \"centerPointCounts\": {\"outside\": "
           << classification.center_outside_count << ", \"inside\": "
           << classification.center_inside_count << ", \"onSurface\": "
           << classification.center_on_surface_count << ", \"conflict\": "
           << classification.center_conflict_count << "},\n"
           << "  \"centerSampleVolumeValid\": "
           << (classification.center_on_surface_count == 0 &&
               classification.center_conflict_count == 0)
           << ",\n"
           << "  \"exactPolyhedralVolume\": " << exact_polyhedral_volume << ",\n"
           << "  \"estimatedCenterSampleVolume\": " << estimated_volume << ",\n"
           << "  \"definitelyInsideVolumeLowerBound\": " << definitely_inside_volume << ",\n"
           << "  \"insidePlusIntersectedVolumeUpperBound\": "
           << inside_plus_intersected_volume << ",\n"
           << "  \"polyhedralVolumeInsideClassificationBounds\": "
           << std::boolalpha << volume_bracket_pass
           << ",\n"
           << "  \"absoluteVolumeError\": " << error << ",\n"
           << "  \"relativeVolumeError\": "
           << (exact_polyhedral_volume > 0.0 ? error / exact_polyhedral_volume : 0.0) << ",\n"
           << "  \"resultHashFnv1a64\": \"" << result_hash << "\",\n"
           << "  \"threads\": 1,\n"
           << "  \"buildType\": \"" << build_type_name() << "\",\n"
           << "  \"compiler\": \"" << json_escape(compiler_name()) << "\",\n"
           << "  \"peakRssBytes\": " << cartmesh::peak_rss_bytes() << ",\n"
           << "  \"timingsSeconds\": {\n"
           << "    \"stlImport\": " << import_seconds << ",\n"
           << "    \"surfaceDiagnostics\": " << diagnostic_seconds << ",\n"
           << "    \"bvhBuild\": " << bvh_seconds << ",\n"
           << "    \"gridConstruction\": " << grid_seconds << ",\n"
           << "    \"classification\": " << classification_seconds << ",\n"
           << "    \"vtkWrite\": " << vtk_seconds << ",\n"
           << "    \"total\": " << total_seconds << "\n"
           << "  },\n"
           << "  \"vtkOutput\": ";
    if (options.output) {
        report << "\"" << json_escape(options.output->string()) << "\"\n";
    } else {
        report << "null\n";
    }
    report << "}\n";
    if (!report) {
        throw std::runtime_error("写入阶段 1 报告失败：" + options.report.string());
    }
}

int run_stl(const Options& options) {
    const auto total_start = Clock::now();
    const auto input_hash = file_hash_fnv1a64(*options.stl);
    const auto import_start = Clock::now();
    const auto surface = cartmesh::read_stl(*options.stl);
    const auto import_end = Clock::now();
    const auto diagnostic_start = Clock::now();
    const auto diagnostics = cartmesh::diagnose_surface(surface);
    const auto diagnostic_end = Clock::now();
    if (!diagnostics.valid_for_stage1_classification()) {
        const auto marker_path = diagnostic_marker_path(options);
        cartmesh::write_surface_diagnostic_vtp(marker_path, surface, diagnostics);
        write_rejected_stl_report(options, surface, diagnostics, marker_path, input_hash,
                                  seconds_between(import_start, import_end),
                                  seconds_between(diagnostic_start, diagnostic_end),
                                  seconds_between(total_start, Clock::now()));
        throw std::runtime_error("STL 未通过阶段 1 封闭、定向、流形诊断；详情见 " +
                                 options.report.string());
    }

    const auto bvh_start = Clock::now();
    const cartmesh::TriangleBvh bvh(surface);
    const auto bvh_end = Clock::now();
    const auto grid_start = Clock::now();
    const cartmesh::UniformCartesianGrid grid(padded_domain(surface.bounds(),
                                                            options.padding_fraction),
                                              options.nx, options.ny, options.nz);
    const auto grid_end = Clock::now();
    const cartmesh::SurfaceClassifier classifier(bvh);
    const auto classification_start = Clock::now();
    const auto classification = cartmesh::classify_uniform_cells(grid, classifier);
    const auto classification_end = Clock::now();

    std::uint64_t hash = 14695981039346656037ULL;
    hash_u64(hash, options.nx);
    hash_u64(hash, options.ny);
    hash_u64(hash, options.nz);
    hash_double(hash, grid.domain().minimum().x);
    hash_double(hash, grid.domain().minimum().y);
    hash_double(hash, grid.domain().minimum().z);
    hash_double(hash, grid.domain().maximum().x);
    hash_double(hash, grid.domain().maximum().y);
    hash_double(hash, grid.domain().maximum().z);
    for (std::size_t index = 0; index < classification.cell_classification.size(); ++index) {
        hash_byte(hash, classification.cell_classification[index]);
        hash_byte(hash, classification.center_point_classification[index]);
    }
    const auto result_hash = hash_hex(hash);
    const double estimated_volume =
        static_cast<double>(classification.center_inside_count) * grid.cell_volume();

    const auto vtk_start = Clock::now();
    if (options.output) {
        std::vector<double> classification_field(classification.cell_classification.size());
        std::vector<double> inside_field(classification.cell_classification.size());
        for (std::size_t index = 0; index < classification.cell_classification.size(); ++index) {
            classification_field[index] =
                static_cast<double>(classification.cell_classification[index]);
            inside_field[index] =
                classification.center_point_classification[index] ==
                        static_cast<std::uint8_t>(cartmesh::PointClassification::inside)
                    ? 1.0
                    : 0.0;
        }
        cartmesh::write_vtu(*options.output, grid,
                            {{"stl_cell_classification", std::move(classification_field)},
                             {"inside_stl_center_sample", std::move(inside_field)}});
    }
    const auto vtk_end = Clock::now();
    write_stage1_report(options, surface, diagnostics, bvh, grid, classification, input_hash,
                        result_hash, estimated_volume, classifier.surface_tolerance(),
                        seconds_between(import_start, import_end),
                        seconds_between(diagnostic_start, diagnostic_end),
                        seconds_between(bvh_start, bvh_end), seconds_between(grid_start, grid_end),
                        seconds_between(classification_start, classification_end),
                        seconds_between(vtk_start, vtk_end),
                        seconds_between(total_start, Clock::now()));

    const auto unresolved = classification.conflict_count;
    const double exact_polyhedral_volume = diagnostics.material_volume;
    const double definitely_inside_volume =
        static_cast<double>(classification.inside_count) * grid.cell_volume();
    const double inside_plus_intersected_volume =
        static_cast<double>(classification.inside_count + classification.intersected_count) *
        grid.cell_volume();
    const bool volume_bracket_pass =
        definitely_inside_volume <= exact_polyhedral_volume &&
        exact_polyhedral_volume <= inside_plus_intersected_volume;
    std::cout << std::setprecision(10) << "STL=" << options.stl->string()
              << " 格式=" << cartmesh::surface_format_name(surface.format())
              << " 三角形=" << surface.triangles().size() << " 单元=" << grid.cell_count()
              << " 内部=" << classification.inside_count << " 外部="
              << classification.outside_count << " 相交=" << classification.intersected_count
              << " 冲突=" << classification.conflict_count << " 中心估算体积="
              << estimated_volume
              << " 多面体体积=" << diagnostics.material_volume << " 哈希=" << result_hash
              << " 峰值RSS字节=" << cartmesh::peak_rss_bytes()
              << " 总秒数=" << seconds_between(total_start, Clock::now()) << '\n';
    if (unresolved != 0 || !volume_bracket_pass) {
        std::cerr << "cartmesh_cli 错误：存在 " << unresolved
                  << " 个分类冲突，体积包络=" << volume_bracket_pass
                  << "；详情已写入报告\n";
        return 2;
    }
    return 0;
}

} // 匿名命名空间

int main(int argc, char** argv) {
    try {
        const auto options = parse_options(argc, argv);
        if (options.stl) {
            return run_stl(options);
        }
        const auto total_start = Clock::now();

        const bool is_sphere = options.case_name == "sphere";
        const cartmesh::AABB domain =
            is_sphere ? cartmesh::AABB({-1.25, -1.25, -1.25}, {1.25, 1.25, 1.25})
                      : cartmesh::AABB({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5});
        const auto grid_start = Clock::now();
        const cartmesh::UniformCartesianGrid grid(domain, options.nx, options.ny, options.nz);
        const auto grid_end = Clock::now();

        const cartmesh::AnalyticSphere sphere({0.0, 0.0, 0.0}, 1.0);
        std::vector<double> region(grid.cell_count(), 1.0);
        std::vector<double> signed_distance;
        if (is_sphere) {
            signed_distance.resize(grid.cell_count());
        }
        std::uint64_t selected_cells = 0;
        std::uint64_t hash = 14695981039346656037ULL;
        hash_u64(hash, options.nx);
        hash_u64(hash, options.ny);
        hash_u64(hash, options.nz);
        hash_byte(hash, static_cast<std::uint8_t>(is_sphere));
        const auto classify_start = Clock::now();
        for (std::uint64_t id = 0; id < grid.cell_count(); ++id) {
            const auto center = grid.cell_center(grid.cell_key(id));
            const bool selected = !is_sphere || sphere.contains(center);
            region[static_cast<std::size_t>(id)] = selected ? 1.0 : 0.0;
            if (is_sphere) {
                signed_distance[static_cast<std::size_t>(id)] = sphere.signed_distance(center);
            }
            selected_cells += static_cast<std::uint64_t>(selected);
            hash_byte(hash, static_cast<std::uint8_t>(selected));
        }
        const auto classify_end = Clock::now();

        const double exact_volume = is_sphere ? sphere.volume() : domain.volume();
        const double estimated_volume = static_cast<double>(selected_cells) * grid.cell_volume();
        const auto vtk_start = Clock::now();
        if (options.output) {
            std::vector<cartmesh::VtkCellData> fields;
            fields.push_back({is_sphere ? "inside_sphere_center_sample" : "inside_analytic_cube",
                              std::move(region)});
            if (is_sphere) {
                fields.push_back({"signed_distance_to_sphere", std::move(signed_distance)});
            }
            cartmesh::write_vtu(*options.output, grid, fields);
        }
        const auto vtk_end = Clock::now();
        const auto result_hash = hash_hex(hash);
        write_report(options, grid, selected_cells, exact_volume, estimated_volume,
                     is_sphere ? "cell_center_sample" : "analytic_domain_exact", result_hash,
                     seconds_between(grid_start, grid_end),
                     seconds_between(classify_start, classify_end),
                     seconds_between(vtk_start, vtk_end), seconds_between(total_start, Clock::now()));

        std::cout << std::setprecision(10) << "案例=" << options.case_name
                  << " 单元数=" << grid.cell_count() << " 选中单元数=" << selected_cells
                  << " 估算体积=" << estimated_volume << " 解析体积=" << exact_volume
                  << " 哈希=" << result_hash << " 峰值RSS字节=" << cartmesh::peak_rss_bytes()
                  << " 总秒数=" << seconds_between(total_start, Clock::now()) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "cartmesh_cli 错误：" << error.what() << '\n';
        return 1;
    }
}
