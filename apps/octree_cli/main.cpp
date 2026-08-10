#include "cartmesh/classify/SurfaceClassifier.hpp"
#include "cartmesh/geometry/SurfaceDiagnostics.hpp"
#include "cartmesh/grid/OctreeRefinement.hpp"
#include "cartmesh/io/DiagnosticVtkWriter.hpp"
#include "cartmesh/io/StlReader.hpp"
#include "cartmesh/io/VtkWriter.hpp"
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
    std::filesystem::path stl;
    std::uint8_t base_level{2};
    std::uint8_t maximum_level{6};
    std::optional<std::uint8_t> surface_level;
    double padding_fraction{0.05};
    cartmesh::OctreeRefinementConfiguration refinement;
    std::optional<std::filesystem::path> output;
    std::optional<std::filesystem::path> report;
    bool strict_gap_resolution{};
};

[[nodiscard]] double elapsed(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
}

[[nodiscard]] std::uint8_t parse_level(std::string_view text, std::string_view option) {
    std::size_t parsed = 0;
    const auto value = std::stoul(std::string(text), &parsed);
    if (parsed != text.size() || value > cartmesh::maximum_octree_level) {
        throw std::invalid_argument(std::string(option) + " 需要 0..21 的整数层级");
    }
    return static_cast<std::uint8_t>(value);
}

[[nodiscard]] std::uint32_t parse_positive_u32(std::string_view text,
                                               std::string_view option) {
    std::size_t parsed = 0;
    const auto value = std::stoull(std::string(text), &parsed);
    if (parsed != text.size() || value == 0 ||
        value > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument(std::string(option) + " 需要 32 位正整数");
    }
    return static_cast<std::uint32_t>(value);
}

[[nodiscard]] double parse_nonnegative(std::string_view text, std::string_view option) {
    std::size_t parsed = 0;
    const double value = std::stod(std::string(text), &parsed);
    if (parsed != text.size() || !std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(std::string(option) + " 需要非负有限数");
    }
    return value;
}

[[nodiscard]] std::vector<double> parse_csv(std::string_view text, std::size_t count,
                                            std::string_view option) {
    std::vector<double> result;
    std::size_t begin = 0;
    while (begin <= text.size()) {
        const auto comma = text.find(',', begin);
        const auto token = text.substr(begin, comma == std::string_view::npos
                                                  ? text.size() - begin
                                                  : comma - begin);
        std::size_t parsed = 0;
        const double value = std::stod(std::string(token), &parsed);
        if (parsed != token.size() || !std::isfinite(value)) {
            throw std::invalid_argument(std::string(option) + " 包含无效数值");
        }
        result.push_back(value);
        if (comma == std::string_view::npos) {
            break;
        }
        begin = comma + 1U;
    }
    if (result.size() != count) {
        throw std::invalid_argument(std::string(option) + " 数值个数不正确");
    }
    return result;
}

[[nodiscard]] std::pair<std::string_view, std::uint8_t>
split_level_suffix(std::string_view text, std::string_view option) {
    const auto colon = text.rfind(':');
    if (colon == std::string_view::npos || colon == 0 || colon + 1U == text.size()) {
        throw std::invalid_argument(std::string(option) + " 需要 values:level 格式");
    }
    return {text.substr(0, colon), parse_level(text.substr(colon + 1U), option)};
}

void print_help() {
    std::cout
        << "用法：cartmesh_octree_cli --stl FILE [选项]\n"
        << "  --base-level L                 基础层级（默认 2）\n"
        << "  --max-level L                  最大层级（默认 6）\n"
        << "  --surface-level L              表面相交细化层级（默认=max）\n"
        << "  --distance D:L                 表面距离 D 内细化至 L，可重复\n"
        << "  --curvature DEG:L              法向夹角至少 DEG 时细化至 L\n"
        << "  --gap D:N                      搜索距离 D 内的相对面，间隙至少 N 单元\n"
        << "  --box xmin,ymin,zmin,xmax,ymax,zmax:L\n"
        << "  --sphere cx,cy,cz,r:L\n"
        << "  --cylinder ax,ay,az,bx,by,bz,r:L\n"
        << "  --padding-fraction X           STL 包围盒留白（默认 0.05）\n"
        << "  --no-balance                   仅用于诊断；不执行面 2:1 平衡\n"
        << "  --strict-gaps                  gap 在最大层仍不足时以非零状态退出\n"
        << "  --output FILE.vtu | --no-vtk   可视化输出\n"
        << "  --report FILE.json             机器可读报告\n";
}

[[nodiscard]] Options parse_options(int argc, char** argv) {
    Options options;
    bool output_set = false;
    for (int argument = 1; argument < argc; ++argument) {
        const std::string_view current(argv[argument]);
        const auto value = [&](std::string_view option) -> std::string_view {
            if (argument + 1 >= argc) {
                throw std::invalid_argument(std::string(option) + " 缺少值");
            }
            return argv[++argument];
        };
        if (current == "--help") {
            print_help();
            std::exit(0);
        } else if (current == "--stl") {
            options.stl = value(current);
        } else if (current == "--base-level") {
            options.base_level = parse_level(value(current), current);
        } else if (current == "--max-level") {
            options.maximum_level = parse_level(value(current), current);
        } else if (current == "--surface-level") {
            options.surface_level = parse_level(value(current), current);
        } else if (current == "--padding-fraction") {
            options.padding_fraction = parse_nonnegative(value(current), current);
        } else if (current == "--distance") {
            const auto [distance, level] = split_level_suffix(value(current), current);
            options.refinement.distance_bands.push_back(
                {parse_nonnegative(distance, current), level});
        } else if (current == "--curvature") {
            const auto [degrees, level] = split_level_suffix(value(current), current);
            options.refinement.curvature = cartmesh::CurvatureRefinementRule{
                parse_nonnegative(degrees, current), 1.5, level};
        } else if (current == "--gap") {
            const auto text = value(current);
            const auto colon = text.rfind(':');
            if (colon == std::string_view::npos) {
                throw std::invalid_argument("--gap 需要 D:N 格式");
            }
            options.refinement.gap = cartmesh::GapRefinementRule{
                parse_nonnegative(text.substr(0, colon), current), -0.75, 0.5,
                parse_positive_u32(text.substr(colon + 1U), current)};
        } else if (current == "--box") {
            const auto [coordinates, level] = split_level_suffix(value(current), current);
            const auto v = parse_csv(coordinates, 6, current);
            options.refinement.boxes.push_back(
                {cartmesh::AABB({v[0], v[1], v[2]}, {v[3], v[4], v[5]}), level});
        } else if (current == "--sphere") {
            const auto [coordinates, level] = split_level_suffix(value(current), current);
            const auto v = parse_csv(coordinates, 4, current);
            options.refinement.spheres.push_back({{v[0], v[1], v[2]}, v[3], level});
        } else if (current == "--cylinder") {
            const auto [coordinates, level] = split_level_suffix(value(current), current);
            const auto v = parse_csv(coordinates, 7, current);
            options.refinement.cylinders.push_back(
                {{v[0], v[1], v[2]}, {v[3], v[4], v[5]}, v[6], level});
        } else if (current == "--no-balance") {
            options.refinement.enforce_face_2_to_1_balance = false;
        } else if (current == "--strict-gaps") {
            options.strict_gap_resolution = true;
        } else if (current == "--output") {
            options.output = std::filesystem::path(value(current));
            output_set = true;
        } else if (current == "--no-vtk") {
            options.output.reset();
            output_set = true;
        } else if (current == "--report") {
            options.report = std::filesystem::path(value(current));
        } else {
            throw std::invalid_argument("未知选项：" + std::string(current));
        }
    }
    if (options.stl.empty()) {
        throw std::invalid_argument("阶段 2 需要 --stl FILE");
    }
    if (options.base_level > options.maximum_level) {
        throw std::invalid_argument("基础层级不得超过最大层级");
    }
    options.refinement.surface_target_level =
        options.surface_level.value_or(options.maximum_level);
    if (!output_set) {
        options.output = std::filesystem::path("artifacts") /
                         (options.stl.stem().string() + "_octree.vtu");
    }
    if (!options.report) {
        options.report = std::filesystem::path("artifacts") /
                         (options.stl.stem().string() + "_octree.json");
    }
    return options;
}

[[nodiscard]] cartmesh::AABB padded_domain(const cartmesh::AABB& bounds,
                                           double fraction) {
    const auto extent = bounds.extent();
    const double reference = std::max({extent.x, extent.y, extent.z});
    if (!(reference > 0.0)) {
        throw std::invalid_argument("STL 必须具有正几何尺度");
    }
    const double amount = std::max(reference * fraction, reference * 1.0e-9);
    const cartmesh::Vec3 padding{amount, amount, amount};
    return cartmesh::AABB(bounds.minimum() - padding, bounds.maximum() + padding);
}

void hash_byte(std::uint64_t& hash, std::uint8_t value) noexcept {
    hash ^= value;
    hash *= 1099511628211ULL;
}

[[nodiscard]] std::string hex_hash(std::uint64_t hash) {
    std::ostringstream output;
    output << std::hex << std::setw(16) << std::setfill('0') << hash;
    return output.str();
}

[[nodiscard]] std::string json_escape(std::string_view value) {
    std::ostringstream output;
    for (const char character : value) {
        switch (character) {
        case '\\':
            output << "\\\\";
            break;
        case '"':
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

[[nodiscard]] std::string file_hash(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("无法打开 STL 计算哈希：" + path.string());
    }
    std::uint64_t hash = 14695981039346656037ULL;
    std::array<char, 65536> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        for (std::streamsize index = 0; index < input.gcount(); ++index) {
            hash_byte(hash, static_cast<std::uint8_t>(buffer[static_cast<std::size_t>(index)]));
        }
    }
    if (!input.eof()) {
        throw std::runtime_error("读取 STL 哈希失败");
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

void write_vec(std::ostream& output, const cartmesh::Vec3& value) {
    output << '[' << value.x << ',' << value.y << ',' << value.z << ']';
}

void ensure_parent(const std::filesystem::path& path) {
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }
}

int run(const Options& options) {
    const auto total_start = Clock::now();
    const auto input_hash = file_hash(options.stl);
    const auto import_start = Clock::now();
    const auto surface = cartmesh::read_stl(options.stl);
    const auto import_end = Clock::now();
    const auto diagnostic_start = Clock::now();
    const auto diagnostics = cartmesh::diagnose_surface(surface);
    const auto diagnostic_end = Clock::now();
    if (!diagnostics.valid_for_stage1_classification()) {
        auto marker = *options.report;
        marker.replace_extension();
        marker += "_geometry_issues.vtp";
        cartmesh::write_surface_diagnostic_vtp(marker, surface, diagnostics);
        ensure_parent(*options.report);
        std::ofstream report(*options.report, std::ios::trunc);
        report << "{\"schemaVersion\":4,\"projectStage\":2,"
                  "\"status\":\"rejected_invalid_surface\","
                  "\"solverReadyCutCellMesh\":false,\"diagnosticMarkerVtp\":\""
               << json_escape(marker.string()) << "\",\"boundaryEdgeCount\":"
               << diagnostics.boundary_edge_count << ",\"nonManifoldEdgeCount\":"
               << diagnostics.non_manifold_edge_count << ",\"nonManifoldVertexCount\":"
               << diagnostics.non_manifold_vertex_count
               << ",\"orientationConflictEdgeCount\":"
               << diagnostics.orientation_conflict_edge_count
               << ",\"componentOrientationMismatchCount\":"
               << diagnostics.component_orientation_mismatch_count
               << ",\"overlappingTrianglePairCount\":"
               << diagnostics.overlapping_triangle_pair_count
               << ",\"selfIntersectionPairCount\":"
               << diagnostics.self_intersection_pair_count
               << ",\"nonAdjacentContactPairCount\":"
               << diagnostics.non_adjacent_contact_pair_count << "}\n";
        throw std::runtime_error("STL 未通过封闭、定向、流形诊断");
    }
    const auto bvh_start = Clock::now();
    const cartmesh::TriangleBvh bvh(surface);
    const auto bvh_end = Clock::now();
    const auto adapt_start = Clock::now();
    cartmesh::LinearOctree tree(padded_domain(surface.bounds(), options.padding_fraction),
                                options.base_level, options.maximum_level);
    const cartmesh::OctreeRefinementEngine engine(options.refinement, &bvh);
    const auto adaptation = engine.apply(tree);
    const auto adapt_end = Clock::now();
    const auto classification_start = Clock::now();
    const cartmesh::SurfaceClassifier classifier(bvh);
    const auto classification = cartmesh::classify_octree_leaves(tree, classifier);
    const auto classification_end = Clock::now();
    const auto balance_report = tree.check_face_balance();
    const auto levels = tree.level_statistics();
    const double exact_polyhedral_volume = diagnostics.material_volume;
    const bool volume_bracket_pass =
        classification.inside_volume <= exact_polyhedral_volume &&
        exact_polyhedral_volume <= classification.inside_plus_intersected_volume;
    const bool gap_constraint_pass = adaptation.gap_resolution_failure_count == 0 ||
                                     !options.strict_gap_resolution;
    const bool invariants_pass =
        classification.conflict_count == 0 && tree.validate_partition() &&
        (!options.refinement.enforce_face_2_to_1_balance || balance_report.balanced) &&
        volume_bracket_pass && gap_constraint_pass;
    std::uint64_t result_hash = tree.result_hash_fnv1a64();
    for (std::size_t index = 0; index < classification.cell_classification.size(); ++index) {
        hash_byte(result_hash, classification.cell_classification[index]);
        hash_byte(result_hash, classification.center_point_classification[index]);
    }
    const auto vtk_start = Clock::now();
    if (options.output) {
        std::vector<double> level(tree.leaf_count());
        std::vector<double> code_low(tree.leaf_count());
        std::vector<double> code_high(tree.leaf_count());
        std::vector<double> cell_classification(tree.leaf_count());
        std::vector<double> center_inside(tree.leaf_count());
        for (std::uint64_t leaf_id = 0; leaf_id < tree.leaf_count(); ++leaf_id) {
            const auto code = tree.leaf_code(leaf_id);
            const auto index = static_cast<std::size_t>(leaf_id);
            level[index] = static_cast<double>(cartmesh::decode_octree_node(code).level);
            code_low[index] = static_cast<double>(code & 0xffffffffULL);
            code_high[index] = static_cast<double>(code >> 32U);
            cell_classification[index] =
                static_cast<double>(classification.cell_classification[index]);
            center_inside[index] =
                classification.center_point_classification[index] ==
                        static_cast<std::uint8_t>(cartmesh::PointClassification::inside)
                    ? 1.0
                    : 0.0;
        }
        cartmesh::write_octree_vtu(
            *options.output, tree,
            {{"octree_level", std::move(level)},
             {"octree_node_code_low32", std::move(code_low)},
             {"octree_node_code_high32", std::move(code_high)},
             {"stl_cell_classification", std::move(cell_classification)},
             {"inside_stl_center_sample", std::move(center_inside)}});
    }
    const auto vtk_end = Clock::now();
    ensure_parent(*options.report);
    std::ofstream report(*options.report, std::ios::trunc);
    if (!report) {
        throw std::runtime_error("无法打开阶段 2 报告：" + options.report->string());
    }
    const auto& bvh_statistics = bvh.statistics();
    report << std::setprecision(17) << "{\n"
           << "  \"schemaVersion\": 4,\n  \"projectStage\": 2,\n"
           << "  \"status\": \""
           << (invariants_pass
                   ? (adaptation.gap_resolution_failure_count == 0
                          ? "pass"
                          : "pass_with_gap_resolution_warning")
                   : (gap_constraint_pass ? "failed_invariant" : "failed_gap_resolution"))
           << "\",\n  \"meshKind\": \"adaptive_cartesian_background_hexahedra\",\n"
           << "  \"solverReadyCutCellMesh\": false,\n"
           << "  \"warnings\": ["
           << (adaptation.gap_resolution_failure_count == 0
                   ? ""
                   : "\"maximumLevelCannotSatisfyRequestedCellsAcrossSomeDetectedGaps\"")
           << "],\n"
           << "  \"strictGapResolution\": " << options.strict_gap_resolution << ",\n"
           << "  \"input\": \"" << json_escape(options.stl.string()) << "\",\n"
           << "  \"inputFormat\": \"" << cartmesh::surface_format_name(surface.format())
           << "\",\n  \"inputHashFnv1a64\": \"" << input_hash << "\",\n"
           << "  \"surfaceDiagnostics\": {\"triangleCount\": "
           << diagnostics.triangle_count << ", \"closed\": " << std::boolalpha
           << diagnostics.closed << ", \"manifold\": " << diagnostics.manifold
           << ", \"consistentlyOriented\": " << diagnostics.consistently_oriented
           << ", \"componentOrientationMismatchCount\": "
           << diagnostics.component_orientation_mismatch_count
           << ", \"signedVolume\": " << diagnostics.signed_volume
           << ", \"materialVolume\": " << diagnostics.material_volume << "},\n"
           << "  \"bvh\": {\"nodeCount\": " << bvh_statistics.node_count
           << ", \"leafCount\": " << bvh_statistics.leaf_count
           << ", \"maximumDepth\": " << bvh_statistics.maximum_depth << "},\n"
           << "  \"domain\": {\"minimum\": ";
    write_vec(report, tree.domain().minimum());
    report << ", \"maximum\": ";
    write_vec(report, tree.domain().maximum());
    report << "},\n  \"baseLevel\": " << static_cast<unsigned>(tree.base_level())
           << ",\n  \"maximumLevel\": " << static_cast<unsigned>(tree.maximum_level())
           << ",\n  \"leafCount\": " << tree.leaf_count()
           << ",\n  \"compactLeafStorageBytes\": " << tree.compact_storage_bytes() << ",\n"
           << "  \"partitionValid\": " << tree.validate_partition() << ",\n"
           << "  \"faceBalance\": {\"enforced\": "
           << options.refinement.enforce_face_2_to_1_balance << ", \"balanced\": "
           << balance_report.balanced << ", \"violatingFacePairCount\": "
           << balance_report.violating_face_pair_count << "},\n"
           << "  \"refinement\": {\"ruleSplitCount\": "
           << adaptation.rule_refinement.split_count << ", \"balanceSplitCount\": "
           << adaptation.balance.split_count << ", \"surfaceHits\": "
           << adaptation.surface_rule_hits << ", \"distanceHits\": "
           << adaptation.distance_rule_hits << ", \"curvatureHits\": "
           << adaptation.curvature_rule_hits << ", \"gapHits\": "
           << adaptation.gap_rule_hits << ", \"gapResolutionFailureCount\": "
           << adaptation.gap_resolution_failure_count
           << ", \"maximumRequiredGapLevel\": "
           << adaptation.maximum_required_gap_level << ", \"userRegionHits\": "
           << adaptation.user_region_rule_hits << "},\n"
           << "  \"refinementRules\": {\"surfaceTargetLevel\": "
           << static_cast<unsigned>(*options.refinement.surface_target_level)
           << ", \"distanceBands\": [";
    for (std::size_t index = 0; index < options.refinement.distance_bands.size(); ++index) {
        const auto& band = options.refinement.distance_bands[index];
        report << (index == 0 ? "" : ",") << "{\"maximumDistance\":"
               << band.maximum_distance << ",\"targetLevel\":"
               << static_cast<unsigned>(band.target_level) << '}';
    }
    report << "], \"curvature\": ";
    if (options.refinement.curvature) {
        report << "{\"minimumNormalAngleDegrees\":"
               << options.refinement.curvature->minimum_normal_angle_degrees
               << ",\"neighborhoodCellDiagonals\":"
               << options.refinement.curvature->neighborhood_cell_diagonals
               << ",\"targetLevel\":"
               << static_cast<unsigned>(options.refinement.curvature->target_level) << '}';
    } else {
        report << "null";
    }
    report << ", \"gap\": ";
    if (options.refinement.gap) {
        report << "{\"maximumSearchDistance\":"
               << options.refinement.gap->maximum_search_distance
               << ",\"minimumCellsAcrossGap\":"
               << options.refinement.gap->minimum_cells_across_gap << '}';
    } else {
        report << "null";
    }
    report << ", \"boxCount\":" << options.refinement.boxes.size()
           << ", \"sphereCount\":" << options.refinement.spheres.size()
           << ", \"cylinderCount\":" << options.refinement.cylinders.size()
           << ", \"boxes\":[";
    for (std::size_t index = 0; index < options.refinement.boxes.size(); ++index) {
        const auto& region = options.refinement.boxes[index];
        report << (index == 0 ? "" : ",") << "{\"minimum\":";
        write_vec(report, region.bounds.minimum());
        report << ",\"maximum\":";
        write_vec(report, region.bounds.maximum());
        report << ",\"targetLevel\":" << static_cast<unsigned>(region.target_level) << '}';
    }
    report << "], \"spheres\":[";
    for (std::size_t index = 0; index < options.refinement.spheres.size(); ++index) {
        const auto& region = options.refinement.spheres[index];
        report << (index == 0 ? "" : ",") << "{\"center\":";
        write_vec(report, region.center);
        report << ",\"radius\":" << region.radius << ",\"targetLevel\":"
               << static_cast<unsigned>(region.target_level) << '}';
    }
    report << "], \"cylinders\":[";
    for (std::size_t index = 0; index < options.refinement.cylinders.size(); ++index) {
        const auto& region = options.refinement.cylinders[index];
        report << (index == 0 ? "" : ",") << "{\"firstAxisPoint\":";
        write_vec(report, region.first_axis_point);
        report << ",\"secondAxisPoint\":";
        write_vec(report, region.second_axis_point);
        report << ",\"radius\":" << region.radius << ",\"targetLevel\":"
               << static_cast<unsigned>(region.target_level) << '}';
    }
    report << "]},\n  \"leafCountByLevel\": [";
    for (std::size_t index = 0; index < levels.leaf_count_by_level.size(); ++index) {
        report << (index == 0 ? "" : ",") << levels.leaf_count_by_level[index];
    }
    report << "],\n  \"classificationCounts\": {\"outside\": "
           << classification.outside_count << ", \"inside\": " << classification.inside_count
           << ", \"intersected\": " << classification.intersected_count
           << ", \"conflict\": " << classification.conflict_count << "},\n"
           << "  \"classificationValueLegend\": {\"outside\":0,\"inside\":1,"
              "\"intersected\":2,\"conflict\":3},\n"
           << "  \"centerPointCounts\": {\"outside\":"
           << classification.center_outside_count << ",\"inside\":"
           << classification.center_inside_count << ",\"onSurface\":"
           << classification.center_on_surface_count << ",\"conflict\":"
           << classification.center_conflict_count << "},\n"
           << "  \"definitelyInsideVolumeLowerBound\": " << classification.inside_volume
           << ",\n  \"insidePlusIntersectedVolumeUpperBound\": "
           << classification.inside_plus_intersected_volume
           << ",\n  \"exactPolyhedralVolume\": " << exact_polyhedral_volume
           << ",\n  \"polyhedralVolumeInsideClassificationBounds\": "
           << volume_bracket_pass
           << ",\n  \"resultHashFnv1a64\": \"" << hex_hash(result_hash) << "\",\n"
           << "  \"threads\": 1,\n  \"buildType\": \"" << build_type()
           << "\",\n  \"compiler\": \"" << json_escape(compiler()) << "\",\n"
           << "  \"peakRssBytes\": " << cartmesh::peak_rss_bytes() << ",\n"
           << "  \"timingsSeconds\": {\"stlImport\": "
           << elapsed(import_start, import_end) << ", \"surfaceDiagnostics\": "
           << elapsed(diagnostic_start, diagnostic_end) << ", \"bvhBuild\": "
           << elapsed(bvh_start, bvh_end) << ", \"adaptation\": "
           << elapsed(adapt_start, adapt_end) << ", \"classification\": "
           << elapsed(classification_start, classification_end) << ", \"vtkWrite\": "
           << elapsed(vtk_start, vtk_end) << ", \"total\": "
           << elapsed(total_start, Clock::now()) << "},\n  \"vtkOutput\": ";
    if (options.output) {
        report << "\"" << json_escape(options.output->string()) << "\"\n";
    } else {
        report << "null\n";
    }
    report << "}\n";
    if (!report) {
        throw std::runtime_error("写入阶段 2 报告失败");
    }
    std::cout << "阶段2 STL=" << options.stl.string() << " 叶数=" << tree.leaf_count()
              << " 内部=" << classification.inside_count
              << " 外部=" << classification.outside_count
              << " 相交=" << classification.intersected_count
              << " 冲突=" << classification.conflict_count
              << " 2:1=" << balance_report.balanced << " 哈希=" << hex_hash(result_hash)
              << " 总秒数=" << elapsed(total_start, Clock::now()) << '\n';
    return invariants_pass ? 0 : 2;
}

} // 匿名命名空间

int main(int argc, char** argv) {
    try {
        return run(parse_options(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "cartmesh_octree_cli 错误：" << error.what() << '\n';
        return 1;
    }
}
