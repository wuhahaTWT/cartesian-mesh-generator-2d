#include "cartmesh/cutcell/CutCellMeshFingerprint.hpp"
#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"
#include "cartmesh/geometry/SurfaceDiagnostics.hpp"
#include "cartmesh/grid/OctreeRefinement.hpp"
#include "cartmesh/incremental/IncrementalRemesher.hpp"
#include "cartmesh/incremental/MeshMapping.hpp"
#include "cartmesh/incremental/SurfaceChangeSet.hpp"
#include "cartmesh/io/CutCellJsonWriter.hpp"
#include "cartmesh/io/CutCellVtkWriter.hpp"
#include "cartmesh/io/StlReader.hpp"
#include "cartmesh/io/VtkWriter.hpp"

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
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#if defined(__APPLE__)
#include <sys/resource.h>
#include <sys/sysctl.h>
#elif defined(__unix__)
#include <sys/resource.h>
#endif

namespace {

struct Options {
    std::filesystem::path old_stl;
    std::filesystem::path new_stl;
    std::filesystem::path old_output{"artifacts/stage5_old.vtu"};
    std::filesystem::path new_output{"artifacts/stage5_incremental.vtu"};
    std::filesystem::path full_output{"artifacts/stage5_full_rebuild.vtu"};
    std::filesystem::path boundary_output{"artifacts/stage5_boundary.vtp"};
    std::filesystem::path geometry_output{"artifacts/stage5_geometry.json"};
    std::filesystem::path mapping_output{"artifacts/stage5_mapping.json"};
    std::filesystem::path report{"artifacts/stage5_report.json"};
    std::uint8_t base_level{2};
    std::uint8_t maximum_level{5};
    std::optional<std::uint8_t> surface_level;
    double padding_fraction{0.2};
    std::optional<double> change_margin;
    double geometric_tolerance{};
    cartmesh::OctreeRefinementConfiguration refinement;
    bool write_vtk{true};
};

[[nodiscard]] double parse_nonnegative(std::string_view text,
                                       std::string_view option) {
    std::size_t used = 0;
    const double value = std::stod(std::string(text), &used);
    if (used != text.size() || !std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(std::string(option) +
                                    " 需要非负有限数");
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

[[nodiscard]] std::pair<double, std::uint8_t> parse_distance(
    std::string_view text) {
    const auto separator = text.rfind(':');
    if (separator == std::string_view::npos) {
        throw std::invalid_argument("--distance 需要 D:L 格式");
    }
    return {parse_nonnegative(text.substr(0, separator), "--distance"),
            parse_level(text.substr(separator + 1U), "--distance")};
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
        if (value == "--old-stl") options.old_stl = next();
        else if (value == "--new-stl") options.new_stl = next();
        else if (value == "--base-level") {
            options.base_level = parse_level(next(), value);
        } else if (value == "--max-level") {
            options.maximum_level = parse_level(next(), value);
        } else if (value == "--surface-level") {
            options.surface_level = parse_level(next(), value);
        } else if (value == "--distance") {
            const auto [distance, level] = parse_distance(next());
            options.refinement.distance_bands.push_back({distance, level});
        } else if (value == "--padding-fraction") {
            options.padding_fraction = parse_nonnegative(next(), value);
        } else if (value == "--change-margin") {
            options.change_margin = parse_nonnegative(next(), value);
        } else if (value == "--geometric-tolerance") {
            options.geometric_tolerance = parse_nonnegative(next(), value);
        } else if (value == "--old-output") options.old_output = next();
        else if (value == "--new-output") options.new_output = next();
        else if (value == "--full-output") options.full_output = next();
        else if (value == "--boundary-output") options.boundary_output = next();
        else if (value == "--geometry-output") options.geometry_output = next();
        else if (value == "--mapping-output") options.mapping_output = next();
        else if (value == "--report") options.report = next();
        else if (value == "--no-vtk") options.write_vtk = false;
        else if (value == "--help") {
            std::cout
                << "用法：cartmesh_incremental_cli --old-stl OLD --new-stl NEW [选项]\n"
                << "  --base-level L --max-level L --surface-level L\n"
                << "  --distance D:L --padding-fraction F\n"
                << "  --change-margin D --geometric-tolerance D\n"
                << "  --old-output FILE.vtu --new-output FILE.vtu --full-output FILE.vtu\n"
                << "  --boundary-output FILE.vtp --geometry-output FILE.json\n"
                << "  --mapping-output FILE.json --report FILE.json --no-vtk\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("未知参数：" + value);
        }
    }
    if (options.old_stl.empty() || options.new_stl.empty()) {
        throw std::invalid_argument("必须同时提供 --old-stl 和 --new-stl");
    }
    if (options.base_level > options.maximum_level) {
        throw std::invalid_argument("基础层级不得超过最大层级");
    }
    options.refinement.surface_target_level =
        options.surface_level.value_or(options.maximum_level);
    return options;
}

void ensure_parent(const std::filesystem::path& path) {
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
}

[[nodiscard]] cartmesh::AABB union_bounds(const cartmesh::AABB& first,
                                          const cartmesh::AABB& second) {
    return cartmesh::AABB(
        {std::min(first.minimum().x, second.minimum().x),
         std::min(first.minimum().y, second.minimum().y),
         std::min(first.minimum().z, second.minimum().z)},
        {std::max(first.maximum().x, second.maximum().x),
         std::max(first.maximum().y, second.maximum().y),
         std::max(first.maximum().z, second.maximum().z)});
}

[[nodiscard]] cartmesh::AABB padded_domain(const cartmesh::SurfaceMesh& old_surface,
                                           const cartmesh::SurfaceMesh& new_surface,
                                           double padding_fraction) {
    const auto combined = union_bounds(old_surface.bounds(), new_surface.bounds());
    const auto extent = combined.extent();
    const double scale = std::max({extent.x, extent.y, extent.z});
    if (!(scale > 0.0)) {
        throw std::invalid_argument("旧、新几何联合包围盒必须具有正尺度");
    }
    const double padding = padding_fraction * scale;
    const cartmesh::Vec3 vector{padding, padding, padding};
    return cartmesh::AABB(combined.minimum() - vector,
                          combined.maximum() + vector);
}

void hash_byte(std::uint64_t& hash, std::uint8_t value) noexcept {
    hash ^= value;
    hash *= 1099511628211ULL;
}

[[nodiscard]] std::string file_hash(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("无法读取输入哈希：" + path.string());
    std::uint64_t hash = 14695981039346656037ULL;
    std::array<char, 65536> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        for (std::streamsize index = 0; index < input.gcount(); ++index) {
            hash_byte(hash, static_cast<std::uint8_t>(
                                buffer[static_cast<std::size_t>(index)]));
        }
    }
    return cartmesh::fingerprint_hex(hash);
}

[[nodiscard]] std::uint64_t peak_rss_bytes() noexcept {
#if defined(__APPLE__) || defined(__unix__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#if defined(__APPLE__)
    return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
#else
    return 0;
#endif
}

[[nodiscard]] std::string sysctl_string(const char* name) {
#if defined(__APPLE__)
    std::size_t size = 0;
    if (sysctlbyname(name, nullptr, &size, nullptr, 0) != 0 || size == 0) {
        return "unknown";
    }
    std::string value(size, '\0');
    if (sysctlbyname(name, value.data(), &size, nullptr, 0) != 0) {
        return "unknown";
    }
    while (!value.empty() && value.back() == '\0') value.pop_back();
    return value;
#else
    static_cast<void>(name);
    return "unknown";
#endif
}

[[nodiscard]] std::uint64_t sysctl_u64(const char* name) noexcept {
#if defined(__APPLE__)
    std::uint64_t value = 0;
    std::size_t size = sizeof(value);
    if (sysctlbyname(name, &value, &size, nullptr, 0) != 0) return 0;
    return value;
#else
    static_cast<void>(name);
    return 0;
#endif
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

[[nodiscard]] std::vector<cartmesh::VtkCellData> vtk_fields(
    const cartmesh::LinearOctree& tree,
    const cartmesh::ConvexCutCellMesh& mesh,
    const std::vector<std::uint8_t>* rebuilt_mask) {
    std::vector<double> fluid_fraction(static_cast<std::size_t>(tree.leaf_count()),
                                       0.0);
    std::vector<double> cut_cell(static_cast<std::size_t>(tree.leaf_count()), 0.0);
    for (const auto& cell : mesh.fluid_cells) {
        const auto index = static_cast<std::size_t>(cell.background_cell_id);
        fluid_fraction[index] = cell.volume_fraction;
        cut_cell[index] = cell.cut ? 1.0 : 0.0;
    }
    std::vector<double> level(static_cast<std::size_t>(tree.leaf_count()));
    std::vector<double> code_low(static_cast<std::size_t>(tree.leaf_count()));
    std::vector<double> code_high(static_cast<std::size_t>(tree.leaf_count()));
    std::vector<double> rebuilt(static_cast<std::size_t>(tree.leaf_count()), 0.0);
    std::vector<double> reused(static_cast<std::size_t>(tree.leaf_count()), 0.0);
    for (std::uint64_t leaf = 0; leaf < tree.leaf_count(); ++leaf) {
        const auto code = tree.leaf_code(leaf);
        const auto index = static_cast<std::size_t>(leaf);
        level[index] = static_cast<double>(cartmesh::decode_octree_node(code).level);
        code_low[index] = static_cast<double>(code & 0xffffffffULL);
        code_high[index] = static_cast<double>(code >> 32U);
        if (rebuilt_mask) {
            rebuilt[index] = (*rebuilt_mask)[index] != 0 ? 1.0 : 0.0;
            reused[index] = (*rebuilt_mask)[index] == 0 ? 1.0 : 0.0;
        }
    }
    return {{"fluid_volume_fraction", std::move(fluid_fraction)},
            {"cut_cell", std::move(cut_cell)},
            {"octree_level", std::move(level)},
            {"octree_node_code_low32", std::move(code_low)},
            {"octree_node_code_high32", std::move(code_high)},
            {"stage5_rebuilt", std::move(rebuilt)},
            {"stage5_reused", std::move(reused)}};
}

void write_vec3(std::ostream& output, const cartmesh::Vec3& value) {
    output << '[' << value.x << ',' << value.y << ',' << value.z << ']';
}

void write_mapping(const std::filesystem::path& path,
                   const cartmesh::IncrementalMeshMapping& mapping) {
    ensure_parent(path);
    std::ofstream output(path, std::ios::trunc);
    if (!output) throw std::runtime_error("无法写入增量映射：" + path.string());
    output << std::setprecision(17)
           << "{\n  \"schema\":\"cartmesh-stage5-mapping-v1\",\n"
           << "  \"oldFluidVolume\":" << mapping.old_fluid_volume << ",\n"
           << "  \"newFluidVolume\":" << mapping.new_fluid_volume << ",\n"
           << "  \"sharedFluidOverlapVolume\":"
           << mapping.total_shared_fluid_overlap_volume << ",\n"
           << "  \"removedFluidVolume\":" << mapping.removed_fluid_volume << ",\n"
           << "  \"createdFluidVolume\":" << mapping.created_fluid_volume << ",\n"
           << "  \"entries\":[\n";
    for (std::size_t index = 0; index < mapping.entries.size(); ++index) {
        const auto& entry = mapping.entries[index];
        output << "    {\"oldCode\":\"" << entry.old_code
               << "\",\"newCode\":\"" << entry.new_code
               << "\",\"relation\":\""
               << cartmesh::incremental_topology_relation_name(entry.relation)
               << "\",\"backgroundOverlapVolume\":"
               << entry.background_overlap_volume
               << ",\"oldFluidVolume\":" << entry.old_fluid_volume
               << ",\"newFluidVolume\":" << entry.new_fluid_volume
               << ",\"sharedFluidOverlapVolume\":"
               << entry.shared_fluid_overlap_volume
               << ",\"oldVolumePreservedFraction\":"
               << entry.old_volume_preserved_fraction
               << ",\"newVolumeFromOldFraction\":"
               << entry.new_volume_from_old_fraction
               << ",\"exactFluidOverlap\":true}";
        output << (index + 1U == mapping.entries.size() ? "\n" : ",\n");
    }
    output << "  ]\n}\n";
    if (!output) throw std::runtime_error("写入增量映射失败：" + path.string());
}

[[nodiscard]] bool mesh_invariants_pass(const cartmesh::ConvexCutCellMesh& mesh) {
    return mesh.nonclosed_cell_count == 0 &&
           mesh.negative_volume_cell_count == 0 &&
           mesh.shared_face_mismatch_count == 0 &&
           mesh.classification_conflict_count == 0 &&
           mesh.component_analysis_pending_cell_count == 0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        const auto old_surface = cartmesh::read_stl(options.old_stl);
        const auto new_surface = cartmesh::read_stl(options.new_stl);
        const auto old_diagnostics = cartmesh::diagnose_surface(old_surface);
        const auto new_diagnostics = cartmesh::diagnose_surface(new_surface);
        if (!old_diagnostics.valid_for_stage1_classification() ||
            !new_diagnostics.valid_for_stage1_classification()) {
            throw std::invalid_argument(
                "旧、新 STL 都必须通过封闭、流形、定向和相交诊断");
        }
        const auto domain = padded_domain(old_surface, new_surface,
                                          options.padding_fraction);
        const auto old_cutter = cartmesh::TriangulatedSurfaceCutter(old_surface, 0);
        const auto new_cutter = cartmesh::TriangulatedSurfaceCutter(new_surface, 0);
        const auto changes = cartmesh::detect_surface_changes(old_surface, new_surface);
        double automatic_margin = std::ldexp(
            cartmesh::norm(domain.extent()),
            -static_cast<int>(options.maximum_level));
        for (const auto& band : options.refinement.distance_bands) {
            automatic_margin = std::max(automatic_margin,
                                        band.maximum_distance);
        }
        const double change_margin =
            options.change_margin.value_or(automatic_margin);
        const auto affected = cartmesh::conservative_affected_bounds(
            changes, domain, change_margin);

        const auto old_start = std::chrono::steady_clock::now();
        cartmesh::LinearOctree old_tree(domain, options.base_level,
                                        options.maximum_level);
        const auto old_adaptation = cartmesh::OctreeRefinementEngine(
            options.refinement, &old_cutter.bvh()).apply(old_tree);
        const auto old_mesh = cartmesh::build_triangulated_cut_cell_mesh(
            old_tree, old_cutter, options.geometric_tolerance);
        const auto old_end = std::chrono::steady_clock::now();

        const cartmesh::OctreeRefinementEngine new_engine(
            options.refinement, &new_cutter.bvh());
        const auto incremental_start = std::chrono::steady_clock::now();
        auto incremental_tree = cartmesh::update_octree_incrementally(
            old_tree, new_engine, affected);
        const auto tree_end = std::chrono::steady_clock::now();
        auto incremental_mesh =
            cartmesh::build_incremental_triangulated_cut_cell_mesh(
                old_tree, old_mesh, incremental_tree.tree, new_cutter,
                affected, options.geometric_tolerance);
        const auto cutcell_end = std::chrono::steady_clock::now();
        const auto mapping = cartmesh::build_incremental_mesh_mapping(
            old_tree, old_mesh, incremental_tree.tree,
            incremental_mesh.mesh, affected, options.geometric_tolerance);
        const auto mapping_end = std::chrono::steady_clock::now();

        const auto full_start = std::chrono::steady_clock::now();
        cartmesh::LinearOctree full_tree(domain, options.base_level,
                                         options.maximum_level);
        const auto full_adaptation = new_engine.apply(full_tree);
        const auto full_mesh = cartmesh::build_triangulated_cut_cell_mesh(
            full_tree, new_cutter, options.geometric_tolerance);
        const auto full_end = std::chrono::steady_clock::now();

        const bool topology_equal =
            incremental_tree.tree.leaf_codes().size() ==
                full_tree.leaf_codes().size() &&
            std::equal(incremental_tree.tree.leaf_codes().begin(),
                       incremental_tree.tree.leaf_codes().end(),
                       full_tree.leaf_codes().begin());
        const auto incremental_fingerprint =
            cartmesh::cut_cell_mesh_fingerprint_fnv1a64(
                incremental_tree.tree, incremental_mesh.mesh);
        const auto full_fingerprint =
            cartmesh::cut_cell_mesh_fingerprint_fnv1a64(full_tree, full_mesh);
        const bool geometry_equal = incremental_fingerprint == full_fingerprint;
        const double volume_tolerance =
            1.0e-11 * std::max(1.0, domain.volume());
        const bool mapping_complete =
            std::abs(mapping.total_background_overlap_volume - domain.volume()) <=
            volume_tolerance;
        const bool internal_pass = topology_equal && geometry_equal &&
            mapping_complete && incremental_tree.tree.validate_partition() &&
            incremental_tree.tree.check_face_balance().balanced &&
            mesh_invariants_pass(incremental_mesh.mesh);

        const auto export_start = std::chrono::steady_clock::now();
        if (options.write_vtk) {
            ensure_parent(options.old_output);
            ensure_parent(options.new_output);
            ensure_parent(options.full_output);
            cartmesh::write_octree_vtu(options.old_output, old_tree,
                                       vtk_fields(old_tree, old_mesh, nullptr));
            cartmesh::write_octree_vtu(
                options.new_output, incremental_tree.tree,
                vtk_fields(incremental_tree.tree, incremental_mesh.mesh,
                           &incremental_mesh.rebuilt_leaf_mask));
            cartmesh::write_octree_vtu(options.full_output, full_tree,
                                       vtk_fields(full_tree, full_mesh, nullptr));
            ensure_parent(options.boundary_output);
            cartmesh::write_embedded_boundary_vtp(options.boundary_output,
                                                   incremental_mesh.mesh);
        }
        ensure_parent(options.geometry_output);
        cartmesh::write_cut_cell_geometry_json(
            options.geometry_output, incremental_tree.tree,
            incremental_mesh.mesh, false);
        write_mapping(options.mapping_output, mapping);
        const auto export_end = std::chrono::steady_clock::now();

        const auto seconds = [](auto begin, auto end) {
            return std::chrono::duration<double>(end - begin).count();
        };
        const double old_seconds = seconds(old_start, old_end);
        const double incremental_tree_seconds =
            seconds(incremental_start, tree_end);
        const double incremental_cutcell_seconds = seconds(tree_end, cutcell_end);
        const double mapping_seconds = seconds(cutcell_end, mapping_end);
        const double incremental_seconds =
            seconds(incremental_start, mapping_end);
        const double full_seconds = seconds(full_start, full_end);
        const double export_seconds = seconds(export_start, export_end);

        ensure_parent(options.report);
        std::ofstream report(options.report, std::ios::trunc);
        if (!report) throw std::runtime_error("无法写入阶段5报告");
        report << std::setprecision(17)
               << "{\n"
               << "  \"schema\":\"cartmesh-stage5-incremental-v1\",\n"
               << "  \"status\":\""
               << (internal_pass ? "internal_pass_external_pending" : "fail")
               << "\",\n"
               << "  \"stage5Complete\":false,\n"
               << "  \"externalIndependentReaderAccepted\":false,\n"
               << "  \"acceptanceBlockers\":[\"independent_external_reader\","
                  "\"three_case_matrix\",\"isolated_performance_runs\"],\n"
               << "  \"oldInput\":\"" << options.old_stl.string() << "\",\n"
               << "  \"newInput\":\"" << options.new_stl.string() << "\",\n"
               << "  \"oldInputHashFnv1a64\":\"" << file_hash(options.old_stl)
               << "\",\n"
               << "  \"newInputHashFnv1a64\":\"" << file_hash(options.new_stl)
               << "\",\n"
               << "  \"changedOldTriangleCount\":"
               << changes.old_triangle_ids.size() << ",\n"
               << "  \"changedNewTriangleCount\":"
               << changes.new_triangle_ids.size() << ",\n"
               << "  \"changeMargin\":" << change_margin << ",\n"
               << "  \"affectedBounds\":";
        if (affected) {
            report << "{\"minimum\":";
            write_vec3(report, affected->minimum());
            report << ",\"maximum\":";
            write_vec3(report, affected->maximum());
            report << "}";
        } else {
            report << "null";
        }
        report << ",\n"
               << "  \"oldLeafCount\":" << old_tree.leaf_count() << ",\n"
               << "  \"newLeafCount\":" << incremental_tree.tree.leaf_count()
               << ",\n"
               << "  \"preservedLeafCount\":"
               << incremental_tree.statistics.preserved_leaf_count << ",\n"
               << "  \"createdLeafCount\":"
               << incremental_tree.statistics.created_leaf_count << ",\n"
               << "  \"removedLeafCount\":"
               << incremental_tree.statistics.removed_leaf_count << ",\n"
               << "  \"coarsenedParentCount\":"
               << incremental_tree.statistics.coarsened_parent_count << ",\n"
               << "  \"ruleSplitCount\":"
               << incremental_tree.statistics.rule_split_count << ",\n"
               << "  \"balanceSplitCount\":"
               << incremental_tree.statistics.balance_split_count << ",\n"
               << "  \"geometryReusedLeafCount\":"
               << incremental_mesh.statistics.reused_leaf_count << ",\n"
               << "  \"geometryRebuiltLeafCount\":"
               << incremental_mesh.statistics.rebuilt_leaf_count << ",\n"
               << "  \"geometryReuseFraction\":"
               << incremental_mesh.statistics.geometry_reuse_fraction << ",\n"
               << "  \"topologyEqualsFullRebuild\":"
               << (topology_equal ? "true" : "false") << ",\n"
               << "  \"geometryEqualsFullRebuild\":"
               << (geometry_equal ? "true" : "false") << ",\n"
               << "  \"incrementalResultHashFnv1a64\":\""
               << cartmesh::fingerprint_hex(incremental_fingerprint) << "\",\n"
               << "  \"fullRebuildResultHashFnv1a64\":\""
               << cartmesh::fingerprint_hex(full_fingerprint) << "\",\n"
               << "  \"partitionValid\":"
               << (incremental_tree.tree.validate_partition() ? "true" : "false")
               << ",\n"
               << "  \"faceBalanceValid\":"
               << (incremental_tree.tree.check_face_balance().balanced ? "true"
                                                                       : "false")
               << ",\n"
               << "  \"cutCellCount\":" << incremental_mesh.mesh.cut_cell_count
               << ",\n"
               << "  \"nonclosedCellCount\":"
               << incremental_mesh.mesh.nonclosed_cell_count << ",\n"
               << "  \"negativeVolumeCellCount\":"
               << incremental_mesh.mesh.negative_volume_cell_count << ",\n"
               << "  \"sharedFaceMismatchCount\":"
               << incremental_mesh.mesh.shared_face_mismatch_count << ",\n"
               << "  \"classificationConflictCount\":"
               << incremental_mesh.mesh.classification_conflict_count << ",\n"
               << "  \"mappingEntryCount\":" << mapping.entries.size() << ",\n"
               << "  \"mappingBackgroundVolume\":"
               << mapping.total_background_overlap_volume << ",\n"
               << "  \"sharedFluidOverlapVolume\":"
               << mapping.total_shared_fluid_overlap_volume << ",\n"
               << "  \"createdFluidVolume\":" << mapping.created_fluid_volume
               << ",\n"
               << "  \"removedFluidVolume\":" << mapping.removed_fluid_volume
               << ",\n"
               << "  \"oldBaselineSeconds\":" << old_seconds << ",\n"
               << "  \"incrementalTreeSeconds\":"
               << incremental_tree_seconds << ",\n"
               << "  \"incrementalCutCellSeconds\":"
               << incremental_cutcell_seconds << ",\n"
               << "  \"mappingSeconds\":" << mapping_seconds << ",\n"
               << "  \"incrementalTotalSeconds\":" << incremental_seconds
               << ",\n"
               << "  \"fullRebuildSeconds\":" << full_seconds << ",\n"
               << "  \"measuredSpeedup\":"
               << (incremental_seconds > 0.0 ? full_seconds / incremental_seconds
                                             : 0.0)
               << ",\n"
               << "  \"exportSeconds\":" << export_seconds << ",\n"
               << "  \"peakRssBytes\":" << peak_rss_bytes() << ",\n"
               << "  \"peakRssScope\":\"old_baseline_plus_incremental_plus_full_validation_plus_export\",\n"
               << "  \"runtimeThreads\":1,\n"
               << "  \"hardwareModel\":\"" << sysctl_string("hw.model")
               << "\",\n"
               << "  \"hardwareMemoryBytes\":" << sysctl_u64("hw.memsize")
               << ",\n"
               << "  \"hardwareLogicalCpu\":" << sysctl_u64("hw.logicalcpu")
               << ",\n"
               << "  \"buildType\":\"" << build_type() << "\",\n"
               << "  \"compiler\":\"" << compiler() << "\",\n"
               << "  \"oldRuleSplitCount\":"
               << old_adaptation.rule_refinement.split_count << ",\n"
               << "  \"fullRuleSplitCount\":"
               << full_adaptation.rule_refinement.split_count << ",\n"
               << "  \"mappingOutput\":\"" << options.mapping_output.string()
               << "\",\n"
               << "  \"geometryOutput\":\"" << options.geometry_output.string()
               << "\",\n"
               << "  \"incrementalVtkOutput\":\""
               << options.new_output.string() << "\",\n"
               << "  \"fullRebuildVtkOutput\":\""
               << options.full_output.string() << "\"\n"
               << "}\n";
        if (!report) throw std::runtime_error("写入阶段5报告失败");

        std::cout << "阶段5增量重构：oldLeaves=" << old_tree.leaf_count()
                  << " newLeaves=" << incremental_tree.tree.leaf_count()
                  << " reused=" << incremental_mesh.statistics.reused_leaf_count
                  << " rebuilt=" << incremental_mesh.statistics.rebuilt_leaf_count
                  << " fullHash=" << cartmesh::fingerprint_hex(full_fingerprint)
                  << " incrementalHash="
                  << cartmesh::fingerprint_hex(incremental_fingerprint)
                  << " speedup="
                  << (incremental_seconds > 0.0 ? full_seconds / incremental_seconds
                                                : 0.0)
                  << '\n';
        return internal_pass ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "错误：" << error.what() << '\n';
        return 1;
    }
}
