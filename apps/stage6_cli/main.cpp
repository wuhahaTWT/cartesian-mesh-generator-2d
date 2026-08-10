#include "cartmesh/cutcell/CutCellMeshFingerprint.hpp"
#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"
#include "cartmesh/geometry/SurfaceDiagnostics.hpp"
#include "cartmesh/io/StlReader.hpp"
#include "cartmesh/io/ScalableOpenFoamWriter.hpp"
#include "cartmesh/io/VtkWriter.hpp"
#include "cartmesh/scalable/CompactUniformCutCellMesh.hpp"
#include "cartmesh/util/Performance.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif

namespace {

struct Options {
    std::filesystem::path stl;
    std::filesystem::path report{"artifacts/stage6_compact.json"};
    std::optional<std::filesystem::path> openfoam_case;
    std::optional<std::filesystem::path> preview;
    std::uint32_t resolution{48};
    double padding_fraction{0.1};
    double geometric_tolerance{};
    double small_cell_threshold{0.01};
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

[[nodiscard]] std::uint32_t parse_resolution(std::string_view text) {
    std::size_t used = 0;
    const auto value = std::stoull(std::string(text), &used);
    if (used != text.size() || value == 0 ||
        value > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument("--resolution 需要正的 32 位整数");
    }
    return static_cast<std::uint32_t>(value);
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
        if (value == "--stl") options.stl = next();
        else if (value == "--resolution") {
            options.resolution = parse_resolution(next());
        } else if (value == "--padding-fraction") {
            options.padding_fraction = parse_nonnegative(next(), value);
        } else if (value == "--geometric-tolerance") {
            options.geometric_tolerance = parse_nonnegative(next(), value);
        } else if (value == "--small-cell-threshold") {
            options.small_cell_threshold = parse_nonnegative(next(), value);
            if (options.small_cell_threshold > 1.0) {
                throw std::invalid_argument("--small-cell-threshold 不得超过 1");
            }
        } else if (value == "--report") options.report = next();
        else if (value == "--openfoam-case") options.openfoam_case = next();
        else if (value == "--preview") options.preview = next();
        else if (value == "--help") {
            std::cout
                << "用法：cartmesh_stage6_cli --stl FILE [选项]\n"
                << "  --resolution N              三向均匀背景分辨率\n"
                << "  --padding-fraction F        域外扩比例（默认 0.1）\n"
                << "  --geometric-tolerance F     显式几何长度容差\n"
                << "  --small-cell-threshold F    小 Cut-cell 阈值\n"
                << "  --openfoam-case DIR        写完整二进制 OpenFOAM case\n"
                << "  --preview FILE.vtu         写局部背景检查预览（非完整网格）\n"
                << "  --report FILE.json          紧凑几何/拓扑报告\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("未知参数：" + value);
        }
    }
    if (options.stl.empty()) throw std::invalid_argument("必须提供 --stl");
    return options;
}

void hash_byte(std::uint64_t& hash, std::uint8_t value) noexcept {
    hash ^= value;
    hash *= 1099511628211ULL;
}

[[nodiscard]] std::string file_hash(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("无法读取输入 hash");
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

[[nodiscard]] std::string hardware_model() {
#if defined(__APPLE__)
    std::size_t size = 0;
    if (sysctlbyname("hw.model", nullptr, &size, nullptr, 0) != 0 ||
        size == 0) {
        return "unknown";
    }
    std::string value(size, '\0');
    if (sysctlbyname("hw.model", value.data(), &size, nullptr, 0) != 0) {
        return "unknown";
    }
    while (!value.empty() && value.back() == '\0') value.pop_back();
    return value;
#else
    return "unknown";
#endif
}

[[nodiscard]] std::uint64_t hardware_memory() noexcept {
#if defined(__APPLE__)
    std::uint64_t value = 0;
    std::size_t size = sizeof(value);
    return sysctlbyname("hw.memsize", &value, &size, nullptr, 0) == 0
               ? value
               : 0;
#else
    return 0;
#endif
}

[[nodiscard]] std::uint32_t logical_cpu() noexcept {
#if defined(__APPLE__)
    std::uint32_t value = 0;
    std::size_t size = sizeof(value);
    return sysctlbyname("hw.logicalcpu", &value, &size, nullptr, 0) == 0
               ? value
               : 0;
#else
    return 0;
#endif
}

void write_vec3(std::ostream& output, const cartmesh::Vec3& value) {
    output << '[' << value.x << ',' << value.y << ',' << value.z << ']';
}

void write_u64_array(std::ostream& output,
                     const std::vector<std::uint64_t>& values) {
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) output << ',';
        output << values[index];
    }
    output << ']';
}

void write_face_failures(
    std::ostream& output,
    const std::vector<cartmesh::CompactFaceFailure>& failures) {
    output << '[';
    for (std::size_t index = 0; index < failures.size(); ++index) {
        if (index != 0) output << ',';
        const auto& failure = failures[index];
        output << "{\"firstBackgroundCellId\":"
               << failure.first_background_cell_id
               << ",\"secondBackgroundCellId\":"
               << failure.second_background_cell_id
               << ",\"firstLocalFace\":"
               << static_cast<unsigned>(failure.first_local_face)
               << ",\"areaMismatch\":" << failure.area_mismatch
               << ",\"firstMomentMismatch\":"
               << failure.first_moment_mismatch << '}';
    }
    output << ']';
}

void write_cell_failures(
    std::ostream& output,
    const std::vector<cartmesh::CompactCellFailure>& failures) {
    output << '[';
    for (std::size_t index = 0; index < failures.size(); ++index) {
        if (index != 0) output << ',';
        const auto& failure = failures[index];
        output << "{\"backgroundCellId\":" << failure.background_cell_id
               << ",\"areaClosureResidual\":"
               << failure.area_closure_residual
               << ",\"boundaryEdgeImbalanceCount\":"
               << failure.boundary_edge_imbalance_count
               << ",\"fluidPieceCount\":" << failure.fluid_piece_count
               << '}';
    }
    output << ']';
}

void ensure_parent(const std::filesystem::path& path) {
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
}

} // 匿名命名空间

int main(int argc, char** argv) {
    try {
        const auto options = parse_options(argc, argv);
        const auto process_start = std::chrono::steady_clock::now();
        const auto surface = cartmesh::read_stl(options.stl);
        const auto diagnostics = cartmesh::diagnose_surface(surface);
        if (!diagnostics.valid_for_stage1_classification()) {
            throw std::invalid_argument(
                "阶段6紧凑路径要求封闭、流形、定向且无自交的 STL");
        }
        const auto extent = surface.bounds().extent();
        const double scale = std::max({extent.x, extent.y, extent.z});
        const double padding = options.padding_fraction * scale;
        const cartmesh::Vec3 padding_vector{padding, padding, padding};
        const cartmesh::AABB domain(
            surface.bounds().minimum() - padding_vector,
            surface.bounds().maximum() + padding_vector);
        const cartmesh::UniformCartesianGrid grid(
            domain, options.resolution, options.resolution,
            options.resolution);
        const cartmesh::TriangulatedSurfaceCutter cutter(surface, 0);
        const auto mesh = cartmesh::build_compact_uniform_cut_cell_mesh(
            grid, cutter,
            {options.geometric_tolerance, options.small_cell_threshold});
        if (options.preview) {
            std::vector<double> fluid_volume_fraction(
                static_cast<std::size_t>(grid.cell_count()), 0.0);
            std::vector<double> cut_cell(fluid_volume_fraction.size(), 0.0);
            std::vector<double> compact_cell_state(
                fluid_volume_fraction.size(), 0.0);
            std::vector<double> small_cut_cell(
                fluid_volume_fraction.size(), 0.0);
            std::size_t explicit_index = 0;
            for (std::uint64_t background = 0;
                 background < grid.cell_count(); ++background) {
                const auto state = mesh.state(background);
                compact_cell_state[static_cast<std::size_t>(background)] =
                    static_cast<double>(state);
                if (state == cartmesh::CompactCellState::full_fluid) {
                    fluid_volume_fraction[static_cast<std::size_t>(background)] =
                        1.0;
                    continue;
                }
                if (state != cartmesh::CompactCellState::explicit_surface) {
                    continue;
                }
                while (explicit_index < mesh.explicit_cells.size() &&
                       mesh.explicit_cells[explicit_index]
                               .background_cell_id < background) {
                    ++explicit_index;
                }
                if (explicit_index == mesh.explicit_cells.size() ||
                    mesh.explicit_cells[explicit_index]
                            .background_cell_id != background) {
                    throw std::runtime_error(
                        "阶段6预览找不到显式单元记录");
                }
                const auto& geometry =
                    mesh.explicit_cells[explicit_index].geometry;
                fluid_volume_fraction[static_cast<std::size_t>(background)] =
                    geometry.volume_fraction;
                cut_cell[static_cast<std::size_t>(background)] =
                    geometry.cut ? 1.0 : 0.0;
                small_cut_cell[static_cast<std::size_t>(background)] =
                    geometry.cut && geometry.volume_fraction <
                                        options.small_cell_threshold
                        ? 1.0
                        : 0.0;
            }
            cartmesh::write_vtu(
                *options.preview, grid,
                {{"fluid_volume_fraction", fluid_volume_fraction},
                 {"cut_cell", cut_cell},
                 {"compact_cell_state", compact_cell_state},
                 {"small_cut_cell", small_cut_cell}});
        }
        std::optional<cartmesh::ScalableOpenFoamWriteStats> openfoam;
        if (options.openfoam_case) {
            openfoam = cartmesh::write_scalable_openfoam_poly_mesh(
                *options.openfoam_case, grid, mesh,
                {{0U, "embedded_wall"}}, options.geometric_tolerance);
        }
        const auto process_end = std::chrono::steady_clock::now();

        ensure_parent(options.report);
        std::ofstream report(options.report, std::ios::trunc);
        if (!report) throw std::runtime_error("无法写入阶段6报告");
        std::vector<std::string_view> run_blockers;
        if (!openfoam) {
            run_blockers.push_back("binary_complete_mesh_export");
        }
        run_blockers.push_back("external_reader");
        if (mesh.background_cell_count < 10000000U) {
            run_blockers.push_back("ten_million_scale");
        }
        if (!options.preview) run_blockers.push_back("viewer_preview");
        report << std::setprecision(17)
               << "{\n"
               << "  \"schema\":\"cartmesh-stage6-compact-v2\",\n"
               << "  \"status\":\""
               << (mesh.invariants_pass()
                       ? (openfoam
                              ? "complete_mesh_exported_external_validation_pending"
                              : "compact_geometry_pass_export_pending")
                       : "fail")
               << "\",\n"
               << "  \"stage6Complete\":false,\n"
               << "  \"solverReadyCutCellMesh\":false,\n"
               << "  \"acceptanceBlockers\":[";
        for (std::size_t index = 0; index < run_blockers.size(); ++index) {
            if (index != 0U) report << ',';
            report << '\"' << run_blockers[index] << '\"';
        }
        report << "],\n"
               << "  \"input\":\"" << options.stl.string() << "\",\n"
               << "  \"inputHashFnv1a64\":\"" << file_hash(options.stl)
               << "\",\n"
               << "  \"triangleCount\":" << surface.triangles().size()
               << ",\n  \"resolution\":" << options.resolution
               << ",\n  \"domain\":{\"minimum\":";
        write_vec3(report, domain.minimum());
        report << ",\"maximum\":";
        write_vec3(report, domain.maximum());
        report << "},\n"
               << "  \"backgroundCellCount\":"
               << mesh.background_cell_count << ",\n"
               << "  \"surfaceCandidateCellCount\":"
               << mesh.surface_candidate_cell_count << ",\n"
               << "  \"fullFluidCellCount\":"
               << mesh.full_fluid_cell_count << ",\n"
               << "  \"fullSolidCellCount\":"
               << mesh.full_solid_cell_count << ",\n"
               << "  \"explicitSurfaceCellCount\":"
               << mesh.explicit_surface_cell_count << ",\n"
               << "  \"cutCellCount\":" << mesh.cut_cell_count << ",\n"
               << "  \"boundaryCellCount\":" << mesh.boundary_cell_count
               << ",\n  \"solverCellCount\":" << mesh.solver_cell_count
               << ",\n  \"globalFluidRegionCount\":"
               << mesh.global_fluid_region_count << ",\n"
               << "  \"totalFluidVolume\":" << mesh.total_fluid_volume
               << ",\n  \"totalEmbeddedBoundaryArea\":"
               << mesh.total_embedded_boundary_area << ",\n"
               << "  \"smallCutCellThreshold\":"
               << mesh.small_cut_cell_threshold << ",\n"
               << "  \"minimumCutCellVolumeFraction\":"
               << mesh.minimum_cut_cell_volume_fraction << ",\n"
               << "  \"smallCutCellCount\":"
               << mesh.small_cut_cell_count << ",\n"
               << "  \"nonclosedCellCount\":"
               << mesh.nonclosed_cell_count << ",\n"
               << "  \"negativeVolumeCellCount\":"
               << mesh.negative_volume_cell_count << ",\n"
               << "  \"sharedFaceMismatchCount\":"
               << mesh.shared_face_mismatch_count << ",\n"
               << "  \"aggregateBoundaryEdgeImbalanceCellCount\":"
               << mesh.aggregate_boundary_edge_imbalance_cell_count << ",\n"
               << "  \"explicitPieceTopologyFailureCount\":"
               << mesh.explicit_piece_topology_failure_count << ",\n"
               << "  \"directFluidSolidFaceCount\":"
               << mesh.direct_fluid_solid_face_count << ",\n"
               << "  \"classificationConflictCount\":"
               << mesh.classification_conflict_count << ",\n"
               << "  \"componentAnalysisPendingCellCount\":"
               << mesh.component_analysis_pending_cell_count << ",\n"
               << "  \"discardedNumericalPieceCount\":"
               << mesh.discarded_numerical_piece_count << ",\n"
               << "  \"discardedNumericalPieceVolume\":"
               << mesh.discarded_numerical_piece_volume << ",\n"
               << "  \"numericallySealedCartesianFaceCount\":"
               << mesh.numerically_sealed_cartesian_face_count << ",\n"
               << "  \"numericallySealedCartesianFaceArea\":"
               << mesh.numerically_sealed_cartesian_face_area << ",\n"
               << "  \"nonclosedCellSamples\":";
        write_cell_failures(report, mesh.nonclosed_cell_samples);
        report << ",\n  \"pendingCellSamples\":";
        write_u64_array(report, mesh.pending_cell_samples);
        report << ",\n  \"sharedFaceFailureSamples\":";
        write_face_failures(report, mesh.shared_face_failure_samples);
        report << ",\n  \"directFluidSolidFaceSamples\":";
        write_face_failures(report,
                            mesh.direct_fluid_solid_face_samples);
        report << ",\n"
               << "  \"maximumCellAreaClosureResidual\":"
               << mesh.maximum_cell_area_closure_residual << ",\n"
               << "  \"maximumSharedFaceAreaMismatch\":"
               << mesh.maximum_shared_face_area_mismatch << ",\n"
               << "  \"maximumSharedFaceFirstMomentMismatch\":"
               << mesh.maximum_shared_face_first_moment_mismatch << ",\n"
               << "  \"compactStorageBytes\":"
               << mesh.compact_storage_bytes << ",\n"
               << "  \"resultHashFnv1a64\":\""
               << cartmesh::fingerprint_hex(mesh.result_hash_fnv1a64)
               << "\",\n"
               << "  \"previewWritten\":"
               << (options.preview ? "true" : "false") << ",\n"
               << "  \"previewScope\":\""
               << (options.preview
                       ? "local_background_cut_cell_preview_only_not_complete_solver_mesh"
                       : "none")
               << "\",\n"
               << "  \"openfoamCellMode\":\"connected_component\",\n";
        if (options.preview) {
            report << "  \"previewVtu\":\""
                   << options.preview->string() << "\",\n";
        }
        if (openfoam) {
            report << "  \"openfoamCase\":\""
                   << options.openfoam_case->string() << "\",\n"
                   << "  \"openfoamBinary\":true,\n"
                   << "  \"openfoamSolverCellCount\":"
                   << openfoam->solver_cell_count << ",\n"
                   << "  \"openfoamPointCount\":" << openfoam->point_count
                   << ",\n  \"openfoamFaceCount\":" << openfoam->face_count
                   << ",\n  \"openfoamInternalFaceCount\":"
                   << openfoam->internal_face_count
                   << ",\n  \"openfoamBoundaryFaceCount\":"
                   << openfoam->boundary_face_count
                   << ",\n  \"openfoamFaceVertexReferenceCount\":"
                   << openfoam->face_vertex_reference_count
                   << ",\n  \"openfoamBackgroundInterfaceCoverageRelaxationCount\":"
                   << openfoam->background_interface_coverage_relaxation_count
                   << ",\n  \"openfoamMaximumBackgroundInterfaceCoverageError\":"
                   << openfoam->maximum_background_interface_coverage_error
                   << ",\n  \"openfoamExplicitPatchCoverageRelaxationCount\":"
                   << openfoam->explicit_patch_coverage_relaxation_count
                   << ",\n  \"openfoamMaximumExplicitPatchCoverageError\":"
                   << openfoam->maximum_explicit_patch_coverage_error
                   << ",\n  \"openfoamTopologyCollapsedFaceCount\":"
                   << openfoam->topology_collapsed_face_count
                   << ",\n  \"openfoamTopologyCollapsedFaceArea\":"
                   << openfoam->topology_collapsed_face_area
                   << ",\n  \"openfoamTopologySealedLoopCount\":"
                   << openfoam->topology_sealed_loop_count
                   << ",\n  \"openfoamTopologySealedEdgeCount\":"
                   << openfoam->topology_sealed_edge_count
                   << ",\n  \"openfoamTopologySealedLoopArea\":"
                   << openfoam->topology_sealed_loop_area
                   << ",\n  \"openfoamMaximumTopologySealedLoopArea\":"
                   << openfoam->maximum_topology_sealed_loop_area
                   << ",\n  \"openfoamWrittenBytes\":"
                   << openfoam->written_bytes
                   << ",\n  \"openfoamTopologyHashFnv1a64\":\""
                   << cartmesh::fingerprint_hex(
                          openfoam->topology_hash_fnv1a64)
                   << "\",\n";
        }
        report
               << "  \"timingsSeconds\":{\"surfaceRasterization\":"
               << mesh.timings.surface_rasterization_seconds
               << ",\"connectedClassification\":"
               << mesh.timings.connected_classification_seconds
               << ",\"localCutCells\":"
               << mesh.timings.local_cut_cell_seconds
               << ",\"topologyAndRegions\":"
               << mesh.timings.topology_and_regions_seconds
               << ",\"compactBuildTotal\":"
               << mesh.timings.total_seconds << ",\"processTotal\":"
               << std::chrono::duration<double>(process_end - process_start)
                      .count();
        if (openfoam) {
            report << ",\"openfoamPreparation\":"
                   << openfoam->preparation_seconds
                   << ",\"openfoamWriting\":"
                   << openfoam->writing_seconds
                   << ",\"openfoamTotal\":"
                   << openfoam->total_seconds;
        }
        report << "},\n"
               << "  \"peakRssBytes\":" << cartmesh::peak_rss_bytes()
               << ",\n  \"runtimeThreads\":1,\n"
               << "  \"hardwareModel\":\"" << hardware_model() << "\",\n"
               << "  \"hardwareMemoryBytes\":" << hardware_memory()
               << ",\n  \"hardwareLogicalCpu\":" << logical_cpu() << ",\n"
               << "  \"buildType\":\"" << build_type() << "\",\n"
               << "  \"compiler\":\"" << compiler() << "\"\n"
               << "}\n";
        report.flush();
        if (!report) throw std::runtime_error("阶段6报告写入失败");
        std::cout << "阶段6紧凑网格：background="
                  << mesh.background_cell_count << " explicit="
                  << mesh.explicit_surface_cell_count << " cut="
                  << mesh.cut_cell_count << " solverCells="
                  << mesh.solver_cell_count << " seconds="
                  << mesh.timings.total_seconds << " hash="
                  << cartmesh::fingerprint_hex(mesh.result_hash_fnv1a64)
                  << '\n';
        if (openfoam) {
            std::cout << "阶段6二进制 OpenFOAM：points="
                      << openfoam->point_count << " faces="
                      << openfoam->face_count << " internal="
                      << openfoam->internal_face_count << " cells="
                      << openfoam->solver_cell_count << " bytes="
                      << openfoam->written_bytes << " seconds="
                      << openfoam->total_seconds << '\n';
        }
        return mesh.invariants_pass() ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "cartmesh_stage6_cli 错误：" << error.what() << '\n';
        return 1;
    }
}
