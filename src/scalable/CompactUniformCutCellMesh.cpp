#include "cartmesh/scalable/CompactUniformCutCellMesh.hpp"

#include "cartmesh/geometry/TriangleBoxIntersection.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace cartmesh {
namespace {

using Clock = std::chrono::steady_clock;
constexpr std::uint32_t no_component =
    std::numeric_limits<std::uint32_t>::max();
constexpr std::uint8_t work_unknown = 0;
constexpr std::uint8_t work_surface = 1;
constexpr std::uint8_t work_fluid = 2;
constexpr std::uint8_t work_solid = 3;
constexpr std::uint8_t work_explicit = 4;
constexpr std::uint8_t work_conflict = 5;
constexpr std::array<std::uint8_t, 6> opposite_face = {1, 0, 3, 2, 5, 4};
constexpr std::size_t maximum_failure_samples = 64;

[[nodiscard]] double seconds(Clock::time_point first,
                             Clock::time_point second) noexcept {
    return std::chrono::duration<double>(second - first).count();
}

[[nodiscard]] AABB expanded_box(const AABB& box, double tolerance) {
    const Vec3 padding{tolerance, tolerance, tolerance};
    return AABB(box.minimum() - padding, box.maximum() + padding);
}

class DisjointSets {
  public:
    explicit DisjointSets(std::size_t size)
        : parent_(size), rank_(size, 0), minimum_key_(size) {
        std::iota(parent_.begin(), parent_.end(), std::size_t{0});
        std::iota(minimum_key_.begin(), minimum_key_.end(), std::uint64_t{0});
    }

    [[nodiscard]] std::size_t find(std::size_t value) {
        if (parent_[value] != value) parent_[value] = find(parent_[value]);
        return parent_[value];
    }

    void set_key(std::size_t value, std::uint64_t key) {
        minimum_key_[value] = key;
    }

    void unite(std::size_t first, std::size_t second) {
        first = find(first);
        second = find(second);
        if (first == second) return;
        if (rank_[first] < rank_[second]) std::swap(first, second);
        parent_[second] = first;
        minimum_key_[first] =
            std::min(minimum_key_[first], minimum_key_[second]);
        if (rank_[first] == rank_[second]) ++rank_[first];
    }

    [[nodiscard]] std::uint64_t minimum_key(std::size_t value) {
        return minimum_key_[find(value)];
    }

  private:
    std::vector<std::size_t> parent_;
    std::vector<std::uint8_t> rank_;
    std::vector<std::uint64_t> minimum_key_;
};

struct Point2 {
    double x{};
    double y{};
};

[[nodiscard]] Point2 project(const Vec3& value, std::size_t normal_axis) {
    if (normal_axis == 0) return {value.y, value.z};
    if (normal_axis == 1) return {value.x, value.z};
    return {value.x, value.y};
}

[[nodiscard]] double signed_area(const std::vector<Point2>& polygon) noexcept {
    double twice_area = 0.0;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const auto& first = polygon[index];
        const auto& second = polygon[(index + 1U) % polygon.size()];
        twice_area += first.x * second.y - first.y * second.x;
    }
    return 0.5 * twice_area;
}

[[nodiscard]] double cross2(const Point2& first, const Point2& second,
                            const Point2& point) noexcept {
    return (second.x - first.x) * (point.y - first.y) -
           (second.y - first.y) * (point.x - first.x);
}

[[nodiscard]] Point2 line_intersection(const Point2& segment_first,
                                       const Point2& segment_second,
                                       const Point2& clip_first,
                                       const Point2& clip_second) noexcept {
    const double segment_x = segment_second.x - segment_first.x;
    const double segment_y = segment_second.y - segment_first.y;
    const double clip_x = clip_second.x - clip_first.x;
    const double clip_y = clip_second.y - clip_first.y;
    const double denominator = segment_x * clip_y - segment_y * clip_x;
    if (denominator == 0.0) return segment_second;
    const double offset_x = clip_first.x - segment_first.x;
    const double offset_y = clip_first.y - segment_first.y;
    const double parameter =
        (offset_x * clip_y - offset_y * clip_x) / denominator;
    return {segment_first.x + parameter * segment_x,
            segment_first.y + parameter * segment_y};
}

[[nodiscard]] double convex_overlap_area(
    const std::vector<Vec3>& first, const std::vector<Vec3>& second,
    std::size_t normal_axis, double length_tolerance) {
    if (first.size() < 3 || second.size() < 3) return 0.0;
    std::vector<Point2> subject;
    std::vector<Point2> clip;
    subject.reserve(first.size());
    clip.reserve(second.size());
    for (const auto& point : first) subject.push_back(project(point, normal_axis));
    for (const auto& point : second) clip.push_back(project(point, normal_axis));
    if (signed_area(subject) < 0.0) std::reverse(subject.begin(), subject.end());
    if (signed_area(clip) < 0.0) std::reverse(clip.begin(), clip.end());
    for (std::size_t edge = 0; edge < clip.size(); ++edge) {
        const Point2 clip_first = clip[edge];
        const Point2 clip_second = clip[(edge + 1U) % clip.size()];
        std::vector<Point2> output;
        if (subject.empty()) break;
        Point2 previous = subject.back();
        bool previous_inside =
            cross2(clip_first, clip_second, previous) >= -length_tolerance;
        for (const auto& current : subject) {
            const bool current_inside =
                cross2(clip_first, clip_second, current) >= -length_tolerance;
            if (current_inside != previous_inside) {
                output.push_back(line_intersection(
                    previous, current, clip_first, clip_second));
            }
            if (current_inside) output.push_back(current);
            previous = current;
            previous_inside = current_inside;
        }
        subject = std::move(output);
    }
    return subject.size() < 3 ? 0.0 : std::abs(signed_area(subject));
}

[[nodiscard]] std::vector<Vec3> piece_face_vertices(
    const FluidPolyhedronPiece& piece, std::size_t face_index) {
    std::vector<Vec3> result;
    const auto& face = piece.polyhedron.faces[face_index];
    result.reserve(face.vertex_indices.size());
    for (const auto vertex : face.vertex_indices) {
        result.push_back(piece.polyhedron.vertices[vertex]);
    }
    return result;
}

template <class Callback>
void for_each_component_on_face(const FluidCellGeometry& cell,
                                std::uint8_t local_face,
                                double area_tolerance,
                                Callback&& callback) {
    if (!cell.cut || cell.fluid_polyhedron_pieces.empty()) {
        if (cell.cartesian_faces[local_face].area > area_tolerance) {
            callback(std::size_t{0}, std::vector<Vec3>{},
                     cell.cartesian_faces[local_face].area);
        }
        return;
    }
    for (const auto& piece : cell.fluid_polyhedron_pieces) {
        for (std::size_t face = 0; face < piece.polyhedron.faces.size(); ++face) {
            const auto& record = piece.polyhedron.faces[face];
            if (record.kind != PolyhedronFaceKind::cartesian ||
                record.source_id != local_face ||
                piece.geometry.faces[face].area <= area_tolerance) {
                continue;
            }
            callback(piece.component_id, piece_face_vertices(piece, face),
                     piece.geometry.faces[face].area);
        }
    }
}

[[nodiscard]] std::uint64_t fluid_cell_storage_bytes(
    const FluidCellGeometry& cell) noexcept {
    std::uint64_t bytes = sizeof(FluidCellGeometry);
    for (const auto& face : cell.cartesian_faces) {
        bytes += face.oriented_boundary_loops.capacity() *
                 sizeof(std::vector<Vec3>);
        for (const auto& loop : face.oriented_boundary_loops) {
            bytes += loop.capacity() * sizeof(Vec3);
        }
    }
    bytes += cell.embedded_boundary_faces.capacity() *
             sizeof(EmbeddedBoundaryFaceGeometry);
    for (const auto& face : cell.embedded_boundary_faces) {
        bytes += face.vertices.capacity() * sizeof(Vec3);
    }
    bytes += cell.fluid_polyhedron_pieces.capacity() *
             sizeof(FluidPolyhedronPiece);
    for (const auto& piece : cell.fluid_polyhedron_pieces) {
        bytes += piece.polyhedron.vertices.capacity() * sizeof(Vec3);
        bytes += piece.polyhedron.faces.capacity() * sizeof(PolyhedronFace);
        for (const auto& face : piece.polyhedron.faces) {
            bytes += face.vertex_indices.capacity() * sizeof(std::uint32_t);
        }
        bytes += piece.geometry.faces.capacity() * sizeof(PolygonGeometry);
    }
    bytes += cell.fluid_component_region_ids.capacity() *
             sizeof(std::uint64_t);
    return bytes;
}

[[nodiscard]] bool explicit_piece_topology_closed(
    const FluidPolyhedronPiece& piece, double area_tolerance) {
    if (!piece.geometry.positive_volume || !(piece.geometry.volume > 0.0) ||
        piece.polyhedron.faces.size() != piece.geometry.faces.size()) {
        return false;
    }
    std::map<std::pair<std::uint32_t, std::uint32_t>,
             std::pair<std::uint64_t, std::int64_t>> edges;
    Vec3 area_vector_sum{};
    for (std::size_t face_index = 0;
         face_index < piece.polyhedron.faces.size(); ++face_index) {
        const auto& face = piece.polyhedron.faces[face_index];
        if (face.vertex_indices.size() < 3) return false;
        area_vector_sum =
            area_vector_sum + piece.geometry.faces[face_index].area_vector;
        for (std::size_t edge = 0; edge < face.vertex_indices.size(); ++edge) {
            const auto first = face.vertex_indices[edge];
            const auto second =
                face.vertex_indices[(edge + 1U) % face.vertex_indices.size()];
            if (first >= piece.polyhedron.vertices.size() ||
                second >= piece.polyhedron.vertices.size() || first == second) {
                return false;
            }
            const auto key = std::minmax(first, second);
            auto& entry = edges[{key.first, key.second}];
            ++entry.first;
            entry.second += first < second ? 1 : -1;
        }
    }
    if (norm(area_vector_sum) > area_tolerance) return false;
    return !edges.empty() &&
           std::all_of(edges.begin(), edges.end(), [](const auto& entry) {
               return entry.second.first == 2U && entry.second.second == 0;
           });
}

void hash_byte(std::uint64_t& hash, std::uint8_t value) noexcept {
    hash ^= value;
    hash *= 1099511628211ULL;
}

void hash_u64(std::uint64_t& hash, std::uint64_t value) noexcept {
    for (std::uint32_t shift = 0; shift < 64U; shift += 8U) {
        hash_byte(hash,
                  static_cast<std::uint8_t>((value >> shift) & 0xffU));
    }
}

void hash_double(std::uint64_t& hash, double value) noexcept {
    hash_u64(hash, std::bit_cast<std::uint64_t>(value));
}

void hash_vec3(std::uint64_t& hash, const Vec3& value) noexcept {
    hash_double(hash, value.x);
    hash_double(hash, value.y);
    hash_double(hash, value.z);
}

[[nodiscard]] std::uint64_t mesh_hash(
    const UniformCartesianGrid& grid,
    const CompactUniformCutCellMesh& mesh) noexcept {
    std::uint64_t hash = 14695981039346656037ULL;
    hash_vec3(hash, grid.domain().minimum());
    hash_vec3(hash, grid.domain().maximum());
    hash_u64(hash, grid.nx());
    hash_u64(hash, grid.ny());
    hash_u64(hash, grid.nz());
    for (const auto state : mesh.cell_states) hash_byte(hash, state);
    for (const auto region : mesh.full_component_region_ids) {
        hash_u64(hash, region);
    }
    for (const auto& record : mesh.explicit_cells) {
        const auto& cell = record.geometry;
        hash_u64(hash, record.background_cell_id);
        hash_double(hash, cell.volume);
        hash_double(hash, cell.volume_fraction);
        hash_vec3(hash, cell.centroid);
        hash_byte(hash, cell.cut ? 1U : 0U);
        for (const auto& face : cell.cartesian_faces) {
            hash_double(hash, face.area);
            hash_vec3(hash, face.centroid);
        }
        for (const auto& face : cell.embedded_boundary_faces) {
            hash_u64(hash, face.boundary_id);
            hash_double(hash, face.area);
            hash_vec3(hash, face.centroid);
            for (const auto& vertex : face.vertices) hash_vec3(hash, vertex);
        }
        for (const auto& piece : cell.fluid_polyhedron_pieces) {
            hash_u64(hash, piece.component_id);
            hash_u64(hash, piece.global_region_id);
            hash_double(hash, piece.geometry.volume);
            for (const auto& vertex : piece.polyhedron.vertices) {
                hash_vec3(hash, vertex);
            }
        }
    }
    for (const auto volume : mesh.global_fluid_region_volumes) {
        hash_double(hash, volume);
    }
    return hash;
}

} // 匿名命名空间

CompactCellState CompactUniformCutCellMesh::state(
    std::uint64_t background_id) const {
    if (background_id >= cell_states.size()) {
        throw std::out_of_range("\u7d27\u51d1\u7f51\u683c\u80cc\u666f ID \u8d8a\u754c");
    }
    return static_cast<CompactCellState>(
        cell_states[static_cast<std::size_t>(background_id)]);
}

const CompactCutCellRecord* CompactUniformCutCellMesh::find_explicit_cell(
    std::uint64_t background_id) const noexcept {
    const auto found = std::lower_bound(
        explicit_cells.begin(), explicit_cells.end(), background_id,
        [](const CompactCutCellRecord& record, std::uint64_t id) {
            return record.background_cell_id < id;
        });
    return found != explicit_cells.end() &&
                   found->background_cell_id == background_id
               ? &*found
               : nullptr;
}

CompactCutCellRecord* CompactUniformCutCellMesh::find_explicit_cell(
    std::uint64_t background_id) noexcept {
    const auto found = std::lower_bound(
        explicit_cells.begin(), explicit_cells.end(), background_id,
        [](const CompactCutCellRecord& record, std::uint64_t id) {
            return record.background_cell_id < id;
        });
    return found != explicit_cells.end() &&
                   found->background_cell_id == background_id
               ? &*found
               : nullptr;
}

bool CompactUniformCutCellMesh::invariants_pass() const noexcept {
    return background_cell_count == cell_states.size() &&
           full_fluid_cell_count + full_solid_cell_count ==
               background_cell_count &&
           explicit_surface_cell_count == explicit_cells.size() &&
           nonclosed_cell_count == 0 && negative_volume_cell_count == 0 &&
           component_analysis_pending_cell_count == 0 &&
           classification_conflict_count == 0 &&
           shared_face_mismatch_count == 0 &&
           direct_fluid_solid_face_count == 0 &&
           global_fluid_region_count == global_fluid_region_volumes.size();
}

const char* compact_cell_state_name(CompactCellState state) noexcept {
    switch (state) {
    case CompactCellState::solid:
        return "solid";
    case CompactCellState::full_fluid:
        return "full_fluid";
    case CompactCellState::explicit_surface:
        return "explicit_surface";
    case CompactCellState::conflict:
        return "conflict";
    }
    return "conflict";
}

CompactUniformCutCellMesh build_compact_uniform_cut_cell_mesh(
    const UniformCartesianGrid& grid,
    const TriangulatedSurfaceCutter& cutter,
    const CompactUniformBuildOptions& options) {
    if (options.geometric_tolerance < 0.0 ||
        !std::isfinite(options.geometric_tolerance)) {
        throw std::invalid_argument("\u7d27\u51d1 Cut-cell \u5bb9\u5dee\u5fc5\u987b\u662f\u975e\u8d1f\u6709\u9650\u6570");
    }
    if (!(options.small_cut_cell_threshold >= 0.0) ||
        !(options.small_cut_cell_threshold <= 1.0) ||
        !std::isfinite(options.small_cut_cell_threshold)) {
        throw std::invalid_argument("\u5c0f Cut-cell \u9608\u503c\u5fc5\u987b\u4f4d\u4e8e 0..1");
    }
    if (grid.cell_count() >
        static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("\u7d27\u51d1\u7f51\u683c\u8d85\u51fa\u5f53\u524d\u5e73\u53f0 size_t \u8303\u56f4");
    }

    const auto total_start = Clock::now();
    CompactUniformCutCellMesh result;
    result.background_cell_count = grid.cell_count();
    result.small_cut_cell_threshold = options.small_cut_cell_threshold;
    result.cell_states.assign(static_cast<std::size_t>(grid.cell_count()),
                              work_unknown);
    result.full_component_labels.assign(
        static_cast<std::size_t>(grid.cell_count()), no_component);
    std::vector<std::uint64_t> candidate_ids;

    const LocalTriangulatedCutCellBuilder local_builder(
        cutter, grid.domain(), grid.spacing(), options.geometric_tolerance,
        true);
    const double length_tolerance = local_builder.length_tolerance();
    const double area_tolerance = local_builder.area_tolerance();
    const double closure_tolerance = local_builder.closure_tolerance();
    const auto raster_start = Clock::now();
    const auto coordinate_range = [&](double minimum, double maximum,
                                      double domain_minimum, double spacing,
                                      std::uint32_t count) {
        const auto first_raw = static_cast<std::int64_t>(std::floor(
            (minimum - length_tolerance - domain_minimum) / spacing));
        const auto last_raw = static_cast<std::int64_t>(std::floor(
            (maximum + length_tolerance - domain_minimum) / spacing));
        const auto first = std::clamp<std::int64_t>(
            first_raw, 0, static_cast<std::int64_t>(count) - 1);
        const auto last = std::clamp<std::int64_t>(
            last_raw, 0, static_cast<std::int64_t>(count) - 1);
        return std::pair<std::uint32_t, std::uint32_t>{
            static_cast<std::uint32_t>(first),
            static_cast<std::uint32_t>(last)};
    };
    const auto& domain_minimum = grid.domain().minimum();
    const auto spacing = grid.spacing();
    for (const auto& triangle : cutter.oriented_surface().triangles()) {
        const auto& bounds = triangle.bounds();
        const auto [first_i, last_i] = coordinate_range(
            bounds.minimum().x, bounds.maximum().x, domain_minimum.x,
            spacing.x, grid.nx());
        const auto [first_j, last_j] = coordinate_range(
            bounds.minimum().y, bounds.maximum().y, domain_minimum.y,
            spacing.y, grid.ny());
        const auto [first_k, last_k] = coordinate_range(
            bounds.minimum().z, bounds.maximum().z, domain_minimum.z,
            spacing.z, grid.nz());
        for (std::uint32_t k = first_k; k <= last_k; ++k) {
            for (std::uint32_t j = first_j; j <= last_j; ++j) {
                for (std::uint32_t i = first_i; i <= last_i; ++i) {
                    const CellKey key{0, i, j, k};
                    const auto id = grid.linear_id(key);
                    auto& state = result.cell_states[static_cast<std::size_t>(id)];
                    if (state == work_surface) continue;
                    if (!triangle_intersects_aabb(
                            triangle,
                            expanded_box(grid.cell_bounds(key),
                                         length_tolerance))) {
                        continue;
                    }
                    state = work_surface;
                    candidate_ids.push_back(id);
                }
            }
        }
    }
    std::sort(candidate_ids.begin(), candidate_ids.end());
    result.surface_candidate_cell_count = candidate_ids.size();
    const auto raster_end = Clock::now();

    const auto classify_start = Clock::now();
    SurfaceClassifier classifier(cutter.bvh(), cutter.length_tolerance());
    std::vector<std::uint64_t> queue;
    std::vector<double> component_volumes;
    std::vector<std::uint64_t> component_minimum_ids;
    queue.reserve(std::min<std::uint64_t>(grid.cell_count(), 1'000'000ULL));
    const std::uint64_t plane =
        static_cast<std::uint64_t>(grid.nx()) * grid.ny();
    const auto push_neighbor = [&](std::uint64_t neighbor,
                                   std::uint8_t target,
                                   std::uint32_t component) {
        auto& state = result.cell_states[static_cast<std::size_t>(neighbor)];
        if (state != work_unknown) return;
        state = target;
        if (target == work_fluid) {
            result.full_component_labels[static_cast<std::size_t>(neighbor)] =
                component;
        }
        queue.push_back(neighbor);
    };
    for (std::uint64_t seed = 0; seed < grid.cell_count(); ++seed) {
        if (result.cell_states[static_cast<std::size_t>(seed)] != work_unknown) {
            continue;
        }
        const auto classification =
            classifier.classify(grid.cell_center(grid.cell_key(seed)))
                .classification;
        std::uint8_t target = work_conflict;
        if (classification == PointClassification::outside) {
            target = work_fluid;
        } else if (classification == PointClassification::inside) {
            target = work_solid;
        } else {
            ++result.classification_conflict_count;
        }
        const std::uint32_t component =
            target == work_fluid
                ? static_cast<std::uint32_t>(component_volumes.size())
                : no_component;
        queue.clear();
        result.cell_states[static_cast<std::size_t>(seed)] = target;
        if (target == work_fluid) {
            result.full_component_labels[static_cast<std::size_t>(seed)] =
                component;
        }
        queue.push_back(seed);
        std::size_t head = 0;
        while (head < queue.size()) {
            const std::uint64_t id = queue[head++];
            const std::uint64_t k = id / plane;
            const std::uint64_t remainder = id - k * plane;
            const std::uint64_t j = remainder / grid.nx();
            const std::uint64_t i = remainder - j * grid.nx();
            if (i > 0) push_neighbor(id - 1U, target, component);
            if (i + 1U < grid.nx()) push_neighbor(id + 1U, target, component);
            if (j > 0) push_neighbor(id - grid.nx(), target, component);
            if (j + 1U < grid.ny()) {
                push_neighbor(id + grid.nx(), target, component);
            }
            if (k > 0) push_neighbor(id - plane, target, component);
            if (k + 1U < grid.nz()) push_neighbor(id + plane, target, component);
        }
        if (target == work_fluid) {
            component_volumes.push_back(
                static_cast<double>(queue.size()) * grid.cell_volume());
            component_minimum_ids.push_back(seed);
            result.full_fluid_cell_count += queue.size();
            result.total_fluid_volume +=
                static_cast<double>(queue.size()) * grid.cell_volume();
        } else if (target == work_solid) {
            result.full_solid_cell_count += queue.size();
        }
    }
    const auto classify_end = Clock::now();

    const auto cut_start = Clock::now();
    result.explicit_cells.reserve(candidate_ids.size());
    for (const auto id : candidate_ids) {
        auto local = local_builder.build(id, grid.cell_bounds(grid.cell_key(id)));
        result.discarded_numerical_piece_count +=
            local.discarded_numerical_piece_count;
        result.discarded_numerical_piece_volume +=
            local.discarded_numerical_piece_volume;
        if (local.component_analysis_pending) {
            ++result.component_analysis_pending_cell_count;
            if (result.pending_cell_samples.size() < maximum_failure_samples) {
                result.pending_cell_samples.push_back(id);
            }
        }
        if (local.classification_conflict) {
            ++result.classification_conflict_count;
        }
        if (!local.has_fluid) {
            result.cell_states[static_cast<std::size_t>(id)] = work_solid;
            ++result.full_solid_cell_count;
            continue;
        }
        ++result.full_fluid_cell_count;
        result.total_fluid_volume += local.cell.volume;
        if (local.cell.cut || !local.cell.embedded_boundary_faces.empty()) {
            result.cell_states[static_cast<std::size_t>(id)] = work_explicit;
            result.total_embedded_boundary_area +=
                local.embedded_boundary_area;
            if (!local.cell.embedded_boundary_faces.empty()) {
                ++result.boundary_cell_count;
            }
            if (local.cell.cut) {
                ++result.cut_cell_count;
                result.minimum_cut_cell_volume_fraction = std::min(
                    result.minimum_cut_cell_volume_fraction,
                    local.cell.volume_fraction);
                if (local.cell.volume_fraction <
                    options.small_cut_cell_threshold) {
                    ++result.small_cut_cell_count;
                }
            }
            result.explicit_cells.push_back(
                {id, std::move(local.cell), 0, {}});
        } else {
            result.cell_states[static_cast<std::size_t>(id)] = work_fluid;
            const auto component =
                static_cast<std::uint32_t>(component_volumes.size());
            result.full_component_labels[static_cast<std::size_t>(id)] =
                component;
            component_volumes.push_back(local.cell.volume);
            component_minimum_ids.push_back(id);
        }
    }
    result.explicit_surface_cell_count = result.explicit_cells.size();
    result.full_fluid_component_count = component_volumes.size();
    const auto cut_end = Clock::now();

    const auto topology_start = Clock::now();
    std::uint64_t component_node_count = component_volumes.size();
    for (auto& record : result.explicit_cells) {
        record.component_node_offset = component_node_count;
        const auto count = std::max<std::size_t>(
            record.geometry.fluid_component_count, 1U);
        component_node_count += count;
    }
    if (component_node_count >
        static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("\u7d27\u51d1\u6d41\u4f53\u5206\u91cf\u8282\u70b9\u8d85\u51fa size_t \u8303\u56f4");
    }
    DisjointSets regions(static_cast<std::size_t>(component_node_count));
    std::vector<double> node_volumes(
        static_cast<std::size_t>(component_node_count), 0.0);
    for (std::size_t component = 0; component < component_volumes.size();
         ++component) {
        node_volumes[component] = component_volumes[component];
        regions.set_key(component, component_minimum_ids[component] << 8U);
    }
    for (auto& record : result.explicit_cells) {
        const auto count = std::max<std::size_t>(
            record.geometry.fluid_component_count, 1U);
        for (std::size_t component = 0; component < count; ++component) {
            const auto node = static_cast<std::size_t>(
                record.component_node_offset + component);
            regions.set_key(node,
                            (record.background_cell_id << 8U) + component);
        }
        if (record.geometry.fluid_polyhedron_pieces.empty()) {
            node_volumes[static_cast<std::size_t>(
                record.component_node_offset)] += record.geometry.volume;
        } else {
            for (const auto& piece : record.geometry.fluid_polyhedron_pieces) {
                if (piece.component_id >= count) {
                    ++result.component_analysis_pending_cell_count;
                    continue;
                }
                node_volumes[static_cast<std::size_t>(
                    record.component_node_offset + piece.component_id)] +=
                    piece.geometry.volume;
            }
        }
        result.maximum_cell_area_closure_residual = std::max(
            result.maximum_cell_area_closure_residual,
            record.geometry.area_vector_closure_residual);
        if (!record.geometry.boundary_edge_closed) {
            ++result.aggregate_boundary_edge_imbalance_cell_count;
        }
        bool pieces_closed = true;
        if (record.geometry.cut) {
            pieces_closed =
                !record.geometry.fluid_polyhedron_pieces.empty() &&
                std::all_of(
                    record.geometry.fluid_polyhedron_pieces.begin(),
                    record.geometry.fluid_polyhedron_pieces.end(),
                    [&](const FluidPolyhedronPiece& piece) {
                        return explicit_piece_topology_closed(
                            piece, closure_tolerance);
                    });
            if (!pieces_closed) {
                ++result.explicit_piece_topology_failure_count;
            }
        }
        const bool topology_closed =
            record.geometry.cut ? pieces_closed
                                : record.geometry.boundary_edge_closed;
        if (record.geometry.area_vector_closure_residual > closure_tolerance ||
            !topology_closed) {
            ++result.nonclosed_cell_count;
            if (result.nonclosed_cell_samples.size() < maximum_failure_samples) {
                result.nonclosed_cell_samples.push_back(
                    {record.background_cell_id,
                     record.geometry.area_vector_closure_residual,
                     record.geometry.boundary_edge_imbalance_count,
                     record.geometry.fluid_piece_count});
            }
        }
        if (!(record.geometry.volume > 0.0)) {
            ++result.negative_volume_cell_count;
        }
    }

    const auto full_node = [&](std::uint64_t id) -> std::size_t {
        const auto label =
            result.full_component_labels[static_cast<std::size_t>(id)];
        if (label == no_component || label >= component_volumes.size()) {
            throw std::runtime_error("\u5168\u6d41\u4f53\u5355\u5143\u7f3a\u5c11\u7d27\u51d1\u5206\u91cf label");
        }
        return label;
    };

    const auto connect_full_to_explicit = [&](std::uint64_t full_id,
                                               const CompactCutCellRecord& record,
                                               std::uint8_t explicit_face) {
        const auto& aperture = record.geometry.cartesian_faces[explicit_face];
        if (aperture.area <= area_tolerance) return;
        if (!record.geometry.cut ||
            record.geometry.fluid_polyhedron_pieces.empty()) {
            regions.unite(full_node(full_id),
                          static_cast<std::size_t>(record.component_node_offset));
            return;
        }
        for_each_component_on_face(
            record.geometry, explicit_face, area_tolerance,
            [&](std::size_t component, const std::vector<Vec3>&, double) {
                regions.unite(
                    full_node(full_id),
                    static_cast<std::size_t>(record.component_node_offset +
                                             component));
            });
    };

    const auto connect_explicit_pair = [&](const CompactCutCellRecord& first,
                                            const CompactCutCellRecord& second,
                                            std::uint8_t first_face) {
        const std::uint8_t second_face = opposite_face[first_face];
        if ((!first.geometry.cut ||
             first.geometry.fluid_polyhedron_pieces.empty()) &&
            (!second.geometry.cut ||
             second.geometry.fluid_polyhedron_pieces.empty())) {
            if (first.geometry.cartesian_faces[first_face].area >
                    area_tolerance &&
                second.geometry.cartesian_faces[second_face].area >
                    area_tolerance) {
                regions.unite(
                    static_cast<std::size_t>(first.component_node_offset),
                    static_cast<std::size_t>(second.component_node_offset));
            }
            return;
        }
        for_each_component_on_face(
            first.geometry, first_face, area_tolerance,
            [&](std::size_t first_component,
                const std::vector<Vec3>& first_polygon, double first_area) {
                for_each_component_on_face(
                    second.geometry, second_face, area_tolerance,
                    [&](std::size_t second_component,
                        const std::vector<Vec3>& second_polygon,
                        double second_area) {
                        bool overlaps = false;
                        if (first_polygon.empty() || second_polygon.empty()) {
                            overlaps = std::min(first_area, second_area) >
                                       area_tolerance;
                        } else {
                            overlaps = convex_overlap_area(
                                first_polygon, second_polygon,
                                first_face / 2U, length_tolerance) >
                                       area_tolerance;
                        }
                        if (overlaps) {
                            regions.unite(
                                static_cast<std::size_t>(
                                    first.component_node_offset +
                                    first_component),
                                static_cast<std::size_t>(
                                    second.component_node_offset +
                                    second_component));
                        }
                    });
            });
    };

    const double numerical_seal_area_tolerance = std::max(
        area_tolerance,
        512.0 * length_tolerance *
            std::max({spacing.x, spacing.y, spacing.z}));
    const auto seal_numerical_cartesian_face = [&](
        CompactCutCellRecord& record, std::uint8_t local_face) {
        std::uint64_t boundary_id = 0;
        if (!record.geometry.embedded_boundary_faces.empty()) {
            boundary_id =
                record.geometry.embedded_boundary_faces.front().boundary_id;
        }
        double sealed_area = 0.0;
        for (std::size_t piece_index = 0;
             piece_index < record.geometry.fluid_polyhedron_pieces.size();
             ++piece_index) {
            auto& piece =
                record.geometry.fluid_polyhedron_pieces[piece_index];
            for (std::size_t face_index = 0;
                 face_index < piece.polyhedron.faces.size(); ++face_index) {
                const auto& face = piece.polyhedron.faces[face_index];
                if (face.kind != PolyhedronFaceKind::cartesian ||
                    face.source_id != local_face ||
                    piece.geometry.faces[face_index].area <= area_tolerance) {
                    continue;
                }
                const auto& geometry = piece.geometry.faces[face_index];
                auto vertices = piece_face_vertices(piece, face_index);
                record.geometry.embedded_boundary_faces.push_back(
                    {boundary_id, geometry.area, geometry.centroid,
                     geometry.outward_normal, std::move(vertices)});
                record.numerically_sealed_faces.push_back(
                    {piece_index, face_index, local_face, boundary_id,
                     geometry.area});
                sealed_area += geometry.area;
                ++result.numerically_sealed_cartesian_face_count;
            }
        }
        auto& aperture = record.geometry.cartesian_faces[local_face];
        aperture.area = 0.0;
        aperture.area_fraction = 0.0;
        result.numerically_sealed_cartesian_face_area += sealed_area;
        result.total_embedded_boundary_area += sealed_area;
        return sealed_area;
    };

    const auto process_pair = [&](std::uint64_t first_id,
                                  std::uint64_t second_id,
                                  std::uint8_t first_face) {
        const auto first_state =
            result.cell_states[static_cast<std::size_t>(first_id)];
        const auto second_state =
            result.cell_states[static_cast<std::size_t>(second_id)];
        const bool first_fluid =
            first_state == work_fluid || first_state == work_explicit;
        const bool second_fluid =
            second_state == work_fluid || second_state == work_explicit;
        if (!first_fluid && !second_fluid) return;
        if (first_state == work_conflict || second_state == work_conflict) {
            ++result.classification_conflict_count;
            return;
        }
        const std::uint8_t second_face = opposite_face[first_face];
        const double full_area = first_face < 2U
                                     ? spacing.y * spacing.z
                                     : (first_face < 4U
                                            ? spacing.x * spacing.z
                                            : spacing.x * spacing.y);
        auto* first_record =
            first_state == work_explicit
                ? result.find_explicit_cell(first_id)
                : nullptr;
        auto* second_record =
            second_state == work_explicit
                ? result.find_explicit_cell(second_id)
                : nullptr;
        const double first_area =
            first_state == work_fluid
                ? full_area
                : (first_record
                       ? first_record->geometry.cartesian_faces[first_face].area
                       : 0.0);
        const double second_area =
            second_state == work_fluid
                ? full_area
                : (second_record
                       ? second_record->geometry.cartesian_faces[second_face].area
                       : 0.0);
        const double area_mismatch = std::abs(first_area - second_area);
        result.maximum_shared_face_area_mismatch = std::max(
            result.maximum_shared_face_area_mismatch, area_mismatch);
        if (!first_fluid || !second_fluid) {
            const double fluid_area = first_fluid ? first_area : second_area;
            if (fluid_area > area_tolerance) {
                CompactCutCellRecord* fluid_record =
                    first_fluid ? first_record : second_record;
                const std::uint8_t fluid_face =
                    first_fluid ? first_face : second_face;
                if (fluid_record &&
                    fluid_area <= numerical_seal_area_tolerance) {
                    const double sealed = seal_numerical_cartesian_face(
                        *fluid_record, fluid_face);
                    if (std::abs(sealed - fluid_area) >
                        numerical_seal_area_tolerance) {
                        ++result.direct_fluid_solid_face_count;
                    }
                } else {
                    ++result.direct_fluid_solid_face_count;
                    if (result.direct_fluid_solid_face_samples.size() <
                        maximum_failure_samples) {
                        result.direct_fluid_solid_face_samples.push_back(
                            {first_id, second_id, first_face, fluid_area, 0.0});
                    }
                }
            }
            return;
        }
        if (area_mismatch > closure_tolerance) {
            ++result.shared_face_mismatch_count;
            if (result.shared_face_failure_samples.size() <
                maximum_failure_samples) {
                result.shared_face_failure_samples.push_back(
                    {first_id, second_id, first_face, area_mismatch, 0.0});
            }
        }
        if (std::max(first_area, second_area) <= area_tolerance) return;
        ++result.internal_background_connection_count;
        if (first_state == work_fluid && second_state == work_fluid) {
            regions.unite(full_node(first_id), full_node(second_id));
        } else if (first_state == work_fluid && second_record) {
            connect_full_to_explicit(first_id, *second_record, second_face);
        } else if (second_state == work_fluid && first_record) {
            connect_full_to_explicit(second_id, *first_record, first_face);
        } else if (first_record && second_record) {
            connect_explicit_pair(*first_record, *second_record, first_face);
        }
        if (first_record && second_record && first_area > area_tolerance &&
            second_area > area_tolerance) {
            const auto first_box = grid.cell_bounds(grid.cell_key(first_id));
            Vec3 reference = first_box.center();
            if (first_face == 0U) reference.x = first_box.minimum().x;
            if (first_face == 1U) reference.x = first_box.maximum().x;
            if (first_face == 2U) reference.y = first_box.minimum().y;
            if (first_face == 3U) reference.y = first_box.maximum().y;
            if (first_face == 4U) reference.z = first_box.minimum().z;
            if (first_face == 5U) reference.z = first_box.maximum().z;
            const double moment_mismatch = norm(
                (first_record->geometry.cartesian_faces[first_face].centroid -
                 reference) * first_area -
                (second_record->geometry.cartesian_faces[second_face].centroid -
                 reference) * second_area);
            result.maximum_shared_face_first_moment_mismatch = std::max(
                result.maximum_shared_face_first_moment_mismatch,
                moment_mismatch);
            const double moment_tolerance =
                area_tolerance *
                    std::max({spacing.x, spacing.y, spacing.z}) +
                local_builder.topology_length_tolerance() *
                    std::max(first_area, second_area);
            if (moment_mismatch > moment_tolerance) {
                ++result.shared_face_mismatch_count;
                if (result.shared_face_failure_samples.size() <
                    maximum_failure_samples) {
                    result.shared_face_failure_samples.push_back(
                        {first_id, second_id, first_face, area_mismatch,
                         moment_mismatch});
                }
            }
        }
    };

    for (std::uint32_t k = 0; k < grid.nz(); ++k) {
        for (std::uint32_t j = 0; j < grid.ny(); ++j) {
            for (std::uint32_t i = 0; i < grid.nx(); ++i) {
                const std::uint64_t id = grid.linear_id({0, i, j, k});
                if (i + 1U < grid.nx()) process_pair(id, id + 1U, 1U);
                if (j + 1U < grid.ny()) {
                    process_pair(id, id + grid.nx(), 3U);
                }
                if (k + 1U < grid.nz()) process_pair(id, id + plane, 5U);
            }
        }
    }

    std::vector<std::size_t> roots;
    roots.reserve(static_cast<std::size_t>(component_node_count));
    for (std::size_t node = 0; node < component_node_count; ++node) {
        const auto root = regions.find(node);
        if (std::find(roots.begin(), roots.end(), root) == roots.end()) {
            roots.push_back(root);
        }
    }
    std::sort(roots.begin(), roots.end(), [&](std::size_t first,
                                               std::size_t second) {
        return regions.minimum_key(first) < regions.minimum_key(second);
    });
    std::vector<std::uint64_t> root_region(component_node_count,
                                           std::uint64_t{0});
    for (std::size_t region = 0; region < roots.size(); ++region) {
        root_region[regions.find(roots[region])] = region;
    }
    result.global_fluid_region_count = roots.size();
    result.global_fluid_region_volumes.assign(roots.size(), 0.0);
    for (std::size_t node = 0; node < component_node_count; ++node) {
        const auto region = root_region[regions.find(node)];
        result.global_fluid_region_volumes[region] += node_volumes[node];
    }
    result.full_component_region_ids.resize(component_volumes.size());
    for (std::size_t component = 0; component < component_volumes.size();
         ++component) {
        result.full_component_region_ids[component] =
            root_region[regions.find(component)];
    }
    for (auto& record : result.explicit_cells) {
        const auto count = std::max<std::size_t>(
            record.geometry.fluid_component_count, 1U);
        record.geometry.fluid_component_region_ids.resize(count);
        for (std::size_t component = 0; component < count; ++component) {
            record.geometry.fluid_component_region_ids[component] =
                root_region[regions.find(static_cast<std::size_t>(
                    record.component_node_offset + component))];
        }
        for (auto& piece : record.geometry.fluid_polyhedron_pieces) {
            if (piece.component_id < count) {
                piece.global_region_id =
                    record.geometry.fluid_component_region_ids[
                        piece.component_id];
            }
        }
    }

    result.solver_cell_count =
        result.full_fluid_cell_count - result.explicit_cells.size();
    for (const auto& record : result.explicit_cells) {
        result.solver_cell_count +=
            record.geometry.cut
                ? record.geometry.fluid_polyhedron_pieces.size()
                : 1U;
    }

    for (auto& state : result.cell_states) {
        if (state == work_solid) {
            state = static_cast<std::uint8_t>(CompactCellState::solid);
        } else if (state == work_fluid) {
            state = static_cast<std::uint8_t>(CompactCellState::full_fluid);
        } else if (state == work_explicit) {
            state =
                static_cast<std::uint8_t>(CompactCellState::explicit_surface);
        } else {
            state = static_cast<std::uint8_t>(CompactCellState::conflict);
        }
    }

    result.compact_storage_bytes =
        result.cell_states.capacity() * sizeof(std::uint8_t) +
        result.full_component_labels.capacity() * sizeof(std::uint32_t) +
        result.full_component_region_ids.capacity() * sizeof(std::uint64_t) +
        result.global_fluid_region_volumes.capacity() * sizeof(double) +
        result.nonclosed_cell_samples.capacity() * sizeof(CompactCellFailure) +
        result.pending_cell_samples.capacity() * sizeof(std::uint64_t) +
        result.shared_face_failure_samples.capacity() *
            sizeof(CompactFaceFailure) +
        result.direct_fluid_solid_face_samples.capacity() *
            sizeof(CompactFaceFailure) +
        result.explicit_cells.capacity() * sizeof(CompactCutCellRecord);
    for (const auto& record : result.explicit_cells) {
        result.compact_storage_bytes += fluid_cell_storage_bytes(record.geometry);
        result.compact_storage_bytes +=
            record.numerically_sealed_faces.capacity() *
            sizeof(NumericallySealedCartesianFace);
    }
    result.result_hash_fnv1a64 = mesh_hash(grid, result);
    const auto topology_end = Clock::now();
    result.timings.surface_rasterization_seconds =
        seconds(raster_start, raster_end);
    result.timings.connected_classification_seconds =
        seconds(classify_start, classify_end);
    result.timings.local_cut_cell_seconds = seconds(cut_start, cut_end);
    result.timings.topology_and_regions_seconds =
        seconds(topology_start, topology_end);
    result.timings.total_seconds = seconds(total_start, topology_end);
    return result;
}

} // 命名空间 cartmesh
