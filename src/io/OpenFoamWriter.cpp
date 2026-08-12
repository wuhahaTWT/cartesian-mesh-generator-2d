#include "cartmesh/io/OpenFoamWriter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace cartmesh {
namespace {

constexpr std::array<std::uint8_t, 6> opposite_face = {1, 0, 3, 2, 5, 4};

struct Point2 {
    double u{};
    double v{};
};

struct ComponentCell {
    std::uint64_t background_cell_id{};
    std::size_t component_id{};
    std::size_t local_piece_id{};
    const FluidPolyhedronPiece* piece{};
    bool full_cartesian{};
};

struct CartesianPatch {
    std::size_t owner{};
    std::uint8_t local_face{};
    std::vector<Vec3> vertices;
    double area{};
    double covered_area{};
};

struct PartitionPatch {
    std::size_t owner{};
    Vec3 outward_normal{};
    std::vector<Vec3> vertices;
    double area{};
    double covered_area{};
};

struct OutputFace {
    std::vector<Vec3> vertices;
    std::size_t owner{};
    std::size_t neighbor{std::numeric_limits<std::size_t>::max()};
    std::uint64_t boundary_id{};
    bool farfield{};
};

struct BackgroundInterface {
    std::uint64_t first{};
    std::vector<std::uint64_t> seconds;
    std::uint8_t first_face{};
};

struct BackgroundGridView {
    AABB domain;
    std::string kind;
    std::vector<std::uint64_t> stable_cell_ids;
    std::vector<AABB> cells;
    std::vector<std::array<bool, 6>> exterior_faces;
    std::vector<BackgroundInterface> interfaces;
};

struct EdgePoint {
    double parameter{};
    Vec3 point{};
};

[[nodiscard]] std::size_t face_axis(std::uint8_t local_face) noexcept {
    return static_cast<std::size_t>(local_face / 2U);
}

[[nodiscard]] Vec3 face_normal(std::uint8_t local_face) noexcept {
    Vec3 result{};
    const double sign = local_face % 2U == 0U ? -1.0 : 1.0;
    if (face_axis(local_face) == 0U) result.x = sign;
    else if (face_axis(local_face) == 1U) result.y = sign;
    else result.z = sign;
    return result;
}

[[nodiscard]] Point2 project(const Vec3& point, std::size_t axis) noexcept {
    if (axis == 0U) return {point.y, point.z};
    if (axis == 1U) return {point.z, point.x};
    return {point.x, point.y};
}

[[nodiscard]] Vec3 lift(const Point2& point, std::size_t axis,
                        double coordinate) noexcept {
    if (axis == 0U) return {coordinate, point.u, point.v};
    if (axis == 1U) return {point.v, coordinate, point.u};
    return {point.u, point.v, coordinate};
}

[[nodiscard]] double signed_area_2d(const std::vector<Point2>& polygon) noexcept {
    double sum = 0.0;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const auto& first = polygon[index];
        const auto& second = polygon[(index + 1U) % polygon.size()];
        sum += first.u * second.v - second.u * first.v;
    }
    return 0.5 * sum;
}

[[nodiscard]] double polygon_area(const std::vector<Vec3>& polygon) noexcept {
    if (polygon.size() < 3U) return 0.0;
    Vec3 area_vector{};
    const Vec3 origin = polygon.front();
    for (std::size_t index = 1; index + 1U < polygon.size(); ++index) {
        area_vector = area_vector +
                      cross(polygon[index] - origin,
                            polygon[index + 1U] - origin) * 0.5;
    }
    return norm(area_vector);
}

[[nodiscard]] double cell_face_area_scale(const AABB& bounds) noexcept {
    const Vec3 extent = bounds.extent();
    return std::max({extent.x * extent.y, extent.x * extent.z,
                     extent.y * extent.z});
}

[[nodiscard]] BackgroundGridView make_background_view(
    const UniformCartesianGrid& grid) {
    BackgroundGridView view{grid.domain(), "uniform_linear_id", {}, {}, {}, {}};
    view.stable_cell_ids.reserve(static_cast<std::size_t>(grid.cell_count()));
    view.cells.reserve(static_cast<std::size_t>(grid.cell_count()));
    view.exterior_faces.reserve(static_cast<std::size_t>(grid.cell_count()));
    view.interfaces.reserve(static_cast<std::size_t>(grid.cell_count()) * 3U);
    for (std::uint64_t background = 0; background < grid.cell_count(); ++background) {
        const CellKey key = grid.cell_key(background);
        view.stable_cell_ids.push_back(background);
        view.cells.push_back(grid.cell_bounds(key));
        view.exterior_faces.push_back({
            key.i == 0U, key.i + 1U == grid.nx(),
            key.j == 0U, key.j + 1U == grid.ny(),
            key.k == 0U, key.k + 1U == grid.nz()});
        if (key.i + 1U < grid.nx()) {
            CellKey neighbor = key;
            ++neighbor.i;
            view.interfaces.push_back(
                {background, {grid.linear_id(neighbor)}, 1U});
        }
        if (key.j + 1U < grid.ny()) {
            CellKey neighbor = key;
            ++neighbor.j;
            view.interfaces.push_back(
                {background, {grid.linear_id(neighbor)}, 3U});
        }
        if (key.k + 1U < grid.nz()) {
            CellKey neighbor = key;
            ++neighbor.k;
            view.interfaces.push_back(
                {background, {grid.linear_id(neighbor)}, 5U});
        }
    }
    return view;
}

[[nodiscard]] BackgroundGridView make_background_view(
    const LinearOctree& tree) {
    BackgroundGridView view{tree.domain(), "octree_node_code", {}, {}, {}, {}};
    view.stable_cell_ids.reserve(static_cast<std::size_t>(tree.leaf_count()));
    view.cells.reserve(static_cast<std::size_t>(tree.leaf_count()));
    view.exterior_faces.reserve(static_cast<std::size_t>(tree.leaf_count()));
    view.interfaces.reserve(static_cast<std::size_t>(tree.leaf_count()) * 3U);
    constexpr std::array<FaceDirection, 6> directions = {
        FaceDirection::negative_x, FaceDirection::positive_x,
        FaceDirection::negative_y, FaceDirection::positive_y,
        FaceDirection::negative_z, FaceDirection::positive_z};
    for (std::uint64_t leaf_id = 0; leaf_id < tree.leaf_count(); ++leaf_id) {
        const OctreeNodeCode code = tree.leaf_code(leaf_id);
        const AABB bounds = tree.cell_bounds(code);
        view.stable_cell_ids.push_back(code);
        view.cells.push_back(bounds);
        std::array<bool, 6> exterior{};
        for (std::size_t face = 0; face < directions.size(); ++face) {
            const auto neighbors = tree.face_neighbors(code, directions[face]);
            exterior[face] = neighbors.empty();
            if (neighbors.empty()) continue;
            const auto node = decode_octree_node(code);
            bool has_coarser_neighbor = false;
            bool has_finer_neighbor = false;
            std::vector<std::uint64_t> neighbor_ids;
            neighbor_ids.reserve(neighbors.size());
            for (const OctreeNodeCode neighbor_code : neighbors) {
                const auto neighbor_node = decode_octree_node(neighbor_code);
                has_coarser_neighbor =
                    has_coarser_neighbor || neighbor_node.level < node.level;
                has_finer_neighbor =
                    has_finer_neighbor || neighbor_node.level > node.level;
                const auto neighbor_id = tree.find_leaf(neighbor_code);
                if (!neighbor_id) {
                    throw std::logic_error("八叉树 OpenFOAM 邻居不在叶数组中");
                }
                neighbor_ids.push_back(*neighbor_id);
            }
            // 跨层接口始终由 coarse face 一次覆盖全部 fine faces；同层接口由较小
            // Morton 叶 ID 处理一次。
            if (has_coarser_neighbor) continue;
            if (!has_finer_neighbor && leaf_id > neighbor_ids.front()) continue;
            view.interfaces.push_back(
                {leaf_id, std::move(neighbor_ids),
                 static_cast<std::uint8_t>(face)});
        }
        view.exterior_faces.push_back(exterior);
    }
    return view;
}

void remove_short_edges(std::vector<Point2>& polygon, double tolerance) {
    std::vector<Point2> compact;
    compact.reserve(polygon.size());
    for (const auto point : polygon) {
        if (compact.empty() ||
            std::hypot(point.u - compact.back().u,
                       point.v - compact.back().v) > tolerance) {
            compact.push_back(point);
        }
    }
    if (compact.size() > 1U &&
        std::hypot(compact.front().u - compact.back().u,
                   compact.front().v - compact.back().v) <= tolerance) {
        compact.pop_back();
    }
    polygon = std::move(compact);
}

[[nodiscard]] std::vector<Vec3> intersect_convex_cartesian_polygons(
    const std::vector<Vec3>& first, const std::vector<Vec3>& second,
    std::uint8_t first_local_face, double length_tolerance) {
    if (first.size() < 3U || second.size() < 3U) return {};
    const std::size_t axis = face_axis(first_local_face);
    std::vector<Point2> subject;
    std::vector<Point2> clip;
    subject.reserve(first.size());
    clip.reserve(second.size());
    for (const auto point : first) subject.push_back(project(point, axis));
    for (const auto point : second) clip.push_back(project(point, axis));
    remove_short_edges(subject, length_tolerance);
    remove_short_edges(clip, length_tolerance);
    if (subject.size() < 3U || clip.size() < 3U) return {};
    if (signed_area_2d(clip) < 0.0) std::reverse(clip.begin(), clip.end());

    for (std::size_t edge = 0; edge < clip.size(); ++edge) {
        const Point2 a = clip[edge];
        const Point2 b = clip[(edge + 1U) % clip.size()];
        const auto side = [&](const Point2& point) {
            return (b.u - a.u) * (point.v - a.v) -
                   (b.v - a.v) * (point.u - a.u);
        };
        std::vector<Point2> output;
        if (subject.empty()) break;
        Point2 previous = subject.back();
        double previous_side = side(previous);
        for (const auto current : subject) {
            const double current_side = side(current);
            const bool previous_inside = previous_side >= -length_tolerance;
            const bool current_inside = current_side >= -length_tolerance;
            if (previous_inside != current_inside) {
                const double denominator = previous_side - current_side;
                if (denominator != 0.0) {
                    const double parameter = previous_side / denominator;
                    output.push_back(
                        {previous.u + (current.u - previous.u) * parameter,
                         previous.v + (current.v - previous.v) * parameter});
                }
            }
            if (current_inside) output.push_back(current);
            previous = current;
            previous_side = current_side;
        }
        remove_short_edges(output, length_tolerance);
        subject = std::move(output);
    }
    if (subject.size() < 3U || std::abs(signed_area_2d(subject)) <=
                                    length_tolerance * length_tolerance) {
        return {};
    }
    const double coordinate = axis == 0U ? first.front().x
                              : axis == 1U ? first.front().y
                                           : first.front().z;
    std::vector<Vec3> result;
    result.reserve(subject.size());
    for (const auto point : subject) result.push_back(lift(point, axis, coordinate));
    return result;
}

[[nodiscard]] std::vector<Vec3> intersect_convex_planar_polygons(
    const std::vector<Vec3>& first, const std::vector<Vec3>& second,
    const Vec3& normal, double length_tolerance) {
    if (first.size() < 3U || second.size() < 3U) return {};
    if (std::abs(dot(second.front() - first.front(), normal)) >
        8.0 * length_tolerance) {
        return {};
    }
    const std::array<Vec3, 3> axes = {
        Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0},
        Vec3{0.0, 0.0, 1.0}};
    const auto reference_axis = *std::min_element(
        axes.begin(), axes.end(), [&](const Vec3& lhs, const Vec3& rhs) {
            return std::abs(dot(lhs, normal)) < std::abs(dot(rhs, normal));
        });
    const Vec3 tangent_u = cross(reference_axis, normal) /
                           norm(cross(reference_axis, normal));
    const Vec3 tangent_v = cross(normal, tangent_u);
    const Vec3 origin = first.front();
    const auto project_local = [&](const Vec3& point) {
        const Vec3 delta = point - origin;
        return Point2{dot(delta, tangent_u), dot(delta, tangent_v)};
    };
    std::vector<Point2> subject;
    std::vector<Point2> clip;
    for (const auto point : first) subject.push_back(project_local(point));
    for (const auto point : second) clip.push_back(project_local(point));
    remove_short_edges(subject, length_tolerance);
    remove_short_edges(clip, length_tolerance);
    if (subject.size() < 3U || clip.size() < 3U) return {};
    if (signed_area_2d(clip) < 0.0) std::reverse(clip.begin(), clip.end());
    for (std::size_t edge = 0; edge < clip.size(); ++edge) {
        const Point2 a = clip[edge];
        const Point2 b = clip[(edge + 1U) % clip.size()];
        const auto side = [&](const Point2& point) {
            return (b.u - a.u) * (point.v - a.v) -
                   (b.v - a.v) * (point.u - a.u);
        };
        std::vector<Point2> output;
        if (subject.empty()) break;
        Point2 previous = subject.back();
        double previous_side = side(previous);
        for (const auto current : subject) {
            const double current_side = side(current);
            const bool previous_inside = previous_side >= -length_tolerance;
            const bool current_inside = current_side >= -length_tolerance;
            if (previous_inside != current_inside) {
                const double denominator = previous_side - current_side;
                if (denominator != 0.0) {
                    const double parameter = previous_side / denominator;
                    output.push_back(
                        {previous.u + (current.u - previous.u) * parameter,
                         previous.v + (current.v - previous.v) * parameter});
                }
            }
            if (current_inside) output.push_back(current);
            previous = current;
            previous_side = current_side;
        }
        remove_short_edges(output, length_tolerance);
        subject = std::move(output);
    }
    if (subject.size() < 3U || std::abs(signed_area_2d(subject)) <=
                                    length_tolerance * length_tolerance) {
        return {};
    }
    std::vector<Vec3> result;
    result.reserve(subject.size());
    for (const auto point : subject) {
        result.push_back(origin + tangent_u * point.u + tangent_v * point.v);
    }
    return result;
}

void orient_face(std::vector<Vec3>& vertices, const Vec3& desired_normal) {
    if (vertices.size() < 3U) return;
    Vec3 area_vector{};
    const Vec3 origin = vertices.front();
    for (std::size_t index = 1; index + 1U < vertices.size(); ++index) {
        area_vector = area_vector +
                      cross(vertices[index] - origin,
                            vertices[index + 1U] - origin);
    }
    if (dot(area_vector, desired_normal) < 0.0) {
        std::reverse(vertices.begin(), vertices.end());
    }
}

[[nodiscard]] Vec3 face_area_vector(
    const std::vector<Vec3>& vertices) noexcept {
    Vec3 result{};
    if (vertices.size() < 3U) return result;
    const Vec3 origin = vertices.front();
    for (std::size_t index = 1; index + 1U < vertices.size(); ++index) {
        result = result +
                 cross(vertices[index] - origin,
                       vertices[index + 1U] - origin) * 0.5;
    }
    return result;
}

[[nodiscard]] std::vector<Vec3> convex_hull_on_plane(
    const std::vector<Vec3>& points, const Vec3& normal,
    double length_tolerance) {
    const std::array<Vec3, 3> axes = {
        Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0},
        Vec3{0.0, 0.0, 1.0}};
    const auto reference_axis = *std::min_element(
        axes.begin(), axes.end(), [&](const Vec3& lhs, const Vec3& rhs) {
            return std::abs(dot(lhs, normal)) < std::abs(dot(rhs, normal));
        });
    const Vec3 tangent_u = cross(reference_axis, normal) /
                           norm(cross(reference_axis, normal));
    const Vec3 tangent_v = cross(normal, tangent_u);
    const Vec3 origin = points.front();
    std::vector<Point2> projected;
    projected.reserve(points.size());
    for (const auto point : points) {
        const Vec3 delta = point - origin;
        projected.push_back({dot(delta, tangent_u), dot(delta, tangent_v)});
    }
    std::sort(projected.begin(), projected.end(), [](const Point2& lhs,
                                                      const Point2& rhs) {
        if (lhs.u != rhs.u) return lhs.u < rhs.u;
        return lhs.v < rhs.v;
    });
    std::vector<Point2> unique;
    for (const auto point : projected) {
        if (unique.empty() ||
            std::hypot(point.u - unique.back().u,
                       point.v - unique.back().v) > length_tolerance) {
            unique.push_back(point);
        }
    }
    if (unique.size() < 3U) return {};
    double scale = 0.0;
    for (const auto point : unique) scale = std::max(scale, std::hypot(point.u, point.v));
    const double turn_tolerance = length_tolerance * std::max(scale, length_tolerance);
    const auto turn = [](const Point2& a, const Point2& b, const Point2& c) {
        return (b.u - a.u) * (c.v - a.v) -
               (b.v - a.v) * (c.u - a.u);
    };
    std::vector<Point2> hull;
    hull.reserve(unique.size() * 2U);
    for (const auto point : unique) {
        while (hull.size() >= 2U &&
               turn(hull[hull.size() - 2U], hull.back(), point) <=
                   turn_tolerance) {
            hull.pop_back();
        }
        hull.push_back(point);
    }
    const std::size_t lower_size = hull.size();
    for (std::size_t index = unique.size() - 1U; index-- > 0U;) {
        const auto point = unique[index];
        while (hull.size() > lower_size &&
               turn(hull[hull.size() - 2U], hull.back(), point) <=
                   turn_tolerance) {
            hull.pop_back();
        }
        hull.push_back(point);
    }
    if (!hull.empty()) hull.pop_back();
    if (hull.size() < 3U) return {};
    std::vector<Vec3> result;
    result.reserve(hull.size());
    for (const auto point : hull) {
        result.push_back(origin + tangent_u * point.u + tangent_v * point.v);
    }
    orient_face(result, normal);
    return result;
}

void merge_coplanar_boundary_faces(std::vector<OutputFace>& boundary_faces,
                                   double length_tolerance,
                                   double area_tolerance) {
    std::vector<OutputFace> merged;
    merged.reserve(boundary_faces.size());
    const double normal_tolerance =
        2048.0 * std::numeric_limits<double>::epsilon();
    for (auto& face : boundary_faces) {
        if (face.farfield) {
            merged.push_back(std::move(face));
            continue;
        }
        const Vec3 area_vector = face_area_vector(face.vertices);
        const Vec3 normal = area_vector / norm(area_vector);
        bool consumed = false;
        for (auto& candidate : merged) {
            if (candidate.farfield || candidate.owner != face.owner ||
                candidate.boundary_id != face.boundary_id) {
                continue;
            }
            const Vec3 candidate_area_vector =
                face_area_vector(candidate.vertices);
            const Vec3 candidate_normal =
                candidate_area_vector / norm(candidate_area_vector);
            if (norm(candidate_normal - normal) > normal_tolerance ||
                std::abs(dot(face.vertices.front() - candidate.vertices.front(),
                             normal)) > 4.0 * length_tolerance) {
                continue;
            }
            std::vector<Vec3> points = candidate.vertices;
            points.insert(points.end(), face.vertices.begin(), face.vertices.end());
            auto hull = convex_hull_on_plane(points, normal, length_tolerance);
            const double hull_area = polygon_area(hull);
            const double separate_area = polygon_area(candidate.vertices) +
                                         polygon_area(face.vertices);
            const double merge_tolerance = std::max(
                area_tolerance,
                8192.0 * std::numeric_limits<double>::epsilon() *
                    std::max(hull_area, separate_area));
            if (hull.size() >= 3U &&
                std::abs(hull_area - separate_area) <= merge_tolerance) {
                candidate.vertices = std::move(hull);
                consumed = true;
                break;
            }
        }
        if (!consumed) merged.push_back(std::move(face));
    }
    boundary_faces = std::move(merged);
}

void split_all_face_edges(std::vector<OutputFace>& internal_faces,
                          std::vector<OutputFace>& boundary_faces,
                          double length_tolerance) {
    std::vector<Vec3> points;
    for (const auto& face : internal_faces) {
        points.insert(points.end(), face.vertices.begin(), face.vertices.end());
    }
    for (const auto& face : boundary_faces) {
        points.insert(points.end(), face.vertices.begin(), face.vertices.end());
    }
    const auto split = [&](OutputFace& face) {
        std::vector<Vec3> result;
        for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
            const Vec3 first = face.vertices[edge];
            const Vec3 second = face.vertices[(edge + 1U) % face.vertices.size()];
            const Vec3 direction = second - first;
            const double squared_length = dot(direction, direction);
            if (!(squared_length > 0.0)) continue;
            std::vector<EdgePoint> edge_points{{0.0, first}};
            for (const auto point : points) {
                const double parameter = dot(point - first, direction) /
                                         squared_length;
                if (parameter <= 0.0 || parameter >= 1.0) continue;
                const Vec3 projection = first + direction * parameter;
                if (norm(point - projection) <= length_tolerance) {
                    edge_points.push_back({parameter, projection});
                }
            }
            std::sort(edge_points.begin(), edge_points.end(),
                      [](const auto& lhs, const auto& rhs) {
                          return lhs.parameter < rhs.parameter;
                      });
            for (const auto& edge_point : edge_points) {
                const Vec3 point = edge_point.point;
                if (result.empty() ||
                    norm(result.back() - point) > length_tolerance) {
                    result.push_back(point);
                }
            }
        }
        face.vertices = std::move(result);
    };
    for (auto& face : internal_faces) split(face);
    for (auto& face : boundary_faces) split(face);
}

[[nodiscard]] std::vector<Vec3> face_vertices(
    const FluidPolyhedronPiece& piece, std::size_t face) {
    std::vector<Vec3> result;
    for (const auto id : piece.polyhedron.faces[face].vertex_indices) {
        result.push_back(piece.polyhedron.vertices[id]);
    }
    return result;
}

[[nodiscard]] std::string foam_name(std::string name) {
    for (char& character : name) {
        const bool valid = (character >= 'a' && character <= 'z') ||
                           (character >= 'A' && character <= 'Z') ||
                           (character >= '0' && character <= '9') ||
                           character == '_';
        if (!valid) character = '_';
    }
    if (name.empty() || (name.front() >= '0' && name.front() <= '9')) {
        name = "patch_" + name;
    }
    return name;
}

void write_header(std::ofstream& output, const std::string& klass,
                  const std::string& object) {
    output << "FoamFile\n{\n"
           << "    version 2.0;\n"
           << "    format ascii;\n"
           << "    class " << klass << ";\n"
           << "    location \"constant/polyMesh\";\n"
           << "    object " << object << ";\n"
           << "}\n\n";
}

class PointWelder {
  public:
    PointWelder(Vec3 origin, double tolerance)
        : origin_(origin), tolerance_(tolerance) {}

    [[nodiscard]] std::size_t insert(const Vec3& point) {
        const auto key = quantize(point);
        for (std::int64_t dx = -1; dx <= 1; ++dx) {
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    const Key neighbor{std::get<0>(key) + dx,
                                       std::get<1>(key) + dy,
                                       std::get<2>(key) + dz};
                    const auto found = buckets_.find(neighbor);
                    if (found == buckets_.end()) continue;
                    for (const auto id : found->second) {
                        if (norm(points_[id] - point) <= tolerance_) return id;
                    }
                }
            }
        }
        const std::size_t id = points_.size();
        points_.push_back(point);
        buckets_[key].push_back(id);
        return id;
    }

    [[nodiscard]] const std::vector<Vec3>& points() const noexcept {
        return points_;
    }

  private:
    using Key = std::tuple<std::int64_t, std::int64_t, std::int64_t>;
    [[nodiscard]] Key quantize(const Vec3& point) const {
        const Vec3 local = (point - origin_) / tolerance_;
        return {static_cast<std::int64_t>(std::llround(local.x)),
                static_cast<std::int64_t>(std::llround(local.y)),
                static_cast<std::int64_t>(std::llround(local.z))};
    }
    Vec3 origin_{};
    double tolerance_{};
    std::vector<Vec3> points_;
    std::map<Key, std::vector<std::size_t>> buckets_;
};

} // 匿名命名空间

static void write_openfoam_poly_mesh_impl(
    const std::filesystem::path& case_directory,
    const BackgroundGridView& grid, const ConvexCutCellMesh& mesh,
    const std::vector<std::pair<std::uint64_t, std::string>>& boundary_names,
    double length_tolerance) {
    if (length_tolerance < 0.0 || !std::isfinite(length_tolerance)) {
        throw std::invalid_argument("OpenFOAM 点合并容差必须是非负有限数");
    }
    const double scale = norm(grid.domain.extent());
    length_tolerance = std::max(
        length_tolerance,
        1024.0 * std::numeric_limits<double>::epsilon() * scale);
    const double area_tolerance = length_tolerance * length_tolerance;

    std::vector<ComponentCell> cells;
    std::vector<const FluidCellGeometry*> fluid_by_background(
        grid.cells.size(), nullptr);
    for (const auto& cell : mesh.fluid_cells) {
        if (cell.background_cell_id >= grid.cells.size()) {
            throw std::runtime_error("OpenFOAM 流体单元 background ID 越界");
        }
        fluid_by_background[static_cast<std::size_t>(cell.background_cell_id)] = &cell;
        if (cell.cut && (cell.fluid_component_count == 0U ||
                         cell.fluid_polyhedron_pieces.empty())) {
            throw std::runtime_error("OpenFOAM 输出遇到未解析的流体分量");
        }
        if (!cell.cut) {
            cells.push_back({cell.background_cell_id, 0U, 0U, nullptr, true});
        } else {
            for (std::size_t piece_id = 0;
                 piece_id < cell.fluid_polyhedron_pieces.size(); ++piece_id) {
                const auto& piece = cell.fluid_polyhedron_pieces[piece_id];
                if (piece.component_id >= cell.fluid_component_count) {
                    throw std::runtime_error("Cut-cell 凸片 component ID 越界");
                }
                cells.push_back({cell.background_cell_id, piece.component_id,
                                 piece_id, &piece, false});
            }
        }
    }

    std::vector<std::array<std::vector<CartesianPatch>, 6>> cartesian(
        grid.cells.size());
    std::vector<std::vector<PartitionPatch>> partitions(
        grid.cells.size());
    for (std::size_t cell_id = 0; cell_id < cells.size(); ++cell_id) {
        const auto& component = cells[cell_id];
        if (component.full_cartesian) {
            const auto polyhedron = make_box_polyhedron(
                grid.cells[static_cast<std::size_t>(
                    component.background_cell_id)]);
            for (std::size_t face = 0; face < polyhedron.faces.size(); ++face) {
                std::vector<Vec3> vertices;
                for (const auto vertex : polyhedron.faces[face].vertex_indices) {
                    vertices.push_back(polyhedron.vertices[vertex]);
                }
                const auto local_face = static_cast<std::uint8_t>(face);
                orient_face(vertices, face_normal(local_face));
                cartesian[static_cast<std::size_t>(component.background_cell_id)]
                         [face]
                    .push_back({cell_id, local_face, vertices,
                                polygon_area(vertices), 0.0});
            }
            continue;
        }
        const auto& piece = *component.piece;
        for (std::size_t face = 0; face < piece.polyhedron.faces.size(); ++face) {
                const auto& record = piece.polyhedron.faces[face];
                auto vertices = face_vertices(piece, face);
                if (record.kind != PolyhedronFaceKind::cartesian) {
                    const auto normal = piece.geometry.faces[face].outward_normal;
                    orient_face(vertices, normal);
                    const double area = polygon_area(vertices);
                    if (area > area_tolerance) {
                        partitions[static_cast<std::size_t>(
                            component.background_cell_id)]
                            .push_back({cell_id, normal, std::move(vertices),
                                        area, 0.0});
                    }
                    continue;
                }
                if (record.source_id >= 6U) {
                    throw std::runtime_error("Cut-cell Cartesian 面编号越界");
                }
                const auto local_face = static_cast<std::uint8_t>(record.source_id);
                orient_face(vertices, face_normal(local_face));
                const double area = polygon_area(vertices);
                if (area > area_tolerance) {
                    cartesian[static_cast<std::size_t>(component.background_cell_id)]
                             [local_face]
                        .push_back({cell_id, local_face, std::move(vertices), area});
                }
        }
    }

    std::vector<OutputFace> internal_faces;
    std::vector<OutputFace> boundary_faces;
    const auto add_interface = [&](const BackgroundInterface& interface) {
        const std::uint64_t first_background = interface.first;
        const std::uint8_t first_face = interface.first_face;
        const std::uint8_t second_face = opposite_face[first_face];
        auto& first = cartesian[static_cast<std::size_t>(first_background)]
                               [first_face];
        std::vector<double> first_covered(first.size(), 0.0);
        double covered_area = 0.0;
        double expected_second_sum = 0.0;
        double face_scale = cell_face_area_scale(
            grid.cells[static_cast<std::size_t>(first_background)]);
        for (const std::uint64_t second_background : interface.seconds) {
            auto& second = cartesian[static_cast<std::size_t>(second_background)]
                                    [second_face];
            std::vector<double> second_covered(second.size(), 0.0);
            for (std::size_t i = 0; i < first.size(); ++i) {
                for (std::size_t j = 0; j < second.size(); ++j) {
                    auto overlap = intersect_convex_cartesian_polygons(
                        first[i].vertices, second[j].vertices, first_face,
                        length_tolerance);
                    const double area = polygon_area(overlap);
                    if (area <= area_tolerance) continue;
                    std::size_t owner = first[i].owner;
                    std::size_t neighbor = second[j].owner;
                    Vec3 normal = face_normal(first_face);
                    if (owner > neighbor) {
                        std::swap(owner, neighbor);
                        normal = normal * -1.0;
                    }
                    orient_face(overlap, normal);
                    internal_faces.push_back(
                        {std::move(overlap), owner, neighbor, 0U, false});
                    first_covered[i] += area;
                    second_covered[j] += area;
                    first[i].covered_area += area;
                    second[j].covered_area += area;
                }
            }
            const auto* second_fluid = fluid_by_background[
                static_cast<std::size_t>(second_background)];
            const double expected_second =
                second_fluid == nullptr
                    ? 0.0
                    : second_fluid->cartesian_faces[second_face].area;
            const double second_total = std::accumulate(
                second_covered.begin(), second_covered.end(), 0.0);
            expected_second_sum += expected_second;
            covered_area += second_total;
            face_scale = std::max(
                face_scale,
                cell_face_area_scale(
                    grid.cells[static_cast<std::size_t>(second_background)]));
            const double local_tolerance = std::max(
                area_tolerance,
                8192.0 * std::numeric_limits<double>::epsilon() * face_scale);
            if (std::abs(second_total - expected_second) > local_tolerance) {
                throw std::runtime_error(
                    "OpenFOAM 粗细接口 fine 面覆盖不守恒：first=" +
                    std::to_string(first_background) +
                    " second=" + std::to_string(second_background) +
                    " face=" + std::to_string(first_face) +
                    " covered=" + std::to_string(second_total) +
                    " expected=" + std::to_string(expected_second));
            }
        }
        const double coverage_tolerance = std::max(
            area_tolerance,
            8192.0 * std::numeric_limits<double>::epsilon() * face_scale);
        const auto* first_fluid = fluid_by_background[
            static_cast<std::size_t>(first_background)];
        const double expected_first =
            first_fluid == nullptr
                ? 0.0
                : first_fluid->cartesian_faces[first_face].area;
        if (std::abs(covered_area - expected_first) > coverage_tolerance ||
            std::abs(covered_area - expected_second_sum) > coverage_tolerance) {
            throw std::runtime_error(
                "OpenFOAM 内部面公共细分与控制体开口面积不守恒：first=" +
                std::to_string(first_background) +
                " face=" + std::to_string(first_face) +
                " covered=" + std::to_string(covered_area) +
                " expectedFirst=" + std::to_string(expected_first) +
                " expectedSecond=" + std::to_string(expected_second_sum));
        }
    };

    for (const auto& interface : grid.interfaces) add_interface(interface);

    // 同一背景单元内的凸片以 arrangement 面成对连接。
    for (std::uint64_t background = 0; background < grid.cells.size(); ++background) {
        auto& patches = partitions[static_cast<std::size_t>(background)];
        for (std::size_t first = 0; first < patches.size(); ++first) {
            for (std::size_t second = first + 1U; second < patches.size(); ++second) {
                if (patches[first].owner == patches[second].owner ||
                    dot(patches[first].outward_normal,
                        patches[second].outward_normal) > -1.0 + 1.0e-10) {
                    continue;
                }
                auto overlap = intersect_convex_planar_polygons(
                    patches[first].vertices, patches[second].vertices,
                    patches[first].outward_normal, length_tolerance);
                const double area = polygon_area(overlap);
                if (area <= area_tolerance) continue;
                std::size_t owner = patches[first].owner;
                std::size_t neighbor = patches[second].owner;
                Vec3 normal = patches[first].outward_normal;
                if (owner > neighbor) {
                    std::swap(owner, neighbor);
                    normal = normal * -1.0;
                }
                orient_face(overlap, normal);
                internal_faces.push_back(
                    {std::move(overlap), owner, neighbor, 0U, false});
                patches[first].covered_area += area;
                patches[second].covered_area += area;
            }
        }
    }

    // 计算域外边界的 Cartesian 面。
    for (std::uint64_t background = 0; background < grid.cells.size(); ++background) {
        const auto& exterior =
            grid.exterior_faces[static_cast<std::size_t>(background)];
        for (std::size_t face = 0; face < 6U; ++face) {
            if (!exterior[face]) continue;
            for (auto& patch :
                 cartesian[static_cast<std::size_t>(background)][face]) {
                auto vertices = patch.vertices;
                orient_face(vertices, face_normal(static_cast<std::uint8_t>(face)));
                boundary_faces.push_back(
                    {std::move(vertices), patch.owner,
                     std::numeric_limits<std::size_t>::max(), 0U, true});
                patch.covered_area += patch.area;
            }
        }
    }

    // 嵌入边界与凸片的外向 arrangement 面取交，得到每个
    // OpenFOAM 凸单元所拥有的 wall 子面。表面与 Cartesian 面共面时
    // arrangement 面可以退化掉，再与 Cartesian 凸片取交。
    for (const auto& fluid_cell : mesh.fluid_cells) {
        auto& local_partitions = partitions[
            static_cast<std::size_t>(fluid_cell.background_cell_id)];
        auto& local_cartesian = cartesian[
            static_cast<std::size_t>(fluid_cell.background_cell_id)];
        const double face_scale = cell_face_area_scale(
            grid.cells[static_cast<std::size_t>(fluid_cell.background_cell_id)]);
        const double coverage_tolerance = std::max(
            area_tolerance,
            8192.0 * std::numeric_limits<double>::epsilon() * face_scale);
        for (const auto& embedded : fluid_cell.embedded_boundary_faces) {
            double covered = 0.0;
            for (auto& patch : local_partitions) {
                if (dot(patch.outward_normal, embedded.outward_normal) <
                    1.0 - 1.0e-10) {
                    continue;
                }
                auto overlap = intersect_convex_planar_polygons(
                    patch.vertices, embedded.vertices,
                    embedded.outward_normal, length_tolerance);
                const double area = polygon_area(overlap);
                if (area <= area_tolerance) continue;
                orient_face(overlap, embedded.outward_normal);
                boundary_faces.push_back(
                    {std::move(overlap), patch.owner,
                     std::numeric_limits<std::size_t>::max(),
                     embedded.boundary_id, false});
                patch.covered_area += area;
                covered += area;
            }
            if (covered < embedded.area - coverage_tolerance) {
                for (auto& face_patches : local_cartesian) {
                    for (auto& patch : face_patches) {
                        const Vec3 normal = face_normal(patch.local_face);
                        if (dot(normal, embedded.outward_normal) <
                            1.0 - 1.0e-10) {
                            continue;
                        }
                        auto overlap = intersect_convex_planar_polygons(
                            patch.vertices, embedded.vertices,
                            embedded.outward_normal, length_tolerance);
                        const double area = polygon_area(overlap);
                        if (area <= area_tolerance) continue;
                        orient_face(overlap, embedded.outward_normal);
                        boundary_faces.push_back(
                            {std::move(overlap), patch.owner,
                             std::numeric_limits<std::size_t>::max(),
                             embedded.boundary_id, false});
                        patch.covered_area += area;
                        covered += area;
                    }
                }
            }
            if (std::abs(covered - embedded.area) > coverage_tolerance) {
                throw std::runtime_error(
                    "OpenFOAM 嵌入边界凸片细分面积不守恒：background=" +
                    std::to_string(fluid_cell.background_cell_id) +
                    " expected=" + std::to_string(embedded.area) +
                    " covered=" + std::to_string(covered));
            }
        }
    }

    for (std::uint64_t background = 0; background < grid.cells.size(); ++background) {
        const double coverage_tolerance = std::max(
            area_tolerance,
            16384.0 * std::numeric_limits<double>::epsilon() *
                cell_face_area_scale(
                    grid.cells[static_cast<std::size_t>(background)]));
        for (const auto& face_patches :
             cartesian[static_cast<std::size_t>(background)]) {
            for (const auto& patch : face_patches) {
                if (std::abs(patch.covered_area - patch.area) >
                    coverage_tolerance) {
                    throw std::runtime_error(
                        "OpenFOAM Cartesian 凸片面未被内部面或边界完整覆盖：background=" +
                        std::to_string(background) +
                        " area=" + std::to_string(patch.area) +
                        " covered=" + std::to_string(patch.covered_area));
                }
            }
        }
        for (const auto& patch : partitions[static_cast<std::size_t>(background)]) {
            if (std::abs(patch.covered_area - patch.area) >
                coverage_tolerance) {
                throw std::runtime_error(
                    "OpenFOAM arrangement 凸片面未被邻居或 wall 完整覆盖：background=" +
                    std::to_string(background) +
                    " area=" + std::to_string(patch.area) +
                    " covered=" + std::to_string(patch.covered_area));
            }
        }
    }

    merge_coplanar_boundary_faces(boundary_faces, length_tolerance,
                                  area_tolerance);
    split_all_face_edges(internal_faces, boundary_faces, length_tolerance);
    std::sort(internal_faces.begin(), internal_faces.end(),
              [](const OutputFace& first, const OutputFace& second) {
                  if (first.owner != second.owner) return first.owner < second.owner;
                  return first.neighbor < second.neighbor;
              });

    std::sort(boundary_faces.begin(), boundary_faces.end(),
              [](const OutputFace& first, const OutputFace& second) {
                  if (first.farfield != second.farfield) return first.farfield;
                  if (first.boundary_id != second.boundary_id) {
                      return first.boundary_id < second.boundary_id;
                  }
                  return first.owner < second.owner;
              });

    PointWelder welder(grid.domain.minimum(), length_tolerance);
    std::vector<std::vector<std::size_t>> face_point_ids;
    face_point_ids.reserve(internal_faces.size() + boundary_faces.size());
    const auto append_face = [&](const OutputFace& face) {
        std::vector<std::size_t> ids;
        ids.reserve(face.vertices.size());
        for (const auto vertex : face.vertices) {
            const std::size_t id = welder.insert(vertex);
            if (ids.empty() || ids.back() != id) ids.push_back(id);
        }
        if (ids.size() > 1U && ids.front() == ids.back()) ids.pop_back();
        if (ids.size() < 3U) {
            throw std::runtime_error("OpenFOAM 面在全局点合并后少于 3 个顶点");
        }
        face_point_ids.push_back(std::move(ids));
    };
    for (const auto& face : internal_faces) append_face(face);
    for (const auto& face : boundary_faces) append_face(face);

    const std::filesystem::path poly_mesh =
        case_directory / "constant" / "polyMesh";
    std::filesystem::create_directories(poly_mesh);
    std::filesystem::create_directories(case_directory / "system");
    {
        std::ofstream output(poly_mesh / "cartmeshCellMapping.json",
                             std::ios::trunc);
        if (!output) {
            throw std::runtime_error("无法写入 OpenFOAM cell mapping");
        }
        output << "{\n  \"schema\": \"cartmesh-openfoam-cell-mapping-v1\",\n"
               << "  \"backgroundStableIdKind\": \"" << grid.kind
               << "\",\n  \"solverCellCount\": " << cells.size()
               << ",\n  \"cells\": [\n";
        for (std::size_t solver_cell_id = 0;
             solver_cell_id < cells.size(); ++solver_cell_id) {
            const auto& cell = cells[solver_cell_id];
            if (solver_cell_id != 0U) output << ",\n";
            output << "    {\"solverCellId\":" << solver_cell_id
                   << ",\"backgroundCellId\":" << cell.background_cell_id
                   << ",\"backgroundStableId\":\""
                   << grid.stable_cell_ids[static_cast<std::size_t>(
                          cell.background_cell_id)]
                   << "\",\"componentId\":" << cell.component_id
                   << ",\"localPieceId\":" << cell.local_piece_id
                   << ",\"fullCartesian\":"
                   << (cell.full_cartesian ? "true" : "false") << '}';
        }
        output << "\n  ]\n}\n";
    }
    {
        std::ofstream output(case_directory / "system" / "controlDict",
                             std::ios::trunc);
        if (!output) {
            throw std::runtime_error("无法写入 OpenFOAM controlDict");
        }
        output << "FoamFile\n{\n"
               << "    version 2.0;\n"
               << "    format ascii;\n"
               << "    class dictionary;\n"
               << "    object controlDict;\n"
               << "}\n\n"
               << "application checkMesh;\n"
               << "startFrom startTime;\n"
               << "startTime 0;\n"
               << "stopAt endTime;\n"
               << "endTime 1;\n"
               << "deltaT 1;\n"
               << "writeControl timeStep;\n"
               << "writeInterval 1;\n"
               << "runTimeModifiable false;\n";
    }
    {
        std::ofstream output(case_directory / "system" / "fvSchemes",
                             std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 OpenFOAM fvSchemes");
        output << "FoamFile\n{\n"
               << "    version 2.0;\n    format ascii;\n"
               << "    class dictionary;\n    object fvSchemes;\n}\n\n"
               << "ddtSchemes { default Euler; }\n"
               << "gradSchemes { default Gauss linear; }\n"
               << "divSchemes { default none; }\n"
               << "laplacianSchemes { default Gauss linear corrected; }\n"
               << "interpolationSchemes { default linear; }\n"
               << "snGradSchemes { default corrected; }\n";
    }
    {
        std::ofstream output(case_directory / "system" / "fvSolution",
                             std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 OpenFOAM fvSolution");
        output << "FoamFile\n{\n"
               << "    version 2.0;\n    format ascii;\n"
               << "    class dictionary;\n    object fvSolution;\n}\n\n"
               << "solvers {}\n";
    }
    {
        std::ofstream output(poly_mesh / "points", std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 OpenFOAM points");
        output << std::setprecision(17);
        write_header(output, "vectorField", "points");
        output << welder.points().size() << "\n(\n";
        for (const auto point : welder.points()) {
            output << '(' << point.x << ' ' << point.y << ' ' << point.z << ")\n";
        }
        output << ")\n";
    }
    {
        std::ofstream output(poly_mesh / "faces", std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 OpenFOAM faces");
        write_header(output, "faceList", "faces");
        output << face_point_ids.size() << "\n(\n";
        for (const auto& face : face_point_ids) {
            output << face.size() << '(';
            for (std::size_t index = 0; index < face.size(); ++index) {
                if (index != 0U) output << ' ';
                output << face[index];
            }
            output << ")\n";
        }
        output << ")\n";
    }
    {
        std::ofstream output(poly_mesh / "owner", std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 OpenFOAM owner");
        write_header(output, "labelList", "owner");
        output << face_point_ids.size() << "\n(\n";
        for (const auto& face : internal_faces) output << face.owner << '\n';
        for (const auto& face : boundary_faces) output << face.owner << '\n';
        output << ")\n";
    }
    {
        std::ofstream output(poly_mesh / "neighbour", std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 OpenFOAM neighbour");
        write_header(output, "labelList", "neighbour");
        output << internal_faces.size() << "\n(\n";
        for (const auto& face : internal_faces) output << face.neighbor << '\n';
        output << ")\n";
    }
    {
        std::map<std::uint64_t, std::string> names;
        for (const auto& [id, name] : boundary_names) names[id] = foam_name(name);
        struct Patch { std::string name; std::string type; std::size_t count{}; };
        std::vector<Patch> patches;
        for (const auto& face : boundary_faces) {
            const std::string name = face.farfield
                                         ? "farfield"
                                         : names.contains(face.boundary_id)
                                               ? names[face.boundary_id]
                                               : "boundary_" +
                                                     std::to_string(face.boundary_id);
            const std::string type = face.farfield ? "patch" : "wall";
            if (patches.empty() || patches.back().name != name) {
                patches.push_back({name, type, 0U});
            }
            ++patches.back().count;
        }
        std::ofstream output(poly_mesh / "boundary", std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 OpenFOAM boundary");
        write_header(output, "polyBoundaryMesh", "boundary");
        output << patches.size() << "\n(\n";
        std::size_t start = internal_faces.size();
        for (const auto& patch : patches) {
            output << patch.name << "\n{\n"
                   << "    type " << patch.type << ";\n"
                   << "    nFaces " << patch.count << ";\n"
                   << "    startFace " << start << ";\n"
                   << "}\n";
            start += patch.count;
        }
        output << ")\n";
    }
}

void write_openfoam_poly_mesh(
    const std::filesystem::path& case_directory,
    const UniformCartesianGrid& grid, const ConvexCutCellMesh& mesh,
    const std::vector<std::pair<std::uint64_t, std::string>>& boundary_names,
    double length_tolerance) {
    write_openfoam_poly_mesh_impl(case_directory, make_background_view(grid),
                                  mesh, boundary_names, length_tolerance);
}

void write_openfoam_poly_mesh(
    const std::filesystem::path& case_directory,
    const LinearOctree& tree, const ConvexCutCellMesh& mesh,
    const std::vector<std::pair<std::uint64_t, std::string>>& boundary_names,
    double length_tolerance) {
    write_openfoam_poly_mesh_impl(case_directory, make_background_view(tree),
                                  mesh, boundary_names, length_tolerance);
}

} // 命名空间 cartmesh
