#include "cartmesh/cutcell/ConvexPolyhedron.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace cartmesh {
namespace {

struct WorkingFace {
    std::vector<Vec3> vertices;
    PolyhedronFaceKind kind{PolyhedronFaceKind::cartesian};
    std::uint64_t source_id{};
};

struct EdgeUse {
    std::uint32_t count{};
    int orientation_balance{};
};

[[nodiscard]] double coordinate_ulp(double value) noexcept {
    const double magnitude = std::abs(value);
    return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
}

[[nodiscard]] double derived_length_tolerance(const ConvexPolyhedron& polyhedron,
                                              const OrientedHalfSpace& half_space,
                                              double requested) noexcept {
    double maximum_ulp = std::max({coordinate_ulp(half_space.point().x),
                                   coordinate_ulp(half_space.point().y),
                                   coordinate_ulp(half_space.point().z)});
    Vec3 minimum = half_space.point();
    Vec3 maximum = half_space.point();
    for (const auto& vertex : polyhedron.vertices) {
        maximum_ulp = std::max(
            {maximum_ulp, coordinate_ulp(vertex.x), coordinate_ulp(vertex.y),
             coordinate_ulp(vertex.z)});
        minimum.x = std::min(minimum.x, vertex.x);
        minimum.y = std::min(minimum.y, vertex.y);
        minimum.z = std::min(minimum.z, vertex.z);
        maximum.x = std::max(maximum.x, vertex.x);
        maximum.y = std::max(maximum.y, vertex.y);
        maximum.z = std::max(maximum.z, vertex.z);
    }
    const double local_scale = std::max(norm(maximum - minimum),
                                        std::numeric_limits<double>::min());
    return std::max({requested, 4.0 * maximum_ulp,
                     64.0 * std::numeric_limits<double>::epsilon() * local_scale});
}

[[nodiscard]] bool near(const Vec3& first, const Vec3& second, double tolerance) noexcept {
    return norm(first - second) <= tolerance;
}

void append_unique(std::vector<Vec3>& vertices, const Vec3& vertex, double tolerance) {
    if (vertices.empty() || !near(vertices.back(), vertex, tolerance)) {
        vertices.push_back(vertex);
    }
}

void close_polygon(std::vector<Vec3>& vertices, double tolerance) {
    if (vertices.size() > 1 && near(vertices.front(), vertices.back(), tolerance)) {
        vertices.pop_back();
    }
}

[[nodiscard]] Vec3 polygon_area_vector(const std::vector<Vec3>& vertices) noexcept {
    if (vertices.size() < 3) {
        return {};
    }
    Vec3 result{};
    const Vec3 reference = vertices.front();
    for (std::size_t index = 1; index + 1 < vertices.size(); ++index) {
        result = result + cross(vertices[index] - reference,
                                vertices[index + 1] - reference) *
                              0.5;
    }
    return result;
}

[[nodiscard]] std::uint32_t find_or_append_vertex(std::vector<Vec3>& vertices,
                                                  const Vec3& candidate,
                                                  double tolerance) {
    for (std::size_t index = 0; index < vertices.size(); ++index) {
        if (near(vertices[index], candidate, tolerance)) {
            return static_cast<std::uint32_t>(index);
        }
    }
    if (vertices.size() >= std::numeric_limits<std::uint32_t>::max()) {
        throw std::overflow_error("Cut-cell 多面体顶点索引超过 32 位范围");
    }
    vertices.push_back(candidate);
    return static_cast<std::uint32_t>(vertices.size() - 1U);
}

[[nodiscard]] std::vector<Vec3> ordered_cut_polygon(std::vector<Vec3> points,
                                                    const Vec3& outward_normal,
                                                    double tolerance) {
    std::vector<Vec3> unique;
    for (const auto& point : points) {
        const bool exists = std::any_of(unique.begin(), unique.end(), [&](const Vec3& value) {
            return near(value, point, tolerance);
        });
        if (!exists) {
            unique.push_back(point);
        }
    }
    if (unique.size() < 3) {
        return {};
    }
    Vec3 center{};
    for (const auto& point : unique) center = center + point;
    center = center / static_cast<double>(unique.size());
    const std::array<Vec3, 3> axes = {
        Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    const auto reference_axis = *std::min_element(
        axes.begin(), axes.end(), [&](const Vec3& first, const Vec3& second) {
            return std::abs(dot(first, outward_normal)) <
                   std::abs(dot(second, outward_normal));
        });
    const Vec3 tangent_u = cross(reference_axis, outward_normal) /
                           norm(cross(reference_axis, outward_normal));
    const Vec3 tangent_v = cross(outward_normal, tangent_u);
    std::sort(unique.begin(), unique.end(), [&](const Vec3& first,
                                                const Vec3& second) {
        const Vec3 first_delta = first - center;
        const Vec3 second_delta = second - center;
        const double first_angle =
            std::atan2(dot(first_delta, tangent_v), dot(first_delta, tangent_u));
        const double second_angle =
            std::atan2(dot(second_delta, tangent_v), dot(second_delta, tangent_u));
        if (first_angle != second_angle) return first_angle < second_angle;
        if (first.x != second.x) return first.x < second.x;
        if (first.y != second.y) return first.y < second.y;
        return first.z < second.z;
    });
    if (dot(polygon_area_vector(unique), outward_normal) < 0.0) {
        std::reverse(unique.begin(), unique.end());
    }
    return unique;
}

} // 匿名命名空间

OrientedHalfSpace::OrientedHalfSpace(Vec3 point, Vec3 outward_normal,
                                     std::uint64_t boundary_id)
    : point_(point), distance_normal_(outward_normal), boundary_id_(boundary_id) {
    if (!is_finite(point) || !is_finite(outward_normal)) {
        throw std::invalid_argument("Cut-cell 半空间的点和法向必须有限");
    }
    distance_scale_ = norm(outward_normal);
    if (!(distance_scale_ > 0.0)) {
        throw std::invalid_argument("Cut-cell 半空间法向不得为零");
    }
    outward_normal_ = outward_normal / distance_scale_;
}

ConvexPolyhedron make_box_polyhedron(const AABB& box) {
    if (!box.has_positive_volume()) {
        throw std::invalid_argument("Cut-cell 背景盒必须具有正体积");
    }
    const auto& minimum = box.minimum();
    const auto& maximum = box.maximum();
    ConvexPolyhedron result;
    result.vertices = {
        {minimum.x, minimum.y, minimum.z}, {maximum.x, minimum.y, minimum.z},
        {maximum.x, maximum.y, minimum.z}, {minimum.x, maximum.y, minimum.z},
        {minimum.x, minimum.y, maximum.z}, {maximum.x, minimum.y, maximum.z},
        {maximum.x, maximum.y, maximum.z}, {minimum.x, maximum.y, maximum.z}};
    result.faces = {
        {{0, 4, 7, 3}, PolyhedronFaceKind::cartesian, 0},
        {{1, 2, 6, 5}, PolyhedronFaceKind::cartesian, 1},
        {{0, 1, 5, 4}, PolyhedronFaceKind::cartesian, 2},
        {{3, 7, 6, 2}, PolyhedronFaceKind::cartesian, 3},
        {{0, 3, 2, 1}, PolyhedronFaceKind::cartesian, 4},
        {{4, 5, 6, 7}, PolyhedronFaceKind::cartesian, 5}};
    return result;
}

ConvexPolyhedron make_tetrahedron_polyhedron(const Vec3& a, const Vec3& b,
                                             const Vec3& c, const Vec3& d) {
    if (!is_finite(a) || !is_finite(b) || !is_finite(c) || !is_finite(d)) {
        throw std::invalid_argument("Cut-cell 四面体顶点必须有限");
    }
    ConvexPolyhedron result;
    result.vertices = {a, b, c, d};
    const double determinant = dot(b - a, cross(c - a, d - a));
    if (determinant == 0.0) {
        throw std::invalid_argument("Cut-cell 四面体不得退化");
    }
    if (determinant < 0.0) {
        std::swap(result.vertices[2], result.vertices[3]);
    }
    result.faces = {
        {{0, 2, 1}, PolyhedronFaceKind::internal_partition, 0},
        {{0, 1, 3}, PolyhedronFaceKind::internal_partition, 0},
        {{0, 3, 2}, PolyhedronFaceKind::internal_partition, 0},
        {{1, 2, 3}, PolyhedronFaceKind::internal_partition, 0}};
    return result;
}

ConvexPolyhedron clip_convex_polyhedron(const ConvexPolyhedron& polyhedron,
                                        const OrientedHalfSpace& half_space,
                                        double length_tolerance,
                                        PolyhedronFaceKind cut_face_kind) {
    if (length_tolerance < 0.0 || !std::isfinite(length_tolerance)) {
        throw std::invalid_argument("Cut-cell 裁剪长度容差必须为非负有限数");
    }
    if (polyhedron.empty()) {
        return {};
    }
    const double tolerance =
        derived_length_tolerance(polyhedron, half_space, length_tolerance);
    const double distance_tolerance = tolerance * half_space.distance_scale();
    bool any_inside = false;
    bool any_strictly_inside = false;
    bool any_outside = false;
    for (const auto& vertex : polyhedron.vertices) {
        const double distance = half_space.signed_distance(vertex);
        any_inside = any_inside || distance <= distance_tolerance;
        any_strictly_inside = any_strictly_inside || distance < -distance_tolerance;
        any_outside = any_outside || distance > distance_tolerance;
    }
    if (!any_inside || !any_strictly_inside) {
        return {};
    }
    if (!any_outside) {
        return polyhedron;
    }

    std::vector<WorkingFace> working_faces;
    std::vector<Vec3> cut_points;
    for (const auto& face : polyhedron.faces) {
        if (face.vertex_indices.size() < 3) {
            continue;
        }
        WorkingFace clipped{{}, face.kind, face.source_id};
        for (std::size_t edge = 0; edge < face.vertex_indices.size(); ++edge) {
            const Vec3 first =
                polyhedron.vertices[face.vertex_indices[edge]];
            const Vec3 second = polyhedron.vertices[
                face.vertex_indices[(edge + 1U) % face.vertex_indices.size()]];
            const double first_distance = half_space.signed_distance(first);
            const double second_distance = half_space.signed_distance(second);
            const bool first_inside = first_distance <= distance_tolerance;
            const bool second_inside = second_distance <= distance_tolerance;
            if (first_inside) {
                append_unique(clipped.vertices, first, tolerance);
            }
            if (first_inside != second_inside) {
                const double denominator = first_distance - second_distance;
                const double fraction = std::clamp(first_distance / denominator, 0.0, 1.0);
                const Vec3 intersection = first + (second - first) * fraction;
                append_unique(clipped.vertices, intersection, tolerance);
                cut_points.push_back(intersection);
            }
        }
        close_polygon(clipped.vertices, tolerance);
        if (clipped.vertices.size() >= 3 &&
            norm(polygon_area_vector(clipped.vertices)) > tolerance * tolerance) {
            working_faces.push_back(std::move(clipped));
        }
    }

    auto cut_polygon =
        ordered_cut_polygon(std::move(cut_points), half_space.outward_normal(), tolerance);
    if (cut_polygon.size() >= 3 &&
        norm(polygon_area_vector(cut_polygon)) > tolerance * tolerance) {
        working_faces.push_back({std::move(cut_polygon),
                                 cut_face_kind,
                                 half_space.boundary_id()});
    }

    ConvexPolyhedron result;
    for (const auto& face : working_faces) {
        PolyhedronFace output_face{{}, face.kind, face.source_id};
        for (const auto& vertex : face.vertices) {
            const auto index = find_or_append_vertex(result.vertices, vertex, tolerance);
            if (output_face.vertex_indices.empty() ||
                output_face.vertex_indices.back() != index) {
                output_face.vertex_indices.push_back(index);
            }
        }
        if (output_face.vertex_indices.size() > 1 &&
            output_face.vertex_indices.front() == output_face.vertex_indices.back()) {
            output_face.vertex_indices.pop_back();
        }
        if (output_face.vertex_indices.size() >= 3) {
            result.faces.push_back(std::move(output_face));
        }
    }
    return result;
}

PolyhedronGeometry measure_polyhedron(const ConvexPolyhedron& polyhedron) {
    PolyhedronGeometry result;
    if (polyhedron.vertices.empty() || polyhedron.faces.empty()) {
        return result;
    }
    Vec3 reference{};
    for (const auto& vertex : polyhedron.vertices) {
        reference = reference + vertex;
    }
    reference = reference / static_cast<double>(polyhedron.vertices.size());

    std::map<std::pair<std::uint32_t, std::uint32_t>, EdgeUse> edge_uses;
    long double signed_volume = 0.0L;
    Vec3 relative_first_moment{};
    result.faces.reserve(polyhedron.faces.size());
    for (const auto& face : polyhedron.faces) {
        PolygonGeometry face_geometry;
        if (face.vertex_indices.size() >= 3) {
            const Vec3 first = polyhedron.vertices[face.vertex_indices.front()];
            Vec3 centroid_sum{};
            for (std::size_t index = 1; index + 1 < face.vertex_indices.size(); ++index) {
                const Vec3 second = polyhedron.vertices[face.vertex_indices[index]];
                const Vec3 third = polyhedron.vertices[face.vertex_indices[index + 1U]];
                const Vec3 triangle_area_vector =
                    cross(second - first, third - first) * 0.5;
                const double triangle_area = norm(triangle_area_vector);
                face_geometry.area_vector =
                    face_geometry.area_vector + triangle_area_vector;
                face_geometry.area += triangle_area;
                centroid_sum = centroid_sum +
                               (first + second + third) * (triangle_area / 3.0);

                const long double six_volume = static_cast<long double>(dot(
                    first - reference,
                    cross(second - reference, third - reference)));
                const long double tetrahedron_volume = six_volume / 6.0L;
                signed_volume += tetrahedron_volume;
                const Vec3 tetrahedron_centroid_from_reference =
                    ((first - reference) + (second - reference) +
                     (third - reference)) /
                    4.0;
                relative_first_moment =
                    relative_first_moment + tetrahedron_centroid_from_reference *
                                                static_cast<double>(tetrahedron_volume);
            }
            if (face_geometry.area > 0.0) {
                face_geometry.centroid = centroid_sum / face_geometry.area;
            }
            const double vector_area = norm(face_geometry.area_vector);
            if (vector_area > 0.0) {
                face_geometry.outward_normal =
                    face_geometry.area_vector / vector_area;
            }
        }
        result.oriented_area_vector_sum =
            result.oriented_area_vector_sum + face_geometry.area_vector;
        result.faces.push_back(face_geometry);

        for (std::size_t edge = 0; edge < face.vertex_indices.size(); ++edge) {
            const std::uint32_t first = face.vertex_indices[edge];
            const std::uint32_t second =
                face.vertex_indices[(edge + 1U) % face.vertex_indices.size()];
            const auto key = std::minmax(first, second);
            auto& use = edge_uses[{key.first, key.second}];
            ++use.count;
            use.orientation_balance += first < second ? 1 : -1;
        }
    }
    result.closed = !edge_uses.empty() &&
                    std::all_of(edge_uses.begin(), edge_uses.end(), [](const auto& entry) {
                        return entry.second.count == 2 &&
                               entry.second.orientation_balance == 0;
                    });
    result.volume = static_cast<double>(signed_volume);
    result.positive_volume = result.volume > 0.0;
    if (result.volume != 0.0) {
        result.centroid = reference + relative_first_moment / result.volume;
    }
    return result;
}

} // 命名空间 cartmesh
