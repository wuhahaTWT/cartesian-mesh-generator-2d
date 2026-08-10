#include "cartmesh/io/ScalableOpenFoamWriter.hpp"

#include "cartmesh/cutcell/ConvexPolyhedron.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace cartmesh {
namespace {

using Clock = std::chrono::steady_clock;
constexpr std::array<std::uint8_t, 6> opposite_face = {1, 0, 3, 2, 5, 4};
constexpr std::uint32_t no_cell = std::numeric_limits<std::uint32_t>::max();

struct Point2 {
    double u{};
    double v{};
};

struct PiecePatch {
    std::uint32_t owner{};
    std::size_t piece_index{};
    std::size_t face_index{};
    std::uint8_t local_face{};
    Vec3 outward_normal{};
    double area{};
    double covered_area{};
    bool sealed{};
};

struct ExplicitCellPatches {
    std::array<std::vector<PiecePatch>, 6> cartesian;
    std::vector<PiecePatch> partitions;
};

enum class BoundaryKind : std::uint8_t { none, farfield, wall };

struct OutputFace {
    std::vector<Vec3> vertices;
    std::uint32_t owner{};
    std::uint32_t neighbor{no_cell};
    std::uint64_t boundary_id{};
    std::uint64_t anchor_background{};
    BoundaryKind boundary_kind{BoundaryKind::none};
};

struct LatticePoint {
    std::uint32_t i{};
    std::uint32_t j{};
    std::uint32_t k{};
    std::uint64_t linear{};
};

[[nodiscard]] double seconds(Clock::time_point first,
                             Clock::time_point second) noexcept {
    return std::chrono::duration<double>(second - first).count();
}

[[nodiscard]] std::string precise(double value) {
    std::ostringstream output;
    output << std::setprecision(17) << value;
    return output.str();
}

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

[[nodiscard]] double signed_area_2d(
    const std::vector<Point2>& polygon) noexcept {
    double result = 0.0;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const auto& first = polygon[index];
        const auto& second = polygon[(index + 1U) % polygon.size()];
        result += first.u * second.v - second.u * first.v;
    }
    return 0.5 * result;
}

[[nodiscard]] Vec3 face_area_vector(
    const std::vector<Vec3>& polygon) noexcept {
    Vec3 result{};
    if (polygon.size() < 3U) return result;
    const Vec3 origin = polygon.front();
    for (std::size_t index = 1; index + 1U < polygon.size(); ++index) {
        result = result +
                 cross(polygon[index] - origin,
                       polygon[index + 1U] - origin) * 0.5;
    }
    return result;
}

[[nodiscard]] double polygon_area(
    const std::vector<Vec3>& polygon) noexcept {
    return norm(face_area_vector(polygon));
}

void orient_face(std::vector<Vec3>& vertices, const Vec3& normal) {
    if (dot(face_area_vector(vertices), normal) < 0.0) {
        std::reverse(vertices.begin(), vertices.end());
    }
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
    double coordinate_scale = length_tolerance;
    const Point2 reference = subject.front();
    for (const auto point : subject) {
        coordinate_scale = std::max(
            coordinate_scale,
            std::hypot(point.u - reference.u, point.v - reference.v));
    }
    for (const auto point : clip) {
        coordinate_scale = std::max(
            coordinate_scale,
            std::hypot(point.u - reference.u, point.v - reference.v));
    }
    const double side_tolerance = length_tolerance * coordinate_scale;
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
            const bool previous_inside = previous_side >= -side_tolerance;
            const bool current_inside = current_side >= -side_tolerance;
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
    if (subject.size() < 3U ||
        std::abs(signed_area_2d(subject)) <=
            length_tolerance * length_tolerance) {
        return {};
    }
    const double coordinate = axis == 0U ? first.front().x
                              : axis == 1U ? first.front().y
                                           : first.front().z;
    std::vector<Vec3> result;
    result.reserve(subject.size());
    for (const auto point : subject) {
        result.push_back(lift(point, axis, coordinate));
    }
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
    const Vec3 cross_reference = cross(reference_axis, normal);
    const Vec3 tangent_u = cross_reference / norm(cross_reference);
    const Vec3 tangent_v = cross(normal, tangent_u);
    const Vec3 origin = first.front();
    const auto local = [&](const Vec3& point) {
        const Vec3 delta = point - origin;
        return Point2{dot(delta, tangent_u), dot(delta, tangent_v)};
    };
    std::vector<Point2> subject;
    std::vector<Point2> clip;
    subject.reserve(first.size());
    clip.reserve(second.size());
    for (const auto point : first) subject.push_back(local(point));
    for (const auto point : second) clip.push_back(local(point));
    remove_short_edges(subject, length_tolerance);
    remove_short_edges(clip, length_tolerance);
    if (subject.size() < 3U || clip.size() < 3U) return {};
    double coordinate_scale = length_tolerance;
    const Point2 reference = subject.front();
    for (const auto point : subject) {
        coordinate_scale = std::max(
            coordinate_scale,
            std::hypot(point.u - reference.u, point.v - reference.v));
    }
    for (const auto point : clip) {
        coordinate_scale = std::max(
            coordinate_scale,
            std::hypot(point.u - reference.u, point.v - reference.v));
    }
    const double side_tolerance = length_tolerance * coordinate_scale;
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
            const bool previous_inside = previous_side >= -side_tolerance;
            const bool current_inside = current_side >= -side_tolerance;
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
    if (subject.size() < 3U ||
        std::abs(signed_area_2d(subject)) <=
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

[[nodiscard]] std::vector<Vec3> convex_hull_on_plane(
    const std::vector<Vec3>& points, const Vec3& normal,
    double length_tolerance) {
    if (points.size() < 3U) return {};
    const std::array<Vec3, 3> axes = {
        Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0},
        Vec3{0.0, 0.0, 1.0}};
    const auto reference_axis = *std::min_element(
        axes.begin(), axes.end(), [&](const Vec3& lhs, const Vec3& rhs) {
            return std::abs(dot(lhs, normal)) < std::abs(dot(rhs, normal));
        });
    const Vec3 reference_cross = cross(reference_axis, normal);
    const Vec3 tangent_u = reference_cross / norm(reference_cross);
    const Vec3 tangent_v = cross(normal, tangent_u);
    const Vec3 origin = points.front();
    std::vector<Point2> projected;
    projected.reserve(points.size());
    for (const auto point : points) {
        const Vec3 delta = point - origin;
        projected.push_back({dot(delta, tangent_u), dot(delta, tangent_v)});
    }
    std::sort(projected.begin(), projected.end(),
              [](const Point2& lhs, const Point2& rhs) {
                  if (lhs.u != rhs.u) return lhs.u < rhs.u;
                  return lhs.v < rhs.v;
              });
    std::vector<Point2> unique;
    unique.reserve(projected.size());
    for (const auto point : projected) {
        if (unique.empty() ||
            std::hypot(point.u - unique.back().u,
                       point.v - unique.back().v) > length_tolerance) {
            unique.push_back(point);
        }
    }
    if (unique.size() < 3U) return {};
    double coordinate_scale = length_tolerance;
    for (const auto point : unique) {
        coordinate_scale =
            std::max(coordinate_scale, std::hypot(point.u, point.v));
    }
    const double turn_tolerance = length_tolerance * coordinate_scale;
    const auto turn = [](const Point2& first, const Point2& second,
                         const Point2& third) {
        return (second.u - first.u) * (third.v - first.v) -
               (second.v - first.v) * (third.u - first.u);
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

[[nodiscard]] std::vector<Point2> clip_planar_half_plane(
    const std::vector<Point2>& polygon, const Point2& first,
    const Point2& second, bool keep_inside, double side_tolerance) {
    std::vector<Point2> result;
    if (polygon.empty()) return result;
    const auto signed_distance = [&](const Point2& point) {
        const double side =
            (second.u - first.u) * (point.v - first.v) -
            (second.v - first.v) * (point.u - first.u);
        return keep_inside ? side + side_tolerance
                           : -(side + side_tolerance);
    };
    Point2 previous = polygon.back();
    double previous_distance = signed_distance(previous);
    for (const auto current : polygon) {
        const double current_distance = signed_distance(current);
        const bool previous_kept = previous_distance >= 0.0;
        const bool current_kept = current_distance >= 0.0;
        if (previous_kept != current_kept) {
            const double denominator =
                previous_distance - current_distance;
            if (denominator != 0.0) {
                const double parameter = previous_distance / denominator;
                result.push_back(
                    {previous.u + (current.u - previous.u) * parameter,
                     previous.v + (current.v - previous.v) * parameter});
            }
        }
        if (current_kept) result.push_back(current);
        previous = current;
        previous_distance = current_distance;
    }
    return result;
}

[[nodiscard]] std::vector<std::vector<Vec3>>
subtract_convex_planar_polygon(const std::vector<Vec3>& subject,
                               const std::vector<Vec3>& clip,
                               const Vec3& normal,
                               double length_tolerance,
                               double area_tolerance) {
    if (subject.size() < 3U || clip.size() < 3U) return {subject};
    const std::array<Vec3, 3> axes = {
        Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0},
        Vec3{0.0, 0.0, 1.0}};
    const auto reference_axis = *std::min_element(
        axes.begin(), axes.end(), [&](const Vec3& lhs, const Vec3& rhs) {
            return std::abs(dot(lhs, normal)) < std::abs(dot(rhs, normal));
        });
    const Vec3 reference_cross = cross(reference_axis, normal);
    const Vec3 tangent_u = reference_cross / norm(reference_cross);
    const Vec3 tangent_v = cross(normal, tangent_u);
    const Vec3 origin = subject.front();
    const auto project_local = [&](const Vec3& point) {
        const Vec3 delta = point - origin;
        return Point2{dot(delta, tangent_u), dot(delta, tangent_v)};
    };
    std::vector<Point2> current;
    std::vector<Point2> clip_local;
    current.reserve(subject.size());
    clip_local.reserve(clip.size());
    for (const auto point : subject) current.push_back(project_local(point));
    for (const auto point : clip) clip_local.push_back(project_local(point));
    if (signed_area_2d(clip_local) < 0.0) {
        std::reverse(clip_local.begin(), clip_local.end());
    }
    double coordinate_scale = length_tolerance;
    for (const auto point : current) {
        coordinate_scale =
            std::max(coordinate_scale, std::hypot(point.u, point.v));
    }
    for (const auto point : clip_local) {
        coordinate_scale =
            std::max(coordinate_scale, std::hypot(point.u, point.v));
    }
    const double side_tolerance = length_tolerance * coordinate_scale;
    std::vector<std::vector<Point2>> outside_parts;
    for (std::size_t edge = 0; edge < clip_local.size(); ++edge) {
        const Point2 first = clip_local[edge];
        const Point2 second = clip_local[(edge + 1U) % clip_local.size()];
        auto outside = clip_planar_half_plane(
            current, first, second, false, side_tolerance);
        remove_short_edges(outside, length_tolerance);
        if (outside.size() >= 3U &&
            std::abs(signed_area_2d(outside)) > area_tolerance) {
            outside_parts.push_back(std::move(outside));
        }
        current = clip_planar_half_plane(
            current, first, second, true, side_tolerance);
        remove_short_edges(current, length_tolerance);
        if (current.size() < 3U) break;
    }
    std::vector<std::vector<Vec3>> result;
    result.reserve(outside_parts.size());
    for (auto& polygon : outside_parts) {
        std::vector<Vec3> lifted;
        lifted.reserve(polygon.size());
        for (const auto point : polygon) {
            lifted.push_back(origin + tangent_u * point.u +
                             tangent_v * point.v);
        }
        orient_face(lifted, normal);
        if (polygon_area(lifted) > area_tolerance) {
            result.push_back(std::move(lifted));
        }
    }
    return result;
}

[[maybe_unused]] void make_coplanar_internal_faces_disjoint(
    std::vector<OutputFace>& faces, double length_tolerance,
    double area_tolerance) {
    std::sort(faces.begin(), faces.end(),
              [](const OutputFace& lhs, const OutputFace& rhs) {
                  if (lhs.owner != rhs.owner) return lhs.owner < rhs.owner;
                  return lhs.neighbor < rhs.neighbor;
              });
    std::vector<OutputFace> disjoint;
    disjoint.reserve(faces.size());
    const double normal_tolerance =
        2048.0 * std::numeric_limits<double>::epsilon();
    std::size_t group_begin = 0;
    while (group_begin < faces.size()) {
        std::size_t group_end = group_begin + 1U;
        while (group_end < faces.size() &&
               faces[group_end].owner == faces[group_begin].owner &&
               faces[group_end].neighbor == faces[group_begin].neighbor) {
            ++group_end;
        }
        const std::size_t output_group_begin = disjoint.size();
        for (std::size_t input = group_begin; input < group_end; ++input) {
            std::vector<OutputFace> remaining;
            remaining.push_back(std::move(faces[input]));
            for (std::size_t prior = output_group_begin;
                 prior < disjoint.size() && !remaining.empty(); ++prior) {
                const Vec3 prior_area =
                    face_area_vector(disjoint[prior].vertices);
                const Vec3 prior_normal = prior_area / norm(prior_area);
                std::vector<OutputFace> next;
                for (auto& candidate : remaining) {
                    const Vec3 candidate_area =
                        face_area_vector(candidate.vertices);
                    const Vec3 candidate_normal =
                        candidate_area / norm(candidate_area);
                    if (norm(candidate_normal - prior_normal) >
                            normal_tolerance ||
                        std::abs(dot(candidate.vertices.front() -
                                         disjoint[prior].vertices.front(),
                                     prior_normal)) >
                            4.0 * length_tolerance) {
                        next.push_back(std::move(candidate));
                        continue;
                    }
                    auto pieces = subtract_convex_planar_polygon(
                        candidate.vertices, disjoint[prior].vertices,
                        candidate_normal, length_tolerance, area_tolerance);
                    for (auto& vertices : pieces) {
                        OutputFace piece = candidate;
                        piece.vertices = std::move(vertices);
                        next.push_back(std::move(piece));
                    }
                }
                remaining = std::move(next);
            }
            for (auto& face : remaining) disjoint.push_back(std::move(face));
        }
        group_begin = group_end;
    }
    faces = std::move(disjoint);
}

[[maybe_unused]] void make_incident_faces_disjoint(std::vector<OutputFace>& faces,
                                  bool group_by_neighbor,
                                  const std::vector<std::uint8_t>& selected,
                                  double length_tolerance,
                                  double area_tolerance) {
    const auto cell = [&](const OutputFace& face) {
        return group_by_neighbor ? face.neighbor : face.owner;
    };
    std::sort(faces.begin(), faces.end(), [&](const OutputFace& lhs,
                                             const OutputFace& rhs) {
        if (cell(lhs) != cell(rhs)) return cell(lhs) < cell(rhs);
        const double lhs_area = polygon_area(lhs.vertices);
        const double rhs_area = polygon_area(rhs.vertices);
        if (lhs_area != rhs_area) return lhs_area < rhs_area;
        if (lhs.owner != rhs.owner) return lhs.owner < rhs.owner;
        if (lhs.neighbor != rhs.neighbor) return lhs.neighbor < rhs.neighbor;
        if (lhs.boundary_kind != rhs.boundary_kind) {
            return lhs.boundary_kind < rhs.boundary_kind;
        }
        return lhs.boundary_id < rhs.boundary_id;
    });
    std::vector<OutputFace> disjoint;
    disjoint.reserve(faces.size());
    const double normal_tolerance =
        2048.0 * std::numeric_limits<double>::epsilon();
    std::size_t group_begin = 0;
    while (group_begin < faces.size()) {
        std::size_t group_end = group_begin + 1U;
        while (group_end < faces.size() &&
               cell(faces[group_end]) == cell(faces[group_begin])) {
            ++group_end;
        }
        const std::uint32_t group_cell = cell(faces[group_begin]);
        if (group_cell == no_cell ||
            selected[static_cast<std::size_t>(group_cell)] == 0U) {
            for (std::size_t input = group_begin; input < group_end; ++input) {
                disjoint.push_back(std::move(faces[input]));
            }
            group_begin = group_end;
            continue;
        }
        const std::size_t output_group_begin = disjoint.size();
        for (std::size_t input = group_begin; input < group_end; ++input) {
            std::vector<OutputFace> remaining;
            remaining.push_back(std::move(faces[input]));
            for (std::size_t prior = output_group_begin;
                 prior < disjoint.size() && !remaining.empty(); ++prior) {
                const Vec3 prior_area =
                    face_area_vector(disjoint[prior].vertices);
                const Vec3 prior_normal = prior_area / norm(prior_area);
                std::vector<OutputFace> next;
                for (auto& candidate : remaining) {
                    const Vec3 candidate_area =
                        face_area_vector(candidate.vertices);
                    const Vec3 candidate_normal =
                        candidate_area / norm(candidate_area);
                    const double aligned = std::min(
                        norm(candidate_normal - prior_normal),
                        norm(candidate_normal + prior_normal));
                    if (aligned > normal_tolerance ||
                        std::abs(dot(candidate.vertices.front() -
                                         disjoint[prior].vertices.front(),
                                     prior_normal)) >
                            4.0 * length_tolerance) {
                        next.push_back(std::move(candidate));
                        continue;
                    }
                    auto pieces = subtract_convex_planar_polygon(
                        candidate.vertices, disjoint[prior].vertices,
                        candidate_normal, length_tolerance, area_tolerance);
                    for (auto& vertices : pieces) {
                        OutputFace piece = candidate;
                        piece.vertices = std::move(vertices);
                        next.push_back(std::move(piece));
                    }
                }
                remaining = std::move(next);
            }
            for (auto& face : remaining) disjoint.push_back(std::move(face));
        }
        group_begin = group_end;
    }
    faces = std::move(disjoint);
}

[[maybe_unused]] void merge_coplanar_boundary_faces(std::vector<OutputFace>& faces,
                                   double length_tolerance,
                                   double area_tolerance) {
    std::sort(faces.begin(), faces.end(),
              [](const OutputFace& lhs, const OutputFace& rhs) {
                  if (lhs.boundary_kind != rhs.boundary_kind) {
                      return lhs.boundary_kind < rhs.boundary_kind;
                  }
                  if (lhs.boundary_id != rhs.boundary_id) {
                      return lhs.boundary_id < rhs.boundary_id;
                  }
                  return lhs.owner < rhs.owner;
              });
    std::vector<OutputFace> merged;
    merged.reserve(faces.size());
    const double normal_tolerance =
        2048.0 * std::numeric_limits<double>::epsilon();
    std::size_t group_begin = 0;
    while (group_begin < faces.size()) {
        std::size_t group_end = group_begin + 1U;
        while (group_end < faces.size() &&
               faces[group_end].boundary_kind ==
                   faces[group_begin].boundary_kind &&
               faces[group_end].boundary_id ==
                   faces[group_begin].boundary_id &&
               faces[group_end].owner == faces[group_begin].owner) {
            ++group_end;
        }
        std::vector<OutputFace> group;
        group.reserve(group_end - group_begin);
        for (std::size_t index = group_begin; index < group_end; ++index) {
            group.push_back(std::move(faces[index]));
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (std::size_t first = 0; first < group.size() && !changed;
                 ++first) {
                const Vec3 first_area =
                    face_area_vector(group[first].vertices);
                const Vec3 normal = first_area / norm(first_area);
                for (std::size_t second = first + 1U;
                     second < group.size(); ++second) {
                    const Vec3 second_area =
                        face_area_vector(group[second].vertices);
                    const Vec3 second_normal =
                        second_area / norm(second_area);
                    if (norm(second_normal - normal) > normal_tolerance ||
                        std::abs(dot(group[second].vertices.front() -
                                         group[first].vertices.front(),
                                     normal)) >
                            4.0 * length_tolerance) {
                        continue;
                    }
                    std::vector<Vec3> points = group[first].vertices;
                    points.insert(points.end(),
                                  group[second].vertices.begin(),
                                  group[second].vertices.end());
                    auto hull = convex_hull_on_plane(
                        points, normal, length_tolerance);
                    const double hull_area = polygon_area(hull);
                    const double separate_area =
                        polygon_area(group[first].vertices) +
                        polygon_area(group[second].vertices);
                    const double merge_tolerance = std::max(
                        area_tolerance,
                        8192.0 * std::numeric_limits<double>::epsilon() *
                            std::max(hull_area, separate_area));
                    if (hull.size() < 3U ||
                        hull_area > separate_area + merge_tolerance) {
                        continue;
                    }
                    group[first].vertices = std::move(hull);
                    group[first].anchor_background = std::min(
                        group[first].anchor_background,
                        group[second].anchor_background);
                    group.erase(group.begin() +
                                static_cast<std::ptrdiff_t>(second));
                    changed = true;
                    break;
                }
            }
        }
        for (auto& face : group) merged.push_back(std::move(face));
        group_begin = group_end;
    }
    faces = std::move(merged);
}

void merge_coplanar_internal_faces(std::vector<OutputFace>& faces,
                                   double length_tolerance,
                                   double area_tolerance) {
    std::sort(faces.begin(), faces.end(),
              [](const OutputFace& lhs, const OutputFace& rhs) {
                  if (lhs.owner != rhs.owner) return lhs.owner < rhs.owner;
                  return lhs.neighbor < rhs.neighbor;
              });
    std::vector<OutputFace> merged;
    merged.reserve(faces.size());
    const double normal_tolerance =
        2048.0 * std::numeric_limits<double>::epsilon();
    std::size_t group_begin = 0;
    while (group_begin < faces.size()) {
        std::size_t group_end = group_begin + 1U;
        while (group_end < faces.size() &&
               faces[group_end].owner == faces[group_begin].owner &&
               faces[group_end].neighbor == faces[group_begin].neighbor) {
            ++group_end;
        }
        std::vector<OutputFace> group;
        group.reserve(group_end - group_begin);
        for (std::size_t index = group_begin; index < group_end; ++index) {
            group.push_back(std::move(faces[index]));
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (std::size_t first = 0; first < group.size() && !changed;
                 ++first) {
                const Vec3 first_area =
                    face_area_vector(group[first].vertices);
                const Vec3 normal = first_area / norm(first_area);
                for (std::size_t second = first + 1U;
                     second < group.size(); ++second) {
                    const Vec3 second_area =
                        face_area_vector(group[second].vertices);
                    const Vec3 second_normal =
                        second_area / norm(second_area);
                    if (norm(second_normal - normal) > normal_tolerance ||
                        std::abs(dot(group[second].vertices.front() -
                                         group[first].vertices.front(),
                                     normal)) >
                            4.0 * length_tolerance) {
                        continue;
                    }
                    std::vector<Vec3> points = group[first].vertices;
                    points.insert(points.end(),
                                  group[second].vertices.begin(),
                                  group[second].vertices.end());
                    auto hull = convex_hull_on_plane(
                        points, normal, length_tolerance);
                    const double hull_area = polygon_area(hull);
                    const double separate_area =
                        polygon_area(group[first].vertices) +
                        polygon_area(group[second].vertices);
                    const double merge_tolerance = std::max(
                        area_tolerance,
                        8192.0 * std::numeric_limits<double>::epsilon() *
                            std::max(hull_area, separate_area));
                    if (hull.size() < 3U ||
                        hull_area > separate_area + merge_tolerance) {
                        continue;
                    }
                    group[first].vertices = std::move(hull);
                    group[first].anchor_background = std::min(
                        group[first].anchor_background,
                        group[second].anchor_background);
                    group.erase(group.begin() +
                                static_cast<std::ptrdiff_t>(second));
                    changed = true;
                    break;
                }
            }
        }
        for (auto& face : group) merged.push_back(std::move(face));
        group_begin = group_end;
    }
    faces = std::move(merged);
}

[[nodiscard]] std::vector<Vec3> piece_face_vertices(
    const CompactCutCellRecord& record, const PiecePatch& patch) {
    const auto& piece = record.geometry.fluid_polyhedron_pieces[patch.piece_index];
    const auto& face = piece.polyhedron.faces[patch.face_index];
    std::vector<Vec3> result;
    result.reserve(face.vertex_indices.size());
    for (const auto vertex : face.vertex_indices) {
        result.push_back(piece.polyhedron.vertices[vertex]);
    }
    return result;
}

[[nodiscard]] bool is_sealed(const CompactCutCellRecord& record,
                             std::size_t piece_index,
                             std::size_t face_index) noexcept {
    return std::any_of(
        record.numerically_sealed_faces.begin(),
        record.numerically_sealed_faces.end(),
        [&](const NumericallySealedCartesianFace& face) {
            return face.fluid_piece_index == piece_index &&
                   face.polyhedron_face_index == face_index;
        });
}

void set_regular_face_vertices(std::vector<Vec3>& result,
                               const UniformCartesianGrid& grid,
                               std::uint64_t background,
                               std::uint8_t local_face) {
    const auto box = grid.cell_bounds(grid.cell_key(background));
    const Vec3 minimum = box.minimum();
    const Vec3 maximum = box.maximum();
    result.clear();
    result.reserve(4U);
    if (local_face == 0U) {
        result = {{minimum.x, minimum.y, minimum.z},
                  {minimum.x, minimum.y, maximum.z},
                  {minimum.x, maximum.y, maximum.z},
                  {minimum.x, maximum.y, minimum.z}};
    } else if (local_face == 1U) {
        result = {{maximum.x, minimum.y, minimum.z},
                  {maximum.x, maximum.y, minimum.z},
                  {maximum.x, maximum.y, maximum.z},
                  {maximum.x, minimum.y, maximum.z}};
    } else if (local_face == 2U) {
        result = {{minimum.x, minimum.y, minimum.z},
                  {maximum.x, minimum.y, minimum.z},
                  {maximum.x, minimum.y, maximum.z},
                  {minimum.x, minimum.y, maximum.z}};
    } else if (local_face == 3U) {
        result = {{minimum.x, maximum.y, minimum.z},
                  {minimum.x, maximum.y, maximum.z},
                  {maximum.x, maximum.y, maximum.z},
                  {maximum.x, maximum.y, minimum.z}};
    } else if (local_face == 4U) {
        result = {{minimum.x, minimum.y, minimum.z},
                  {minimum.x, maximum.y, minimum.z},
                  {maximum.x, maximum.y, minimum.z},
                  {maximum.x, minimum.y, minimum.z}};
    } else {
        result = {{minimum.x, minimum.y, maximum.z},
                  {maximum.x, minimum.y, maximum.z},
                  {maximum.x, maximum.y, maximum.z},
                  {minimum.x, maximum.y, maximum.z}};
    }
}

[[nodiscard]] std::vector<Vec3> regular_face_vertices(
    const UniformCartesianGrid& grid, std::uint64_t background,
    std::uint8_t local_face) {
    std::vector<Vec3> result;
    set_regular_face_vertices(result, grid, background, local_face);
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
                  const std::string& object, bool binary) {
    output << "FoamFile\n{\n"
           << "    version 2.0;\n"
           << "    format " << (binary ? "binary" : "ascii") << ";\n";
    if (binary) {
        output << "    arch \"LSB;label=32;scalar=64\";\n";
    }
    output << "    class " << klass << ";\n"
           << "    location \"constant/polyMesh\";\n"
           << "    object " << object << ";\n"
           << "}\n\n";
}

void write_label(std::ofstream& output, std::uint64_t value) {
    if (value > static_cast<std::uint64_t>(
                    std::numeric_limits<std::int32_t>::max())) {
        throw std::overflow_error("OpenFOAM label=32 范围不足");
    }
    const auto label = static_cast<std::int32_t>(value);
    output.write(reinterpret_cast<const char*>(&label), sizeof(label));
}

class ActiveLatticePoints {
  public:
    explicit ActiveLatticePoints(std::uint64_t count)
        : count_(count), bits_(static_cast<std::size_t>((count + 63U) / 64U),
                              0U) {}

    void mark(std::uint64_t point) {
        if (point >= count_) throw std::out_of_range("结构点 ID 越界");
        bits_[static_cast<std::size_t>(point / 64U)] |=
            std::uint64_t{1} << (point % 64U);
    }

    void finalize() {
        offsets_.resize(bits_.size() + 1U, 0U);
        for (std::size_t word = 0; word < bits_.size(); ++word) {
            offsets_[word + 1U] =
                offsets_[word] + std::popcount(bits_[word]);
        }
    }

    [[nodiscard]] bool contains(std::uint64_t point) const noexcept {
        return (bits_[static_cast<std::size_t>(point / 64U)] &
                (std::uint64_t{1} << (point % 64U))) != 0U;
    }

    [[nodiscard]] std::uint64_t id(std::uint64_t point) const {
        if (offsets_.empty() || !contains(point)) {
            throw std::runtime_error("结构点尚未激活");
        }
        const std::size_t word = static_cast<std::size_t>(point / 64U);
        const std::uint32_t bit = static_cast<std::uint32_t>(point % 64U);
        const std::uint64_t mask =
            bit == 0U ? 0U : ((std::uint64_t{1} << bit) - 1U);
        return offsets_[word] + std::popcount(bits_[word] & mask);
    }

    [[nodiscard]] std::uint64_t active_count() const noexcept {
        return offsets_.empty() ? 0U : offsets_.back();
    }

    [[nodiscard]] const std::vector<std::uint64_t>& bits() const noexcept {
        return bits_;
    }

  private:
    std::uint64_t count_{};
    std::vector<std::uint64_t> bits_;
    std::vector<std::uint64_t> offsets_;
};

class OffLatticePointWelder {
  public:
    OffLatticePointWelder(Vec3 origin, double tolerance)
        : origin_(origin), tolerance_(tolerance) {}

    [[nodiscard]] std::uint32_t insert(const Vec3& point) {
        const Key key = quantize(point);
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
        if (points_.size() >= no_cell) {
            throw std::overflow_error("非结构点数量超出 32 位范围");
        }
        const auto id = static_cast<std::uint32_t>(points_.size());
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
    std::map<Key, std::vector<std::uint32_t>> buckets_;
};

[[nodiscard]] std::optional<LatticePoint> lattice_point(
    const UniformCartesianGrid& grid, const Vec3& point, double tolerance) {
    const Vec3 local = point - grid.domain().minimum();
    const Vec3 spacing = grid.spacing();
    const std::array<double, 3> coordinate = {
        local.x / spacing.x, local.y / spacing.y, local.z / spacing.z};
    const std::array<std::uint32_t, 3> limit = {
        grid.nx(), grid.ny(), grid.nz()};
    std::array<std::uint32_t, 3> index{};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
        const auto rounded = static_cast<std::int64_t>(
            std::llround(coordinate[axis]));
        if (rounded < 0 ||
            rounded > static_cast<std::int64_t>(limit[axis])) {
            return std::nullopt;
        }
        const double physical_error =
            std::abs(coordinate[axis] - static_cast<double>(rounded)) *
            (axis == 0U ? spacing.x : axis == 1U ? spacing.y : spacing.z);
        if (physical_error > tolerance) return std::nullopt;
        index[axis] = static_cast<std::uint32_t>(rounded);
    }
    const std::uint64_t nx = static_cast<std::uint64_t>(grid.nx()) + 1U;
    const std::uint64_t ny = static_cast<std::uint64_t>(grid.ny()) + 1U;
    return LatticePoint{
        index[0], index[1], index[2],
        static_cast<std::uint64_t>(index[0]) +
            nx * (static_cast<std::uint64_t>(index[1]) +
                  ny * static_cast<std::uint64_t>(index[2]))};
}

[[nodiscard]] Vec3 lattice_coordinate(const UniformCartesianGrid& grid,
                                      std::uint64_t linear) noexcept {
    const std::uint64_t nx = static_cast<std::uint64_t>(grid.nx()) + 1U;
    const std::uint64_t ny = static_cast<std::uint64_t>(grid.ny()) + 1U;
    const std::uint64_t plane = nx * ny;
    const std::uint64_t k = linear / plane;
    const std::uint64_t remainder = linear - k * plane;
    const std::uint64_t j = remainder / nx;
    const std::uint64_t i = remainder - j * nx;
    const Vec3 minimum = grid.domain().minimum();
    const Vec3 spacing = grid.spacing();
    return {minimum.x + spacing.x * static_cast<double>(i),
            minimum.y + spacing.y * static_cast<double>(j),
            minimum.z + spacing.z * static_cast<double>(k)};
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

struct LocalPointIndex {
    std::array<std::vector<Vec3>, 3> by_axis;
};

[[nodiscard]] double coordinate(const Vec3& point,
                                std::size_t axis) noexcept {
    if (axis == 0U) return point.x;
    if (axis == 1U) return point.y;
    return point.z;
}

void split_face_edges(OutputFace& face,
                      const std::vector<const LocalPointIndex*>& candidates,
                      double tolerance) {
    std::vector<Vec3> result;
    for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
        const Vec3 first = face.vertices[edge];
        const Vec3 second = face.vertices[(edge + 1U) % face.vertices.size()];
        const Vec3 direction = second - first;
        const double squared_length = dot(direction, direction);
        if (!(squared_length > 0.0)) continue;
        const std::array<double, 3> magnitude = {
            std::abs(direction.x), std::abs(direction.y),
            std::abs(direction.z)};
        const std::size_t dominant_axis = static_cast<std::size_t>(
            std::max_element(magnitude.begin(), magnitude.end()) -
            magnitude.begin());
        const double coordinate_minimum =
            std::min(coordinate(first, dominant_axis),
                     coordinate(second, dominant_axis)) - tolerance;
        const double coordinate_maximum =
            std::max(coordinate(first, dominant_axis),
                     coordinate(second, dominant_axis)) + tolerance;
        std::vector<std::pair<double, Vec3>> edge_points{{0.0, first}};
        for (const auto* index : candidates) {
            const auto& points = index->by_axis[dominant_axis];
            const auto begin = std::lower_bound(
                points.begin(), points.end(), coordinate_minimum,
                [&](const Vec3& point, double value) {
                    return coordinate(point, dominant_axis) < value;
                });
            const auto end = std::upper_bound(
                begin, points.end(), coordinate_maximum,
                [&](double value, const Vec3& point) {
                    return value < coordinate(point, dominant_axis);
                });
            for (auto candidate = begin; candidate != end; ++candidate) {
                const Vec3 point = *candidate;
                const double parameter =
                    dot(point - first, direction) / squared_length;
                if (parameter <= 0.0 || parameter >= 1.0) continue;
                const Vec3 projection = first + direction * parameter;
                if (norm(point - projection) <= tolerance) {
                    const bool already_a_face_vertex = std::any_of(
                        face.vertices.begin(), face.vertices.end(),
                        [&](const Vec3& vertex) {
                            return norm(vertex - point) <= 2.0 * tolerance;
                    });
                    if (already_a_face_vertex) continue;
                    edge_points.push_back({parameter, point});
                }
            }
        }
        std::sort(edge_points.begin(), edge_points.end(),
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.first < rhs.first;
                  });
        for (const auto& entry : edge_points) {
            if (result.empty() || norm(result.back() - entry.second) > tolerance) {
                result.push_back(entry.second);
            }
        }
    }
    if (result.size() > 1U && norm(result.front() - result.back()) <= tolerance) {
        result.pop_back();
    }
    face.vertices = std::move(result);
}

[[nodiscard]] std::vector<std::uint64_t> neighboring_cells(
    const UniformCartesianGrid& grid, std::uint64_t background) {
    const CellKey key = grid.cell_key(background);
    std::vector<std::uint64_t> result;
    result.reserve(27);
    for (std::int32_t dk = -1; dk <= 1; ++dk) {
        for (std::int32_t dj = -1; dj <= 1; ++dj) {
            for (std::int32_t di = -1; di <= 1; ++di) {
                const std::int64_t i = static_cast<std::int64_t>(key.i) + di;
                const std::int64_t j = static_cast<std::int64_t>(key.j) + dj;
                const std::int64_t k = static_cast<std::int64_t>(key.k) + dk;
                if (i < 0 || j < 0 || k < 0 ||
                    i >= static_cast<std::int64_t>(grid.nx()) ||
                    j >= static_cast<std::int64_t>(grid.ny()) ||
                    k >= static_cast<std::int64_t>(grid.nz())) {
                    continue;
                }
                result.push_back(grid.linear_id(
                    {0, static_cast<std::uint32_t>(i),
                     static_cast<std::uint32_t>(j),
                     static_cast<std::uint32_t>(k)}));
            }
        }
    }
    return result;
}

} // 匿名命名空间

ScalableOpenFoamWriteStats write_scalable_openfoam_poly_mesh(
    const std::filesystem::path& case_directory,
    const UniformCartesianGrid& grid,
    const CompactUniformCutCellMesh& mesh,
    const std::vector<std::pair<std::uint64_t, std::string>>& boundary_names,
    double length_tolerance) {
    static_assert(std::endian::native == std::endian::little,
                  "阶段6二进制 OpenFOAM 写入器当前要求小端主机");
    if (!mesh.invariants_pass()) {
        throw std::invalid_argument("阶段6导出拒绝未通过完整不变量的紧凑网格");
    }
    if (mesh.background_cell_count != grid.cell_count() ||
        mesh.cell_states.size() != grid.cell_count()) {
        throw std::invalid_argument("紧凑网格与 Cartesian 背景域不匹配");
    }
    if (length_tolerance < 0.0 || !std::isfinite(length_tolerance)) {
        throw std::invalid_argument("OpenFOAM 点合并容差必须是非负有限数");
    }
    const double scale = norm(grid.domain().extent());
    const double requested_length_tolerance = length_tolerance;
    length_tolerance = std::max(
        length_tolerance,
        1024.0 * std::numeric_limits<double>::epsilon() * scale);
    const double point_merge_tolerance = std::max(
        requested_length_tolerance,
        64.0 * std::numeric_limits<double>::epsilon() * scale);
    const double topology_point_tolerance = std::max(
        {requested_length_tolerance,
         32768.0 * std::numeric_limits<double>::epsilon() * scale,
         1.0e-8 * std::min({grid.spacing().x, grid.spacing().y,
                            grid.spacing().z})});
    const double face_scale = std::max(
        {grid.spacing().x * grid.spacing().y,
         grid.spacing().x * grid.spacing().z,
         grid.spacing().y * grid.spacing().z});
    const double area_tolerance = length_tolerance * length_tolerance;
    const double overlap_roundoff_area = std::max(
        area_tolerance,
        4096.0 * std::numeric_limits<double>::epsilon() * face_scale);
    const double coverage_tolerance = std::max(
        {area_tolerance,
         32768.0 * std::numeric_limits<double>::epsilon() * face_scale,
         16.0 * length_tolerance *
             std::max({grid.spacing().x, grid.spacing().y,
                       grid.spacing().z}),
         16.0 * topology_point_tolerance *
             std::max({grid.spacing().x, grid.spacing().y,
                       grid.spacing().z})});

    const auto total_start = Clock::now();
    ScalableOpenFoamWriteStats stats;
    struct RecordExportPlan {
        std::vector<std::uint32_t> piece_group;
        std::vector<Vec3> group_centroid;
    };
    std::vector<RecordExportPlan> export_plans(mesh.explicit_cells.size());
    for (std::size_t record_index = 0;
         record_index < mesh.explicit_cells.size(); ++record_index) {
        const auto& record = mesh.explicit_cells[record_index];
        if (!record.geometry.cut) continue;
        auto& plan = export_plans[record_index];
        plan.piece_group.resize(
            record.geometry.fluid_polyhedron_pieces.size(), 0U);
        const std::size_t component_count = std::max<std::size_t>(
            record.geometry.fluid_component_count, 1U);
        for (std::size_t component = 0; component < component_count;
             ++component) {
            std::vector<std::size_t> pieces;
            double volume = 0.0;
            Vec3 first_moment{};
            for (std::size_t piece_index = 0;
                 piece_index <
                 record.geometry.fluid_polyhedron_pieces.size();
                 ++piece_index) {
                const auto& piece =
                    record.geometry.fluid_polyhedron_pieces[piece_index];
                if (piece.component_id != component) continue;
                pieces.push_back(piece_index);
                volume += piece.geometry.volume;
                first_moment = first_moment +
                               piece.geometry.centroid * piece.geometry.volume;
            }
            if (pieces.empty()) continue;
            const Vec3 component_centroid = first_moment / volume;
            const auto group = static_cast<std::uint32_t>(
                plan.group_centroid.size());
            for (const auto piece_index : pieces) {
                plan.piece_group[piece_index] = group;
            }
            plan.group_centroid.push_back(component_centroid);
        }
    }

    std::vector<std::uint32_t> cell_start(
        static_cast<std::size_t>(grid.cell_count()), no_cell);
    std::vector<std::uint8_t> solver_cell_is_cut;
    std::vector<std::uint64_t> solver_cell_background;
    solver_cell_is_cut.reserve(
        static_cast<std::size_t>(mesh.solver_cell_count));
    solver_cell_background.reserve(
        static_cast<std::size_t>(mesh.solver_cell_count));
    const auto append_solver_cell = [&](bool cut,
                                        std::uint64_t background) {
        solver_cell_is_cut.push_back(cut ? 1U : 0U);
        solver_cell_background.push_back(background);
    };
    std::uint64_t next_cell = 0;
    for (std::uint64_t background = 0; background < grid.cell_count();
         ++background) {
        const auto state = mesh.state(background);
        if (state == CompactCellState::solid) continue;
        if (state == CompactCellState::conflict) {
            throw std::runtime_error("导出遇到分类冲突单元");
        }
        cell_start[static_cast<std::size_t>(background)] =
            static_cast<std::uint32_t>(next_cell);
        if (state == CompactCellState::full_fluid) {
            const auto component = mesh.full_component_labels[
                static_cast<std::size_t>(background)];
            if (component >= mesh.full_component_region_ids.size()) {
                throw std::runtime_error(
                    "完整流体单元缺少全局区域 ID");
            }
            append_solver_cell(false, background);
            ++next_cell;
            continue;
        }
        const auto* record = mesh.find_explicit_cell(background);
        if (record == nullptr) {
            throw std::runtime_error("显式表面单元缺少几何记录");
        }
        if (!record->geometry.cut) {
            if (record->geometry.fluid_component_region_ids.empty()) {
                throw std::runtime_error(
                    "显式完整流体单元缺少全局区域 ID");
            }
            append_solver_cell(false, background);
            ++next_cell;
        } else {
            const std::size_t record_index = static_cast<std::size_t>(
                record - mesh.explicit_cells.data());
            const auto& plan = export_plans[record_index];
            std::vector<std::uint64_t> group_region(
                plan.group_centroid.size(),
                std::numeric_limits<std::uint64_t>::max());
            for (std::size_t piece_index = 0;
                 piece_index <
                 record->geometry.fluid_polyhedron_pieces.size();
                 ++piece_index) {
                const auto group = plan.piece_group[piece_index];
                const auto region =
                    record->geometry.fluid_polyhedron_pieces[piece_index]
                        .global_region_id;
                auto& stored = group_region[group];
                if (stored ==
                    std::numeric_limits<std::uint64_t>::max()) {
                    stored = region;
                } else if (stored != region) {
                    throw std::runtime_error(
                        "凸片组跨越了全局流体区域");
                }
            }
            for (std::size_t group = 0;
                 group < plan.group_centroid.size(); ++group) {
                if (group_region[group] ==
                    std::numeric_limits<std::uint64_t>::max()) {
                    throw std::runtime_error(
                        "凸片组缺少全局流体区域 ID");
                }
                append_solver_cell(true, background);
            }
            next_cell += plan.group_centroid.size();
        }
    }
    if (next_cell > no_cell) {
        throw std::overflow_error("OpenFOAM label=32 无法表示求解器单元数");
    }
    stats.solver_cell_count = next_cell;

    std::vector<ExplicitCellPatches> explicit_patches(mesh.explicit_cells.size());
    for (std::size_t record_index = 0;
         record_index < mesh.explicit_cells.size(); ++record_index) {
        const auto& record = mesh.explicit_cells[record_index];
        if (!record.geometry.cut) continue;
        auto& output = explicit_patches[record_index];
        for (std::size_t piece_index = 0;
             piece_index < record.geometry.fluid_polyhedron_pieces.size();
             ++piece_index) {
            const auto& piece =
                record.geometry.fluid_polyhedron_pieces[piece_index];
            if (piece.component_id >=
                std::max<std::size_t>(record.geometry.fluid_component_count,
                                      1U)) {
                throw std::runtime_error("Cut-cell 凸片 component ID 越界");
            }
            const auto owner = static_cast<std::uint32_t>(
                static_cast<std::uint64_t>(cell_start[static_cast<std::size_t>(
                    record.background_cell_id)]) +
                export_plans[record_index].piece_group[piece_index]);
            for (std::size_t face_index = 0;
                 face_index < piece.polyhedron.faces.size(); ++face_index) {
                const auto& face = piece.polyhedron.faces[face_index];
                const auto& geometry = piece.geometry.faces[face_index];
                if (geometry.area <= area_tolerance) continue;
                PiecePatch patch{owner, piece_index, face_index, 0,
                                 geometry.outward_normal, geometry.area, 0.0,
                                 false};
                if (face.kind == PolyhedronFaceKind::cartesian) {
                    if (face.source_id >= 6U) {
                        throw std::runtime_error("Cut-cell Cartesian 面编号越界");
                    }
                    patch.local_face =
                        static_cast<std::uint8_t>(face.source_id);
                    patch.sealed = is_sealed(record, piece_index, face_index);
                    output.cartesian[patch.local_face].push_back(patch);
                } else {
                    output.partitions.push_back(patch);
                }
            }
        }
    }

    const auto explicit_index = [&](std::uint64_t background) -> std::size_t {
        const auto found = std::lower_bound(
            mesh.explicit_cells.begin(), mesh.explicit_cells.end(), background,
            [](const CompactCutCellRecord& record, std::uint64_t id) {
                return record.background_cell_id < id;
            });
        if (found == mesh.explicit_cells.end() ||
            found->background_cell_id != background) {
            throw std::runtime_error("找不到显式单元 patch 索引");
        }
        return static_cast<std::size_t>(found - mesh.explicit_cells.begin());
    };

    std::vector<OutputFace> special_internal;
    std::vector<OutputFace> special_boundary;
    special_internal.reserve(mesh.explicit_surface_cell_count * 8U);
    special_boundary.reserve(mesh.boundary_cell_count * 8U);

    const auto cartesian_face_area = [&](std::uint8_t local_face) noexcept {
        if (local_face < 2U) {
            return grid.spacing().y * grid.spacing().z;
        }
        if (local_face < 4U) {
            return grid.spacing().x * grid.spacing().z;
        }
        return grid.spacing().x * grid.spacing().y;
    };

    const auto process_background_interface = [&](std::uint64_t first_id,
                                                   std::uint64_t second_id,
                                                   std::uint8_t first_face) {
        const auto first_state = mesh.state(first_id);
        const auto second_state = mesh.state(second_id);
        const bool first_fluid = first_state == CompactCellState::full_fluid ||
                                 first_state == CompactCellState::explicit_surface;
        const bool second_fluid = second_state == CompactCellState::full_fluid ||
                                  second_state == CompactCellState::explicit_surface;
        if (!first_fluid || !second_fluid) return;
        const bool first_explicit =
            first_state == CompactCellState::explicit_surface;
        const bool second_explicit =
            second_state == CompactCellState::explicit_surface;
        if (!first_explicit && !second_explicit) return;
        const auto* first_record =
            first_explicit ? mesh.find_explicit_cell(first_id) : nullptr;
        const auto* second_record =
            second_explicit ? mesh.find_explicit_cell(second_id) : nullptr;
        const bool first_cut = first_record != nullptr && first_record->geometry.cut;
        const bool second_cut =
            second_record != nullptr && second_record->geometry.cut;
        const std::uint8_t second_face = opposite_face[first_face];
        const double expected_first =
            first_record != nullptr
                ? first_record->geometry.cartesian_faces[first_face].area
                : cartesian_face_area(first_face);
        const double expected_second =
            second_record != nullptr
                ? second_record->geometry.cartesian_faces[second_face].area
                : cartesian_face_area(second_face);
        if (!first_cut && !second_cut) {
            if (std::max(expected_first, expected_second) <= area_tolerance) {
                return;
            }
            const double full_area = cartesian_face_area(first_face);
            if (std::abs(expected_first - full_area) > coverage_tolerance ||
                std::abs(expected_second - full_area) > coverage_tolerance) {
                throw std::runtime_error(
                    "未切割共面 wall 产生了尚不支持的部分 Cartesian 开口");
            }
            auto vertices = regular_face_vertices(grid, first_id, first_face);
            special_internal.push_back(
                {std::move(vertices),
                 cell_start[static_cast<std::size_t>(first_id)],
                 cell_start[static_cast<std::size_t>(second_id)], 0U,
                 first_id, BoundaryKind::none});
            return;
        }
        double first_covered = 0.0;
        double second_covered = 0.0;
        if (first_cut && !second_cut) {
            auto& patches = explicit_patches[explicit_index(first_id)]
                                .cartesian[first_face];
            for (auto& patch : patches) {
                if (patch.sealed) continue;
                auto vertices = piece_face_vertices(*first_record, patch);
                orient_face(vertices, face_normal(first_face));
                special_internal.push_back(
                    {std::move(vertices), patch.owner,
                     cell_start[static_cast<std::size_t>(second_id)], 0U,
                     first_id, BoundaryKind::none});
                patch.covered_area += patch.area;
                first_covered += patch.area;
                second_covered += patch.area;
            }
        } else if (!first_cut && second_cut) {
            auto& patches = explicit_patches[explicit_index(second_id)]
                                .cartesian[second_face];
            for (auto& patch : patches) {
                if (patch.sealed) continue;
                auto vertices = piece_face_vertices(*second_record, patch);
                orient_face(vertices, face_normal(first_face));
                special_internal.push_back(
                    {std::move(vertices),
                     cell_start[static_cast<std::size_t>(first_id)],
                     patch.owner, 0U, first_id, BoundaryKind::none});
                patch.covered_area += patch.area;
                first_covered += patch.area;
                second_covered += patch.area;
            }
        } else {
            auto& first_patches = explicit_patches[explicit_index(first_id)]
                                      .cartesian[first_face];
            auto& second_patches = explicit_patches[explicit_index(second_id)]
                                       .cartesian[second_face];
            const auto unique_open_owner = [](const auto& patches) {
                std::vector<std::uint32_t> owners;
                for (const auto& patch : patches) {
                    if (patch.sealed) continue;
                    if (std::find(owners.begin(), owners.end(), patch.owner) ==
                        owners.end()) {
                        owners.push_back(patch.owner);
                    }
                }
                return owners;
            };
            const auto first_owners = unique_open_owner(first_patches);
            const auto second_owners = unique_open_owner(second_patches);
            if (first_owners.size() == 1U &&
                second_owners.size() == 1U) {
                for (auto& patch : first_patches) {
                    if (patch.sealed) continue;
                    auto vertices =
                        piece_face_vertices(*first_record, patch);
                    orient_face(vertices, face_normal(first_face));
                    special_internal.push_back(
                        {std::move(vertices), first_owners.front(),
                         second_owners.front(), 0U, first_id,
                         BoundaryKind::none});
                    patch.covered_area += patch.area;
                    first_covered += patch.area;
                }
                for (auto& patch : second_patches) {
                    if (patch.sealed) continue;
                    patch.covered_area += patch.area;
                    second_covered += patch.area;
                }
            } else {
                for (auto& first_patch : first_patches) {
                    if (first_patch.sealed) continue;
                    const auto first_vertices =
                        piece_face_vertices(*first_record, first_patch);
                    for (auto& second_patch : second_patches) {
                        if (second_patch.sealed) continue;
                        const auto second_vertices =
                            piece_face_vertices(*second_record, second_patch);
                        auto overlap = intersect_convex_cartesian_polygons(
                            first_vertices, second_vertices, first_face,
                            length_tolerance);
                        double area = polygon_area(overlap);
                        if (area <= area_tolerance) continue;
                        if (area <= overlap_roundoff_area &&
                            area < 1.0e-6 *
                                       std::min(first_patch.area,
                                                second_patch.area)) {
                            continue;
                        }
                        if (std::abs(area - first_patch.area) <=
                                coverage_tolerance &&
                            std::abs(area - second_patch.area) <=
                                coverage_tolerance) {
                            overlap = first_vertices;
                            area = first_patch.area;
                        }
                        orient_face(overlap, face_normal(first_face));
                        special_internal.push_back(
                            {std::move(overlap), first_patch.owner,
                             second_patch.owner, 0U, first_id,
                             BoundaryKind::none});
                        first_patch.covered_area += area;
                        second_patch.covered_area += area;
                        first_covered += area;
                        second_covered += area;
                    }
                }
            }
        }
        const double first_error =
            std::abs(first_covered - expected_first);
        const double second_error =
            std::abs(second_covered - expected_second);
        const double interface_error = std::max(first_error, second_error);
        const double diagnostic_coverage_limit =
            std::max(coverage_tolerance, 1.0e-5 * face_scale);
        if (interface_error > coverage_tolerance) {
            ++stats.background_interface_coverage_relaxation_count;
            stats.maximum_background_interface_coverage_error = std::max(
                stats.maximum_background_interface_coverage_error,
                interface_error);
        }
        if (interface_error > diagnostic_coverage_limit) {
            throw std::runtime_error(
                "阶段6 OpenFOAM 背景接口公共细分面积不守恒：first=" +
                std::to_string(first_id) + " second=" +
                std::to_string(second_id) + " firstCovered=" +
                precise(first_covered) + " firstExpected=" +
                precise(expected_first) + " firstError=" +
                precise(first_covered - expected_first) +
                " secondCovered=" + precise(second_covered) +
                " secondExpected=" + precise(expected_second) +
                " secondError=" +
                precise(second_covered - expected_second) + " tolerance=" +
                precise(coverage_tolerance) + " firstState=" +
                std::to_string(static_cast<unsigned>(first_state)) +
                " secondState=" +
                std::to_string(static_cast<unsigned>(second_state)));
        }
    };

    for (std::uint32_t k = 0; k < grid.nz(); ++k) {
        for (std::uint32_t j = 0; j < grid.ny(); ++j) {
            for (std::uint32_t i = 0; i < grid.nx(); ++i) {
                const auto id = grid.linear_id({0, i, j, k});
                if (i + 1U < grid.nx()) {
                    process_background_interface(id, id + 1U, 1U);
                }
                if (j + 1U < grid.ny()) {
                    process_background_interface(id, id + grid.nx(), 3U);
                }
                if (k + 1U < grid.nz()) {
                    process_background_interface(
                        id, id + static_cast<std::uint64_t>(grid.nx()) *
                                     grid.ny(),
                        5U);
                }
            }
        }
    }

    for (std::size_t record_index = 0;
         record_index < mesh.explicit_cells.size(); ++record_index) {
        const auto& record = mesh.explicit_cells[record_index];
        if (!record.geometry.cut) {
            const CellKey key = grid.cell_key(record.background_cell_id);
            const std::array<bool, 6> exterior = {
                key.i == 0U, key.i + 1U == grid.nx(),
                key.j == 0U, key.j + 1U == grid.ny(),
                key.k == 0U, key.k + 1U == grid.nz()};
            for (std::uint8_t face_id = 0; face_id < 6U; ++face_id) {
                if (!exterior[face_id] ||
                    record.geometry.cartesian_faces[face_id].area <=
                        area_tolerance) {
                    continue;
                }
                auto vertices =
                    regular_face_vertices(grid, record.background_cell_id,
                                          face_id);
                special_boundary.push_back(
                    {std::move(vertices),
                     cell_start[static_cast<std::size_t>(
                         record.background_cell_id)],
                     no_cell, 0U, record.background_cell_id,
                     BoundaryKind::farfield});
            }
            for (const auto& embedded :
                 record.geometry.embedded_boundary_faces) {
                auto vertices = embedded.vertices;
                orient_face(vertices, embedded.outward_normal);
                special_boundary.push_back(
                    {std::move(vertices),
                     cell_start[static_cast<std::size_t>(
                         record.background_cell_id)],
                     no_cell, embedded.boundary_id,
                     record.background_cell_id, BoundaryKind::wall});
            }
            continue;
        }
        auto& data = explicit_patches[record_index];
        for (std::size_t first = 0; first < data.partitions.size(); ++first) {
            for (std::size_t second = first + 1U;
                 second < data.partitions.size(); ++second) {
                auto& first_patch = data.partitions[first];
                auto& second_patch = data.partitions[second];
                if (dot(first_patch.outward_normal,
                        second_patch.outward_normal) > -1.0 + 1.0e-10) {
                    continue;
                }
                const auto first_vertices =
                    piece_face_vertices(record, first_patch);
                const auto second_vertices =
                    piece_face_vertices(record, second_patch);
                auto overlap = intersect_convex_planar_polygons(
                    first_vertices, second_vertices,
                    first_patch.outward_normal, length_tolerance);
                double area = polygon_area(overlap);
                if (area <= area_tolerance) continue;
                if (area <= overlap_roundoff_area &&
                    area < 1.0e-6 *
                               std::min(first_patch.area,
                                        second_patch.area)) {
                    continue;
                }
                if (std::abs(area - first_patch.area) <=
                        coverage_tolerance &&
                    std::abs(area - second_patch.area) <=
                        coverage_tolerance) {
                    overlap = first_vertices;
                    area = first_patch.area;
                }
                if (first_patch.owner == second_patch.owner) {
                    first_patch.covered_area += area;
                    second_patch.covered_area += area;
                    continue;
                }
                std::uint32_t owner = first_patch.owner;
                std::uint32_t neighbor = second_patch.owner;
                Vec3 normal = first_patch.outward_normal;
                if (owner > neighbor) {
                    std::swap(owner, neighbor);
                    normal = normal * -1.0;
                }
                orient_face(overlap, normal);
                special_internal.push_back(
                    {std::move(overlap), owner, neighbor, 0U,
                     record.background_cell_id, BoundaryKind::none});
                first_patch.covered_area += area;
                second_patch.covered_area += area;
            }
        }

        const CellKey key = grid.cell_key(record.background_cell_id);
        const std::array<bool, 6> exterior = {
            key.i == 0U, key.i + 1U == grid.nx(),
            key.j == 0U, key.j + 1U == grid.ny(),
            key.k == 0U, key.k + 1U == grid.nz()};
        for (std::size_t face_index = 0; face_index < 6U; ++face_index) {
            if (!exterior[face_index]) continue;
            for (auto& patch : data.cartesian[face_index]) {
                if (patch.sealed) continue;
                auto vertices = piece_face_vertices(record, patch);
                orient_face(vertices,
                            face_normal(static_cast<std::uint8_t>(face_index)));
                special_boundary.push_back(
                    {std::move(vertices), patch.owner, no_cell, 0U,
                     record.background_cell_id, BoundaryKind::farfield});
                patch.covered_area += patch.area;
            }
        }

        for (const auto& embedded : record.geometry.embedded_boundary_faces) {
            double covered = 0.0;
            std::vector<OutputFace> embedded_fragments;
            for (auto& patch : data.partitions) {
                if (dot(patch.outward_normal, embedded.outward_normal) <
                    1.0 - 1.0e-10) {
                    continue;
                }
                const auto patch_vertices = piece_face_vertices(record, patch);
                auto overlap = intersect_convex_planar_polygons(
                    patch_vertices, embedded.vertices,
                    embedded.outward_normal, length_tolerance);
                const double area = polygon_area(overlap);
                if (area <= area_tolerance) continue;
                orient_face(overlap, embedded.outward_normal);
                embedded_fragments.push_back(
                    {std::move(overlap), patch.owner, no_cell,
                     embedded.boundary_id, record.background_cell_id,
                     BoundaryKind::wall});
                patch.covered_area += area;
                covered += area;
            }
            if (covered < embedded.area - coverage_tolerance) {
                for (auto& face_patches : data.cartesian) {
                    for (auto& patch : face_patches) {
                        if (dot(face_normal(patch.local_face),
                                embedded.outward_normal) <
                            1.0 - 1.0e-10) {
                            continue;
                        }
                        const auto patch_vertices =
                            piece_face_vertices(record, patch);
                        auto overlap = intersect_convex_planar_polygons(
                            patch_vertices, embedded.vertices,
                            embedded.outward_normal, length_tolerance);
                        const double area = polygon_area(overlap);
                        if (area <= area_tolerance) continue;
                        orient_face(overlap, embedded.outward_normal);
                        embedded_fragments.push_back(
                            {std::move(overlap), patch.owner, no_cell,
                             embedded.boundary_id,
                             record.background_cell_id,
                             BoundaryKind::wall});
                        patch.covered_area += area;
                        covered += area;
                    }
                }
            }
            std::vector<std::uint32_t> fragment_owners;
            for (const auto& fragment : embedded_fragments) {
                if (std::find(fragment_owners.begin(), fragment_owners.end(),
                              fragment.owner) == fragment_owners.end()) {
                    fragment_owners.push_back(fragment.owner);
                }
            }
            if (fragment_owners.size() == 1U) {
                auto vertices = embedded.vertices;
                orient_face(vertices, embedded.outward_normal);
                special_boundary.push_back(
                    {std::move(vertices), fragment_owners.front(), no_cell,
                     embedded.boundary_id, record.background_cell_id,
                     BoundaryKind::wall});
            } else {
                for (auto& fragment : embedded_fragments) {
                    special_boundary.push_back(std::move(fragment));
                }
            }
            const double embedded_error =
                std::abs(covered - embedded.area);
            if (embedded_error > coverage_tolerance) {
                ++stats.explicit_patch_coverage_relaxation_count;
                stats.maximum_explicit_patch_coverage_error = std::max(
                    stats.maximum_explicit_patch_coverage_error,
                    embedded_error);
            }
            if (embedded_error > std::max(coverage_tolerance,
                                          1.0e-5 * face_scale)) {
                throw std::runtime_error(
                    "阶段6 OpenFOAM 嵌入边界细分面积不守恒：background=" +
                    std::to_string(record.background_cell_id));
            }
        }
    }

    for (std::size_t record_index = 0;
         record_index < mesh.explicit_cells.size(); ++record_index) {
        const auto& record = mesh.explicit_cells[record_index];
        const auto& data = explicit_patches[record_index];
        for (const auto& patches : data.cartesian) {
            for (const auto& patch : patches) {
                const double error =
                    std::abs(patch.covered_area - patch.area);
                if (error > coverage_tolerance) {
                    ++stats.explicit_patch_coverage_relaxation_count;
                    stats.maximum_explicit_patch_coverage_error = std::max(
                        stats.maximum_explicit_patch_coverage_error, error);
                }
                if (error > std::max(coverage_tolerance,
                                     1.0e-5 * face_scale)) {
                    throw std::runtime_error(
                        "阶段6 OpenFOAM Cartesian patch 未完整覆盖：background=" +
                        std::to_string(record.background_cell_id));
                }
            }
        }
        for (const auto& patch : data.partitions) {
            const double error =
                std::abs(patch.covered_area - patch.area);
            if (error > coverage_tolerance) {
                ++stats.explicit_patch_coverage_relaxation_count;
                stats.maximum_explicit_patch_coverage_error = std::max(
                    stats.maximum_explicit_patch_coverage_error, error);
            }
            if (error > std::max(coverage_tolerance,
                                 1.0e-5 * face_scale)) {
                throw std::runtime_error(
                    "阶段6 OpenFOAM arrangement patch 未完整覆盖：background=" +
                    std::to_string(record.background_cell_id));
            }
        }
    }

    merge_coplanar_internal_faces(special_internal, length_tolerance,
                                  area_tolerance);

    std::map<std::uint64_t, std::vector<Vec3>> local_special_points;
    std::map<std::uint64_t, LocalPointIndex> local_special_point_indices;
    const auto collect_points = [&](const OutputFace& face) {
        auto& bucket = local_special_points[face.anchor_background];
        bucket.insert(bucket.end(), face.vertices.begin(), face.vertices.end());
    };
    const auto rebuild_local_point_indices = [&]() {
        local_special_point_indices.clear();
        for (const auto& [background, points] : local_special_points) {
            auto& index = local_special_point_indices[background];
            for (std::size_t axis = 0; axis < 3U; ++axis) {
                index.by_axis[axis] = points;
                std::sort(
                    index.by_axis[axis].begin(),
                    index.by_axis[axis].end(),
                    [axis](const Vec3& lhs, const Vec3& rhs) {
                        const double lhs_coordinate = coordinate(lhs, axis);
                        const double rhs_coordinate = coordinate(rhs, axis);
                        if (lhs_coordinate != rhs_coordinate) {
                            return lhs_coordinate < rhs_coordinate;
                        }
                        if (lhs.x != rhs.x) return lhs.x < rhs.x;
                        if (lhs.y != rhs.y) return lhs.y < rhs.y;
                        return lhs.z < rhs.z;
                    });
            }
        }
    };
    for (const auto& face : special_internal) collect_points(face);
    for (const auto& face : special_boundary) collect_points(face);
    for (auto& [background, points] : local_special_points) {
        static_cast<void>(background);
        std::sort(points.begin(), points.end(), [](const Vec3& lhs,
                                                   const Vec3& rhs) {
            if (lhs.x != rhs.x) return lhs.x < rhs.x;
            if (lhs.y != rhs.y) return lhs.y < rhs.y;
            return lhs.z < rhs.z;
        });
        std::vector<Vec3> unique;
        unique.reserve(points.size());
        for (const auto point : points) {
            if (unique.empty() ||
                norm(unique.back() - point) > length_tolerance) {
                unique.push_back(point);
            }
        }
        points = std::move(unique);
    }
    rebuild_local_point_indices();
    const auto split_special = [&](OutputFace& face) {
        const std::size_t original_vertex_count = face.vertices.size();
        const double original_area = polygon_area(face.vertices);
        double minimum_edge_length =
            std::numeric_limits<double>::infinity();
        double maximum_edge_length = 0.0;
        for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
            const double edge_length =
                norm(face.vertices[(edge + 1U) % face.vertices.size()] -
                     face.vertices[edge]);
            minimum_edge_length =
                std::min(minimum_edge_length, edge_length);
            maximum_edge_length =
                std::max(maximum_edge_length, edge_length);
        }
        if (original_area <=
            topology_point_tolerance * maximum_edge_length) return;
        const double local_split_tolerance = std::min(
            topology_point_tolerance, 0.125 * minimum_edge_length);
        std::vector<const LocalPointIndex*> candidates;
        for (const auto cell :
             neighboring_cells(grid, face.anchor_background)) {
            const auto found = local_special_point_indices.find(cell);
            if (found != local_special_point_indices.end()) {
                candidates.push_back(&found->second);
            }
        }
        split_face_edges(face, candidates, local_split_tolerance);
        if (face.vertices.size() < 3U ||
            polygon_area(face.vertices) <= area_tolerance) {
            throw std::runtime_error(
                "特殊面局部边分裂后退化：owner=" +
                std::to_string(face.owner) + " neighbour=" +
                std::to_string(face.neighbor) + " background=" +
                std::to_string(face.anchor_background) + " kind=" +
                std::to_string(static_cast<unsigned>(face.boundary_kind)) +
                " beforeVertices=" +
                std::to_string(original_vertex_count) + " afterVertices=" +
                std::to_string(face.vertices.size()) + " beforeArea=" +
                precise(original_area) + " afterArea=" +
                precise(polygon_area(face.vertices)) + " areaTolerance=" +
                precise(area_tolerance));
        }
    };
    for (auto& face : special_internal) split_special(face);
    for (auto& face : special_boundary) split_special(face);

    // 普通尺度的公共顶点使用较宽的拓扑容差统一到同一坐标；若统一后
    // 一个极小正面积面会退化，则该面保留原始双精度坐标，不删除小片。
    OffLatticePointWelder topology_points(grid.domain().minimum(),
                                          topology_point_tolerance);
    const auto normalize_topology_points = [&](OutputFace& face) {
        const double original_area = polygon_area(face.vertices);
        double maximum_edge_length = 0.0;
        for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
            maximum_edge_length = std::max(
                maximum_edge_length,
                norm(face.vertices[(edge + 1U) % face.vertices.size()] -
                     face.vertices[edge]));
        }
        if (original_area <=
            topology_point_tolerance * maximum_edge_length) {
            ++stats.topology_collapsed_face_count;
            stats.topology_collapsed_face_area += original_area;
            return false;
        }
        std::vector<Vec3> normalized;
        std::vector<std::uint32_t> ids;
        normalized.reserve(face.vertices.size());
        ids.reserve(face.vertices.size());
        for (const auto point : face.vertices) {
            const auto id = topology_points.insert(point);
            if (std::find(ids.begin(), ids.end(), id) != ids.end()) continue;
            ids.push_back(id);
            normalized.push_back(topology_points.points()[id]);
        }
        if (normalized.size() >= 3U &&
            polygon_area(normalized) > area_tolerance) {
            face.vertices = std::move(normalized);
            return true;
        }
        ++stats.topology_collapsed_face_count;
        stats.topology_collapsed_face_area += original_area;
        return false;
    };
    std::erase_if(special_internal, [&](OutputFace& face) {
        return !normalize_topology_points(face);
    });
    std::erase_if(special_boundary, [&](OutputFace& face) {
        return !normalize_topology_points(face);
    });

    // 后续规则面的边分裂必须查询规范化后的点，而不是规范化前的副本。
    local_special_points.clear();
    for (const auto& face : special_internal) collect_points(face);
    for (const auto& face : special_boundary) collect_points(face);
    for (auto& [background, points] : local_special_points) {
        static_cast<void>(background);
        std::sort(points.begin(), points.end(), [](const Vec3& lhs,
                                                   const Vec3& rhs) {
            if (lhs.x != rhs.x) return lhs.x < rhs.x;
            if (lhs.y != rhs.y) return lhs.y < rhs.y;
            return lhs.z < rhs.z;
        });
        std::vector<Vec3> unique;
        unique.reserve(points.size());
        for (const auto point : points) {
            if (unique.empty() ||
                norm(unique.back() - point) > topology_point_tolerance) {
                unique.push_back(point);
            }
        }
        points = std::move(unique);
    }
    rebuild_local_point_indices();

    std::vector<std::uint8_t> near_special(
        static_cast<std::size_t>(grid.cell_count()), 0U);
    for (const auto& [background, points] : local_special_points) {
        static_cast<void>(points);
        for (const auto neighbor : neighboring_cells(grid, background)) {
            near_special[static_cast<std::size_t>(neighbor)] = 1U;
        }
    }

    ActiveLatticePoints active(grid.point_count());
    std::vector<Vec3> regular_vertices;
    regular_vertices.reserve(4U);
    for (std::uint64_t background = 0; background < grid.cell_count();
         ++background) {
        const auto state = mesh.state(background);
        if (state != CompactCellState::full_fluid) continue;
        const CellKey key = grid.cell_key(background);
        const std::array<bool, 6> exterior = {
            key.i == 0U, key.i + 1U == grid.nx(),
            key.j == 0U, key.j + 1U == grid.ny(),
            key.k == 0U, key.k + 1U == grid.nz()};
        const std::array<std::uint64_t, 6> neighbor = {
            background - static_cast<std::uint64_t>(key.i != 0U),
            background + static_cast<std::uint64_t>(key.i + 1U < grid.nx()),
            background - static_cast<std::uint64_t>(grid.nx()) *
                             static_cast<std::uint64_t>(key.j != 0U),
            background + static_cast<std::uint64_t>(grid.nx()) *
                             static_cast<std::uint64_t>(key.j + 1U < grid.ny()),
            background - static_cast<std::uint64_t>(grid.nx()) * grid.ny() *
                             static_cast<std::uint64_t>(key.k != 0U),
            background + static_cast<std::uint64_t>(grid.nx()) * grid.ny() *
                             static_cast<std::uint64_t>(key.k + 1U < grid.nz())};
        for (std::uint8_t face = 0; face < 6U; ++face) {
            if (!exterior[face] &&
                mesh.state(neighbor[face]) != CompactCellState::full_fluid) {
                // Interfaces involving an explicit cell were already inserted
                // into special_internal and are registered below.
                continue;
            }
            set_regular_face_vertices(regular_vertices, grid, background,
                                      face);
            for (const auto point : regular_vertices) {
                const auto lattice =
                    lattice_point(grid, point, length_tolerance);
                if (!lattice) {
                    throw std::runtime_error(
                        "规则 Cartesian 面顶点未落在背景格点上");
                }
                active.mark(lattice->linear);
            }
        }
    }
    OffLatticePointWelder off_lattice(grid.domain().minimum(),
                                      point_merge_tolerance);
    const auto register_vertex = [&](const Vec3& point) {
        const auto lattice = lattice_point(grid, point, length_tolerance);
        if (lattice) active.mark(lattice->linear);
        else static_cast<void>(off_lattice.insert(point));
    };
    for (const auto& face : special_internal) {
        for (const auto point : face.vertices) register_vertex(point);
    }
    for (const auto& face : special_boundary) {
        for (const auto point : face.vertices) register_vertex(point);
    }
    active.finalize();
    stats.point_count = active.active_count() + off_lattice.points().size();
    if (stats.point_count > no_cell) {
        throw std::overflow_error("OpenFOAM label=32 无法表示点数量");
    }

    const auto resolve_vertex = [&](const Vec3& point) -> std::uint32_t {
        const auto lattice = lattice_point(grid, point, length_tolerance);
        std::uint64_t id = 0;
        if (lattice) id = active.id(lattice->linear);
        else {
            id = active.active_count() + off_lattice.insert(point);
        }
        if (id > no_cell) throw std::overflow_error("OpenFOAM 点 ID 越界");
        return static_cast<std::uint32_t>(id);
    };

    const auto canonicalize_face = [&](OutputFace& face) {
        std::vector<Vec3> unique;
        std::vector<std::uint32_t> ids;
        unique.reserve(face.vertices.size());
        ids.reserve(face.vertices.size());
        for (const auto point : face.vertices) {
            const auto id = resolve_vertex(point);
            if (std::find(ids.begin(), ids.end(), id) != ids.end()) continue;
            ids.push_back(id);
            unique.push_back(point);
        }
        face.vertices = std::move(unique);
        if (face.vertices.size() < 3U ||
            polygon_area(face.vertices) <= area_tolerance) {
            throw std::runtime_error(
                "OpenFOAM 面在最终点编号规范化后退化：owner=" +
                std::to_string(face.owner) + " neighbour=" +
                std::to_string(face.neighbor));
        }
    };
    for (auto& face : special_internal) canonicalize_face(face);
    for (auto& face : special_boundary) canonicalize_face(face);

    struct EdgeUse {
        std::uint32_t count{};
        std::uint32_t outward_start{};
        std::uint32_t outward_end{};
        Vec3 outward_start_point{};
        Vec3 outward_end_point{};
    };
    using EdgeKey = std::pair<std::uint32_t, std::uint32_t>;
    using CellEdges = std::map<EdgeKey, EdgeUse>;
    const auto collect_cut_cell_edges = [&]() {
        std::map<std::uint32_t, CellEdges> result;
        const auto add_for_cell = [&](const OutputFace& face,
                                      std::uint32_t cell,
                                      bool reverse) {
            if (cell == no_cell ||
                solver_cell_is_cut[static_cast<std::size_t>(cell)] == 0U) {
                return;
            }
            auto& edges = result[cell];
            for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
                const Vec3 first_point = face.vertices[edge];
                const Vec3 second_point =
                    face.vertices[(edge + 1U) % face.vertices.size()];
                const std::uint32_t first = resolve_vertex(first_point);
                const std::uint32_t second = resolve_vertex(second_point);
                const EdgeKey key = first < second
                                        ? EdgeKey{first, second}
                                        : EdgeKey{second, first};
                auto& use = edges[key];
                if (use.count == 0U) {
                    use.outward_start = reverse ? second : first;
                    use.outward_end = reverse ? first : second;
                    use.outward_start_point =
                        reverse ? second_point : first_point;
                    use.outward_end_point =
                        reverse ? first_point : second_point;
                }
                ++use.count;
            }
        };
        for (const auto& face : special_internal) {
            add_for_cell(face, face.owner, false);
            add_for_cell(face, face.neighbor, true);
        }
        for (const auto& face : special_boundary) {
            add_for_cell(face, face.owner, false);
        }
        return result;
    };

    auto cut_cell_edges = collect_cut_cell_edges();
    const auto describe_cut_cell_edge =
        [&](std::uint32_t cell, const EdgeKey& target) {
            std::ostringstream description;
            const auto inspect = [&](const OutputFace& face) {
                if (face.owner != cell && face.neighbor != cell) return;
                bool contains = false;
                std::vector<std::uint32_t> ids;
                ids.reserve(face.vertices.size());
                for (const auto point : face.vertices) {
                    ids.push_back(resolve_vertex(point));
                }
                for (std::size_t edge = 0; edge < ids.size(); ++edge) {
                    const std::uint32_t first = ids[edge];
                    const std::uint32_t second =
                        ids[(edge + 1U) % ids.size()];
                    const EdgeKey key = first < second
                                            ? EdgeKey{first, second}
                                            : EdgeKey{second, first};
                    if (key == target) {
                        contains = true;
                        break;
                    }
                }
                if (!contains) return;
                description << " [owner=" << face.owner
                            << " ownerBackground="
                            << solver_cell_background[
                                   static_cast<std::size_t>(face.owner)]
                            << " neighbor=" << face.neighbor
                            << " neighborBackground=";
                if (face.neighbor == no_cell) {
                    description << "boundary";
                } else {
                    description << solver_cell_background[
                                       static_cast<std::size_t>(
                                           face.neighbor)];
                }
                description
                            << " kind="
                            << static_cast<unsigned>(face.boundary_kind)
                            << " anchor=" << face.anchor_background
                            << " area=" << precise(polygon_area(face.vertices))
                            << " ids=";
                for (const auto id : ids) description << id << ',';
                description << ']';
            };
            for (const auto& face : special_internal) inspect(face);
            for (const auto& face : special_boundary) inspect(face);
            return description.str();
        };
    std::vector<OutputFace> topology_seals;
    for (const auto& [cell, edges] : cut_cell_edges) {
        struct DirectedGapEdge {
            std::uint32_t end{};
            Vec3 start_point{};
            Vec3 end_point{};
        };
        std::map<std::uint32_t, DirectedGapEdge> outgoing;
        std::map<std::uint32_t, std::uint32_t> incoming_count;
        std::uint64_t gap_edge_count = 0U;
        for (const auto& [key, use] : edges) {
            if (use.count > 2U) {
                throw std::runtime_error(
                    "Cut-cell 导出边关联超过 2：cell=" +
                    std::to_string(cell) + " background=" +
                    std::to_string(solver_cell_background[
                        static_cast<std::size_t>(cell)]) +
                    " edge=" + std::to_string(key.first) + "," +
                    std::to_string(key.second) + " count=" +
                    std::to_string(use.count) +
                    describe_cut_cell_edge(cell, key));
            }
            if (use.count == 2U) continue;
            ++gap_edge_count;
            const std::uint32_t start = use.outward_end;
            const std::uint32_t end = use.outward_start;
            const auto [position, inserted] = outgoing.emplace(
                start, DirectedGapEdge{end, use.outward_end_point,
                                       use.outward_start_point});
            static_cast<void>(position);
            if (!inserted) {
                throw std::runtime_error(
                    "Cut-cell 缺口边不能构成简单有向环：cell=" +
                    std::to_string(cell));
            }
            ++incoming_count[end];
        }
        if (outgoing.empty()) continue;
        for (const auto& [start, edge] : outgoing) {
            static_cast<void>(edge);
            if (incoming_count[start] != 1U) {
                throw std::runtime_error(
                    "Cut-cell 缺口顶点入度不是 1：cell=" +
                    std::to_string(cell));
            }
        }
        std::map<std::uint32_t, bool> used;
        for (const auto& [loop_start, first_edge] : outgoing) {
            static_cast<void>(first_edge);
            if (used[loop_start]) continue;
            std::vector<Vec3> vertices;
            std::uint32_t current = loop_start;
            for (std::size_t step = 0; step <= outgoing.size(); ++step) {
                const auto found = outgoing.find(current);
                if (found == outgoing.end() || used[current]) {
                    throw std::runtime_error(
                        "Cut-cell 缺口边未形成闭环：cell=" +
                        std::to_string(cell));
                }
                used[current] = true;
                vertices.push_back(found->second.start_point);
                current = found->second.end;
                if (current == loop_start) break;
                if (step == outgoing.size()) {
                    throw std::runtime_error(
                        "Cut-cell 缺口环追踪超限：cell=" +
                        std::to_string(cell));
                }
            }
            if (vertices.size() < 3U ||
                polygon_area(vertices) <= area_tolerance) {
                throw std::runtime_error(
                    "Cut-cell 缺口环无法形成正面积封口：cell=" +
                    std::to_string(cell));
            }
            OutputFace seal{
                std::move(vertices), cell, no_cell, 0U,
                solver_cell_background[static_cast<std::size_t>(cell)],
                BoundaryKind::wall};
            canonicalize_face(seal);
            const double seal_area = polygon_area(seal.vertices);
            ++stats.topology_sealed_loop_count;
            stats.topology_sealed_edge_count += gap_edge_count;
            stats.topology_sealed_loop_area += seal_area;
            stats.maximum_topology_sealed_loop_area =
                std::max(stats.maximum_topology_sealed_loop_area, seal_area);
            topology_seals.push_back(std::move(seal));
            gap_edge_count = 0U;
        }
    }
    for (auto& seal : topology_seals) {
        special_boundary.push_back(std::move(seal));
    }
    cut_cell_edges = collect_cut_cell_edges();
    for (const auto& [cell, edges] : cut_cell_edges) {
        for (const auto& [key, use] : edges) {
            static_cast<void>(key);
            if (use.count != 2U) {
                throw std::runtime_error(
                    "Cut-cell 导出拓扑封口后仍存在非二关联边：cell=" +
                    std::to_string(cell));
            }
        }
    }


    const auto regular_fluid = [&](std::uint64_t background) {
        const auto state = mesh.state(background);
        return state == CompactCellState::full_fluid;
    };

    const auto split_regular = [&](OutputFace& face) {
        if (near_special[static_cast<std::size_t>(face.anchor_background)] == 0U) {
            return;
        }
        std::vector<const LocalPointIndex*> candidates;
        for (const auto cell :
             neighboring_cells(grid, face.anchor_background)) {
            const auto found = local_special_point_indices.find(cell);
            if (found != local_special_point_indices.end()) {
                candidates.push_back(&found->second);
            }
        }
        split_face_edges(face, candidates, topology_point_tolerance);
        canonicalize_face(face);
    };

    std::uint64_t regular_internal_count = 0;
    std::uint64_t regular_boundary_count = 0;
    const auto for_each_regular_internal = [&](auto&& callback) {
        OutputFace face;
        face.boundary_kind = BoundaryKind::none;
        face.vertices.reserve(16U);
        for (std::uint32_t k = 0; k < grid.nz(); ++k) {
            for (std::uint32_t j = 0; j < grid.ny(); ++j) {
                for (std::uint32_t i = 0; i < grid.nx(); ++i) {
                    const auto id = grid.linear_id({0, i, j, k});
                    const std::array<std::pair<bool, std::pair<std::uint64_t,
                                                               std::uint8_t>>, 3>
                        neighbors = {{
                            {i + 1U < grid.nx(), {id + 1U, 1U}},
                            {j + 1U < grid.ny(),
                             {id + grid.nx(), 3U}},
                            {k + 1U < grid.nz(),
                             {id + static_cast<std::uint64_t>(grid.nx()) *
                                       grid.ny(),
                              5U}},
                        }};
                    if (!regular_fluid(id)) continue;
                    for (const auto& neighbor : neighbors) {
                        if (!neighbor.first ||
                            !regular_fluid(neighbor.second.first)) {
                            continue;
                        }
                        set_regular_face_vertices(
                            face.vertices, grid, id, neighbor.second.second);
                        face.owner = cell_start[static_cast<std::size_t>(id)];
                        face.neighbor = cell_start[static_cast<std::size_t>(
                            neighbor.second.first)];
                        face.boundary_id = 0U;
                        face.anchor_background = id;
                        split_regular(face);
                        callback(face);
                    }
                }
            }
        }
    };
    for_each_regular_internal(
        [&](const OutputFace&) { ++regular_internal_count; });

    const auto for_each_regular_boundary = [&](auto&& callback) {
        OutputFace face;
        face.neighbor = no_cell;
        face.boundary_id = 0U;
        face.boundary_kind = BoundaryKind::farfield;
        face.vertices.reserve(16U);
        for (std::uint64_t id = 0; id < grid.cell_count(); ++id) {
            if (!regular_fluid(id)) continue;
            const CellKey key = grid.cell_key(id);
            const std::array<bool, 6> exterior = {
                key.i == 0U, key.i + 1U == grid.nx(),
                key.j == 0U, key.j + 1U == grid.ny(),
                key.k == 0U, key.k + 1U == grid.nz()};
            for (std::uint8_t face_id = 0; face_id < 6U; ++face_id) {
                if (!exterior[face_id]) continue;
                set_regular_face_vertices(face.vertices, grid, id, face_id);
                face.owner = cell_start[static_cast<std::size_t>(id)];
                face.anchor_background = id;
                split_regular(face);
                callback(face);
            }
        }
    };
    for_each_regular_boundary(
        [&](const OutputFace&) { ++regular_boundary_count; });

    std::vector<std::uint8_t> secondary_topology_cell(
        static_cast<std::size_t>(next_cell), 0U);
    const auto mark_secondary_cell = [&](std::uint32_t cell) {
        if (cell != no_cell &&
            solver_cell_is_cut[static_cast<std::size_t>(cell)] == 0U) {
            secondary_topology_cell[static_cast<std::size_t>(cell)] = 1U;
        }
    };
    for (const auto& face : special_internal) {
        mark_secondary_cell(face.owner);
        mark_secondary_cell(face.neighbor);
    }
    for (const auto& face : special_boundary) {
        mark_secondary_cell(face.owner);
    }
    const auto collect_secondary_edges = [&]() {
        std::map<std::uint32_t, CellEdges> result;
        const auto add_for_cell = [&](const OutputFace& face,
                                      std::uint32_t cell,
                                      bool reverse) {
            if (cell == no_cell ||
                secondary_topology_cell[static_cast<std::size_t>(cell)] ==
                    0U) {
                return;
            }
            auto& edges = result[cell];
            for (std::size_t edge = 0; edge < face.vertices.size(); ++edge) {
                const Vec3 first_point = face.vertices[edge];
                const Vec3 second_point =
                    face.vertices[(edge + 1U) % face.vertices.size()];
                const std::uint32_t first = resolve_vertex(first_point);
                const std::uint32_t second = resolve_vertex(second_point);
                const EdgeKey key = first < second
                                        ? EdgeKey{first, second}
                                        : EdgeKey{second, first};
                auto& use = edges[key];
                if (use.count == 0U) {
                    use.outward_start = reverse ? second : first;
                    use.outward_end = reverse ? first : second;
                    use.outward_start_point =
                        reverse ? second_point : first_point;
                    use.outward_end_point =
                        reverse ? first_point : second_point;
                }
                ++use.count;
            }
        };
        for (const auto& face : special_internal) {
            add_for_cell(face, face.owner, false);
            add_for_cell(face, face.neighbor, true);
        }
        for (const auto& face : special_boundary) {
            add_for_cell(face, face.owner, false);
        }
        for_each_regular_internal([&](const OutputFace& face) {
            add_for_cell(face, face.owner, false);
            add_for_cell(face, face.neighbor, true);
        });
        for_each_regular_boundary([&](const OutputFace& face) {
            add_for_cell(face, face.owner, false);
        });
        return result;
    };

    auto secondary_edges = collect_secondary_edges();
    std::vector<OutputFace> secondary_seals;
    for (const auto& [cell, edges] : secondary_edges) {
        struct DirectedGapEdge {
            std::uint32_t end{};
            Vec3 start_point{};
            Vec3 end_point{};
        };
        std::map<std::uint32_t, DirectedGapEdge> outgoing;
        std::map<std::uint32_t, std::uint32_t> incoming_count;
        std::uint64_t gap_edge_count = 0U;
        for (const auto& [key, use] : edges) {
            static_cast<void>(key);
            if (use.count > 2U) {
                throw std::runtime_error(
                    "Cut-cell 邻域规则单元边关联超过 2：cell=" +
                    std::to_string(cell));
            }
            if (use.count == 2U) continue;
            ++gap_edge_count;
            const std::uint32_t start = use.outward_end;
            const std::uint32_t end = use.outward_start;
            const auto [position, inserted] = outgoing.emplace(
                start, DirectedGapEdge{end, use.outward_end_point,
                                       use.outward_start_point});
            static_cast<void>(position);
            if (!inserted) {
                throw std::runtime_error(
                    "Cut-cell 邻域规则单元缺口不是简单环：cell=" +
                    std::to_string(cell));
            }
            ++incoming_count[end];
        }
        if (outgoing.empty()) continue;
        for (const auto& [start, edge] : outgoing) {
            static_cast<void>(edge);
            if (incoming_count[start] != 1U) {
                throw std::runtime_error(
                    "Cut-cell 邻域规则单元缺口入度不是 1：cell=" +
                    std::to_string(cell));
            }
        }
        std::map<std::uint32_t, bool> used;
        for (const auto& [loop_start, first_edge] : outgoing) {
            static_cast<void>(first_edge);
            if (used[loop_start]) continue;
            std::vector<Vec3> vertices;
            std::uint32_t current = loop_start;
            for (std::size_t step = 0; step <= outgoing.size(); ++step) {
                const auto found = outgoing.find(current);
                if (found == outgoing.end() || used[current]) {
                    throw std::runtime_error(
                        "Cut-cell 邻域规则单元缺口未闭合：cell=" +
                        std::to_string(cell));
                }
                used[current] = true;
                vertices.push_back(found->second.start_point);
                current = found->second.end;
                if (current == loop_start) break;
                if (step == outgoing.size()) {
                    throw std::runtime_error(
                        "Cut-cell 邻域规则单元缺口追踪超限：cell=" +
                        std::to_string(cell));
                }
            }
            if (vertices.size() < 3U ||
                polygon_area(vertices) <= area_tolerance) {
                throw std::runtime_error(
                    "Cut-cell 邻域规则单元缺口封口非正面积：cell=" +
                    std::to_string(cell));
            }
            OutputFace seal{
                std::move(vertices), cell, no_cell, 0U,
                solver_cell_background[static_cast<std::size_t>(cell)],
                BoundaryKind::wall};
            canonicalize_face(seal);
            const double seal_area = polygon_area(seal.vertices);
            ++stats.topology_sealed_loop_count;
            stats.topology_sealed_edge_count += gap_edge_count;
            stats.topology_sealed_loop_area += seal_area;
            stats.maximum_topology_sealed_loop_area =
                std::max(stats.maximum_topology_sealed_loop_area, seal_area);
            secondary_seals.push_back(std::move(seal));
            gap_edge_count = 0U;
        }
    }
    for (auto& seal : secondary_seals) {
        special_boundary.push_back(std::move(seal));
    }
    secondary_edges = collect_secondary_edges();
    for (const auto& [cell, edges] : secondary_edges) {
        for (const auto& [key, use] : edges) {
            static_cast<void>(key);
            if (use.count != 2U) {
                throw std::runtime_error(
                    "Cut-cell 邻域规则单元封口后仍有非二关联边：cell=" +
                    std::to_string(cell));
            }
        }
    }

    const auto for_each_internal = [&](auto&& callback) {
        std::size_t special = 0;
        const auto precedes = [](const OutputFace& lhs,
                                 const OutputFace& rhs) {
            if (lhs.owner != rhs.owner) return lhs.owner < rhs.owner;
            return lhs.neighbor < rhs.neighbor;
        };
        for_each_regular_internal([&](const OutputFace& regular) {
            while (special < special_internal.size() &&
                   precedes(special_internal[special], regular)) {
                callback(special_internal[special++]);
            }
            callback(regular);
        });
        while (special < special_internal.size()) {
            callback(special_internal[special++]);
        }
    };
    const auto for_each_boundary = [&](auto&& callback) {
        for_each_regular_boundary(callback);
        for (const auto& face : special_boundary) {
            if (face.boundary_kind == BoundaryKind::farfield) callback(face);
        }
        for (const auto& face : special_boundary) {
            if (face.boundary_kind == BoundaryKind::wall) callback(face);
        }
    };
    const auto for_each_face = [&](auto&& callback) {
        for_each_internal(callback);
        for_each_boundary(callback);
    };

    stats.internal_face_count =
        regular_internal_count + special_internal.size();
    stats.boundary_face_count =
        regular_boundary_count + special_boundary.size();
    stats.face_count = stats.internal_face_count + stats.boundary_face_count;
    if (stats.face_count >= no_cell) {
        throw std::overflow_error("OpenFOAM label=32 无法表示面数量");
    }
    for_each_face([&](const OutputFace& face) {
        stats.face_vertex_reference_count += face.vertices.size();
    });
    if (stats.face_vertex_reference_count > static_cast<std::uint64_t>(
                                                std::numeric_limits<std::int32_t>::max())) {
        throw std::overflow_error("faceCompactList 顶点引用超过 label=32 范围");
    }


    std::uint64_t topology_hash = 14695981039346656037ULL;
    for_each_face([&](const OutputFace& face) {
        hash_u64(topology_hash, face.owner);
        hash_u64(topology_hash, face.neighbor);
        hash_u64(topology_hash, face.boundary_id);
        hash_u64(topology_hash, static_cast<std::uint64_t>(face.boundary_kind));
        hash_u64(topology_hash, face.vertices.size());
        for (const auto point : face.vertices) {
            hash_u64(topology_hash, resolve_vertex(point));
        }
    });
    stats.topology_hash_fnv1a64 = topology_hash;
    const auto preparation_end = Clock::now();
    stats.preparation_seconds = seconds(total_start, preparation_end);

    const auto poly_mesh = case_directory / "constant" / "polyMesh";
    std::filesystem::create_directories(poly_mesh);
    std::filesystem::create_directories(case_directory / "system");
    {
        std::ofstream output(case_directory / "system" / "controlDict",
                             std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 controlDict");
        output << "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
               << "    class dictionary;\n    object controlDict;\n}\n\n"
               << "application checkMesh;\nstartFrom startTime;\nstartTime 0;\n"
               << "stopAt endTime;\nendTime 1;\ndeltaT 1;\n"
               << "writeControl timeStep;\nwriteInterval 1;\n"
               << "writeFormat binary;\nrunTimeModifiable false;\n";
    }
    {
        std::ofstream output(case_directory / "system" / "fvSchemes",
                             std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 fvSchemes");
        output << "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
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
        if (!output) throw std::runtime_error("无法写入 fvSolution");
        output << "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
               << "    class dictionary;\n    object fvSolution;\n}\n\n"
               << "solvers {}\n";
    }
    {
        std::ofstream output(poly_mesh / "points",
                             std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 points");
        write_header(output, "vectorField", "points", true);
        output << stats.point_count << "\n(";
        const auto& bits = active.bits();
        for (std::size_t word = 0; word < bits.size(); ++word) {
            std::uint64_t value = bits[word];
            while (value != 0U) {
                const auto bit = static_cast<std::uint32_t>(std::countr_zero(value));
                const std::uint64_t linear =
                    static_cast<std::uint64_t>(word) * 64U + bit;
                const Vec3 point = lattice_coordinate(grid, linear);
                output.write(reinterpret_cast<const char*>(&point.x),
                             sizeof(double));
                output.write(reinterpret_cast<const char*>(&point.y),
                             sizeof(double));
                output.write(reinterpret_cast<const char*>(&point.z),
                             sizeof(double));
                value &= value - 1U;
            }
        }
        for (const auto& point : off_lattice.points()) {
            output.write(reinterpret_cast<const char*>(&point.x), sizeof(double));
            output.write(reinterpret_cast<const char*>(&point.y), sizeof(double));
            output.write(reinterpret_cast<const char*>(&point.z), sizeof(double));
        }
        output << ")\n";
        if (!output) throw std::runtime_error("points 二进制写入失败");
    }
    {
        std::ofstream output(poly_mesh / "faces",
                             std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 faces");
        write_header(output, "faceCompactList", "faces", true);
        output << stats.face_count + 1U << "\n(";
        std::uint64_t offset = 0;
        write_label(output, offset);
        for_each_face([&](const OutputFace& face) {
            offset += face.vertices.size();
            write_label(output, offset);
        });
        output << ")\n" << stats.face_vertex_reference_count << "\n(";
        for_each_face([&](const OutputFace& face) {
            for (const auto point : face.vertices) {
                write_label(output, resolve_vertex(point));
            }
        });
        output << ")\n";
        if (!output) throw std::runtime_error("faces 二进制写入失败");
    }
    {
        std::ofstream output(poly_mesh / "owner",
                             std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 owner");
        write_header(output, "labelList", "owner", true);
        output << stats.face_count << "\n(";
        for_each_face([&](const OutputFace& face) {
            write_label(output, face.owner);
        });
        output << ")\n";
        if (!output) throw std::runtime_error("owner 二进制写入失败");
    }
    {
        std::ofstream output(poly_mesh / "neighbour",
                             std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 neighbour");
        write_header(output, "labelList", "neighbour", true);
        output << stats.internal_face_count << "\n(";
        for_each_internal([&](const OutputFace& face) {
            write_label(output, face.neighbor);
        });
        output << ")\n";
        if (!output) throw std::runtime_error("neighbour 二进制写入失败");
    }
    {
        struct Patch {
            std::string name;
            std::string type;
            std::uint64_t count{};
        };
        std::map<std::uint64_t, std::string> names;
        for (const auto& [id, name] : boundary_names) {
            names[id] = foam_name(name);
        }
        std::vector<Patch> patches;
        const std::uint64_t farfield_count = regular_boundary_count +
            static_cast<std::uint64_t>(std::count_if(
                special_boundary.begin(), special_boundary.end(),
                [](const OutputFace& face) {
                    return face.boundary_kind == BoundaryKind::farfield;
                }));
        if (farfield_count != 0U) {
            patches.push_back({"farfield", "patch", farfield_count});
        }
        for (const auto& face : special_boundary) {
            if (face.boundary_kind != BoundaryKind::wall) continue;
            const std::string name = names.contains(face.boundary_id)
                                         ? names[face.boundary_id]
                                         : "boundary_" +
                                               std::to_string(face.boundary_id);
            if (patches.empty() || patches.back().name != name) {
                patches.push_back({name, "wall", 0U});
            }
            ++patches.back().count;
        }
        std::ofstream output(poly_mesh / "boundary", std::ios::trunc);
        if (!output) throw std::runtime_error("无法写入 boundary");
        write_header(output, "polyBoundaryMesh", "boundary", false);
        output << patches.size() << "\n(\n";
        std::uint64_t start = stats.internal_face_count;
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

    const std::array<std::filesystem::path, 8> files = {
        case_directory / "system" / "controlDict",
        case_directory / "system" / "fvSchemes",
        case_directory / "system" / "fvSolution",
        poly_mesh / "points", poly_mesh / "faces", poly_mesh / "owner",
        poly_mesh / "neighbour", poly_mesh / "boundary"};
    for (const auto& file : files) stats.written_bytes += std::filesystem::file_size(file);
    const auto write_end = Clock::now();
    stats.writing_seconds = seconds(preparation_end, write_end);
    stats.total_seconds = seconds(total_start, write_end);
    return stats;
}

} // 命名空间 cartmesh
