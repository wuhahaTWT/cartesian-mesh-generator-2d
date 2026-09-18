#pragma once

#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace cartmesh2d::fv::detail {

// Internal operators: the caller validates the mesh and field extents once.
// A single reconstructed value is stored per shared face, never per incidence.
inline std::vector<double> pressureFaceValues(
    const FvMesh2D& mesh, const std::vector<double>& p,
    const std::vector<Vector2D>& gradient, const std::vector<double>& boundary,
    const std::vector<bool>& fixed) {
    std::vector<double> result(mesh.faces.size());
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& f = mesh.faces[id];
        const auto i = f.owner;
        if (f.neighbour) {
            const auto j = *f.neighbour;
            const double w = f.neighbourWeight;
            const Point2D point{(1-w)*mesh.cells[i].centre.x+w*mesh.cells[j].centre.x,
                                (1-w)*mesh.cells[i].centre.y+w*mesh.cells[j].centre.y};
            const Vector2D g{(1-w)*gradient[i].x+w*gradient[j].x,
                             (1-w)*gradient[i].y+w*gradient[j].y};
            const auto d = f.centre-point;
            result[id] = (1-w)*p[i]+w*p[j]+g.x*d.x+g.y*d.y;
        } else if (fixed[id]) {
            result[id] = boundary[id];
        } else {
            const auto d = f.centre-mesh.cells[i].centre;
            result[id] = p[i]+gradient[i].x*d.x+gradient[i].y*d.y;
        }
    }
    return result;
}

inline std::vector<Vector2D> conservativePressureGradient(
    const FvMesh2D& mesh, const std::vector<double>& facePressure) {
    std::vector<Vector2D> result(mesh.cells.size());
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& f = mesh.faces[id];
        const Vector2D force{facePressure[id]*f.areaVector.x,
                             facePressure[id]*f.areaVector.y};
        result[f.owner].x += force.x;
        result[f.owner].y += force.y;
        if (f.neighbour) {
            result[*f.neighbour].x -= force.x;
            result[*f.neighbour].y -= force.y;
        }
    }
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i].x /= mesh.cells[i].area;
        result[i].y /= mesh.cells[i].area;
    }
    return result;
}

// Barth-Jespersen-style cell limiter, applied at every actual face centre.
// It bounds reconstructed face values by the cell's neighbour/Dirichlet stencil;
// it does not assert a maximum principle for the coupled Navier-Stokes solution.
inline std::vector<double> faceReconstructionLimiter(
    const FvMesh2D& mesh, const std::vector<double>& field,
    const std::vector<Vector2D>& gradient, const std::vector<double>& boundary,
    const std::vector<bool>& fixed) {
    std::vector<double> limiter(mesh.cells.size(), 1.);
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        double lo = field[i], hi = field[i];
        for (const auto id : mesh.cells[i].faces) {
            const auto& f = mesh.faces[id];
            if (f.neighbour) {
                const auto j = f.owner == i ? *f.neighbour : f.owner;
                lo = std::min(lo, field[j]); hi = std::max(hi, field[j]);
            } else if (fixed[id]) {
                lo = std::min(lo, boundary[id]); hi = std::max(hi, boundary[id]);
            }
        }
        for (const auto id : mesh.cells[i].faces) {
            const auto d = mesh.faces[id].centre-mesh.cells[i].centre;
            const double delta = gradient[i].x*d.x+gradient[i].y*d.y;
            if (delta > 0) limiter[i] = std::min(limiter[i], (hi-field[i])/delta);
            else if (delta < 0) limiter[i] = std::min(limiter[i], (lo-field[i])/delta);
        }
        limiter[i] = std::clamp(limiter[i], 0., 1.);
    }
    return limiter;
}

inline double upwindFaceValue(
    const FvMesh2D& mesh, std::size_t id, double flux,
    const std::vector<double>& field, const std::vector<Vector2D>& gradient,
    const std::vector<double>& limiter) {
    const auto& f = mesh.faces[id];
    const auto up = (f.neighbour && flux < 0) ? *f.neighbour : f.owner;
    if (limiter.empty()) return field[up];
    const auto d = f.centre-mesh.cells[up].centre;
    return field[up]+limiter[up]*(gradient[up].x*d.x+gradient[up].y*d.y);
}

} // namespace cartmesh2d::fv::detail
