#pragma once
#include "cartmesh2d/topology/Topology2D.hpp"
#include <vector>

namespace cartmesh2d::fv {
// All areas/lengths are the final 2D geometry, never background-grid values.
struct Cell { Point2D centre; double area = 0; std::vector<std::size_t> faces; };
struct Face {
    std::size_t owner = 0;
    std::optional<std::size_t> neighbour;
    BoundaryPatch2D patch = BoundaryPatch2D::None;
    Point2D centre;
    Vector2D areaVector; // edge length times outward unit normal of owner
    Vector2D correction; // S - transmissibility * (C_neighbour - C_owner)
    double transmissibility = 0;
    double neighbourWeight = 0;
};
struct FvMesh2D {
    std::vector<Cell> cells;
    std::vector<Face> faces;
    double maxClosureError = 0; // |sum outward S| / cell perimeter
};
// Throws for inconsistent incidence/geometry or the unchanged Solver quality gate.
// The serialized AUDIT record alone is not trusted.
[[nodiscard]] FvMesh2D makeFvMesh2D(const TopologyMesh2D&, const TolerancePolicy& = {});
// Safety checks on factory-created caches; not a replacement for makeFvMesh2D.
void validateFvMesh2D(const FvMesh2D&);
[[nodiscard]] std::vector<Vector2D> reconstructGradient(
    const FvMesh2D&, const std::vector<double>& cellValues,
    const std::vector<double>& boundaryValues);
}
