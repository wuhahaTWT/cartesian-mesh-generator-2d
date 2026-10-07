#pragma once
#include "cartmesh2d/topology/Topology2D.hpp"
#include <algorithm>
#include <map>
#include <utility>

namespace cartmesh2d::test {
inline TopologyMesh2D fromPolygons(const std::vector<Polygon2D>& polygons) {
    TopologyMesh2D topology;
    std::map<std::pair<double, double>, std::size_t> vertices;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edges;
    for (const auto& polygon : polygons) {
        TopologyCell2D cell;
        cell.id = topology.cells.size();
        cell.geometryArea = polygon.area();
        for (const auto point : polygon.vertices) {
            const auto key = std::make_pair(point.x, point.y);
            const auto [it, inserted] = vertices.emplace(key, topology.vertices.size());
            if (inserted) topology.vertices.push_back({it->second, point});
            cell.vertices.push_back(it->second);
        }
        for (std::size_t i = 0; i < cell.vertices.size(); ++i) {
            const auto a = cell.vertices[i];
            const auto b = cell.vertices[(i + 1) % cell.vertices.size()];
            const auto key = std::minmax(a, b);
            const auto [it, inserted] = edges.emplace(key, topology.edges.size());
            if (inserted) {
                topology.edges.push_back({it->second, a, b, cell.id, {},
                                          BoundaryPatch2D::DomainBoundary});
            } else {
                auto& edge = topology.edges[it->second];
                edge.neighbour = cell.id;
                edge.patch = BoundaryPatch2D::None;
            }
            cell.edges.push_back(it->second);
        }
        topology.cells.push_back(std::move(cell));
    }
    return topology;
}
}
