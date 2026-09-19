#include "cartmesh2d/fv/WallDistance2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cstddef>
#include <limits>
#include <map>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

TopologyMesh2D fromPolygons(const std::vector<Polygon2D>& polygons) {
    TopologyMesh2D result;
    std::map<std::pair<double, double>, std::size_t> vertices;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edges;
    for (const auto& polygon : polygons) {
        TopologyCell2D cell;
        cell.id = result.cells.size();
        cell.geometryArea = polygon.area();
        for (const auto& point : polygon.vertices) {
            const auto key = std::make_pair(point.x, point.y);
            auto [it, inserted] = vertices.emplace(key, result.vertices.size());
            if (inserted) result.vertices.push_back({it->second, point});
            cell.vertices.push_back(it->second);
        }
        for (std::size_t i = 0; i < cell.vertices.size(); ++i) {
            const auto a = cell.vertices[i], b = cell.vertices[(i + 1) % cell.vertices.size()];
            auto [it, inserted] = edges.emplace(std::minmax(a, b), result.edges.size());
            if (inserted) result.edges.push_back({it->second, a, b, cell.id, {}, BoundaryPatch2D::DomainBoundary});
            else { result.edges[it->second].neighbour = cell.id; result.edges[it->second].patch = BoundaryPatch2D::None; }
            cell.edges.push_back(it->second);
        }
        result.cells.push_back(cell);
    }
    return result;
}

FvMesh2D meshFromPolygons(const std::vector<Polygon2D>& polygons) {
    return makeFvMesh2D(fromPolygons(polygons));
}

FvMesh2D rectangleGrid(int nx, int ny) {
    std::vector<Polygon2D> polygons;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const double x0 = static_cast<double>(i) / nx, x1 = static_cast<double>(i + 1) / nx;
        const double y0 = static_cast<double>(j) / ny, y1 = static_cast<double>(j + 1) / ny;
        polygons.push_back({{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}}});
    }
    return meshFromPolygons(polygons);
}

std::pair<Point2D, Point2D> endpoints(const Face& face) {
    const Vector2D half{-face.areaVector.y * .5, face.areaVector.x * .5};
    return {{face.centre.x - half.x, face.centre.y - half.y},
            {face.centre.x + half.x, face.centre.y + half.y}};
}

double segmentDistance(Point2D point, Point2D a, Point2D b) {
    const Vector2D ab{b.x - a.x, b.y - a.y};
    const Vector2D ap{point.x - a.x, point.y - a.y};
    const double length2 = ab.x * ab.x + ab.y * ab.y;
    const double t = length2 == 0 ? 0 : std::clamp((ap.x * ab.x + ap.y * ab.y) / length2, 0., 1.);
    const Point2D closest{a.x + t * ab.x, a.y + t * ab.y};
    return std::hypot(point.x - closest.x, point.y - closest.y);
}

std::pair<double, std::size_t> brute(Point2D point, const FvMesh2D& mesh,
                                     const std::vector<bool>& wallFaces) {
    double best = std::numeric_limits<double>::infinity();
    std::size_t bestId = mesh.faces.size();
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) if (wallFaces[id]) {
        const auto [a, b] = endpoints(mesh.faces[id]);
        const double value = segmentDistance(point, a, b);
        if (value < best || (value == best && id < bestId)) { best = value; bestId = id; }
    }
    return {best, bestId};
}

void bottomWallAndMultipleWalls() {
    const auto mesh = rectangleGrid(1, 2);
    std::vector<bool> bottom(mesh.faces.size(), false);
    bottom[0] = true;
    const auto result = computeWallDistance2D(mesh, bottom);
    require(result.distance.size() == 2 && result.nearestFace.size() == 2, "bottom wall result dimensions");
    require(std::abs(result.distance[0] - .25) < 1e-14 && std::abs(result.distance[1] - .75) < 1e-14,
            "bottom wall distances use cell centres");
    std::vector<bool> sides(mesh.faces.size(), false);
    sides[0] = sides[3] = sides[4] = true;
    const auto multiple = computeWallDistance2D(mesh, sides);
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        const auto expected = brute(mesh.cells[i].centre, mesh, sides);
        require(std::abs(multiple.distance[i] - expected.first) < 1e-14 &&
                multiple.nearestFace[i] == expected.second, "multiple wall brute-force match");
    }
}

void slantedEndpointAndTie() {
    const auto slanted = meshFromPolygons({{{{0., 0.}, {.2, 0.}, {1.4, 1.}, {0., 1.}}}});
    std::vector<bool> endpointMask(slanted.faces.size(), false);
    endpointMask[0] = true;
    const auto endpointResult = computeWallDistance2D(slanted, endpointMask);
    const auto endpointExpected = brute(slanted.cells[0].centre, slanted, endpointMask);
    require(std::abs(endpointResult.distance[0] - endpointExpected.first) < 1e-14 && endpointResult.nearestFace[0] == 0,
            "slanted wall endpoint projection matches independent segment distance");
    std::vector<bool> slantedMask(slanted.faces.size(), false);
    slantedMask[1] = true;
    const auto slantedResult = computeWallDistance2D(slanted, slantedMask);
    const double knownDistance = segmentDistance(slanted.cells[0].centre, {.2, 0.}, {1.4, 1.});
    require(std::abs(slantedResult.distance[0] - knownDistance) < 1e-14 && slantedResult.nearestFace[0] == 1,
            "slanted wall interior projection matches known original endpoints");

    const auto square = rectangleGrid(1, 1);
    std::vector<bool> tie(square.faces.size(), false);
    tie[0] = tie[2] = true;
    const auto tied = computeWallDistance2D(square, tie);
    require(std::abs(tied.distance[0] - .5) < 1e-14 && tied.nearestFace[0] == 0,
            "equal-distance walls choose the smaller face ID");
}

void bvhMatchesBruteForce() {
    const auto mesh = rectangleGrid(4, 3);
    std::vector<bool> walls(mesh.faces.size(), false);
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        walls[id] = !mesh.faces[id].neighbour;
    const auto result = computeWallDistance2D(mesh, walls);
    const auto wallCount = std::count(walls.begin(), walls.end(), true);
    require(wallCount > 8 && result.segmentTests > 0, "large wall set exercises BVH segment testing");
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        const auto expected = brute(mesh.cells[i].centre, mesh, walls);
        require(std::abs(result.distance[i] - expected.first) < 1e-14 &&
                result.nearestFace[i] == expected.second, "BVH result matches brute force");
    }
}

void rejectsInvalidMasks() {
    const auto mesh = rectangleGrid(2, 1);
    std::vector<bool> empty(mesh.faces.size(), false);
    bool rejected = false;
    try { (void)computeWallDistance2D(mesh, empty); } catch (const std::exception&) { rejected = true; }
    require(rejected, "all-false wall mask is rejected");
    rejected = false;
    try { (void)computeWallDistance2D(mesh, std::vector<bool>(mesh.faces.size() - 1, true)); }
    catch (const std::exception&) { rejected = true; }
    require(rejected, "wall mask size mismatch is rejected");
    const auto internal = std::find_if(mesh.faces.begin(), mesh.faces.end(), [](const Face& face) {
        return face.neighbour.has_value();
    });
    require(internal != mesh.faces.end(), "test mesh has an internal face");
    std::vector<bool> internalMask(mesh.faces.size(), false);
    internalMask[static_cast<std::size_t>(internal - mesh.faces.begin())] = true;
    rejected = false;
    try { (void)computeWallDistance2D(mesh, internalMask); } catch (const std::exception&) { rejected = true; }
    require(rejected, "internal face marked as wall is rejected");
}

} // namespace

int main(int argc,char** argv) {
    try {
        bottomWallAndMultipleWalls();
        slantedEndpointAndTie();
        bvhMatchesBruteForce();
        rejectsInvalidMasks();
        if(argc!=1) {
            require(argc==3,"usage: wall_distance_tests [mesh.cm2d outputPrefix]");
            const auto read=readCm2dTopology(argv[1]);require(read.valid(),"invalid native mesh");
            const auto mesh=makeFvMesh2D(read.topology);
            std::vector<bool> walls(mesh.faces.size());
            for(std::size_t id=0;id<walls.size();++id)walls[id]=!mesh.faces[id].neighbour;
            const auto start=std::chrono::steady_clock::now();
            const auto result=computeWallDistance2D(mesh,walls);
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            std::ofstream out(std::string(argv[2])+".csv"),meta(std::string(argv[2])+".json");
            require(bool(out)&&bool(meta),"cannot open wall distance output");
            out<<std::setprecision(17)<<"cell,distance,nearestFace\n";
            for(std::size_t i=0;i<mesh.cells.size();++i)out<<i<<','<<result.distance[i]<<','<<result.nearestFace[i]<<'\n';
            meta<<std::setprecision(17)<<"{\"cells\":"<<mesh.cells.size()<<",\"wallSegments\":"<<std::count(walls.begin(),walls.end(),true)
                <<",\"segmentTests\":"<<result.segmentTests<<",\"distanceSeconds\":"<<seconds<<"}\n";
        }
    } catch (const std::exception& error) {
        std::cerr << "wall distance test failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
