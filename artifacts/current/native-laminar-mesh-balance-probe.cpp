#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// Native mesh-repair sensitivity probe.  It changes no product default and
// solves no surrogate equation: an accepted output is produced only through
// the existing exact-union/split repair and the authoritative Solver gate.
int main(int argc, char** argv) {
    using namespace cartmesh2d;
    if (argc < 5) {
        throw std::runtime_error(
            "usage: mesh-balance-probe input.solver.cm2d boundary.xy "
            "minimum-volume-ratio output.solver.cm2d [mutable-cell ...]");
    }

    const auto input = readCm2dTopology(argv[1]);
    if (!input.valid()) throw std::runtime_error(input.error);

    std::ifstream boundaryInput(argv[2]);
    if (!boundaryInput) throw std::runtime_error("cannot read boundary input");
    std::vector<BoundaryLoop> loops;
    std::vector<Point2D> points;
    const auto finishLoop = [&]() {
        if (points.empty()) return;
        if (points.size() < 3) throw std::runtime_error("boundary loop has fewer than three points");
        loops.emplace_back(points);
        points.clear();
    };
    std::string line;
    while (std::getline(boundaryInput, line)) {
        std::istringstream row(line);
        Point2D point;
        if (!(row >> point.x >> point.y)) {
            finishLoop();
            continue;
        }
        points.push_back(point);
    }
    finishLoop();
    BoundaryRegion2D boundary(std::move(loops));
    if (!boundary.diagnose().valid() || !boundary.normalizeAlternating())
        throw std::runtime_error("invalid boundary region");

    auto bounds = boundary.bounds();
    bool hasDomainBoundary = false;
    Point2D domainMinimum{std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::infinity()};
    Point2D domainMaximum{-std::numeric_limits<double>::infinity(),
                          -std::numeric_limits<double>::infinity()};
    for (const auto& edge : input.topology.edges) {
        if (edge.patch != BoundaryPatch2D::DomainBoundary) continue;
        hasDomainBoundary = true;
        for (const auto vertex : {edge.v0, edge.v1}) {
            const auto point = input.topology.vertices.at(vertex).point;
            domainMinimum.x = std::min(domainMinimum.x, point.x);
            domainMinimum.y = std::min(domainMinimum.y, point.y);
            domainMaximum.x = std::max(domainMaximum.x, point.x);
            domainMaximum.y = std::max(domainMaximum.y, point.y);
        }
    }
    if (!hasDomainBoundary) {
        const double span = std::max(bounds.max.x - bounds.min.x,
                                     bounds.max.y - bounds.min.y);
        const double padding = 0.15 * span;
        domainMinimum = {bounds.min.x - padding, bounds.min.y - padding};
        domainMaximum = {bounds.max.x + padding, bounds.max.y + padding};
    }
    const Domain2D domain{{domainMinimum, domainMaximum}};

    SolverQualityPolicy2D target;
    target.minVolumeRatio = std::stod(argv[3]);
    std::vector<bool> immutable;
    if (argc > 5) {
        immutable.assign(input.topology.cells.size(), true);
        for (int arg = 5; arg < argc; ++arg) {
            const auto cell = static_cast<std::size_t>(std::stoull(argv[arg]));
            if (cell >= immutable.size()) throw std::runtime_error("mutable cell is out of range");
            immutable[cell] = false;
        }
    }
    const auto repaired = improveSolverForTargetPolicy2D(
        input.topology, domain, boundary, immutable, target);
    if (!repaired.valid()) {
        std::cerr << "repair failed";
        for (const auto& issue : repaired.issues) std::cerr << ": " << issue;
        std::cerr << '\n';
        return 2;
    }

    const auto originalQuality = evaluateSolverQuality2D(input.topology);
    const auto defaultQuality = evaluateSolverQuality2D(repaired.topology);
    const auto targetQuality = evaluateSolverQuality2D(repaired.topology, target);
    std::string error;
    if (!writeCm2dTopology(repaired.topology, argv[4], &error))
        throw std::runtime_error(error);

    const auto totalArea = [](const TopologyMesh2D& mesh) {
        double area = 0.0;
        for (const auto& cell : mesh.cells) area += cell.geometryArea;
        return area;
    };
    const auto boundaryEdges = [](const TopologyMesh2D& mesh) {
        return static_cast<std::size_t>(std::count_if(
            mesh.edges.begin(), mesh.edges.end(),
            [](const Edge2D& edge) { return !edge.neighbour.has_value(); }));
    };
    const auto boundaryAtoms = [](const TopologyMesh2D& mesh) {
        using Atom = std::tuple<double, double, double, double, int>;
        std::vector<Atom> atoms;
        for (const auto& edge : mesh.edges) {
            if (edge.neighbour) continue;
            auto a = mesh.vertices.at(edge.v0).point;
            auto b = mesh.vertices.at(edge.v1).point;
            if (std::tie(b.x, b.y) < std::tie(a.x, a.y)) std::swap(a, b);
            atoms.emplace_back(a.x, a.y, b.x, b.y, static_cast<int>(edge.patch));
        }
        std::sort(atoms.begin(), atoms.end());
        return atoms;
    };
    const auto minimumArea = [](const TopologyMesh2D& mesh) {
        double area = std::numeric_limits<double>::infinity();
        for (const auto& cell : mesh.cells) area = std::min(area, cell.geometryArea);
        return area;
    };

    std::cout << std::setprecision(17)
              << "{\"inputCells\":" << input.topology.cells.size()
              << ",\"outputCells\":" << repaired.topology.cells.size()
              << ",\"repartitionCount\":" << repaired.repartitionCount
              << ",\"inputArea\":" << totalArea(input.topology)
              << ",\"outputArea\":" << totalArea(repaired.topology)
              << ",\"inputMinimumCellArea\":" << minimumArea(input.topology)
              << ",\"outputMinimumCellArea\":" << minimumArea(repaired.topology)
              << ",\"inputBoundaryEdges\":" << boundaryEdges(input.topology)
              << ",\"outputBoundaryEdges\":" << boundaryEdges(repaired.topology)
              << ",\"boundaryAtomsIdentical\":"
              << (boundaryAtoms(input.topology) == boundaryAtoms(repaired.topology)
                      ? "true" : "false")
              << ",\"originalDefaultPass\":" << (originalQuality.valid() ? "true" : "false")
              << ",\"outputDefaultPass\":" << (defaultQuality.valid() ? "true" : "false")
              << ",\"outputTargetPass\":" << (targetQuality.valid() ? "true" : "false")
              << ",\"inputMinimumVolumeRatio\":" << originalQuality.minVolumeRatio
              << ",\"outputMinimumVolumeRatio\":" << defaultQuality.minVolumeRatio
              << "}\n";
}
