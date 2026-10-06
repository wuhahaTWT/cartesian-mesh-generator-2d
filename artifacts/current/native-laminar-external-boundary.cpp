#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

// Reproduce the native `external` preset as an explicit custom-boundary file.
// This research adapter exists only so the product's custom-only Stokes mode
// can be compared on exactly the same final FvMesh2D.  It does not infer or
// repair geometry and deliberately leaves the public preset surface unchanged.
int main(int argc, char** argv) {
    using namespace cartmesh2d;
    using namespace cartmesh2d::fv;
    if (argc != 3) throw std::runtime_error(
        "usage: external-boundary mesh.solver.cm2d output.boundaries");
    const auto read = readCm2dTopology(argv[1]);
    if (!read.valid()) throw std::runtime_error(read.error);
    const auto mesh = makeFvMesh2D(read.topology);

    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -xmin;
    double ymin = xmin;
    double ymax = -xmin;
    for (const auto& face : mesh.faces) {
        if (face.neighbour) continue;
        xmin = std::min(xmin, face.centre.x);
        xmax = std::max(xmax, face.centre.x);
        ymin = std::min(ymin, face.centre.y);
        ymax = std::max(ymax, face.centre.y);
    }
    const double epsilon = TolerancePolicy{}.scale(std::max(xmax-xmin, ymax-ymin));
    const auto equal = [&](double a, double b) { return std::abs(a-b) <= epsilon; };

    FlowControls2D controls;
    controls.scenario = "custom";
    controls.speed = 1;
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        if (face.neighbour) continue;
        FlowBoundaryCondition2D condition;
        condition.face = id;
        if (face.patch == BoundaryPatch2D::EmbeddedBoundary) {
            condition.kind = FlowBoundaryKind2D::Wall;
            condition.name = "body";
        } else if (equal(face.centre.x, xmin)) {
            condition.kind = FlowBoundaryKind2D::VelocityInlet;
            condition.velocity = {1, 0};
            condition.name = "inlet";
        } else if (equal(face.centre.x, xmax)) {
            condition.kind = FlowBoundaryKind2D::PressureOutlet;
            condition.name = "outlet";
        } else if (equal(face.centre.y, ymin) || equal(face.centre.y, ymax)) {
            condition.kind = FlowBoundaryKind2D::Symmetry;
            condition.name = "farfield";
        } else {
            throw std::runtime_error("external boundary is not rectangular or embedded");
        }
        controls.boundaryConditions.push_back(std::move(condition));
    }
    std::ofstream output(argv[2]);
    if (!output) throw std::runtime_error("cannot open boundary output");
    writeFlowBoundaryConditions2D(output, mesh, controls);
}
