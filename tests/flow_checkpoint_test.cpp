#include "cartmesh2d/fv/FlowCheckpoint2D.hpp"

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {

TopologyMesh2D topology() {
    TopologyMesh2D mesh;
    mesh.vertices = {{0, {0, 0}}, {1, {1, 0}}, {2, {1, 1}}, {3, {0, 1}}};
    mesh.edges = {{0, 0, 1, 0, {}, BoundaryPatch2D::DomainBoundary},
                 {1, 1, 2, 0, {}, BoundaryPatch2D::DomainBoundary},
                 {2, 2, 3, 0, {}, BoundaryPatch2D::DomainBoundary},
                 {3, 3, 0, 0, {}, BoundaryPatch2D::DomainBoundary}};
    TopologyCell2D cell;
    cell.id = 0;
    cell.geometryArea = 1.0;
    cell.vertices = {0, 1, 2, 3};
    cell.edges = {0, 1, 2, 3};
    mesh.cells.push_back(cell);
    return mesh;
}

FvMesh2D mesh() { return makeFvMesh2D(topology()); }

FlowState2D state() {
    FlowState2D result;
    result.time = 1.0 / 3.0;
    result.u = {1.25}; result.v = {-2.5}; result.p = {3.75}; result.flux = {0.0, 1.0, 0.0, -1.0};
    return result;
}

std::string serialized(const FvMesh2D& m, const FlowControls2D& c, const FlowState2D& s) {
    std::ostringstream out;
    writeFlowCheckpoint2D(out, m, c, s);
    return out.str();
}

template<class Function>
void rejects(Function&& function) {
    bool rejected = false;
    try { function(); } catch (const std::exception&) { rejected = true; }
    if (!rejected) throw std::runtime_error("expected checkpoint input to be rejected");
}

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main() {
    const auto m = mesh();
    FlowControls2D controls;
    controls.scenario = "manufactured";
    controls.nu = .125;
    controls.speed = 2.0;
    controls.convection = ConvectionScheme2D::LimitedLinearUpwind;
    controls.viscousStress = ViscousStress2D::Laplacian;
    controls.manufacturedPressureSlope = -0.75;
    const auto original = state();
    const auto text = serialized(m, controls, original);
    std::ostringstream fixed;
    fixed << std::fixed << std::setprecision(2);
    writeFlowCheckpoint2D(fixed, m, controls, original);
    std::istringstream fixedInput(fixed.str());
    const auto fixedRestored = readFlowCheckpoint2D(fixedInput, m, controls);
    require(fixedRestored.u == original.u && fixedRestored.flux == original.flux,
            "fixed-format writer did not roundtrip");

    // Restart policy ignores iteration/relaxation controls while retaining the
    // physical/discretization identity that affects the state.
    auto relaxed = controls;
    relaxed.tolerance = 1e-10;
    relaxed.maxIterations = 17;
    relaxed.velocityRelaxation = .3;
    relaxed.pressureRelaxation = .8;
    std::istringstream roundtrip(text);
    const auto restored = readFlowCheckpoint2D(roundtrip, m, relaxed);
    require(restored.time == original.time && restored.u == original.u && restored.v == original.v &&
            restored.p == original.p && restored.flux == original.flux, "checkpoint roundtrip mismatch");

    auto incompatible = controls;
    incompatible.nu = .126;
    rejects([&] { std::istringstream in(text); (void)readFlowCheckpoint2D(in, m, incompatible); });
    rejects([&] { std::istringstream in(text.substr(0, text.size() - 5)); (void)readFlowCheckpoint2D(in, m, controls); });
    rejects([&] { std::istringstream in(text + " EXTRA\n"); (void)readFlowCheckpoint2D(in, m, controls); });

    auto reordered = text;
    const auto cell = reordered.find("CELL 0");
    require(cell != std::string::npos, "missing cell record");
    reordered.replace(cell, 6, "CELL 1");
    rejects([&] { std::istringstream in(reordered); (void)readFlowCheckpoint2D(in, m, controls); });

    auto invalidTime = text;
    const auto time = invalidTime.find("TIME ");
    require(time != std::string::npos, "missing time record");
    invalidTime.replace(time, 11, "TIME -1.0\n");
    rejects([&] { std::istringstream in(invalidTime); (void)readFlowCheckpoint2D(in, m, controls); });

    auto nonfinite = text;
    const auto vector = nonfinite.find("U 1 ");
    require(vector != std::string::npos, "missing vector record");
    nonfinite.replace(vector, 8, "U 1 nan\n");
    rejects([&] { std::istringstream in(nonfinite); (void)readFlowCheckpoint2D(in, m, controls); });

    auto badSize = text;
    const auto flux = badSize.find("FLUX 4 ");
    require(flux != std::string::npos, "missing flux record");
    badSize.replace(flux, 7, "FLUX 999999999999999999999 ");
    rejects([&] { std::istringstream in(badSize); (void)readFlowCheckpoint2D(in, m, controls); });
}
