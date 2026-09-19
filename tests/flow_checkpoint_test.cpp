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
    require(text.find("CARTMESH2D_FLOW_CHECKPOINT 2\n") == 0, "writer did not emit checkpoint v2");
    std::ostringstream fixed;
    fixed << std::fixed << std::setprecision(2);
    writeFlowCheckpoint2D(fixed, m, controls, original);
    std::istringstream fixedInput(fixed.str());
    const auto fixedRestored = readFlowCheckpoint2D(fixedInput, m, controls);
    require(fixedRestored.u == original.u && fixedRestored.flux == original.flux,
            "fixed-format writer did not roundtrip");

    auto variable = controls;
    variable.faceViscosity = {.11, .17, .23, .31};
    const auto variableText = serialized(m, variable, original);
    require(variableText.find("CARTMESH2D_FLOW_CHECKPOINT 3\n") == 0,
            "face viscosity writer did not emit checkpoint v3");
    const auto viscosityLine = variableText.find("FACE_VISCOSITY 4 ");
    require(viscosityLine != std::string::npos, "v3 checkpoint omitted face viscosity field");
    std::istringstream variableInput(variableText);
    const auto variableRestored = readFlowCheckpoint2D(variableInput, m, variable);
    require(variableRestored.flux == original.flux && variableRestored.u == original.u,
            "v3 face viscosity checkpoint did not roundtrip");
    auto mismatchViscosity = variable;
    mismatchViscosity.faceViscosity[2] += .01;
    rejects([&] { std::istringstream in(variableText); (void)readFlowCheckpoint2D(in, m, mismatchViscosity); });
    auto missingViscosity = controls;
    rejects([&] { std::istringstream in(variableText); (void)readFlowCheckpoint2D(in, m, missingViscosity); });
    auto badViscosity = variableText;
    const auto valueStart = badViscosity.find("FACE_VISCOSITY 4 ") + std::string("FACE_VISCOSITY 4 ").size();
    badViscosity.replace(valueStart, 4, "nan");
    rejects([&] { std::istringstream in(badViscosity); (void)readFlowCheckpoint2D(in, m, variable); });
    auto badViscosityCount = variableText;
    badViscosityCount.replace(badViscosityCount.find("FACE_VISCOSITY 4"), 17, "FACE_VISCOSITY 3");
    rejects([&] { std::istringstream in(badViscosityCount); (void)readFlowCheckpoint2D(in, m, variable); });
    auto badWriter = controls;
    badWriter.faceViscosity = {.1};
    rejects([&] { std::ostringstream out; writeFlowCheckpoint2D(out, m, badWriter, original); });
    auto legacyWithField = controls;
    legacyWithField.faceViscosity = variable.faceViscosity;
    rejects([&] { std::istringstream in(text); (void)readFlowCheckpoint2D(in, m, legacyWithField); });

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

    // Version 1 remains readable and carries the historical reject-default policy.
    auto legacy = text;
    const auto header = legacy.find("CARTMESH2D_FLOW_CHECKPOINT 2");
    require(header != std::string::npos, "missing v2 header");
    legacy.replace(header, 28, "CARTMESH2D_FLOW_CHECKPOINT 1");
    const auto configEnd = legacy.find('\n', legacy.find("CONFIG "));
    require(configEnd != std::string::npos, "missing config line");
    const auto model = legacy.rfind(" ", configEnd - 1);
    require(model != std::string::npos, "missing backflow model");
    legacy.erase(model, configEnd - model);
    auto legacyControls = controls;
    legacyControls.outletBackflow = OutletBackflow2D::Reject;
    std::istringstream legacyInput(legacy);
    (void)readFlowCheckpoint2D(legacyInput, m, legacyControls);
    auto legacyNormalMismatch = controls;
    legacyNormalMismatch.outletBackflow = OutletBackflow2D::NormalInlet;
    rejects([&] { std::istringstream in(legacy); (void)readFlowCheckpoint2D(in, m, legacyNormalMismatch); });
    auto normalMismatch = controls;
    normalMismatch.outletBackflow = OutletBackflow2D::NormalInlet;
    rejects([&] { std::istringstream in(text); (void)readFlowCheckpoint2D(in, m, normalMismatch); });
    auto unknownModel = text;
    const auto modelText = unknownModel.find(" normal-inlet", unknownModel.find("CONFIG "));
    if (modelText != std::string::npos) unknownModel.replace(modelText + 1, 12, "unknown-model");
    else {
        const auto rejectText = unknownModel.find(" reject", unknownModel.find("CONFIG "));
        require(rejectText != std::string::npos, "missing backflow model token");
        unknownModel.replace(rejectText + 1, 6, "unknown-model");
    }
    rejects([&] { std::istringstream in(unknownModel); (void)readFlowCheckpoint2D(in, m, controls); });

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
