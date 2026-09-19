#include "cartmesh2d/fv/ThermalCheckpoint2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {
int failures = 0;

void check(bool ok, const std::string& message) {
    if (!ok) { ++failures; std::cerr << "FAIL: " << message << '\n'; }
}

template<class F>
void rejects(F&& action, const std::string& expected, const std::string& message) {
    try {
        action();
        check(false, message + " (accepted)");
    } catch (const std::exception& error) {
        check(std::string(error.what()).find(expected) != std::string::npos,
              message + " (wrong reason: " + error.what() + ")");
    }
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
            const auto a = cell.vertices[i];
            const auto b = cell.vertices[(i + 1) % cell.vertices.size()];
            auto [it, inserted] = edges.emplace(std::minmax(a, b), result.edges.size());
            if (inserted) result.edges.push_back({it->second, a, b, cell.id, {}, BoundaryPatch2D::DomainBoundary});
            else { result.edges[it->second].neighbour = cell.id; result.edges[it->second].patch = BoundaryPatch2D::None; }
            cell.edges.push_back(it->second);
        }
        result.cells.push_back(std::move(cell));
    }
    return result;
}

FvMesh2D cavityMesh(int n = 3) {
    auto point = [n](int i, int j) { return Point2D{static_cast<double>(i) / n, static_cast<double>(j) / n}; };
    std::vector<Polygon2D> polygons;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            polygons.push_back({{point(i,j), point(i+1,j), point(i+1,j+1), point(i,j+1)}});
    return makeFvMesh2D(fromPolygons(polygons));
}

FlowControls2D flowControls() {
    FlowControls2D c;
    c.scenario = "cavity";
    c.nu = .1;
    c.speed = 1.;
    c.tolerance = 1e-8;
    c.maxIterations = 250;
    return c;
}

ThermalSetup2D setup(const FvMesh2D& mesh, double source = 2.) {
    ThermalSetup2D s;
    s.diffusivity = .1;
    s.sourceDensity.assign(mesh.cells.size(), source);
    s.boundary.resize(mesh.faces.size());
    for (auto& b : s.boundary) b = {ScalarBoundaryKind2D::DiffusiveFlux, 0., std::nullopt};
    return s;
}

ScalarTransportControls2D scalarControls() {
    ScalarTransportControls2D c;
    c.maxCorrections = 100;
    c.relativeTolerance = 1e-9;
    c.absoluteTolerance = 1e-12;
    c.cellTolerance = 1e-9;
    return c;
}

ThermalFlowState2D initial(const FvMesh2D& mesh, double value = 5.) {
    const auto fc = flowControls();
    auto state = ThermalFlowState2D{initialIncompressibleState2D(mesh, fc), {}};
    state.scalar.assign(mesh.cells.size(), value);
    return state;
}

bool sameBits(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    return std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

void compareState(const ThermalFlowState2D& a, const ThermalFlowState2D& b, const std::string& label) {
    check(a.flow.time == b.flow.time, label + " time is bitwise identical");
    check(sameBits(a.flow.u, b.flow.u), label + " u is bitwise identical");
    check(sameBits(a.flow.v, b.flow.v), label + " v is bitwise identical");
    check(sameBits(a.flow.p, b.flow.p), label + " p is bitwise identical");
    check(sameBits(a.flow.flux, b.flow.flux), label + " flux is bitwise identical");
    check(sameBits(a.scalar, b.scalar), label + " scalar is bitwise identical");
}

ThermalFlowState2D advanceAccepted(const FvMesh2D& mesh, const FlowControls2D& fc,
                                   const ThermalSetup2D& setup, const ScalarTransportControls2D& sc,
                                   ThermalFlowState2D state, int steps, double dt) {
    for (int i = 0; i < steps; ++i) {
        const auto result = advanceThermalFlow2D(mesh, fc, setup, sc, state, dt);
        check(result.flow.converged, "thermal flow step converges");
        check(result.scalar.converged, "thermal scalar step converges");
        check(result.accepted.has_value(), "thermal step has accepted coupled state");
        if (!result.accepted) break;
        state = *result.accepted;
    }
    return state;
}

void uniformSourceAndEvolution() {
    const auto mesh = cavityMesh();
    const auto fc = flowControls();
    const auto setupData = setup(mesh);
    const auto sc = scalarControls();
    const auto before = initial(mesh);
    auto state = before;
    for (int step = 0; step < 4; ++step) {
        const auto result = advanceThermalFlow2D(mesh, fc, setupData, sc, state, .01);
        check(result.accepted.has_value(), "uniform source cavity step accepted");
        if (!result.accepted) return;
        state = *result.accepted;
        for (double value : state.scalar) check(std::abs(value - (5. + 2. * .01 * (step + 1))) < 1e-8,
                                                 "uniform source gives exact scalar rise");
    }
    check(state.flow.time == .04, "thermal clock reaches four physical steps");
    check(!sameBits(state.flow.u, before.flow.u) || !sameBits(state.flow.v, before.flow.v), "cavity velocity evolves");
    check(sameBits(before.scalar, std::vector<double>(mesh.cells.size(), 5.)), "input scalar remains unchanged");
}

void restartMatchesContinuous() {
    const auto mesh = cavityMesh();
    const auto fc = flowControls(); const auto setupData = setup(mesh); const auto sc = scalarControls();
    const auto start = initial(mesh);
    const auto continuous = advanceAccepted(mesh, fc, setupData, sc, start, 4, .01);
    const auto two = advanceAccepted(mesh, fc, setupData, sc, start, 2, .01);
    std::stringstream checkpoint;
    writeThermalCheckpoint2D(checkpoint, mesh, fc, setupData, sc, two);
    const auto restored = readThermalCheckpoint2D(checkpoint, mesh, fc, setupData, sc);
    const auto restarted = advanceAccepted(mesh, fc, setupData, sc, restored, 2, .01);
    compareState(continuous, restarted, "2+checkpoint+2 versus continuous 4");
}

void transientOutletInflow() {
    const auto mesh=cavityMesh(8);
    auto fc=flowControls();fc.scenario="channel";fc.maxIterations=600;
    fc.outletBackflow=OutletBackflow2D::NormalInlet;
    fc.convection=ConvectionScheme2D::LimitedLinearUpwind;
    auto data=setup(mesh,0.);const auto sc=scalarControls();
    ThermalFlowState2D start{initialIncompressibleState2D(mesh,fc),{}};
    start.scalar.assign(mesh.cells.size(),5.);
    const auto profile=[](double y) {return 1.+2.*std::cos(2.*std::acos(-1.)*y);};
    // Divergence-free initial parallel flow, distinct from the new-time
    // parabolic inlet. The right outlet initially has both signs of flux.
    for(std::size_t i=0;i<mesh.cells.size();++i) start.flow.u[i]=profile(mesh.cells[i].centre.y);
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id];start.flow.flux[id]=profile(f.centre.y)*f.areaVector.x;
        data.boundary[id].inflowValue=f.centre.x==1.?7.:5.;
    }
    auto reject=fc;reject.outletBackflow=OutletBackflow2D::Reject;
    rejects([&] {(void)advanceThermalFlow2D(mesh,reject,data,sc,start,.001);},
            "backflow unsupported","default mode rejects true reverse outlet");
    const auto step=advanceThermalFlow2D(mesh,fc,data,sc,start,.001);
    check(step.accepted.has_value(),"transient reverse outlet step accepted");
    if(!step.accepted) return;
    check(step.flow.outletBackflowFaces>0 && step.flow.outletInflow>0,"accepted step retains outlet inflow");
    check(*std::max_element(step.scalar.values.begin(),step.scalar.values.end())>5.,"incoming outlet scalar enters domain");
    for(double value:step.scalar.values) check(value>=5.-1e-9 && value<=7.+1e-9,"reverse-flow thermal field remains bounded");
    // The global defect is a sum of cell residuals, bounded by sqrt(N)||r||2;
    // a separate fixed dimensional threshold would ignore dt-dependent mass.
    check(std::abs(step.scalar.globalBalance)<=std::sqrt(static_cast<double>(mesh.cells.size()))*
          step.scalar.history.back().residualNorm+1e-12,"reverse-flow thermal storage and boundary flux balance");
    std::cout<<"reverse thermal: faces="<<step.flow.outletBackflowFaces
             <<" inflow="<<step.flow.outletInflow<<" globalBalance="<<step.scalar.globalBalance
             <<" cellScaled="<<step.scalar.history.back().maxDiagonalScaledImbalance<<'\n';
    auto missing=data;
    for(std::size_t id=0;id<mesh.faces.size();++id)
        if(mesh.faces[id].centre.x==1.) missing.boundary[id].inflowValue.reset();
    rejects([&] {(void)advanceThermalFlow2D(mesh,fc,missing,sc,start,.001);},
            "inflow","reverse-flow scalar requires explicit inflow value");
    const auto continuous=advanceAccepted(mesh,fc,data,sc,start,2,.001);
    std::stringstream saved;writeThermalCheckpoint2D(saved,mesh,fc,data,sc,*step.accepted);
    const auto restored=readThermalCheckpoint2D(saved,mesh,fc,data,sc);
    compareState(continuous,advanceAccepted(mesh,fc,data,sc,restored,1,.001),"reverse-flow joint restart");
}

void failureDoesNotMutateInputs() {
    const auto mesh = cavityMesh(); const auto fc = flowControls(); const auto setupData = setup(mesh); const auto sc = scalarControls();
    const auto before = initial(mesh);
    auto shortFlow = fc; shortFlow.maxIterations = 1;
    const auto flowFailure = advanceThermalFlow2D(mesh, shortFlow, setupData, sc, before, .01);
    check(!flowFailure.flow.converged && flowFailure.scalar.history.empty() && !flowFailure.accepted,
          "flow failure stops before scalar and has no accepted state");
    compareState(before, initial(mesh), "flow-failure input");

    auto heated = setupData;
    for (std::size_t id = 0; id < mesh.faces.size(); ++id)
        if (!mesh.faces[id].neighbour) heated.boundary[id] = {ScalarBoundaryKind2D::Value, 1., std::nullopt};
    auto one = sc; one.maxCorrections = 1;
    const auto scalarFailure = advanceThermalFlow2D(mesh, fc, heated, one, before, .01);
    check(scalarFailure.flow.converged && !scalarFailure.scalar.converged && !scalarFailure.accepted,
          "scalar correction limit rejects nonuniform heated boundary");
    check(sameBits(before.scalar, std::vector<double>(mesh.cells.size(), 5.)), "scalar failure leaves previous scalar unchanged");

    bool called = false;
    rejects([&] { (void)advanceThermalFlow2D(mesh, fc, setupData, sc, before, .01,
        [&](const FlowIteration2D&) { called = true; throw std::runtime_error("cancelled"); }); },
        "cancelled", "progress cancellation propagates");
    check(called, "progress callback was invoked");
    check(sameBits(before.flow.u, initial(mesh).flow.u) && sameBits(before.scalar, initial(mesh).scalar),
          "cancellation leaves input unchanged");
}

std::string checkpointText(const FvMesh2D& mesh, const FlowControls2D& fc,
                           const ThermalSetup2D& setupData, const ScalarTransportControls2D& sc,
                           const ThermalFlowState2D& state) {
    std::stringstream s; writeThermalCheckpoint2D(s, mesh, fc, setupData, sc, state); return s.str();
}

void strictCheckpointValidation() {
    const auto mesh = cavityMesh(); const auto fc = flowControls(); const auto setupData = setup(mesh); const auto sc = scalarControls();
    const auto state = advanceAccepted(mesh, fc, setupData, sc, initial(mesh), 1, .01);
    const auto valid = checkpointText(mesh, fc, setupData, sc, state);
    auto expect = [&](std::string text, const FvMesh2D& m, const FlowControls2D& f,
                      const ThermalSetup2D& s, const ScalarTransportControls2D& c,
                      const std::string& reason, const std::string& label) {
        rejects([&] { std::stringstream in(text); (void)readThermalCheckpoint2D(in, m, f, s, c); }, reason, label);
    };
    auto changedD = setupData; changedD.diffusivity = .2;
    expect(valid, mesh, fc, changedD, sc, "thermal diffusivity", "changed diffusivity rejected");
    auto changedSource = setupData; changedSource.sourceDensity[0] += 1.;
    expect(valid, mesh, fc, changedSource, sc, "thermal source", "changed source rejected");
    auto changedBoundary = setupData; changedBoundary.boundary[0].value = 1.;
    expect(valid, mesh, fc, changedBoundary, sc, "thermal boundary value", "changed boundary value rejected");
    auto changedKind = setupData; changedKind.boundary[0].kind = ScalarBoundaryKind2D::Value;
    expect(valid, mesh, fc, changedKind, sc, "thermal boundary type/inflow mismatch", "changed boundary kind rejected");
    auto changedInflow = setupData; changedInflow.boundary[0].inflowValue = 3.;
    expect(valid, mesh, fc, changedInflow, sc, "thermal boundary type/inflow mismatch", "changed inflow rejected");
    auto changedConvection = sc; changedConvection.convection = ConvectionScheme2D::LimitedLinearUpwind;
    expect(valid, mesh, fc, setupData, changedConvection, "thermal convection mismatch", "changed scalar convection rejected");
    auto changedFlow = fc; changedFlow.nu = .2;
    expect(valid, mesh, changedFlow, setupData, sc, "incompatible configuration", "changed flow configuration rejected");
    auto alteredGeometry = valid; const auto marker = alteredGeometry.find("CELL 0 ");
    alteredGeometry.replace(marker + 7, 4, "0.125");
    expect(alteredGeometry, mesh, fc, setupData, sc, "cell centre", "changed geometry rejected");
    expect(valid.substr(0, valid.size() - 8), mesh, fc, setupData, sc, "truncated", "truncated checkpoint rejected");
    auto nonfinite = valid; const auto sourceMarker = nonfinite.find("SOURCES ");
    const auto nonfiniteEnd = nonfinite.find('\n', sourceMarker);
    nonfinite.replace(sourceMarker, nonfiniteEnd - sourceMarker, "SOURCES " + std::to_string(mesh.cells.size()) + " nan");
    expect(nonfinite, mesh, fc, setupData, sc, "thermal source", "nonfinite checkpoint rejected");
    expect(valid + " EXTRA", mesh, fc, setupData, sc, "trailing data", "trailing checkpoint data rejected");
    std::stringstream dtStream(valid);
    const auto restored = readThermalCheckpoint2D(dtStream, mesh, fc, setupData, sc);
    auto changedIterationLimits = fc; changedIterationLimits.maxIterations = 17;
    auto changedCorrections = sc; changedCorrections.maxCorrections = 7;
    std::stringstream controlsStream(valid);
    const auto controlsRestored = readThermalCheckpoint2D(controlsStream, mesh, changedIterationLimits,
                                                            setupData, changedCorrections);
    check(controlsRestored.flow.time == state.flow.time, "checkpoint permits changed iteration controls");
    const auto changedDtResult = advanceThermalFlow2D(mesh, fc, setupData, sc, restored, .005);
    check(changedDtResult.accepted.has_value(), "restart permits changed physical time step");
}

void indexedAndCallbackRepresentations() {
    const auto mesh = cavityMesh();
    ScalarTransportProblem2D indexed;
    indexed.diffusivity = .1; indexed.volumeFlux.assign(mesh.faces.size(), 0.); indexed.sourceDensity.assign(mesh.cells.size(), 2.); indexed.boundaryData.resize(mesh.faces.size());
    for (auto& b : indexed.boundaryData) b = {ScalarBoundaryKind2D::Value, 0., std::nullopt};
    auto callback = indexed; callback.sourceDensity.clear(); callback.boundaryData.clear();
    callback.source = [](Point2D) { return 2.; };
    callback.boundary = [](std::size_t, const Face&) { return ScalarBoundary2D{ScalarBoundaryKind2D::Value, 0., std::nullopt}; };
    const auto a = solveScalarTransport2D(mesh, indexed); const auto b = solveScalarTransport2D(mesh, callback);
    check(a.converged && b.converged && sameBits(a.values, b.values), "indexed and callback scalar data are equivalent");
    auto ambiguous = indexed; ambiguous.source = [](Point2D) { return 2.; };
    rejects([&] { (void)solveScalarTransport2D(mesh, ambiguous); }, "exactly one valid source representation", "ambiguous source representation rejected");
    auto wrong = indexed; wrong.sourceDensity.pop_back();
    rejects([&] { (void)solveScalarTransport2D(mesh, wrong); }, "exactly one valid source representation", "wrong source dimension rejected");
    auto ambiguousBoundary = indexed; ambiguousBoundary.boundary = [](std::size_t, const Face&) { return ScalarBoundary2D{}; };
    rejects([&] { (void)solveScalarTransport2D(mesh, ambiguousBoundary); }, "exactly one valid boundary representation", "ambiguous boundary representation rejected");
}
}

int main() {
    try {
        uniformSourceAndEvolution(); restartMatchesContinuous(); failureDoesNotMutateInputs();
        transientOutletInflow();
        strictCheckpointValidation(); indexedAndCallbackRepresentations();
    } catch (const std::exception& error) {
        std::cerr << "unexpected thermal test exception: " << error.what() << '\n';
        return 1;
    }
    std::cout << "thermal flow failures: " << failures << '\n';
    return failures == 0 ? 0 : 1;
}
