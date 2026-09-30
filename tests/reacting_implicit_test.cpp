#include "cartmesh2d/fv/ReactingImplicit2D.hpp"
#include "FvTestMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::chemistry;
using namespace cartmesh2d::fv;
namespace {
void check(bool value, const char* text) { if (!value) throw std::runtime_error(text); }
std::vector<double> mixture(const DetailedGas& gas) {
    std::vector<double> x(gas.mechanism().species.size());
    for (const auto& [name, value] : {std::pair{"H2", 2.}, {"O2", 1.}, {"N2", 3.76}}) {
        const auto k = std::find(gas.mechanism().species.begin(), gas.mechanism().species.end(), name);
        check(k != gas.mechanism().species.end(), "missing species"); x[static_cast<std::size_t>(k - gas.mechanism().species.begin())] = value;
    }
    return x;
}
std::vector<ReactingBoundary2D> sealed(const FvMesh2D& mesh) {
    std::vector<ReactingBoundary2D> b;
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) if (!mesh.faces[i].neighbour)
        b.push_back({i, ReactingBoundaryKind2D::SlipWall, {}, {}, 0});
    return b;
}
std::vector<double> integral(const FvMesh2D& mesh, const ReactingState2D& state) {
    std::vector<double> result(state.cells[0].size());
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) for (std::size_t k = 0; k < result.size(); ++k)
        result[k] += mesh.cells[i].area * state.cells[i][k];
    return result;
}
void budget(const DetailedGas& gas, const FvMesh2D& mesh, const ReactingState2D& initial,
            const ReactingImplicitProgress2D& result) {
    const auto before = integral(mesh, initial), after = integral(mesh, result.lastAccepted);
    double speciesDefect = 0;
    for (std::size_t k = 4; k < before.size(); ++k)
        speciesDefect = std::max(speciesDefect, std::abs(after[k] - before[k] + result.boundaryImpulse[k]
            - result.chemistryChange[k] - result.constraintChange[k]) / before[0]);
    check(speciesDefect < 1e-8, "implicit species integral/quadrature mismatch");
    check(std::abs(after[0] - before[0] + result.boundaryImpulse[0]) / before[0] < 1e-8, "implicit mass budget failed");
    check(std::abs(after[3] - before[3] + result.boundaryImpulse[3]) / std::max(std::abs(before[3]), before[0]) < 1e-8,
          "implicit energy budget failed");
    const auto& mechanism = gas.mechanism();
    for (std::size_t m = 0; m < mechanism.elements.size(); ++m) {
        double physical = 0;
        for (std::size_t k = 0; k < mechanism.species.size(); ++k)
            physical += (after[k + 4] - before[k + 4] + result.boundaryImpulse[k + 4])
                * mechanism.atomCounts[k * mechanism.elements.size() + m] * mechanism.atomicWeights[m] / mechanism.molecularWeights[k];
        check(std::abs(physical / before[0]) < 1e-8, "implicit physical elemental budget failed");
    }
    std::cout << "quadrature species defect per initial mass=" << speciesDefect << '\n';
}
}
int main(int argc, char** argv) {
    try {
        check(argc == 2, "expected detailed H2 mechanism");
        DetailedGas gas(argv[1]); const auto x = mixture(gas);
        const auto mesh = fv_test::rectangle(2, 2, 1, true);
        ReactingFlowStepper2D solver(gas, mesh, sealed(mesh), {true, false});
        const auto q = gas.fromMoleAmounts(1100, 101325, x);
        const auto initial = solver.initialState(std::vector<ReactingConservative2D>(mesh.cells.size(), reactingConservative2D(q)));
        ReactingImplicitControls2D c; c.maximumStep = 1e-4; c.relativeTolerance = 1e-8; c.maximumBdfOrder = 5;
        ReactingImplicitIntegrator2D integrator(gas, solver, initial, c);
        const auto canceled = integrator.advance(1e-3, [] { return true; });
        check(canceled.canceled && canceled.lastAccepted.cells == initial.cells && canceled.lastAccepted.steps == 0,
              "cancel advanced the accepted state");
        const auto result = integrator.advance(1e-3);
        if (!result.reachedEnd) throw std::runtime_error(result.failure);
        check(result.lastAccepted.time == 1e-3 && result.lastAccepted.steps > 0, "wrong implicit endpoint");
        const auto reference = gas.advanceConstantVolume(q, 1e-3); check(reference.accepted.has_value(), "reference reactor failed");
        const auto referenceProperties = gas.properties(*reference.accepted);
        double difference = 0;
        for (const auto& u : result.lastAccepted.cells) {
            const auto p = reactingPrimitive2D(gas, u);
            difference = std::max(difference, std::abs(p.properties.temperature / referenceProperties.temperature - 1));
            for (std::size_t k = 0; k < x.size(); ++k)
                difference = std::max(difference, std::abs(p.properties.massFractions[k] - referenceProperties.massFractions[k]));
        }
        // Same homogeneous-source comparison used by the split integrator;
        // local integration accuracy only, not a mechanism/flame qualification.
        check(difference < 2e-6, "fully coupled BDF disagrees with homogeneous reactor");
        budget(gas, mesh, initial, result);
        std::cout << "full coupled BDF: " << result.internalSteps << " accepted, " << result.rhsCalls << " RHS, "
                  << result.rejectedRhsCalls << " rejected RHS, difference=" << difference << '\n';
        auto invalid = initial; invalid.cells[0][4] = -1e-30;
        bool rejected = false;
        try { ReactingImplicitIntegrator2D bad(gas, solver, invalid, c); } catch (const std::exception&) { rejected = true; }
        check(rejected && invalid.cells[0][4] == -1e-30, "invalid initial species repaired");
        c.maximumAcceptedSteps = 1;
        ReactingImplicitIntegrator2D limited(gas, solver, initial, c);
        const auto stopped = limited.advance(1e-3);
        check(!stopped.reachedEnd && !stopped.failure.empty() && stopped.lastAccepted.steps == 1, "step budget failure lost accepted state");
        const auto again = limited.advance(1e-3);
        check(again.lastAccepted.cells == stopped.lastAccepted.cells && again.lastAccepted.time == stopped.lastAccepted.time,
              "failed integration history was reused");
        check(initial.steps == 0 && initial.time == 0 && initial.cells[0] == reactingConservative2D(q), "implicit integration modified caller input");
        auto rhsLimit = c; rhsLimit.maximumResidualEvaluations = 1;
        ReactingImplicitIntegrator2D rhsLimited(gas, solver, initial, rhsLimit);
        std::size_t observed = 0;
        const auto rhsStopped = rhsLimited.advance(1e-3, {}, {}, [&](const auto& p, double) {
            observed = p.rhsCalls;
            check(p.lastAccepted.cells == initial.cells, "trial observer exposed unaccepted values");
        });
        check(observed == 1 && !rhsStopped.reachedEnd && rhsStopped.lastAccepted.cells == initial.cells
            && rhsStopped.failure.find("residual-evaluation budget") != std::string::npos,
            "startup evaluation budget did not retain the initial state");
        ReactingImplicitIntegrator2D interrupted(gas, solver, initial, c);
        std::size_t evaluations = 0;
        const auto trialCancel = interrupted.advance(1e-3, [&] { return evaluations >= 3; }, {},
            [&](const auto& p, double) { evaluations = p.rhsCalls; });
        check(trialCancel.canceled && !trialCancel.reachedEnd && trialCancel.lastAccepted.cells == initial.cells
            && trialCancel.failure.find("canceled during a trial") != std::string::npos,
            "trial cancellation lost the last accepted state");

        auto spatialMesh = fv_test::rectangle(4, 3, 1, true);
        constexpr double size = .002;
        for (auto& cell : spatialMesh.cells) { cell.centre.x *= size; cell.centre.y *= size; cell.area *= size * size; }
        for (auto& face : spatialMesh.faces) { face.centre.x *= size; face.centre.y *= size; face.areaVector = face.areaVector * size; face.correction = face.correction * size; }
        ReactingFlowStepper2D spatial(gas, spatialMesh, sealed(spatialMesh));
        std::vector<ReactingConservative2D> field;
        for (const auto& cell : spatialMesh.cells)
            field.push_back(reactingConservative2D(gas.fromMoleAmounts(1100 + 75 * std::sin(3.141592653589793 * cell.centre.x / size)
                * std::sin(3.141592653589793 * cell.centre.y / size), 101325, x)));
        const auto nonuniform = spatial.initialState(field);
        const auto withRate = spatial.evaluateResidual(nonuniform), withoutRate = spatial.evaluateResidual(nonuniform, 2, false);
        if (withRate.derivative != withoutRate.derivative || withRate.faceFlux != withoutRate.faceFlux) {
            for (std::size_t k = 0; k < field[0].size(); ++k) {
                double change = 0, magnitude = 0, faceChange = 0;
                for (std::size_t i = 0; i < field.size(); ++i) {
                    change = std::max(change, std::abs(withRate.derivative[i][k] - withoutRate.derivative[i][k]));
                    magnitude = std::max(magnitude, std::abs(withRate.derivative[i][k]));
                }
                for (std::size_t i = 0; i < withRate.faceFlux.size(); ++i)
                    faceChange = std::max(faceChange, std::abs(withRate.faceFlux[i][k] - withoutRate.faceFlux[i][k]));
                std::cerr << "rate toggle component=" << k << " derivative delta=" << change << " magnitude=" << magnitude
                          << " face delta=" << faceChange << '\n';
            }
        }
        check(withRate.derivative == withoutRate.derivative && withRate.faceFlux == withoutRate.faceFlux
            && withoutRate.transportRate.empty(), "optional rate estimation changed physical residual");
        c.maximumAcceptedSteps = 20000; c.maximumStep = 2e-7;
        ReactingImplicitIntegrator2D spatialIntegrator(gas, spatial, nonuniform, c);
        const auto spatialResult = spatialIntegrator.advance(2e-7);
        if (!spatialResult.reachedEnd) throw std::runtime_error(spatialResult.failure);
        budget(gas, spatialMesh, nonuniform, spatialResult);
        auto explicitState = nonuniform;
        ReactingStepControls2D explicitControls; explicitControls.maximumStep = 1e-8; explicitControls.endTime = 2e-7;
        while (explicitState.time < *explicitControls.endTime) {
            const auto step = spatial.advance(explicitState, explicitControls);
            if (!step.accepted) throw std::runtime_error(step.failure);
            explicitState = *step.accepted;
        }
        double temperatureDifference = 0;
        for (std::size_t i = 0; i < field.size(); ++i)
            temperatureDifference = std::max(temperatureDifference, std::abs(reactingPrimitive2D(gas, spatialResult.lastAccepted.cells[i]).properties.temperature
                - reactingPrimitive2D(gas, explicitState.cells[i]).properties.temperature));
        // Report the full spatial response separately; mesh/temporal physical
        // accuracy needs a refinement study, not a fitted one-case threshold.
        std::cout << "nonuniform coupled BDF: " << spatialResult.internalSteps << " steps, explicit " << explicitState.steps
                  << ", temperature difference K=" << temperatureDifference << '\n';
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
