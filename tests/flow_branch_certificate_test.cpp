#include "cartmesh2d/fv/FlowBranchCertificate2D.hpp"

#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

FvMesh2D rectangle(int nx = 10, int ny = 4) {
    TopologyMesh2D topology;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
            topology.vertices.push_back({topology.vertices.size(), {2. * i / nx, double(j) / ny}});
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edges;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        TopologyCell2D cell;
        cell.id = topology.cells.size();
        cell.geometryArea = 2. / nx / ny;
        const auto a = static_cast<std::size_t>(j * (nx + 1) + i);
        cell.vertices = {a, a + 1, a + static_cast<std::size_t>(nx) + 2,
                         a + static_cast<std::size_t>(nx) + 1};
        for (std::size_t k = 0; k < 4; ++k) {
            const auto x = cell.vertices[k], y = cell.vertices[(k + 1) % 4];
            const auto [entry, inserted] = edges.emplace(std::minmax(x, y), topology.edges.size());
            if (inserted)
                topology.edges.push_back({entry->second, x, y, cell.id, {}, BoundaryPatch2D::DomainBoundary});
            else {
                topology.edges[entry->second].neighbour = cell.id;
                topology.edges[entry->second].patch = BoundaryPatch2D::None;
            }
            cell.edges.push_back(entry->second);
        }
        topology.cells.push_back(cell);
    }
    return makeFvMesh2D(topology);
}

FlowControls2D controls(const FvMesh2D& mesh) {
    FlowControls2D result;
    result.scenario = "custom";
    result.nu = .2;
    result.tolerance = 1e-6;
    result.maxIterations = 1000;
    result.steadyAcceleration = SteadyAcceleration2D::NewtonKrylov;
    result.convection = ConvectionScheme2D::FaceLimitedLinearUpwind;
    result.pressurePreconditioner = PressurePreconditioner2D::Aggregation;
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& face = mesh.faces[id];
        if (face.neighbour) continue;
        FlowBoundaryCondition2D boundary{id, FlowBoundaryKind2D::Wall, {}, 0, "wall"};
        if (face.areaVector.x < 0)
            boundary = {id, FlowBoundaryKind2D::VelocityInlet, {1, 0}, 0, "inlet"};
        else if (face.areaVector.x > 0)
            boundary = {id, FlowBoundaryKind2D::PressureOutlet, {}, 0, "outlet"};
        result.boundaryConditions.push_back(boundary);
    }
    return result;
}

void sameState(const FlowResult2D& expected, const FlowBranchCandidateArchive2D& actual) {
    require(actual.state.time == 0 && expected.u == actual.state.u && expected.v == actual.state.v &&
            expected.p == actual.state.p && expected.flux == actual.state.flux &&
            expected.performance.coupledEvaluations == actual.coupledEvaluations &&
            expected.converged == actual.converged && expected.stopped == actual.stopped,
            "branch candidate archive changed field or cost");
}

} // namespace

int main() {
    const auto mesh = rectangle();
    const auto target = controls(mesh);
    FlowBranchCertificateControls2D certificate;
    certificate.maximumVelocityRmsDifference = 1e-4;
    certificate.maximumPressureRmsDifference = 1e-4;

    const auto completed = certifyIncompressibleBranch2D(mesh, target, certificate);
    require(completed.outcome == FlowBranchCertificateOutcome2D::Completed &&
            completed.stage == FlowBranchCertificateStage2D::GuidedTarget &&
            completed.directTarget && completed.guide && completed.guidedTarget &&
            completed.pathsCompared && completed.consistent,
            "normal branch certificate did not preserve three consistent stages");
    require(completed.difference.velocityRms > 0 && completed.difference.pressureRms > 0,
            "branch certificate comparison did not observe independent paths");

    std::stringstream archive;
    writeFlowBranchCertificate2D(archive, mesh, target, certificate, completed);
    const auto restored = readFlowBranchCertificate2D(archive, mesh, target, certificate);
    require(restored.outcome == completed.outcome && restored.stage == completed.stage &&
            restored.reason == completed.reason && restored.pathsCompared && restored.consistent &&
            restored.directTarget && restored.guide && restored.guidedTarget,
            "branch certificate archive lost classification or candidates");
    sameState(*completed.directTarget, *restored.directTarget);
    sameState(*completed.guide, *restored.guide);
    sameState(*completed.guidedTarget, *restored.guidedTarget);

    auto strict = certificate;
    strict.maximumVelocityRmsDifference = 1e-30;
    strict.maximumPressureRmsDifference = 1e-30;
    const auto divergent = certifyIncompressibleBranch2D(mesh, target, strict);
    require(divergent.outcome == FlowBranchCertificateOutcome2D::Completed &&
            divergent.pathsCompared && !divergent.consistent && divergent.directTarget &&
            divergent.guidedTarget,
            "caller-declared branch disagreement was silently selected");

    auto guideFailure = certificate;
    guideFailure.guideMaximumIterations = 1;
    const auto failedGuide = certifyIncompressibleBranch2D(mesh, target, guideFailure);
    require(failedGuide.outcome == FlowBranchCertificateOutcome2D::Unconverged &&
            failedGuide.stage == FlowBranchCertificateStage2D::Guide &&
            failedGuide.directTarget && failedGuide.directTarget->converged && failedGuide.guide &&
            !failedGuide.guide->converged && !failedGuide.guidedTarget,
            "guide failure replaced or hid the direct accepted candidate");

    auto targetFailure = certificate;
    targetFailure.guidedTargetMaximumIterations = 1;
    const auto failedTarget = certifyIncompressibleBranch2D(mesh, target, targetFailure);
    require(failedTarget.outcome == FlowBranchCertificateOutcome2D::Unconverged &&
            failedTarget.stage == FlowBranchCertificateStage2D::GuidedTarget &&
            failedTarget.directTarget && failedTarget.directTarget->converged &&
            failedTarget.guide && failedTarget.guide->converged && failedTarget.guidedTarget &&
            !failedTarget.guidedTarget->converged,
            "guided-target failure replaced an earlier accepted candidate");

    auto stopped = certificate;
    stopped.stopRequested = [](FlowBranchCertificateStage2D stage) {
        return stage == FlowBranchCertificateStage2D::Guide;
    };
    const auto cancelled = certifyIncompressibleBranch2D(mesh, target, stopped);
    require(cancelled.outcome == FlowBranchCertificateOutcome2D::Stopped &&
            cancelled.stage == FlowBranchCertificateStage2D::Guide &&
            cancelled.directTarget && cancelled.directTarget->converged &&
            !cancelled.guide && !cancelled.guidedTarget,
            "stage cancellation replaced the direct accepted candidate");
    std::stringstream cancelledArchive;
    writeFlowBranchCertificate2D(cancelledArchive, mesh, target, stopped, cancelled);
    const auto restoredCancellation = readFlowBranchCertificate2D(
        cancelledArchive, mesh, target, stopped);
    require(restoredCancellation.outcome == FlowBranchCertificateOutcome2D::Stopped &&
            restoredCancellation.stage == FlowBranchCertificateStage2D::Guide &&
            restoredCancellation.directTarget && !restoredCancellation.guide &&
            !restoredCancellation.guidedTarget && !restoredCancellation.pathsCompared,
            "cancelled branch archive invented or lost a candidate");

    auto throwing = certificate;
    throwing.stopRequested = [](FlowBranchCertificateStage2D stage) -> bool {
        if (stage == FlowBranchCertificateStage2D::Guide) throw std::runtime_error("caller cancel failure");
        return false;
    };
    bool propagated = false;
    try { (void)certifyIncompressibleBranch2D(mesh, target, throwing); }
    catch (const std::runtime_error& error) { propagated = std::string(error.what()) == "caller cancel failure"; }
    require(propagated, "branch certificate swallowed callback exception");

    auto incompatible = certificate;
    incompatible.maximumVelocityRmsDifference *= 2;
    bool rejected = false;
    std::stringstream copy(archive.str());
    try { (void)readFlowBranchCertificate2D(copy, mesh, target, incompatible); }
    catch (const std::runtime_error&) { rejected = true; }
    require(rejected, "branch archive accepted different caller limits");

    rejected = false;
    auto truncatedText = archive.str();
    truncatedText.resize(truncatedText.size() - 17);
    std::stringstream truncated(std::move(truncatedText));
    try { (void)readFlowBranchCertificate2D(truncated, mesh, target, certificate); }
    catch (const std::runtime_error&) { rejected = true; }
    require(rejected, "branch archive accepted truncated candidate data");

    auto invalidOutcomeText = archive.str();
    const auto outcome = invalidOutcomeText.find("OUTCOME completed");
    require(outcome != std::string::npos, "branch archive omitted completed outcome");
    invalidOutcomeText.replace(outcome, std::string("OUTCOME completed").size(), "OUTCOME stopped  ");
    std::stringstream invalidOutcome(std::move(invalidOutcomeText));
    rejected = false;
    try { (void)readFlowBranchCertificate2D(invalidOutcome, mesh, target, certificate); }
    catch (const std::runtime_error&) { rejected = true; }
    require(rejected, "branch archive accepted outcome/classification contradiction");

    std::cout << "Flow branch certificate: independent candidates, classification, failure, cancellation and archive verified\n";
}
