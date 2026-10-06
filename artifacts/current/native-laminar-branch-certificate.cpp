// Product branch-certificate transaction on an existing final FvMesh2D.
// This driver supplies reporting limits explicitly; they are not solver gates.
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/fv/FlowBranchCertificate2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

const char* stage(cartmesh2d::fv::FlowBranchCertificateStage2D value) {
    using S = cartmesh2d::fv::FlowBranchCertificateStage2D;
    switch (value) {
    case S::DirectTarget: return "direct-target";
    case S::Guide: return "guide";
    case S::GuidedTarget: return "guided-target";
    }
    return "invalid";
}

const char* outcome(cartmesh2d::fv::FlowBranchCertificateOutcome2D value) {
    using O = cartmesh2d::fv::FlowBranchCertificateOutcome2D;
    switch (value) {
    case O::Completed: return "completed";
    case O::Unconverged: return "unconverged";
    case O::Stopped: return "stopped";
    case O::Failed: return "failed";
    }
    return "invalid";
}

void candidate(std::ostream& output, const char* name,
               const std::optional<cartmesh2d::fv::FlowResult2D>& value) {
    output << std::quoted(name) << ':';
    if (!value) { output << "null"; return; }
    output << "{\"converged\":" << (value->converged ? "true" : "false")
           << ",\"stopped\":" << (value->stopped ? "true" : "false")
           << ",\"evaluations\":" << value->performance.coupledEvaluations
           << ",\"solveSeconds\":" << value->performance.solveSeconds
           << ",\"globalRelativeImbalance\":" << value->globalRelativeImbalance
           << ",\"maximumSpeedRatio\":" << value->fieldAmplitude.maximumSpeedRatio
           << ",\"pressureRangeRatio\":" << value->fieldAmplitude.pressureRangeRatio
           << ",\"pressureForceX\":" << value->pressureForceX
           << ",\"viscousForceX\":" << value->forceX - value->pressureForceX
           << ",\"forceX\":" << value->forceX << '}';
}

} // namespace

int main(int argc, char** argv) {
    using namespace cartmesh2d;
    using namespace cartmesh2d::fv;
    if (argc != 4)
        throw std::runtime_error("usage: branch-certificate mesh.solver.cm2d boundary-file output-prefix");
    const auto read = readCm2dTopology(argv[1]);
    if (!read.valid()) throw std::runtime_error(read.error);
    const auto mesh = makeFvMesh2D(read.topology);
    FlowControls2D target;
    target.scenario = "custom";
    target.nu = .1;
    target.speed = 1;
    target.tolerance = 1e-6;
    target.maxIterations = 3000;
    target.profile = true;
    target.steadyAcceleration = SteadyAcceleration2D::NewtonKrylov;
    target.convection = ConvectionScheme2D::FaceLimitedLinearUpwind;
    std::ifstream boundaries(argv[2]);
    if (!boundaries) throw std::runtime_error("cannot open boundary file");
    target.boundaryConditions = readFlowBoundaryConditions2D(boundaries, mesh, target);

    FlowBranchCertificateControls2D controls;
    // Research reporting resolution for normalized all-cell RMS differences.
    // Classification never selects or overwrites either target candidate.
    controls.maximumVelocityRmsDifference = 1e-4;
    controls.maximumPressureRmsDifference = 1e-4;
    const auto result = certifyIncompressibleBranch2D(mesh, target, controls);

    const std::string prefix = argv[3];
    std::ofstream archive(prefix + ".certificate");
    writeFlowBranchCertificate2D(archive, mesh, target, controls, result);
    archive.close();
    std::ifstream archiveInput(prefix + ".certificate");
    const auto restored = readFlowBranchCertificate2D(archiveInput, mesh, target, controls);

    std::ofstream output(prefix + ".json");
    output << std::setprecision(17)
           << "{\"schema\":1,\"cells\":" << mesh.cells.size()
           << ",\"stage\":" << std::quoted(stage(result.stage))
           << ",\"outcome\":" << std::quoted(outcome(result.outcome))
           << ",\"reason\":" << std::quoted(result.reason)
           << ",\"pathsCompared\":" << (result.pathsCompared ? "true" : "false")
           << ",\"consistent\":" << (result.consistent ? "true" : "false")
           << ",\"limits\":{\"velocityRms\":" << controls.maximumVelocityRmsDifference
           << ",\"pressureRms\":" << controls.maximumPressureRmsDifference << "}"
           << ",\"difference\":{\"velocityRms\":" << result.difference.velocityRms
           << ",\"maximumVelocity\":" << result.difference.maximumVelocity
           << ",\"pressureRms\":" << result.difference.pressureRms
           << ",\"maximumPressure\":" << result.difference.maximumPressure << "},\"candidates\":{";
    candidate(output, "directTarget", result.directTarget); output << ',';
    candidate(output, "guide", result.guide); output << ',';
    candidate(output, "guidedTarget", result.guidedTarget);
    output << "},\"archiveReadback\":{\"pathsCompared\":"
           << (restored.pathsCompared ? "true" : "false")
           << ",\"consistent\":" << (restored.consistent ? "true" : "false")
           << ",\"directTarget\":" << (restored.directTarget ? "true" : "false")
           << ",\"guide\":" << (restored.guide ? "true" : "false")
           << ",\"guidedTarget\":" << (restored.guidedTarget ? "true" : "false") << "}}\n";
    return result.outcome == FlowBranchCertificateOutcome2D::Completed &&
           result.pathsCompared && result.directTarget && result.guide && result.guidedTarget ? 0 : 2;
}
