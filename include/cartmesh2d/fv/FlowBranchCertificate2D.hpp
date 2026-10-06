#pragma once

#include "cartmesh2d/fv/Incompressible2D.hpp"

#include <functional>
#include <iosfwd>
#include <optional>
#include <string>

namespace cartmesh2d::fv {

// Explicit steady-solver branch probe. It never replaces the ordinary result:
// the direct target, high-viscosity guide and guide-seeded target are retained
// as separate candidates, in that order. A divergent certificate reports risk;
// it does not decide which converged branch is physical.
enum class FlowBranchCertificateStage2D { DirectTarget, Guide, GuidedTarget };
enum class FlowBranchCertificateOutcome2D { Completed, Unconverged, Stopped, Failed };

struct FlowBranchCertificateControls2D {
    double guideViscosityMultiplier = 10;
    double guideToleranceExponent = .5;
    // These dimensionless RMS limits are supplied by the caller. They are not
    // flow acceptance thresholds and do not inspect solution amplitudes.
    double maximumVelocityRmsDifference = 0;
    double maximumPressureRmsDifference = 0;
    // Zero inherits the target solver budget independently for each stage.
    std::size_t guideMaximumIterations = 0;
    std::size_t guidedTargetMaximumIterations = 0;
    // Checked before a stage and cooperatively after complete solver updates.
    // Throwing propagates exactly like FlowControls2D::stopRequested.
    std::function<bool(FlowBranchCertificateStage2D)> stopRequested;
};

struct FlowBranchDifference2D {
    double velocityRms = 0;       // area weighted, divided by target Uref
    double maximumVelocity = 0;   // vector difference / target Uref
    double pressureRms = 0;       // area weighted after gauge removal / Uref^2
    double maximumPressure = 0;   // after the same gauge removal / Uref^2
};

struct FlowBranchCertificate2D {
    FlowBranchCertificateStage2D stage = FlowBranchCertificateStage2D::DirectTarget;
    FlowBranchCertificateOutcome2D outcome = FlowBranchCertificateOutcome2D::Failed;
    std::string reason;
    std::optional<FlowResult2D> directTarget;
    std::optional<FlowResult2D> guide;
    std::optional<FlowResult2D> guidedTarget;
    FlowBranchDifference2D difference;
    bool pathsCompared = false;
    bool consistent = false;
};

[[nodiscard]] FlowBranchCertificate2D certifyIncompressibleBranch2D(
    const FvMesh2D&, const FlowControls2D&, const FlowBranchCertificateControls2D&);

// Durable certificate representation. Candidate fields are ordinary versioned
// FlowState2D checkpoints, kept separately so either can seed a later explicit
// solve. The archive records cost and reason but deliberately omits a selected
// or preferred candidate.
struct FlowBranchCandidateArchive2D {
    FlowState2D state;
    std::size_t coupledEvaluations = 0;
    bool converged = false;
    bool stopped = false;
};

struct FlowBranchCertificateArchive2D {
    FlowBranchCertificateStage2D stage = FlowBranchCertificateStage2D::DirectTarget;
    FlowBranchCertificateOutcome2D outcome = FlowBranchCertificateOutcome2D::Failed;
    std::string reason;
    double maximumVelocityRmsDifference = 0;
    double maximumPressureRmsDifference = 0;
    FlowBranchDifference2D difference;
    bool pathsCompared = false;
    bool consistent = false;
    std::optional<FlowBranchCandidateArchive2D> directTarget;
    std::optional<FlowBranchCandidateArchive2D> guide;
    std::optional<FlowBranchCandidateArchive2D> guidedTarget;
};

void writeFlowBranchCertificate2D(
    std::ostream&, const FvMesh2D&, const FlowControls2D&,
    const FlowBranchCertificateControls2D&, const FlowBranchCertificate2D&);

[[nodiscard]] FlowBranchCertificateArchive2D readFlowBranchCertificate2D(
    std::istream&, const FvMesh2D&, const FlowControls2D&,
    const FlowBranchCertificateControls2D&);

} // namespace cartmesh2d::fv
