#include "cartmesh2d/fv/FlowBranchCertificate2D.hpp"

#include "cartmesh2d/fv/FlowCheckpoint2D.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cartmesh2d::fv {
namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

const char* stageName(FlowBranchCertificateStage2D stage) {
    switch (stage) {
    case FlowBranchCertificateStage2D::DirectTarget: return "direct-target";
    case FlowBranchCertificateStage2D::Guide: return "guide";
    case FlowBranchCertificateStage2D::GuidedTarget: return "guided-target";
    }
    throw std::runtime_error("Invalid flow branch certificate stage");
}

FlowBranchCertificateStage2D stageValue(const std::string& value) {
    if (value == "direct-target") return FlowBranchCertificateStage2D::DirectTarget;
    if (value == "guide") return FlowBranchCertificateStage2D::Guide;
    if (value == "guided-target") return FlowBranchCertificateStage2D::GuidedTarget;
    throw std::runtime_error("Flow branch certificate: invalid stage");
}

const char* outcomeName(FlowBranchCertificateOutcome2D outcome) {
    switch (outcome) {
    case FlowBranchCertificateOutcome2D::Completed: return "completed";
    case FlowBranchCertificateOutcome2D::Unconverged: return "unconverged";
    case FlowBranchCertificateOutcome2D::Stopped: return "stopped";
    case FlowBranchCertificateOutcome2D::Failed: return "failed";
    }
    throw std::runtime_error("Invalid flow branch certificate outcome");
}

FlowBranchCertificateOutcome2D outcomeValue(const std::string& value) {
    if (value == "completed") return FlowBranchCertificateOutcome2D::Completed;
    if (value == "unconverged") return FlowBranchCertificateOutcome2D::Unconverged;
    if (value == "stopped") return FlowBranchCertificateOutcome2D::Stopped;
    if (value == "failed") return FlowBranchCertificateOutcome2D::Failed;
    throw std::runtime_error("Flow branch certificate: invalid outcome");
}

FlowControls2D explicitTargetControls(const FlowControls2D& input) {
    auto result = input;
    result.steadyAcceleration = resolveSteadyAcceleration2D(input);
    result.convection = resolveSteadyConvection2D(input);
    return result;
}

FlowControls2D guideControls(const FlowControls2D& target,
                             const FlowBranchCertificateControls2D& certificate) {
    auto result = explicitTargetControls(target);
    result.nu = target.nu * certificate.guideViscosityMultiplier;
    result.tolerance = std::pow(target.tolerance, certificate.guideToleranceExponent);
    if (certificate.guideMaximumIterations != 0)
        result.maxIterations = certificate.guideMaximumIterations;
    return result;
}

void validateCertificateControls(const FvMesh2D& mesh, const FlowControls2D& target,
                                 const FlowBranchCertificateControls2D& certificate) {
    validateFvMesh2D(mesh);
    require(std::isfinite(target.nu) && target.nu > 0 &&
            std::isfinite(target.speed) && target.speed > 0 &&
            std::isfinite(target.tolerance) && target.tolerance > 0 && target.tolerance < 1,
            "Branch certificate requires positive steady controls and target tolerance below one");
    require(resolveSteadyAcceleration2D(target) == SteadyAcceleration2D::NewtonKrylov &&
            target.convergence == FlowConvergence2D::Strict && !target.adaptiveLinear &&
            target.momentumInertia == 1,
            "Branch certificate requires strict constant-property Newton-Krylov flow");
    require(target.faceViscosity.empty() && target.manufacturedViscositySlope == 0,
            "Branch certificate currently requires uniform viscosity");
    require(std::isfinite(certificate.guideViscosityMultiplier) &&
            certificate.guideViscosityMultiplier > 1 &&
            std::isfinite(target.nu * certificate.guideViscosityMultiplier),
            "Branch certificate guide viscosity multiplier must be finite and above one");
    require(std::isfinite(certificate.guideToleranceExponent) &&
            certificate.guideToleranceExponent > 0 && certificate.guideToleranceExponent < 1,
            "Branch certificate guide tolerance exponent must be between zero and one");
    require(std::isfinite(certificate.maximumVelocityRmsDifference) &&
            certificate.maximumVelocityRmsDifference > 0 &&
            std::isfinite(certificate.maximumPressureRmsDifference) &&
            certificate.maximumPressureRmsDifference > 0,
            "Branch certificate comparison tolerances must be positive and finite");
}

FlowInitialGuess2D guess(const FlowResult2D& result) {
    return {result.u, result.v, result.p, result.flux};
}

FlowState2D state(const FlowResult2D& result) {
    return {0, result.u, result.v, result.p, result.flux};
}

FlowBranchDifference2D compare(const FvMesh2D& mesh, const FlowControls2D& target,
                               const FlowResult2D& direct, const FlowResult2D& guided) {
    require(direct.u.size() == mesh.cells.size() && direct.v.size() == mesh.cells.size() &&
            direct.p.size() == mesh.cells.size() && guided.u.size() == mesh.cells.size() &&
            guided.v.size() == mesh.cells.size() && guided.p.size() == mesh.cells.size(),
            "Branch certificate candidate size differs from mesh");
    double area = 0, pressureMean = 0;
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        area += mesh.cells[i].area;
        pressureMean += mesh.cells[i].area * (guided.p[i] - direct.p[i]);
    }
    require(std::isfinite(area) && area > 0, "Branch certificate mesh has invalid area");
    pressureMean /= area;
    double velocitySquared = 0, pressureSquared = 0;
    FlowBranchDifference2D result;
    const double pressureScale = target.speed * target.speed;
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        const double du = (guided.u[i] - direct.u[i]) / target.speed;
        const double dv = (guided.v[i] - direct.v[i]) / target.speed;
        const double dp = (guided.p[i] - direct.p[i] - pressureMean) / pressureScale;
        const double velocity = std::hypot(du, dv);
        velocitySquared += mesh.cells[i].area * velocity * velocity;
        pressureSquared += mesh.cells[i].area * dp * dp;
        result.maximumVelocity = std::max(result.maximumVelocity, velocity);
        result.maximumPressure = std::max(result.maximumPressure, std::abs(dp));
    }
    result.velocityRms = std::sqrt(velocitySquared / area);
    result.pressureRms = std::sqrt(pressureSquared / area);
    require(std::isfinite(result.velocityRms) && std::isfinite(result.maximumVelocity) &&
            std::isfinite(result.pressureRms) && std::isfinite(result.maximumPressure),
            "Branch certificate comparison is nonfinite");
    return result;
}

struct StageRun {
    std::optional<FlowResult2D> result;
    std::string failure;
};

template<class Solve>
StageRun runStage(FlowBranchCertificateStage2D stage, FlowControls2D controls,
                  const FlowBranchCertificateControls2D& certificate, Solve&& solve) {
    if (certificate.stopRequested && certificate.stopRequested(stage))
        return {{}, "stopped before " + std::string(stageName(stage))};
    std::exception_ptr callbackError;
    const auto originalStop = controls.stopRequested;
    controls.stopRequested = [&, originalStop] {
        try {
            return (originalStop && originalStop()) ||
                   (certificate.stopRequested && certificate.stopRequested(stage));
        } catch (...) {
            callbackError = std::current_exception();
            return true;
        }
    };
    try {
        auto result = solve(controls);
        if (callbackError) std::rethrow_exception(callbackError);
        return {std::move(result), {}};
    } catch (const std::exception& error) {
        if (callbackError) std::rethrow_exception(callbackError);
        return {{}, error.what()};
    }
}

void token(std::istream& input, const char* expected) {
    std::string actual;
    if (!(input >> actual) || actual != expected)
        throw std::runtime_error("Flow branch certificate: expected '" + std::string(expected) + "'");
}

void finite(double value, const char* what) {
    if (!std::isfinite(value))
        throw std::runtime_error("Flow branch certificate: nonfinite " + std::string(what));
}

FlowControls2D archivalTarget(const FlowControls2D& target) {
    auto result = explicitTargetControls(target);
    result.stopRequested = {};
    return result;
}

void writeCandidate(std::ostream& output, const char* name,
                    const std::optional<FlowResult2D>& candidate,
                    const FvMesh2D& mesh, const FlowControls2D& controls) {
    std::string checkpoint;
    if (candidate) {
        std::ostringstream buffer;
        writeFlowCheckpoint2D(buffer, mesh, controls, state(*candidate));
        checkpoint = std::move(buffer).str();
    }
    output << "CANDIDATE " << name << ' ' << (candidate ? 1 : 0) << ' '
           << (candidate ? candidate->performance.coupledEvaluations : 0) << ' '
           << (candidate && candidate->converged ? 1 : 0) << ' '
           << (candidate && candidate->stopped ? 1 : 0) << ' '
           << checkpoint.size() << '\n' << checkpoint;
}

std::optional<FlowBranchCandidateArchive2D> readCandidate(
    std::istream& input, const char* expectedName, const FvMesh2D& mesh,
    const FlowControls2D& controls) {
    token(input, "CANDIDATE");
    std::string name;
    int present = 0, converged = 0, stopped = 0;
    unsigned long long evaluations = 0, checkpointSize = 0;
    if (!(input >> name >> present >> evaluations >> converged >> stopped >> checkpointSize) ||
        name != expectedName ||
        (present != 0 && present != 1) ||
        (converged != 0 && converged != 1) || (stopped != 0 && stopped != 1) ||
        evaluations > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()) ||
        checkpointSize > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()))
        throw std::runtime_error("Flow branch certificate: invalid candidate header");
    if (input.get() != '\n')
        throw std::runtime_error("Flow branch certificate: invalid candidate separator");
    if (!present) {
        if (evaluations != 0 || converged || stopped || checkpointSize != 0)
            throw std::runtime_error("Flow branch certificate: absent candidate has work or state");
        return {};
    }
    if (evaluations == 0 || checkpointSize == 0 || (converged && stopped))
        throw std::runtime_error("Flow branch certificate: candidate state flags or data are invalid");
    std::string checkpoint(static_cast<std::size_t>(checkpointSize), '\0');
    if (!input.read(checkpoint.data(), static_cast<std::streamsize>(checkpoint.size())))
        throw std::runtime_error("Flow branch certificate: truncated candidate checkpoint");
    FlowBranchCandidateArchive2D result;
    result.coupledEvaluations = static_cast<std::size_t>(evaluations);
    result.converged = converged != 0;
    result.stopped = stopped != 0;
    std::istringstream buffer(std::move(checkpoint));
    result.state = readFlowCheckpoint2D(buffer, mesh, controls);
    if (result.state.time != 0)
        throw std::runtime_error("Flow branch certificate: steady candidate has nonzero time");
    return result;
}

} // namespace

FlowBranchCertificate2D certifyIncompressibleBranch2D(
    const FvMesh2D& mesh, const FlowControls2D& input,
    const FlowBranchCertificateControls2D& certificate) {
    validateCertificateControls(mesh, input, certificate);
    auto target = explicitTargetControls(input);
    FlowBranchCertificate2D result;

    auto direct = runStage(FlowBranchCertificateStage2D::DirectTarget, target, certificate,
        [&](const FlowControls2D& controls) { return solveIncompressible2D(mesh, controls); });
    result.stage = FlowBranchCertificateStage2D::DirectTarget;
    if (!direct.result) {
        result.outcome = direct.failure.starts_with("stopped before")
            ? FlowBranchCertificateOutcome2D::Stopped : FlowBranchCertificateOutcome2D::Failed;
        result.reason = std::move(direct.failure);
        return result;
    }
    result.directTarget = std::move(direct.result);
    if (result.directTarget->stopped) {
        result.outcome = FlowBranchCertificateOutcome2D::Stopped;
        result.reason = "direct target stopped";
        return result;
    }
    if (!result.directTarget->converged) {
        result.outcome = FlowBranchCertificateOutcome2D::Unconverged;
        result.reason = "direct target did not converge";
        return result;
    }

    auto guideControl = guideControls(target, certificate);
    auto guide = runStage(FlowBranchCertificateStage2D::Guide, guideControl, certificate,
        [&](const FlowControls2D& controls) { return solveIncompressible2D(mesh, controls); });
    result.stage = FlowBranchCertificateStage2D::Guide;
    if (!guide.result) {
        result.outcome = guide.failure.starts_with("stopped before")
            ? FlowBranchCertificateOutcome2D::Stopped : FlowBranchCertificateOutcome2D::Failed;
        result.reason = std::move(guide.failure);
        return result;
    }
    result.guide = std::move(guide.result);
    if (result.guide->stopped) {
        result.outcome = FlowBranchCertificateOutcome2D::Stopped;
        result.reason = "guide stopped";
        return result;
    }
    if (!result.guide->converged) {
        result.outcome = FlowBranchCertificateOutcome2D::Unconverged;
        result.reason = "guide did not converge";
        return result;
    }

    auto guidedControl = target;
    if (certificate.guidedTargetMaximumIterations != 0)
        guidedControl.maxIterations = certificate.guidedTargetMaximumIterations;
    const auto guideGuess = guess(*result.guide);
    auto guided = runStage(FlowBranchCertificateStage2D::GuidedTarget, guidedControl, certificate,
        [&](const FlowControls2D& controls) {
            return solveIncompressibleFromGuess2D(mesh, controls, guideGuess);
        });
    result.stage = FlowBranchCertificateStage2D::GuidedTarget;
    if (!guided.result) {
        result.outcome = guided.failure.starts_with("stopped before")
            ? FlowBranchCertificateOutcome2D::Stopped : FlowBranchCertificateOutcome2D::Failed;
        result.reason = std::move(guided.failure);
        return result;
    }
    result.guidedTarget = std::move(guided.result);
    if (result.guidedTarget->stopped) {
        result.outcome = FlowBranchCertificateOutcome2D::Stopped;
        result.reason = "guided target stopped";
        return result;
    }
    if (!result.guidedTarget->converged) {
        result.outcome = FlowBranchCertificateOutcome2D::Unconverged;
        result.reason = "guided target did not converge";
        return result;
    }

    result.difference = compare(mesh, target, *result.directTarget, *result.guidedTarget);
    result.pathsCompared = true;
    result.consistent = result.difference.velocityRms <= certificate.maximumVelocityRmsDifference &&
                        result.difference.pressureRms <= certificate.maximumPressureRmsDifference;
    result.outcome = FlowBranchCertificateOutcome2D::Completed;
    result.reason = result.consistent ? "independent target paths agree within caller limits"
                                      : "independent target paths differ beyond caller limits";
    return result;
}

void writeFlowBranchCertificate2D(
    std::ostream& output, const FvMesh2D& mesh, const FlowControls2D& input,
    const FlowBranchCertificateControls2D& controls, const FlowBranchCertificate2D& certificate) {
    validateCertificateControls(mesh, input, controls);
    if (certificate.pathsCompared)
        require(certificate.directTarget && certificate.guidedTarget,
                "Compared branch certificate is missing a target candidate");
    struct RestoreFormat {
        std::ostream& stream;
        std::ios_base::fmtflags flags;
        std::streamsize precision;
        ~RestoreFormat() { stream.flags(flags); stream.precision(precision); }
    } restore{output, output.flags(), output.precision()};
    output << std::defaultfloat << std::dec << std::noshowpos << std::noshowbase
           << std::setprecision(17);
    output << "CARTMESH2D_FLOW_BRANCH_CERTIFICATE 1\n"
           << "STAGE " << stageName(certificate.stage) << '\n'
           << "OUTCOME " << outcomeName(certificate.outcome) << '\n'
           << "REASON " << std::quoted(certificate.reason) << '\n'
           << "GUIDE_SETTINGS " << controls.guideViscosityMultiplier << ' '
           << controls.guideToleranceExponent << ' ' << controls.guideMaximumIterations << ' '
           << controls.guidedTargetMaximumIterations << '\n'
           << "LIMITS " << controls.maximumVelocityRmsDifference << ' '
           << controls.maximumPressureRmsDifference << '\n'
           << "DIFFERENCE " << certificate.difference.velocityRms << ' '
           << certificate.difference.maximumVelocity << ' '
           << certificate.difference.pressureRms << ' '
           << certificate.difference.maximumPressure << '\n'
           << "CLASSIFICATION " << (certificate.pathsCompared ? 1 : 0) << ' '
           << (certificate.consistent ? 1 : 0) << '\n';
    const auto target = archivalTarget(input);
    auto guide = guideControls(target, controls);
    guide.stopRequested = {};
    writeCandidate(output, "DIRECT_TARGET", certificate.directTarget, mesh, target);
    writeCandidate(output, "GUIDE", certificate.guide, mesh, guide);
    writeCandidate(output, "GUIDED_TARGET", certificate.guidedTarget, mesh, target);
    output << "END_BRANCH_CERTIFICATE\n";
    if (!output) throw std::runtime_error("Flow branch certificate: write failed");
}

FlowBranchCertificateArchive2D readFlowBranchCertificate2D(
    std::istream& input, const FvMesh2D& mesh, const FlowControls2D& source,
    const FlowBranchCertificateControls2D& controls) {
    validateCertificateControls(mesh, source, controls);
    token(input, "CARTMESH2D_FLOW_BRANCH_CERTIFICATE");
    int version = 0;
    if (!(input >> version) || version != 1)
        throw std::runtime_error("Flow branch certificate: unsupported version");
    FlowBranchCertificateArchive2D result;
    token(input, "STAGE");
    std::string value;
    if (!(input >> value)) throw std::runtime_error("Flow branch certificate: truncated stage");
    result.stage = stageValue(value);
    token(input, "OUTCOME");
    if (!(input >> value)) throw std::runtime_error("Flow branch certificate: truncated outcome");
    result.outcome = outcomeValue(value);
    token(input, "REASON");
    if (!(input >> std::quoted(result.reason)))
        throw std::runtime_error("Flow branch certificate: truncated reason");
    token(input, "GUIDE_SETTINGS");
    double guideMultiplier = 0, guideExponent = 0;
    unsigned long long guideBudget = 0, guidedBudget = 0;
    if (!(input >> guideMultiplier >> guideExponent >> guideBudget >> guidedBudget) ||
        guideBudget > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()) ||
        guidedBudget > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()) ||
        guideMultiplier != controls.guideViscosityMultiplier ||
        guideExponent != controls.guideToleranceExponent ||
        guideBudget != controls.guideMaximumIterations ||
        guidedBudget != controls.guidedTargetMaximumIterations)
        throw std::runtime_error("Flow branch certificate: guide controls differ from caller");
    token(input, "LIMITS");
    if (!(input >> result.maximumVelocityRmsDifference >> result.maximumPressureRmsDifference))
        throw std::runtime_error("Flow branch certificate: truncated limits");
    finite(result.maximumVelocityRmsDifference, "velocity limit");
    finite(result.maximumPressureRmsDifference, "pressure limit");
    if (result.maximumVelocityRmsDifference != controls.maximumVelocityRmsDifference ||
        result.maximumPressureRmsDifference != controls.maximumPressureRmsDifference)
        throw std::runtime_error("Flow branch certificate: comparison limits differ from caller");
    token(input, "DIFFERENCE");
    if (!(input >> result.difference.velocityRms >> result.difference.maximumVelocity >>
          result.difference.pressureRms >> result.difference.maximumPressure))
        throw std::runtime_error("Flow branch certificate: truncated difference");
    finite(result.difference.velocityRms, "velocity RMS difference");
    finite(result.difference.maximumVelocity, "maximum velocity difference");
    finite(result.difference.pressureRms, "pressure RMS difference");
    finite(result.difference.maximumPressure, "maximum pressure difference");
    if (result.difference.velocityRms < 0 || result.difference.maximumVelocity < 0 ||
        result.difference.pressureRms < 0 || result.difference.maximumPressure < 0)
        throw std::runtime_error("Flow branch certificate: negative difference");
    token(input, "CLASSIFICATION");
    int compared = 0, consistent = 0;
    if (!(input >> compared >> consistent) || (compared != 0 && compared != 1) ||
        (consistent != 0 && consistent != 1) || (!compared && consistent))
        throw std::runtime_error("Flow branch certificate: invalid classification");
    result.pathsCompared = compared != 0;
    result.consistent = consistent != 0;
    const auto target = archivalTarget(source);
    auto guide = guideControls(target, controls);
    guide.stopRequested = {};
    result.directTarget = readCandidate(input, "DIRECT_TARGET", mesh, target);
    result.guide = readCandidate(input, "GUIDE", mesh, guide);
    result.guidedTarget = readCandidate(input, "GUIDED_TARGET", mesh, target);
    token(input, "END_BRANCH_CERTIFICATE");
    input >> std::ws;
    if (!input.eof()) throw std::runtime_error("Flow branch certificate: trailing data");
    const bool classifiedConsistent = result.difference.velocityRms <= result.maximumVelocityRmsDifference &&
                                      result.difference.pressureRms <= result.maximumPressureRmsDifference;
    const FlowBranchCandidateArchive2D* current = nullptr;
    if (result.stage == FlowBranchCertificateStage2D::DirectTarget && result.directTarget)
        current = &*result.directTarget;
    else if (result.stage == FlowBranchCertificateStage2D::Guide && result.guide)
        current = &*result.guide;
    else if (result.stage == FlowBranchCertificateStage2D::GuidedTarget && result.guidedTarget)
        current = &*result.guidedTarget;
    const bool completed = result.outcome == FlowBranchCertificateOutcome2D::Completed;
    const bool validOutcome =
        (completed && result.stage == FlowBranchCertificateStage2D::GuidedTarget &&
            result.directTarget && result.guide && result.guidedTarget &&
            result.directTarget->converged && result.guide->converged &&
            result.guidedTarget->converged && !result.directTarget->stopped &&
            !result.guide->stopped && !result.guidedTarget->stopped) ||
        (result.outcome == FlowBranchCertificateOutcome2D::Unconverged && current &&
            !current->converged && !current->stopped) ||
        (result.outcome == FlowBranchCertificateOutcome2D::Stopped &&
            (!current || current->stopped)) ||
        (result.outcome == FlowBranchCertificateOutcome2D::Failed && !current);
    if ((result.pathsCompared != completed) ||
        (result.pathsCompared && result.consistent != classifiedConsistent) ||
        (!result.pathsCompared && result.consistent) || !validOutcome ||
        (result.stage == FlowBranchCertificateStage2D::Guide && !result.directTarget) ||
        (result.stage == FlowBranchCertificateStage2D::GuidedTarget &&
            (!result.directTarget || !result.guide)) ||
        (result.directTarget && result.stage != FlowBranchCertificateStage2D::DirectTarget &&
            !result.directTarget->converged) ||
        (result.guide && result.stage == FlowBranchCertificateStage2D::GuidedTarget &&
            !result.guide->converged))
        throw std::runtime_error("Flow branch certificate: inconsistent stage, candidates or classification");
    return result;
}

} // namespace cartmesh2d::fv
