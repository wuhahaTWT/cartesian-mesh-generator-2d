#include "cartmesh2d/geometry/ConstructionRecovery2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iterator>
#include <map>
#include <stdexcept>
#include <tuple>

namespace cartmesh2d {
namespace {

[[nodiscard]] int decisionPriority(ConstructionDecision2D decision) noexcept {
    switch (decision) {
    case ConstructionDecision2D::Refine: return 0;
    case ConstructionDecision2D::Rephase: return 1;
    case ConstructionDecision2D::Resample: return 2;
    case ConstructionDecision2D::ReusedExactKey: return 3;
    case ConstructionDecision2D::Accepted: return 4;
    case ConstructionDecision2D::Rejected: return 5;
    case ConstructionDecision2D::ShadowOnly: return 6;
    }
    return 6;
}

} // namespace

std::vector<ConstructionRecoveryCandidate2D> planConstructionRecoveryCandidates(
    const ConstructionRecoveryRequest2D& request,
    std::size_t maximumGridPhaseCandidates) {
    const auto d=request.sourceSegment.b-request.sourceSegment.a;
    const double length=std::sqrt(squaredNorm(d));
    if (!std::isfinite(request.localH) || !(request.localH>0.0) ||
        !std::isfinite(request.sourceParameter) || request.sourceParameter<0.0 ||
        request.sourceParameter>1.0 || !std::isfinite(length) || !(length>0.0) ||
        maximumGridPhaseCandidates>16U) {
        throw std::invalid_argument("invalid typed construction recovery request");
    }

    std::vector<ConstructionRecoveryCandidate2D> result;
    result.push_back({0U,ConstructionDecision2D::Refine,std::nullopt,std::nullopt,
                      std::nullopt,false,true,"refine_affected_leaf_with_2_to_1_closure"});

    // Bounded deterministic offsets only. These are proposals, not geometry
    // mutations; the caller must rebuild the local patch and evaluate it.
    static constexpr double phaseFractions[]{-0.125,0.125,-0.25,0.25};
    const auto phaseCount=std::min(maximumGridPhaseCandidates,
                                   std::size(phaseFractions));
    for (std::size_t i=0;i<phaseCount;++i) {
        result.push_back({result.size(),ConstructionDecision2D::Rephase,
                          std::nullopt,std::nullopt,
                          phaseFractions[i]*request.localH,false,true,
                          "bounded_grid_phase_candidate"});
    }

    // Resampling stays strictly on the original source segment and strictly
    // inside its feature endpoints. No normal displacement or cross-feature
    // search is representable by this candidate type.
    const double delta=std::min(0.25,0.25*request.localH/length);
    std::vector<double> parameters;
    for (const double sign:{-1.0,1.0}) {
        const double parameter=request.sourceParameter+sign*delta;
        if (parameter>0.0 && parameter<1.0) parameters.push_back(parameter);
    }
    std::sort(parameters.begin(),parameters.end());
    parameters.erase(std::unique(parameters.begin(),parameters.end()),parameters.end());
    for (const double parameter:parameters) {
        result.push_back({result.size(),ConstructionDecision2D::Resample,
                          request.sourceSegment.a+d*parameter,parameter,std::nullopt,
                          true,true,"source_parameter_resample_candidate"});
    }
    return result;
}

std::optional<ConstructionRecoveryCandidate2D> selectConstructionRecoveryCandidate(
    const std::vector<ConstructionRecoveryCandidate2D>& candidates,
    const std::vector<ConstructionRecoveryCandidateEvaluation2D>& evaluations) {
    std::map<std::size_t,ConstructionRecoveryCandidateEvaluation2D> byOrdinal;
    for (const auto& evaluation:evaluations) {
        if (!std::isfinite(evaluation.qualityRank)) {
            throw std::invalid_argument("construction recovery quality rank must be finite");
        }
        if (!byOrdinal.emplace(evaluation.ordinal,evaluation).second) {
            throw std::invalid_argument("duplicate construction recovery evaluation ordinal");
        }
    }
    std::optional<ConstructionRecoveryCandidate2D> best;
    double bestRank=std::numeric_limits<double>::infinity();
    for (const auto& candidate:candidates) {
        const auto found=byOrdinal.find(candidate.ordinal);
        if (found==byOrdinal.end()) continue;
        const auto& evaluation=found->second;
        if (!evaluation.areaConserved || !evaluation.featureCompatible ||
            !evaluation.hardQualityPass) continue;
        if (!best || std::tuple(evaluation.qualityRank,
                                decisionPriority(candidate.decision),candidate.ordinal)<
                     std::tuple(bestRank,decisionPriority(best->decision),best->ordinal)) {
            best=candidate;bestRank=evaluation.qualityRank;
        }
    }
    return best;
}

} // namespace cartmesh2d
