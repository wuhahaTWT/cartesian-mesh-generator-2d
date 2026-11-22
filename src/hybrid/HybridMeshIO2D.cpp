#include "cartmesh2d/hybrid/HybridMesh2D.hpp"

#include <fstream>
#include <iomanip>
#include <ostream>
#include <string>

namespace cartmesh2d {
namespace {

void setError(std::string* error, const std::string& message) {
    if (error != nullptr) *error = message;
}

void writeJsonString(std::ostream& out, const std::string& value) {
    out << '"';
    for (const char character : value) {
        switch (character) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default: out << character; break;
        }
    }
    out << '"';
}

} // namespace

const char* h4MeshModeName(H4MeshMode2D mode) noexcept {
    switch (mode) {
    case H4MeshMode2D::Failed: return "failed";
    case H4MeshMode2D::Hybrid: return "hybrid";
    case H4MeshMode2D::PureCutCellFallback: return "pure_cutcell_fallback";
    }
    return "unknown";
}

const char* h4FallbackStageName(H4FallbackStage2D stage) noexcept {
    switch (stage) {
    case H4FallbackStage2D::None: return "none";
    case H4FallbackStage2D::RequestedLayers: return "requested_layers";
    case H4FallbackStage2D::LocalReduction: return "local_reduction";
    case H4FallbackStage2D::HybridCandidate: return "hybrid_candidate";
    }
    return "unknown";
}

const char* hybridMeshFailureReasonName(HybridMeshFailureReason2D reason) noexcept {
    switch (reason) {
    case HybridMeshFailureReason2D::None: return "none";
    case HybridMeshFailureReason2D::InvalidInput: return "invalid_input";
    case HybridMeshFailureReason2D::UnsupportedWallSemantics: return "unsupported_wall_semantics";
    case HybridMeshFailureReason2D::OuterEnvelopeOutsideDomain: return "outer_envelope_outside_domain";
    case HybridMeshFailureReason2D::InvalidOuterEnvelope: return "invalid_outer_envelope";
    case HybridMeshFailureReason2D::RemainderRefinementFailed: return "remainder_refinement_failed";
    case HybridMeshFailureReason2D::RemainderCutCellFailed: return "remainder_cutcell_failed";
    case HybridMeshFailureReason2D::RemainderStabilizationFailed: return "remainder_stabilization_failed";
    case HybridMeshFailureReason2D::LayerConversionFailed: return "layer_conversion_failed";
    case HybridMeshFailureReason2D::UnifiedTopologyFailed: return "unified_topology_failed";
    case HybridMeshFailureReason2D::NonConformalInterface: return "nonconformal_interface";
    case HybridMeshFailureReason2D::AreaConservationFailed: return "area_conservation_failed";
    case HybridMeshFailureReason2D::RegionClassificationConflict: return "region_classification_conflict";
    case HybridMeshFailureReason2D::QualityFailed: return "quality_failed";
    case HybridMeshFailureReason2D::SolverTopologyFailed: return "solver_topology_failed";
    case HybridMeshFailureReason2D::SolverQualityFailed: return "solver_quality_failed";
    case HybridMeshFailureReason2D::IoFailure: return "io_failure";
    }
    return "unknown";
}

const char* hybridCellKindName(HybridCellKind2D kind) noexcept {
    switch (kind) {
    case HybridCellKind2D::BoundaryLayer: return "boundary_layer";
    case HybridCellKind2D::RemainderCut: return "remainder_cut";
    case HybridCellKind2D::RemainderCartesian: return "remainder_cartesian";
    case HybridCellKind2D::Termination: return "termination";
    case HybridCellKind2D::Transition: return "transition";
    }
    return "unknown";
}

bool writeHybridLegacyVtk2D(const HybridMeshBuildResult2D& result,
                            const std::filesystem::path& path,
                            std::string* error) {
    if (!result.success() || !result.topology.valid() ||
        result.cellRecords.size() != result.topology.cells.size()) {
        setError(error, "cannot write failed or inconsistent H4-2 hybrid candidate");
        return false;
    }
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open H4-2 VTK output");
        return false;
    }
    std::size_t cellListSize = 0U;
    for (const auto& cell : result.topology.cells) {
        cellListSize += cell.vertices.size() + 1U;
    }
    out << std::setprecision(17);
    out << "# vtk DataFile Version 3.0\n";
    out << "cartmesh2d H4-2 conformal hybrid mesh\nASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << result.topology.vertices.size() << " double\n";
    for (const auto& vertex : result.topology.vertices) {
        out << vertex.point.x << ' ' << vertex.point.y << " 0\n";
    }
    out << "CELLS " << result.topology.cells.size() << ' ' << cellListSize << "\n";
    for (const auto& cell : result.topology.cells) {
        out << cell.vertices.size();
        for (const auto id : cell.vertices) out << ' ' << id;
        out << '\n';
    }
    out << "CELL_TYPES " << result.topology.cells.size() << "\n";
    for (std::size_t i = 0; i < result.topology.cells.size(); ++i) out << "7\n";
    out << "CELL_DATA " << result.topology.cells.size() << "\n";
    out << "SCALARS hybrid_kind int 1\nLOOKUP_TABLE default\n";
    for (const auto& record : result.cellRecords) {
        out << static_cast<int>(record.kind) << '\n';
    }
    out << "SCALARS layer_index int 1\nLOOKUP_TABLE default\n";
    for (const auto& record : result.cellRecords) {
        out << (record.layerIndex ? static_cast<long long>(*record.layerIndex) : -1LL)
            << '\n';
    }
    out << "SCALARS geometry_area double 1\nLOOKUP_TABLE default\n";
    for (const auto& cell : result.topology.cells) out << cell.geometryArea << '\n';
    if (!out.good()) {
        setError(error, "failed while writing H4-2 VTK output");
        return false;
    }
    return true;
}

bool writeHybridReportJson2D(const HybridMeshBuildResult2D& result,
                             const std::filesystem::path& path,
                             std::string* error) {
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open H4-2 JSON report");
        return false;
    }
    out << std::setprecision(17);
    out << "{\n  \"hybrid_status\": \""
        << (result.success() ? "success" : "failed") << "\",\n";
    if (!result.success()) {
        out << "  \"failure_reason\": \""
            << hybridMeshFailureReasonName(result.failure.reason) << "\",\n";
        out << "  \"message\": ";
        writeJsonString(out, result.failure.message);
        out << ",\n  \"strip_id\": ";
        if (result.failure.stripId) out << *result.failure.stripId; else out << "null";
        out << ",\n  \"leaf_id\": ";
        if (result.failure.leafId) out << *result.failure.leafId; else out << "null";
        out << ",\n  \"cell_id\": ";
        if (result.failure.cellId) out << *result.failure.cellId; else out << "null";
        out << ",\n  \"edge_id\": ";
        if (result.failure.edgeId) out << *result.failure.edgeId; else out << "null";
        out << "\n}\n";
    } else {
        const auto& metrics = result.metrics;
        const auto& interface = result.interfaceAudit;
        const auto& profile = result.solverTopologyReport.profile;
        std::size_t recoveryRequestedRefined=0U,recoveryClosureRefined=0U;
        std::size_t recoveryClosureIterations=0U,recoveryLineageSplits=0U;
        for (const auto& refinement:result.constructionRecoveryRefinements) {
            recoveryRequestedRefined+=refinement.requestedRefinedLeaves;
            recoveryClosureRefined+=refinement.closureRefinedLeaves;
            recoveryClosureIterations+=refinement.closureIterations;
            recoveryLineageSplits+=refinement.lineage.size();
        }
        out << "  \"failure_reason\": \"none\",\n";
        out << "  \"directional_repair\": {\"candidate_count\": "
            <<metrics.directionalRepairCandidateCount<<", \"accepted_count\": "
            <<metrics.directionalRepairAcceptedCount<<", \"minimum_before\": ";
        if (metrics.directionalMinimumBefore) out<<*metrics.directionalMinimumBefore;else out<<"null";
        out<<", \"minimum_after\": ";
        if (metrics.directionalMinimumAfter) out<<*metrics.directionalMinimumAfter;else out<<"null";
        out<<", \"failed_cell_count\": "<<metrics.directionalFailedCellCount<<"},\n";
        out << "  \"target_solver_repair\": {\"max_non_orthogonality_deg\": 65, "
            << "\"attempted\": "
            << (metrics.targetPolicyRepairAttempted?"true":"false")
            << ", \"committed\": "
            << (metrics.targetPolicyRepairCommitted?"true":"false")
            << ", \"structural_failure\": "
            << (metrics.targetPolicyRepairStructuralFailure?"true":"false")
            << ", \"issue_count_before\": "
            << metrics.targetPolicyRepairIssueCountBefore
            << ", \"issue_count_after\": "
            << metrics.targetPolicyRepairIssueCountAfter
            << ", \"repartition_count\": "
            << metrics.targetPolicyRepairRepartitionCount << "},\n";
        out << "  \"extruded_determinant_repair_repartition_count\": "
            << metrics.extrudedDeterminantRepairRepartitionCount << ",\n";
        out << "  \"quadtree_leaf_count\": " << metrics.quadtreeLeafCount << ",\n";
        out << "  \"boundary_layer_cell_count\": " << metrics.boundaryLayerCellCount << ",\n";
        out << "  \"requested_boundary_layer_cell_count\": "
            << metrics.requestedBoundaryLayerCellCount << ",\n";
        out << "  \"zero_layer_column_count\": "
            << metrics.zeroLayerColumnCount << ",\n";
        out << "  \"termination_cell_count\": "
            << metrics.terminationCellCount << ",\n";
        out << "  \"termination_edge_count\": "
            << metrics.terminationEdgeCount << ",\n";
        out << "  \"remainder_cut_cell_count\": " << metrics.remainderCutCellCount << ",\n";
        out << "  \"remainder_cartesian_cell_count\": " << metrics.remainderCartesianCellCount << ",\n";
        out << "  \"transition_polygon_count\": " << metrics.transitionPolygonCount << ",\n";
        out << "  \"termination_growth_ratio\": "
            << metrics.terminationGrowthRatio << ",\n";
        out << "  \"vertex_count\": " << metrics.unifiedVertexCount << ",\n";
        out << "  \"edge_count\": " << metrics.unifiedEdgeCount << ",\n";
        out << "  \"cell_count\": " << metrics.unifiedCellCount << ",\n";
        out << "  \"solver_cell_count\": " << metrics.solverCellCount << ",\n";
        out << "  \"r1_short_face_candidate_count\": "
            << metrics.r1ShortFaceCandidates << ",\n";
        out << "  \"r1_local_candidate_count\": "
            << metrics.r1LocalCandidates << ",\n";
        out << "  \"r1_local_quality_evaluation_count\": "
            << metrics.r1LocalQualityEvaluations << ",\n";
        out << "  \"r1_accepted_transaction_count\": "
            << metrics.r1AcceptedTransactions << ",\n";
        out << "  \"r1_candidate_global_topology_build_count\": "
            << metrics.r1CandidateGlobalTopologyBuilds << ",\n";
        out << "  \"r1_candidate_full_global_quality_evaluation_count\": "
            << metrics.r1CandidateFullGlobalQualityEvaluations << ",\n";
        out << "  \"r1_global_oracle_build_count\": "
            << metrics.r1GlobalOracleBuilds << ",\n";
        out << "  \"r1_maximum_winner_global_oracle_builds_per_transaction\": "
            << metrics.r1MaximumWinnerGlobalOracleBuilds << ",\n";
        out << "  \"r1_authoritative_full_quality_evaluation_count\": "
            << metrics.r1AuthoritativeFullQualityEvaluations << ",\n";
        out << "  \"q3_termination_candidate_count\": "
            << metrics.q3TerminationCandidates << ",\n";
        out << "  \"q3_local_candidate_count\": "
            << metrics.q3LocalCandidates << ",\n";
        out << "  \"q3_local_quality_evaluation_count\": "
            << metrics.q3LocalQualityEvaluations << ",\n";
        out << "  \"q3_accepted_transaction_count\": "
            << metrics.q3AcceptedTransactions << ",\n";
        out << "  \"q3_candidate_global_topology_build_count\": "
            << metrics.q3CandidateGlobalTopologyBuilds << ",\n";
        out << "  \"q3_candidate_full_global_quality_evaluation_count\": "
            << metrics.q3CandidateFullGlobalQualityEvaluations << ",\n";
        out << "  \"q3_global_oracle_build_count\": "
            << metrics.q3GlobalOracleBuilds << ",\n";
        out << "  \"q3_maximum_winner_global_oracle_builds_per_transaction\": "
            << metrics.q3MaximumWinnerGlobalOracleBuilds << ",\n";
        out << "  \"q3_authoritative_full_quality_evaluation_count\": "
            << metrics.q3AuthoritativeFullQualityEvaluations << ",\n";
        out << "  \"q3_hard_volume_ratio_count_before\": "
            << metrics.q3HardVolumeRatioBefore << ",\n";
        out << "  \"q3_hard_volume_ratio_count_after\": "
            << metrics.q3HardVolumeRatioAfter << ",\n";
        out << "  \"q3_hard_face_weight_count_before\": "
            << metrics.q3HardFaceWeightBefore << ",\n";
        out << "  \"q3_hard_face_weight_count_after\": "
            << metrics.q3HardFaceWeightAfter << ",\n";
        out << "  \"q3_hard_short_face_count_before\": "
            << metrics.q3HardShortFaceBefore << ",\n";
        out << "  \"q3_hard_short_face_count_after\": "
            << metrics.q3HardShortFaceAfter << ",\n";
        out << "  \"q3_maximum_volume_ratio_severity_before\": "
            << metrics.q3MaximumVolumeRatioSeverityBefore << ",\n";
        out << "  \"q3_maximum_volume_ratio_severity_after\": "
            << metrics.q3MaximumVolumeRatioSeverityAfter << ",\n";
        out << "  \"q3_maximum_face_weight_severity_before\": "
            << metrics.q3MaximumFaceWeightSeverityBefore << ",\n";
        out << "  \"q3_maximum_face_weight_severity_after\": "
            << metrics.q3MaximumFaceWeightSeverityAfter << ",\n";
        out << "  \"q3_patch_outside_stable_ids_unchanged\": "
            << (metrics.q3PatchOutsideStableIdsUnchanged?"true":"false") << ",\n";
        out << "  \"q3_local_winner_matches_global_authority\": "
            << (metrics.q3LocalWinnerMatchesGlobalAuthority?"true":"false") << ",\n";
        out << "  \"q3_local_delta_matches_global_oracle\": "
            << (metrics.q3LocalDeltaMatchesGlobalOracle?"true":"false") << ",\n";
        out << "  \"q3_transaction_bound_reached\": "
            << (metrics.q3TransactionBoundReached?"true":"false") << ",\n";
        out << "  \"q32_repartition_candidate_count\": "
            << metrics.q32RepartitionCandidates << ",\n";
        out << "  \"q32_local_candidate_count\": "
            << metrics.q32LocalCandidates << ",\n";
        out << "  \"q32_local_quality_evaluation_count\": "
            << metrics.q32LocalQualityEvaluations << ",\n";
        out << "  \"q32_accepted_transaction_count\": "
            << metrics.q32AcceptedTransactions << ",\n";
        out << "  \"q32_candidate_global_topology_build_count\": "
            << metrics.q32CandidateGlobalTopologyBuilds << ",\n";
        out << "  \"q32_candidate_full_global_quality_evaluation_count\": "
            << metrics.q32CandidateFullGlobalQualityEvaluations << ",\n";
        out << "  \"q32_global_oracle_build_count\": "
            << metrics.q32GlobalOracleBuilds << ",\n";
        out << "  \"q32_maximum_winner_global_oracle_builds_per_transaction\": "
            << metrics.q32MaximumWinnerGlobalOracleBuilds << ",\n";
        out << "  \"q32_authoritative_full_quality_evaluation_count\": "
            << metrics.q32AuthoritativeFullQualityEvaluations << ",\n";
        out << "  \"q32_hard_volume_ratio_count_before\": "
            << metrics.q32HardVolumeRatioBefore << ",\n";
        out << "  \"q32_hard_volume_ratio_count_after\": "
            << metrics.q32HardVolumeRatioAfter << ",\n";
        out << "  \"q32_hard_face_weight_count_before\": "
            << metrics.q32HardFaceWeightBefore << ",\n";
        out << "  \"q32_hard_face_weight_count_after\": "
            << metrics.q32HardFaceWeightAfter << ",\n";
        out << "  \"q32_hard_short_face_count_before\": "
            << metrics.q32HardShortFaceBefore << ",\n";
        out << "  \"q32_hard_short_face_count_after\": "
            << metrics.q32HardShortFaceAfter << ",\n";
        out << "  \"q32_patch_outside_stable_ids_unchanged\": "
            << (metrics.q32PatchOutsideStableIdsUnchanged?"true":"false") << ",\n";
        out << "  \"q32_local_winner_matches_global_authority\": "
            << (metrics.q32LocalWinnerMatchesGlobalAuthority?"true":"false") << ",\n";
        out << "  \"q32_local_delta_matches_global_oracle\": "
            << (metrics.q32LocalDeltaMatchesGlobalOracle?"true":"false") << ",\n";
        out << "  \"q32_transaction_bound_reached\": "
            << (metrics.q32TransactionBoundReached?"true":"false") << ",\n";
        out << "  \"q33_grouped_candidate_count\": "
            << metrics.q33GroupedCandidates << ",\n";
        out << "  \"q33_local_candidate_count\": "
            << metrics.q33LocalCandidates << ",\n";
        out << "  \"q33_local_quality_evaluation_count\": "
            << metrics.q33LocalQualityEvaluations << ",\n";
        out << "  \"q33_accepted_transaction_count\": "
            << metrics.q33AcceptedTransactions << ",\n";
        out << "  \"q33_candidate_global_topology_build_count\": "
            << metrics.q33CandidateGlobalTopologyBuilds << ",\n";
        out << "  \"q33_candidate_full_global_quality_evaluation_count\": "
            << metrics.q33CandidateFullGlobalQualityEvaluations << ",\n";
        out << "  \"q33_global_oracle_build_count\": "
            << metrics.q33GlobalOracleBuilds << ",\n";
        out << "  \"q33_maximum_winner_global_oracle_builds_per_transaction\": "
            << metrics.q33MaximumWinnerGlobalOracleBuilds << ",\n";
        out << "  \"q33_authoritative_full_quality_evaluation_count\": "
            << metrics.q33AuthoritativeFullQualityEvaluations << ",\n";
        out << "  \"q33_hard_volume_ratio_count_before\": "
            << metrics.q33HardVolumeRatioBefore << ",\n";
        out << "  \"q33_hard_volume_ratio_count_after\": "
            << metrics.q33HardVolumeRatioAfter << ",\n";
        out << "  \"q33_hard_face_weight_count_before\": "
            << metrics.q33HardFaceWeightBefore << ",\n";
        out << "  \"q33_hard_face_weight_count_after\": "
            << metrics.q33HardFaceWeightAfter << ",\n";
        out << "  \"q33_hard_short_face_count_before\": "
            << metrics.q33HardShortFaceBefore << ",\n";
        out << "  \"q33_hard_short_face_count_after\": "
            << metrics.q33HardShortFaceAfter << ",\n";
        out << "  \"q33_patch_outside_stable_ids_unchanged\": "
            << (metrics.q33PatchOutsideStableIdsUnchanged?"true":"false") << ",\n";
        out << "  \"q33_local_winner_matches_global_authority\": "
            << (metrics.q33LocalWinnerMatchesGlobalAuthority?"true":"false") << ",\n";
        out << "  \"q33_local_delta_matches_global_oracle\": "
            << (metrics.q33LocalDeltaMatchesGlobalOracle?"true":"false") << ",\n";
        out << "  \"q33_transaction_bound_reached\": "
            << (metrics.q33TransactionBoundReached?"true":"false") << ",\n";
        out << "  \"q41_construction_candidate_count\": "
            << metrics.q41ConstructionCandidates << ",\n";
        out << "  \"q41_local_candidate_count\": "
            << metrics.q41LocalCandidates << ",\n";
        out << "  \"q41_local_quality_evaluation_count\": "
            << metrics.q41LocalQualityEvaluations << ",\n";
        out << "  \"q41_accepted_construction_count\": "
            << metrics.q41AcceptedConstructions << ",\n";
        out << "  \"q41_candidate_global_topology_build_count\": "
            << metrics.q41CandidateGlobalTopologyBuilds << ",\n";
        out << "  \"q41_candidate_full_global_quality_evaluation_count\": "
            << metrics.q41CandidateFullGlobalQualityEvaluations << ",\n";
        out << "  \"q41_global_oracle_build_count\": "
            << metrics.q41GlobalOracleBuilds << ",\n";
        out << "  \"q41_maximum_winner_global_oracle_builds_per_construction\": "
            << metrics.q41MaximumWinnerGlobalOracleBuilds << ",\n";
        out << "  \"q41_authoritative_full_quality_evaluation_count\": "
            << metrics.q41AuthoritativeFullQualityEvaluations << ",\n";
        out << "  \"q41_hard_volume_ratio_count_before\": "
            << metrics.q41HardVolumeRatioBefore << ",\n";
        out << "  \"q41_hard_volume_ratio_count_after\": "
            << metrics.q41HardVolumeRatioAfter << ",\n";
        out << "  \"q41_hard_face_weight_count_before\": "
            << metrics.q41HardFaceWeightBefore << ",\n";
        out << "  \"q41_hard_face_weight_count_after\": "
            << metrics.q41HardFaceWeightAfter << ",\n";
        out << "  \"q41_hard_short_face_count_before\": "
            << metrics.q41HardShortFaceBefore << ",\n";
        out << "  \"q41_hard_short_face_count_after\": "
            << metrics.q41HardShortFaceAfter << ",\n";
        out << "  \"q41_patch_outside_stable_ids_unchanged\": "
            << (metrics.q41PatchOutsideStableIdsUnchanged?"true":"false") << ",\n";
        out << "  \"q41_local_winner_matches_global_authority\": "
            << (metrics.q41LocalWinnerMatchesGlobalAuthority?"true":"false") << ",\n";
        out << "  \"q41_local_delta_matches_global_oracle\": "
            << (metrics.q41LocalDeltaMatchesGlobalOracle?"true":"false") << ",\n";
        out << "  \"q41_construction_bound_reached\": "
            << (metrics.q41ConstructionBoundReached?"true":"false") << ",\n";
        out << "  \"q41_construction_selection_declined\": "
            << (metrics.q41ConstructionSelectionDeclined?"true":"false") << ",\n";
        out << "  \"q51_outer_transition_radial_subdivision\": "
            << metrics.q51OuterTransitionRadialSubdivision << ",\n";
        out << "  \"q51_outer_transition_radial_declined\": "
            << (metrics.q51OuterTransitionRadialDeclined?"true":"false") << ",\n";
        out << "  \"q51_outer_transition_radial_target_reachable\": "
            << (metrics.q51OuterTransitionRadialTargetReachable?"true":"false")
            << ",\n";
        out << "  \"q52_termination_buffer_radial_committed\": "
            << (metrics.q52TerminationBufferRadialCommitted?"true":"false") << ",\n";
        out << "  \"q52_termination_buffer_target_issues_with_matching\": "
            << metrics.q52TerminationBufferTargetIssuesWithMatching << ",\n";
        out << "  \"q52_termination_buffer_target_issues_with_historical_march\": "
            << metrics.q52TerminationBufferTargetIssuesWithHistoricalMarch << ",\n";
        out << "  \"q52_termination_buffer_rows\": "
            << metrics.q52TerminationBufferRows << ",\n";
        out << "  \"q52_termination_buffer_growth_ratio\": "
            << metrics.q52TerminationBufferGrowthRatio << ",\n";
        out << "  \"q52_termination_buffer_outer_row_thickness\": "
            << metrics.q52TerminationBufferOuterRowThickness << ",\n";
        out << "  \"q52_termination_buffer_row_cap\": "
            << metrics.q52TerminationBufferRowCap << ",\n";
        out << "  \"q52_termination_buffer_row_cap_reachable\": "
            << (metrics.q52TerminationBufferRowCapReachable?"true":"false") << ",\n";
        out << "  \"build_global_topology_call_count\": "
            << metrics.buildGlobalTopologyCalls << ",\n";
        out << "  \"build_global_topology_input_cell_total\": "
            << metrics.buildGlobalTopologyInputCells << ",\n";
        out << "  \"full_solver_quality_evaluation_count\": "
            << metrics.fullSolverQualityEvaluations << ",\n";
        // Attribution of the two totals above. The legacy agglomeration and
        // repartition passes build small patch topologies through the same
        // entry point; those are counted here so the R1 zeros above cannot be
        // mistaken for the whole pipeline being rebuild-free.
        out << "  \"solver_profile_global_topology_rebuild_calls\": "
            << profile.globalTopologyRebuildCalls << ",\n";
        out << "  \"solver_profile_candidate_topology_count\": "
            << profile.candidateTopologyCount << ",\n";
        out << "  \"solver_profile_full_quality_calls\": "
            << profile.fullQualityCalls << ",\n";
        out << "  \"solver_profile_candidate_quality_evaluation_count\": "
            << profile.candidateQualityEvaluationCount << ",\n";
        out << "  \"r1_minimum_face_over_local_h_before\": "
            << metrics.r1MinimumFaceOverLocalHBefore << ",\n";
        out << "  \"r1_minimum_face_over_local_h_after\": "
            << metrics.r1MinimumFaceOverLocalHAfter << ",\n";
        out << "  \"r1_patch_outside_stable_ids_unchanged\": "
            << (metrics.r1PatchOutsideStableIdsUnchanged?"true":"false") << ",\n";
        out << "  \"r1_local_winner_matches_global_authority\": "
            << (metrics.r1LocalWinnerMatchesGlobalAuthority?"true":"false") << ",\n";
        out << "  \"r1_local_delta_matches_global_oracle\": "
            << (metrics.r1LocalDeltaMatchesGlobalOracle?"true":"false") << ",\n";
        out << "  \"construction_half_edge_count\": "
            << result.constructionIncidence.audit.halfEdges << ",\n";
        out << "  \"construction_twin_pair_count\": "
            << result.constructionIncidence.audit.twinPairs << ",\n";
        out << "  \"solver_half_edge_count\": "
            << result.solverIncidence.audit.halfEdges << ",\n";
        out << "  \"solver_twin_pair_count\": "
            << result.solverIncidence.audit.twinPairs << ",\n";
        out << "  \"solid_area\": " << metrics.solidArea << ",\n";
        out << "  \"outer_envelope_area\": " << metrics.outerEnvelopeArea << ",\n";
        out << "  \"layer_area\": " << metrics.layerArea << ",\n";
        out << "  \"remainder_area\": " << metrics.remainderArea << ",\n";
        out << "  \"expected_fluid_area\": " << metrics.expectedFluidArea << ",\n";
        out << "  \"actual_fluid_area\": " << metrics.actualFluidArea << ",\n";
        out << "  \"area_error\": " << metrics.areaError << ",\n";
        out << "  \"source_lineage_solver_cell_count\": "
            << result.sourceLineageAudit.solverCellCount << ",\n";
        out << "  \"source_lineage_candidate_checks\": "
            << result.sourceLineageAudit.lineageCandidateChecks << ",\n";
        out << "  \"source_lineage_oracle_candidate_checks\": "
            << result.sourceLineageAudit.oracleCandidateChecks << ",\n";
        out << "  \"source_lineage_mismatched_cells\": "
            << result.sourceLineageAudit.mismatchedCells << ",\n";
        out << "  \"source_lineage_oracle_verified\": "
            << (result.sourceLineageAudit.oracleVerified?"true":"false") << ",\n";
        out << "  \"construction_recovery_request_count\": "
            << result.constructionRecoveryRequests.size() << ",\n";
        out << "  \"construction_recovery_refinement_passes\": "
            << result.constructionRecoveryRefinements.size() << ",\n";
        out << "  \"construction_recovery_requested_refined_leaves\": "
            << recoveryRequestedRefined << ",\n";
        out << "  \"construction_recovery_closure_refined_leaves\": "
            << recoveryClosureRefined << ",\n";
        out << "  \"construction_recovery_closure_iterations\": "
            << recoveryClosureIterations << ",\n";
        out << "  \"construction_recovery_lineage_splits\": "
            << recoveryLineageSplits << ",\n";
        out << "  \"interface_edge_count\": " << interface.interfaceEdgeCount << ",\n";
        out << "  \"interface_vertex_count\": " << interface.interfaceVertexCount << ",\n";
        out << "  \"single_owner_interface_edges\": " << interface.singleOwnerInterfaceEdges << ",\n";
        out << "  \"wrong_cell_pair_interface_edges\": " << interface.wrongCellPairInterfaceEdges << ",\n";
        out << "  \"non_two_valent_interface_vertices\": " << interface.nonTwoValentInterfaceVertices << ",\n";
        out << "  \"expected_interface_length\": " << interface.expectedLength << ",\n";
        out << "  \"actual_interface_length\": " << interface.actualLength << ",\n";
        out << "  \"interface_length_error\": " << interface.lengthError << ",\n";
        out << "  \"topology_valid\": " << (result.topology.valid() ? "true" : "false") << ",\n";
        out << "  \"mesh_quality_valid\": " << (result.meshQuality.valid() ? "true" : "false") << ",\n";
        out << "  \"solver_quality_valid\": " << (result.solverQuality.valid() ? "true" : "false") << ",\n";
        out << "  \"solver_quality_issue_count\": " << result.solverQuality.issues.size() << ",\n";
        out << "  \"max_non_orthogonality_deg\": "
            << result.solverQuality.maxNonOrthogonalityDeg << ",\n";
        out << "  \"min_face_weight\": "
            << result.solverQuality.minFaceWeight << ",\n";
        out << "  \"min_volume_ratio\": "
            << result.solverQuality.minVolumeRatio << ",\n";
        out << "  \"balance_violations_after\": " << result.balance.violationsAfter << "\n}\n";
    }
    if (!out.good()) {
        setError(error, "failed while writing H4-2 JSON report");
        return false;
    }
    return true;
}

bool writeHybridProfileJson2D(const HybridMeshBuildResult2D& result,
                              const std::filesystem::path& path,
                              std::string* error,
                              const RobustH4Profile2D* robustProfile) {
    std::ofstream out(path);
    if (!out) {
        setError(error, "failed to open hybrid profile JSON");
        return false;
    }
    const auto& metrics=result.metrics;
    const auto& profile=result.solverTopologyReport.profile;
    out << std::setprecision(17);
    out << "{\n";
    out << "  \"measurement_class\": \"wall_time\",\n";
    out << "  \"reproducible\": false,\n";
    out << "  \"r1_repair_seconds\": " << metrics.r1RepairSeconds << ",\n";
    out << "  \"directional_repair_seconds\": " << metrics.directionalRepairSeconds << ",\n";
    // The optional termination-repair variants report their wall time here, not
    // in the byte-compared report: two runs of the same input legitimately
    // differ in seconds, which would otherwise make the determinism comparison
    // fail on an unchanged mesh.
    out << "  \"q3_repair_seconds\": " << metrics.q3RepairSeconds << ",\n";
    out << "  \"q32_repair_seconds\": " << metrics.q32RepairSeconds << ",\n";
    out << "  \"q33_repair_seconds\": " << metrics.q33RepairSeconds << ",\n";
    out << "  \"q41_construction_selection_seconds\": "
        << metrics.q41ConstructionSelectionSeconds << ",\n";
    out << "  \"build_global_topology_seconds\": "
        << metrics.buildGlobalTopologySeconds << ",\n";
    out << "  \"solver_repair_seconds\": " << metrics.solverRepairSeconds << ",\n";
    out << "  \"solver_profile_initial_partition_seconds\": "
        << profile.initialPartitionSeconds << ",\n";
    out << "  \"solver_profile_source_repair_seconds\": "
        << profile.sourceRepairSeconds << ",\n";
    out << "  \"solver_profile_final_repartition_seconds\": "
        << profile.finalRepartitionSeconds << ",\n";
    out << "  \"solver_profile_candidate_polygon_work_seconds\": "
        << profile.candidatePolygonWorkSeconds << ",\n";
    out << "  \"solver_profile_candidate_global_rebuild_seconds\": "
        << profile.candidateGlobalRebuildSeconds << ",\n";
    out << "  \"solver_profile_candidate_quality_seconds\": "
        << profile.candidateQualitySeconds << ",\n";
    out << "  \"solver_profile_full_quality_seconds\": "
        << profile.fullQualitySeconds;
    if (robustProfile) {
        out << ",\n";
        out << "  \"h4_total_seconds\": " << robustProfile->totalSeconds << ",\n";
        out << "  \"h4_requested_layer_seconds\": "
            << robustProfile->requestedLayerSeconds << ",\n";
        out << "  \"h4_requested_hybrid_seconds\": "
            << robustProfile->requestedHybridSeconds << ",\n";
        out << "  \"h4_local_layer_seconds\": "
            << robustProfile->localLayerSeconds << ",\n";
        out << "  \"h4_local_hybrid_seconds\": "
            << robustProfile->localHybridSeconds << ",\n";
        out << "  \"h4_pure_cutcell_fallback_seconds\": "
            << robustProfile->pureCutCellFallbackSeconds << ",\n";
        out << "  \"h4_unattributed_seconds\": "
            << robustProfile->unattributedSeconds() << ",\n";
        out << "  \"h4_requested_hybrid_attempts\": "
            << robustProfile->requestedHybridAttempts << ",\n";
        out << "  \"h4_local_hybrid_attempts\": "
            << robustProfile->localHybridAttempts << ",\n";
        out << "  \"h4_pure_cutcell_fallback_attempts\": "
            << robustProfile->pureCutCellFallbackAttempts << ",\n";
        out << "  \"h4_conformal_hybrid_build_calls\": "
            << robustProfile->conformalHybridBuildCalls << "\n";
    } else {
        out << "\n";
    }
    out << "}\n";
    if (!out.good()) {
        setError(error, "failed while writing hybrid profile JSON");
        return false;
    }
    return true;
}

} // namespace cartmesh2d
