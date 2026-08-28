# One-shot branch applicator; retained only while the quality hardening PR is under validation.
from pathlib import Path

path = Path("src/quality/SolverTopology2D.cpp")
text = path.read_text(encoding="utf-8")


def replace_once(source: str, old: str, new: str, label: str) -> str:
    count = source.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    return source.replace(old, new, 1)


rank_struct = """struct LocalQualityRank2D {
    QualityScore2D issues;
    double maxNonOrthogonality = 0.0;
    double maxInternalSkewness = 0.0;
    double maxBoundarySkewness = 0.0;
    double maxCellAspect = 0.0;
    double negativeMinInteriorAngle = 0.0;
    double negativeMinFaceWeight = 0.0;
    double negativeMinVolumeRatio = 0.0;
};
"""
preferred_policy = rank_struct + """
[[nodiscard]] SolverQualityPolicy2D preferredSolverQualityPolicy() noexcept {
    SolverQualityPolicy2D policy;
    // cfMesh treats 65 degrees as a low-quality non-orthogonality boundary;
    // the remaining targets deliberately retain margin over the hard safety
    // gate without pretending that every geometry can reach ideal values.
    policy.maxNonOrthogonalityDeg = 65.0;
    policy.maxInternalSkewness = 3.5;
    policy.maxBoundarySkewness = 3.0;
    policy.maxConcavityDeg = 60.0;
    policy.maxCellAspect = 200.0;
    policy.minInteriorAngleDeg = 1.0;
    policy.minFaceWeight = 0.08;
    policy.minVolumeRatio = 0.02;
    return policy;
}
"""
text = replace_once(text, rank_struct, preferred_policy, "preferred policy insertion")

old_local_sig = """[[nodiscard]] std::optional<LocalQualityRank2D> localRepartitionQualityRank(
    const TopologyMesh2D& topology,const std::vector<std::size_t>& halo,
    std::size_t first,std::size_t second,
    const Polygon2D& firstPiece,const Polygon2D& secondPiece,
    const TolerancePolicy& tol) {
"""
new_local_sig = """[[nodiscard]] std::optional<LocalQualityRank2D> localRepartitionQualityRank(
    const TopologyMesh2D& topology,const std::vector<std::size_t>& halo,
    std::size_t first,std::size_t second,
    const Polygon2D& firstPiece,const Polygon2D& secondPiece,
    const TolerancePolicy& tol,const SolverQualityPolicy2D& qualityPolicy) {
"""
text = replace_once(text, old_local_sig, new_local_sig, "local rank signature")

local_start = text.index(new_local_sig)
local_end = text.index("template<class Proposal>", local_start)
local_block = text[local_start:local_end]
local_block = replace_once(
    local_block,
    "const auto quality=evaluateSolverQuality2D(patch,{},tol);",
    "const auto quality=evaluateSolverQuality2D(patch,qualityPolicy,tol);",
    "local rank quality policy",
)
text = text[:local_start] + local_block + text[local_end:]

old_timed = """[[nodiscard]] SolverQualityReport2D timedFullQuality(
    const TopologyMesh2D& topology,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool candidate) {
    const auto start=ProfileClock::now();
    auto quality=evaluateSolverQuality2D(topology,{},tol);
"""
new_timed = """[[nodiscard]] SolverQualityReport2D timedFullQuality(
    const TopologyMesh2D& topology,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool candidate,
    const SolverQualityPolicy2D& qualityPolicy = {}) {
    const auto start=ProfileClock::now();
    auto quality=evaluateSolverQuality2D(topology,qualityPolicy,tol);
"""
text = replace_once(text, old_timed, new_timed, "timed quality policy")

old_impl_sig = """SolverLocalRepartitionResult2D repartitionSolverTopologyByQualityImpl(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool useBatch,
    const std::vector<bool>& initialImmutableCells={}) {
"""
new_impl_sig = """SolverLocalRepartitionResult2D repartitionSolverTopologyByQualityImpl(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol,
    SolverTopologyProfile2D* profile,bool useBatch,
    const std::vector<bool>& initialImmutableCells={},
    const SolverQualityPolicy2D& qualityPolicy = {},
    std::size_t maximumIterations = 128U) {
"""
text = replace_once(text, old_impl_sig, new_impl_sig, "repartition implementation signature")

impl_start = text.index(new_impl_sig)
impl_end = text.index("\n} // namespace", impl_start)
impl = text[impl_start:impl_end]
impl = replace_once(
    impl,
    "for (std::size_t iteration=0;iteration<128;++iteration) {",
    "for (std::size_t iteration=0;iteration<maximumIterations;++iteration) {",
    "repartition iteration limit",
)
impl = impl.replace(
    "timedFullQuality(result.topology,tol,profile,false)",
    "timedFullQuality(result.topology,tol,profile,false,qualityPolicy)",
)
impl = impl.replace(
    "timedFullQuality(candidate.topology,tol,profile,true)",
    "timedFullQuality(candidate.topology,tol,profile,true,qualityPolicy)",
)
old_rank_call = """const auto rank=localRepartitionQualityRank(
                        result.topology,halo,first,second,firstPiece,secondPiece,tol);"""
new_rank_call = """const auto rank=localRepartitionQualityRank(
                        result.topology,halo,first,second,firstPiece,secondPiece,
                        tol,qualityPolicy);"""
impl = replace_once(impl, old_rank_call, new_rank_call, "preferred local rank call")
text = text[:impl_start] + impl + text[impl_end:]

old_final = """    const auto repartitionStart=ProfileClock::now();
    auto repartitioned=repartitionSolverTopologyByQualityImpl(
        partition.topology,domain,boundary,tol,&result.profile,true,
        partition.immutableForCell);
    result.profile.finalRepartitionSeconds=profileSeconds(repartitionStart);
    if (!repartitioned.valid()) {
        result.issues.insert(result.issues.end(),repartitioned.issues.begin(),
                             repartitioned.issues.end());
        return result;
    }
    partition.topology=std::move(repartitioned.topology);
    result.qualityRepartitionCount=repartitioned.repartitionCount;
"""
new_final = """    const auto repartitionStart=ProfileClock::now();
    auto repartitioned=repartitionSolverTopologyByQualityImpl(
        partition.topology,domain,boundary,tol,&result.profile,true,
        partition.immutableForCell);
    if (!repartitioned.valid()) {
        result.issues.insert(result.issues.end(),repartitioned.issues.begin(),
                             repartitioned.issues.end());
        return result;
    }

    // A hard-valid mesh is only the beginning. Run a bounded second pass with
    // preferred quality margins, analogous to a mature mesher's low-quality
    // face optimization stage. Every accepted candidate must improve the
    // stricter score; immutable H4 layer cells remain protected. If the
    // preferred target is geometrically unattainable, retain the best hard-
    // valid topology rather than failing or spinning indefinitely.
    std::size_t totalRepartitionCount=repartitioned.repartitionCount;
    const auto preferredPolicy=preferredSolverQualityPolicy();
    auto marginOptimized=repartitionSolverTopologyByQualityImpl(
        repartitioned.topology,domain,boundary,tol,&result.profile,true,
        repartitioned.immutableCells,preferredPolicy,16U);
    if (marginOptimized.valid()) {
        const auto hardQuality=evaluateSolverQuality2D(
            marginOptimized.topology,SolverQualityPolicy2D{},tol);
        if (hardQuality.valid()) {
            totalRepartitionCount+=marginOptimized.repartitionCount;
            repartitioned=std::move(marginOptimized);
        }
    }
    result.profile.finalRepartitionSeconds=profileSeconds(repartitionStart);
    partition.topology=std::move(repartitioned.topology);
    result.qualityRepartitionCount=totalRepartitionCount;
"""
text = replace_once(text, old_final, new_final, "preferred post-repair pass")

path.write_text(text, encoding="utf-8")
print("Applied solver quality-margin optimization patch")
