#pragma once
#include "cartmesh2d/fv/Sst2003m2D.hpp"

namespace cartmesh2d::fv {
struct SstRansControls2D {
    FlowControls2D flow;
    SstTransportControls2D turbulence;
    double inletK=.001, inletOmega=2;
    // Interleave a bounded number of nonlinear SST updates with each SIMPLE
    // iteration. Final original-equation stopping gates are unchanged. Set to
    // turbulence.maxIterations to recover the fully nested solve strategy.
    std::size_t turbulenceUpdatesPerIteration=1;
};
struct SstRansIteration2D {
    std::size_t iteration=0, turbulenceIterations=0;
    double kNorm=0, omegaNorm=0, kCellResidual=0, omegaCellResidual=0;
};
// Enabled by flow.profile. updateSeconds includes the other SST timers;
// never sum them as disjoint phases. Failed calls do not return a result.
struct SstRansPerformance2D {
    ScalarTransportPerformance2D scalarSolves, scalarEvaluations;
    ScalarTransportPerformance2D kSolves, omegaSolves; // disjoint parts of scalarSolves
    std::size_t updates=0;
    double updateSeconds=0, transportSeconds=0, gradientSeconds=0, wallDistanceSeconds=0;
};
struct SstRansResult2D {
    SstRansPerformance2D performance;
    bool converged=false;
    FlowResult2D flow;
    SstTransportResult2D turbulence;
    std::vector<double> faceViscosity, wallDistance;
    std::vector<bool> resolvedWalls;
    std::vector<ScalarBoundary2D> boundaryK,boundaryOmega;
    std::vector<SstVelocityBoundary2D> velocityBoundary;
    std::vector<SstRansIteration2D> history;
};
// Experimental steady, incompressible SST-2003m. TMR's m variant omits the
// isotropic 2/3 k stress: pressure is p/rho, NOT p/rho+2k/3. Complete symmetric
// stress uses nu+nu_t on each shared face, and molecular nu at resolved walls.
// No wall function, heat feedback, URANS or arbitrary boundary patches here.
// Inlet values and optional initial fields are SI; no silent positivity floor.
[[nodiscard]] SstRansResult2D solveSstRans2D(const FvMesh2D&,const SstRansControls2D&,
    const std::vector<double>& initialK={},const std::vector<double>& initialOmega={},
    const std::function<void(const FlowIteration2D&)>& progress={});
}
