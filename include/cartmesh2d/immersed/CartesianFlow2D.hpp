#pragma once
#include "cartmesh2d/geometry/Geometry2D.hpp"
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace cartmesh2d::immersed {
// A separate computational lattice, including auxiliary solid degrees of freedom.
// This is deliberately NOT FvMesh2D or the product background-grid file format.
struct Controls {
    std::size_t nx = 128, ny = 32, maxSteps = 40000;
    double length = 4, height = 1, viscosity = .01, drive = .12;
    double penaltyTime = 1e-4, maskHalfWidthCells = .5;
    double steadyTolerance = 1e-4, continuityTolerance = 1e-6;
    double linearTolerance = 1e-8, maxTimeStep = .02;
    bool systemCholesky = false;
};
struct Grid {
    Controls controls;
    std::vector<BoundaryLoop> solids;
    std::vector<double> maskU, maskV, maskCell;
    std::vector<unsigned> classification; // 0 outside, 1 inside, 2 intersected
    double solidArea = 0, gridSeconds = 0, boundarySeconds = 0;
    double dx() const { return controls.length / static_cast<double>(controls.nx); }
    double dy() const { return controls.height / static_cast<double>(controls.ny); }
    // u: nx*ny, x-periodic. v: nx*(ny+1), with the two wall rows retained.
    std::size_t index(std::size_t i, std::size_t j) const { return j*controls.nx+i; }
};
struct State {
    std::vector<double> u, v, p; // p is periodic kinematic pressure fluctuation
    std::size_t steps = 0;
    double pseudoTime = 0;
};
struct Metrics {
    double momentum = 0, continuity = 0, fieldChange = 0;
    double pressureLinearResidual = 0, pressureLinearRhsNorm = 0;
    double meanVelocity = 0, maxSpeed = 0, fluxSpread = 0;
    double wallSpeed = 0, wallNormalSpeed = 0, wallTangentialSpeed = 0;
    double deepSolidSpeed = 0, penaltyDrag = 0, penaltyLift = 0;
    double channelRelativeL2 = 0, dt = 0;
    std::size_t pressureIterations = 0;
};
struct Result {
    State state;
    Metrics metrics;
    bool converged = false;
    std::string stopReason = "iteration-limit", error;
    double solveSeconds = 0, pressureSeconds = 0;
    std::size_t pressureIterations = 0;
};
struct WallSample { Point2D point; Vector2D normal, velocity; double length = 0; };
Grid makeGrid(const Controls&, const std::vector<BoundaryLoop>& solids = {});
State zeroState(const Grid&);
std::vector<WallSample> sampleWalls(const Grid&, const State&);
Metrics evaluate(const Grid&, const State&);
// Each callback receives an accepted state only. A rejected candidate never
// replaces the state returned to callers. The initial state is not a checkpoint.
Result solve(const Grid&, const std::function<bool(const State&, const Metrics&)>& progress = {});
}
