#pragma once

#include "cartmesh2d/fv/FvMesh2D.hpp"

#include <functional>
#include <string>

namespace cartmesh2d::fv {

enum class ConvectionScheme2D { Upwind, LimitedLinearUpwind };

enum class PressurePreconditioner2D { Jacobi, IncompleteCholesky0 };

struct FlowControls2D {
    std::string scenario = "external";
    double nu = .01;
    double speed = 1;
    double tolerance = 1e-6;
    std::size_t maxIterations = 1500;
    double velocityRelaxation = .6;
    double pressureRelaxation = .25;
    bool profile = false;
    double manufacturedPressureSlope = 0; // verification-only linear pressure addition
    ConvectionScheme2D convection = ConvectionScheme2D::Upwind;
    PressurePreconditioner2D pressurePreconditioner = PressurePreconditioner2D::IncompleteCholesky0;
};

struct FlowPerformance2D {
    std::size_t momentumSolves = 0;
    std::size_t momentumIterations = 0;
    std::size_t maxMomentumIterations = 0;
    std::size_t pressureSolves = 0;
    std::size_t pressureIterations = 0;
    std::size_t maxPressureIterations = 0;
    double momentumLinearSolveSeconds = 0;
    double pressureLinearSolveSeconds = 0;
    double solveSeconds = 0;
};

struct FlowIteration2D {
    std::size_t iteration = 0;
    double momentumResidual = 0;
    double continuity = 0;
    double velocityChange = 0;
    double pressureChange = 0;
};

struct FaceMomentum2D {
    double pressure = 0;
    Vector2D advection{};
    Vector2D diffusion{}; // -nu grad(U) dot S; same Laplacian flux as momentum
};

struct FlowResult2D {
    bool converged = false;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> p;
    std::vector<double> flux;
    std::vector<FaceMomentum2D> faceMomentum;
    std::vector<Vector2D> sourceIntegrals; // populated only for manufactured verification
    std::vector<FlowIteration2D> history;
    double globalImbalance = 0;
    double globalRelativeImbalance = 0;
    double forceX = 0;
    double forceY = 0;
    double pressureForceX = 0;
    double pressureForceY = 0;
    double discreteForceX = 0;
    double discreteForceY = 0;
    double domainHeight = 0;
    FlowPerformance2D performance;
};

// Fixed-grid, constant-property laminar SIMPLE, kinematic pressure.
// No turbulence, heat transport, moving mesh or compressibility.
[[nodiscard]] FlowResult2D solveIncompressible2D(
    const FvMesh2D&,
    const FlowControls2D&,
    const std::function<void(const FlowIteration2D&)>& progress = {});

}
