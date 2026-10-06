#pragma once
#include "cartmesh2d/fv/HeatConduction2D.hpp"
#include <functional>
namespace cartmesh2d::fv {
// Research HMM operator, selected only by optional Euler research controls.
// Cell temperatures and one actual-face trace. Constant isotropic k, 2D per depth.
struct HybridHeatDiagnostics2D {
    double closureRoundoffRatio=0,momentRoundoffRatio=0,minimumNormalDistance=0;
    double maximumLocalAffineResidual=0;
};
struct HybridHeatResult2D {
    std::vector<double> temperature,stageTemperature,trace,faceHeatFlux,cellResidual;
    // Dissipation/trace/flux refer to the implicit stage; temperature is final.
    // All cell balances and flux jumps: W/m per depth; temperature gaps: K.
    double dissipation=0,maximumFaceFluxJump=0,maximumCellBalance=0;
    double maximumLocalCellBalance=0,maximumSharedStageCellBalance=0;
    double maximumStageToFinalDifferenceK=0,maximumStageToFinalBudgetRatio=0;
    double traceResidualNorm=0,traceRhsNorm=0;
    std::size_t iterations=0,traceUnknowns=0;
};
class HybridHeatOperator2D {
public:
    HybridHeatOperator2D(const FvMesh2D&,const std::vector<HeatBoundary2D>&,double conductivity);
    [[nodiscard]] HybridHeatResult2D solveSteady(double relativeTolerance=1e-12) const;
    // Uniform physical dt, first order. Positive heat capacity required.
    // Stage trace/flux retained; final cell values use the unique shared flux.
    // Stage-to-final difference must fit the actual linear/roundoff budget.
    // Candidate returned by value; caller keeps previous state on solve failure.
    // Coercivity is not an M-matrix/point-positivity guarantee. No clipping.
    [[nodiscard]] HybridHeatResult2D backwardEuler(const std::vector<double>& previous,
        const std::vector<double>& volumetricHeatCapacity,double dt,double relativeTolerance=1e-12) const;
    // Static global trace equilibration for supplied cell values. The optional
    // Euler research path pays this global assembly/solve cost per residual.
    [[nodiscard]] HybridHeatResult2D evaluateAtCells(const std::vector<double>& temperature,double relativeTolerance=1e-12) const;
    [[nodiscard]] const HybridHeatDiagnostics2D& diagnostics() const {return diagnostics_;}
    // Native assembled condensed matrix for linear-algebra diagnostics only.
    using MatrixVisitor=std::function<void(std::size_t,std::size_t,double)>;
    void visitSteadyTraceMatrix(const MatrixVisitor&) const;
    void visitLocalMatrices(const std::function<void(std::size_t,std::size_t,std::size_t,double)>&) const;
private:
    struct Local {std::vector<std::size_t> faces;std::vector<double> a,b;double sum=0,area=0;};
    struct Boundary {bool fixed=false;double temperature=0,flux=0;};
    std::vector<Local> local_;
    std::vector<Boundary> boundary_;
    std::vector<std::optional<std::size_t>> neighbour_;
    std::vector<std::size_t> owner_,unknown_;
    std::size_t count_=0;
    HybridHeatDiagnostics2D diagnostics_;
    HybridHeatResult2D solve(const std::vector<double>& previous,const std::vector<double>& mass,
        double relativeTolerance,const MatrixVisitor* visitor=nullptr,bool fixedCells=false) const;
};
}
