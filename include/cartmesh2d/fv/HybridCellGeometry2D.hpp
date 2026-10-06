#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"
namespace cartmesh2d::fv {
// Shared native 2D HMM geometry. Caller must validate the full FvMesh first.
// No material coefficient or unknown/state; eligible for heat and vector stress.
struct HybridCellGeometry2D {
    std::vector<std::size_t> faces;
    std::vector<Vector2D> areaVectors,offsets,gradientColumns;
    std::vector<double> remainder,stabilizationWeights; // R row-major; |S|/d
    double area=0,minimumNormalDistance=0;
    double closureRoundoffRatio=0,momentRoundoffRatio=0;
};
[[nodiscard]] HybridCellGeometry2D hybridCellGeometry2D(const FvMesh2D&,std::size_t cell);
}
