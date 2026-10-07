#include "cartmesh2d/fv/CompatibleFlowBoundary2D.hpp"
#include <stdexcept>
#include <cmath>
namespace cartmesh2d::fv {
std::vector<CompatibleBoundary2D> compatibleFlowBoundaries2D(const FvMesh2D& mesh,const std::vector<FlowBoundaryCondition2D>& conditions,double referenceSpeed) {
    validateFvMesh2D(mesh);FlowControls2D checked;checked.scenario="custom";checked.speed=referenceSpeed;checked.boundaryConditions=conditions;validateFlowBoundaryConditions2D(mesh,checked);
    std::vector<CompatibleBoundary2D> result;
    for(const auto& condition:conditions){CompatibleBoundary2D b;b.face=condition.face;
        switch(condition.kind){
        case FlowBoundaryKind2D::VelocityInlet:case FlowBoundaryKind2D::VelocityOutlet:case FlowBoundaryKind2D::Wall:case FlowBoundaryKind2D::MovingWall:
            b.value=[u=condition.velocity](Point2D){return u;};break;
        case FlowBoundaryKind2D::Symmetry:b.kind=CompatibleBoundaryKind2D::NormalVelocity;break;
        case FlowBoundaryKind2D::PressureOutlet:{
            b.kind=CompatibleBoundaryKind2D::PseudoTraction;b.rejectBackflow=true;
            const auto s=mesh.faces[condition.face].areaVector;const double length=std::hypot(s.x,s.y);const Vector2D traction{-condition.pressure*s.x/length,-condition.pressure*s.y/length};
            b.value=[traction](Point2D){return traction;};break;}
        case FlowBoundaryKind2D::PressureOpening:throw std::invalid_argument("Compatible pressure-opening inflow policy is not yet integrated");
        case FlowBoundaryKind2D::SmoothMovingWall:throw std::invalid_argument("Compatible smooth-moving-wall trace reconstruction is not yet integrated");
        }
        result.push_back(std::move(b));
    }
    return result;
}
}
