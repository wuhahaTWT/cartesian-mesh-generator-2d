#include "cartmesh2d/fv/HybridViscousCell2D.hpp"
#include "cartmesh2d/fv/HybridCellGeometry2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(std::string("hybrid viscous cell: ")+message);}
double checked(double value){require(std::isfinite(value),"nonfinite arithmetic");return value;}
bool finite(Vector2D v){return std::isfinite(v.x)&&std::isfinite(v.y);}
}
HybridViscousCell2D::HybridViscousCell2D(const FvMesh2D& mesh,std::size_t index,double viscosity)
    :viscosity_(viscosity) {
    require(std::isfinite(viscosity)&&viscosity>0,"viscosity must be positive");
    auto geometry=hybridCellGeometry2D(mesh,index);
    area_=geometry.area;faces_=std::move(geometry.faces);gradientWeights_=std::move(geometry.gradientColumns);
    areaVectors_=std::move(geometry.areaVectors);stabilization_=std::move(geometry.stabilizationWeights);residual_=std::move(geometry.remainder);
    const auto count=faces_.size();
    const auto dimension=2*count;matrix_.resize(dimension*dimension);
    // Build each entry once and mirror exactly. The form is
    // mu*V*(2 sym(g):sym(h)-2/3 div(g)div(h)) + mu*R' diag(|S|/d) R.
    // Its strain part equals the sum of three nonnegative squares used below.
    for(std::size_t row=0;row<dimension;++row)for(std::size_t column=0;column<=row;++column) {
        const auto i=row/2,j=column/2;const auto a=gradientWeights_[i],b=gradientWeights_[j];
        double value;
        if(row%2==0&&column%2==0)value=4./3*a.x*b.x+a.y*b.y;
        else if(row%2==1&&column%2==1)value=a.x*b.x+4./3*a.y*b.y;
        else if(row%2==0)value=-2./3*a.x*b.y+a.y*b.x;
        else value=a.x*b.y-2./3*a.y*b.x;
        value*=area_;
        if(row%2==column%2)for(std::size_t k=0;k<count;++k)value+=stabilization_[k]*residual_[k*count+i]*residual_[k*count+j];
        value=checked(viscosity_*value);matrix_[row*dimension+column]=value;matrix_[column*dimension+row]=value;
    }
}

HybridViscousCellResult2D HybridViscousCell2D::evaluate(Vector2D owner,const std::vector<Vector2D>& trace) const {
    const auto count=faces_.size();require(trace.size()==count&&finite(owner),"invalid owner/trace dimensions");
    std::vector<Vector2D> delta(count),r(count);Vector2D gu{},gv{};
    for(std::size_t i=0;i<count;++i) {
        require(finite(trace[i]),"nonfinite trace velocity");delta[i]={trace[i].x-owner.x,trace[i].y-owner.y};
        gu.x+=gradientWeights_[i].x*delta[i].x;gu.y+=gradientWeights_[i].y*delta[i].x;
        gv.x+=gradientWeights_[i].x*delta[i].y;gv.y+=gradientWeights_[i].y*delta[i].y;
    }
    for(std::size_t i=0;i<count;++i)for(std::size_t j=0;j<count;++j){r[i].x+=residual_[i*count+j]*delta[j].x;r[i].y+=residual_[i*count+j]*delta[j].y;}
    const double difference=gu.x-gv.y,shear=gu.y+gv.x,divergence=gu.x+gv.y;
    HybridViscousCellResult2D result;result.outwardFlux.resize(count);
    result.dissipation=area_*(difference*difference+shear*shear+divergence*divergence/3);
    for(std::size_t i=0;i<count;++i)result.dissipation+=stabilization_[i]*dot(r[i],r[i]);
    result.dissipation=checked(viscosity_*result.dissipation);
    const double xx=viscosity_*(2*gu.x-2./3*divergence),yy=viscosity_*(2*gv.y-2./3*divergence),xy=viscosity_*shear;
    for(std::size_t i=0;i<count;++i) {
        const auto s=areaVectors_[i];Vector2D traction{xx*s.x+xy*s.y,xy*s.x+yy*s.y};
        for(std::size_t j=0;j<count;++j){const double weight=viscosity_*residual_[j*count+i]*stabilization_[j];traction.x+=weight*r[j].x;traction.y+=weight*r[j].y;}
        result.outwardFlux[i]={checked(-traction.x),checked(-traction.y),checked(-dot(trace[i],traction))};
    }
    return result;
}
} // namespace cartmesh2d::fv
