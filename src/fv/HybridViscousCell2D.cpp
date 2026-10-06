#include "cartmesh2d/fv/HybridViscousCell2D.hpp"
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
    require(index<mesh.cells.size(),"invalid cell index");
    require(std::isfinite(viscosity)&&viscosity>0,"viscosity must be positive");
    const auto& cell=mesh.cells[index];area_=cell.area;faces_=cell.faces;
    require(std::isfinite(area_)&&area_>0&&faces_.size()>=3,"invalid cell geometry");
    require(std::isfinite(cell.centre.x)&&std::isfinite(cell.centre.y),"invalid cell centre");
    const auto count=faces_.size();std::vector<Vector2D> offsets(count);
    Vector2D closure{};double perimeter=0,magnitude=area_;
    std::array<double,4> moment{};
    for(std::size_t i=0;i<count;++i) {
        require(faces_[i]<mesh.faces.size(),"invalid face index");
        require(std::find(faces_.begin(),faces_.begin()+static_cast<std::ptrdiff_t>(i),faces_[i])==faces_.begin()+static_cast<std::ptrdiff_t>(i),"repeated cell face");
        const auto& face=mesh.faces[faces_[i]];
        require(face.owner==index||(face.neighbour&&*face.neighbour==index),"unrelated face");
        const double sign=face.owner==index?1:-1;
        const Vector2D s{sign*face.areaVector.x,sign*face.areaVector.y};
        const Vector2D d{face.centre.x-cell.centre.x,face.centre.y-cell.centre.y};
        const double length=std::hypot(s.x,s.y),normal=dot(s,d)/length;
        require(finite(s)&&finite(d)&&std::isfinite(normal)&&length>0&&normal>0,"nonpositive/nonfinite face normal distance");
        offsets[i]=d;areaVectors_.push_back(s);gradientWeights_.push_back({s.x/area_,s.y/area_});
        stabilization_.push_back(checked(length/normal));
        closure.x+=s.x;closure.y+=s.y;perimeter+=length;
        moment[0]+=s.x*d.x;moment[1]+=s.x*d.y;moment[2]+=s.y*d.x;moment[3]+=s.y*d.y;
        magnitude+=length*std::hypot(d.x,d.y);
    }
    // Dimensionless roundoff budget for geometric identities, scaled by the
    // actual lengths/moments accumulated above. No physical-error tolerance and
    // no area floor: a resolvable small cell must satisfy the same identities.
    const double budget=1024*std::numeric_limits<double>::epsilon()*static_cast<double>(count);
    require(std::hypot(closure.x,closure.y)<=budget*perimeter,"cell normals do not close");
    require(std::max({std::abs(moment[0]-area_),std::abs(moment[1]),std::abs(moment[2]),std::abs(moment[3]-area_)})<=budget*magnitude,"inconsistent cell moment geometry");
    residual_.resize(count*count);
    for(std::size_t i=0;i<count;++i)for(std::size_t j=0;j<count;++j)
        residual_[i*count+j]=checked((i==j?1.:0.)-dot(offsets[i],gradientWeights_[j]));
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
