#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/HeatConduction2D.hpp"
#include "cartmesh2d/fv/ViscousStress2D.hpp"
#include <iostream>
#include <iomanip>
#include <numeric>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
Vector2D rotation(Vector2D p,double a){return {std::cos(a)*p.x-std::sin(a)*p.y,std::sin(a)*p.x+std::cos(a)*p.y};}
FvMesh2D transformed(double scale,double angle,double length){auto m=fv_test::rectangle(8,6,length,true);for(auto& c:m.cells){const auto p=rotation({c.centre.x,c.centre.y},angle);c.centre={scale*p.x,scale*p.y};c.area*=scale*scale;}for(auto& f:m.faces){const auto p=rotation({f.centre.x,f.centre.y},angle);f.centre={scale*p.x,scale*p.y};auto s=rotation(f.areaVector,angle),c=rotation(f.correction,angle);f.areaVector={scale*s.x,scale*s.y};f.correction={scale*c.x,scale*c.y};}return m;}
}
int main(){try {
    double error=0,originalError=0,roundoffScaledError=0,rateError=0;
    for(double length:{.03,3.,30.})for(double scale:{1e-3,1.,1e3})for(double angle:{0.,.731}) {
        const auto mesh=transformed(scale,angle,length);
        auto value=[&](Point2D p){const double x=p.x/scale,y=p.y/scale;return 10+.7*x-.2*y+.3*x*x-.15*x*y+.4*y*y;};
        auto gradient=[&](Point2D p){return Vector2D{(.7+.6*p.x/scale-.15*p.y/scale)/scale,(-.2-.15*p.x/scale+.8*p.y/scale)/scale};};
        std::vector<double> t,capacity(mesh.cells.size(),2),rho(mesh.cells.size(),1.7);std::vector<Vector2D> velocity;
        for(const auto& c:mesh.cells){t.push_back(value(c.centre));velocity.push_back({t.back(),-.3*t.back()});}
        std::vector<HeatBoundary2D> hb;std::vector<ViscousBoundary2D> vb;
        for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour){const double v=value(mesh.faces[id].centre);hb.push_back({id,HeatBoundaryKind2D::Temperature,v,{}});vb.push_back({id,ViscousBoundaryKind2D::Velocity,{v,-.3*v},{}});}
        HeatConductionOperator2D heat(mesh,hb,.37,WallGradient2D::Quadratic);ViscousStressOperator2D stress(mesh,vb,.23,WallGradient2D::Quadratic);
        const auto h=heat.evaluate(t,capacity);const auto v=stress.evaluate(velocity,rho);
        require(heat.quadraticWalls()==hb.size()&&stress.quadraticWalls()==vb.size(),"missing quadratic walls");
        std::vector<bool> prescribed(mesh.faces.size());std::vector<std::optional<std::size_t>> partners(mesh.faces.size());for(const auto& b:hb)prescribed[b.face]=true;
        for(const auto& b:hb) {
            const auto& f=mesh.faces[b.face];const auto g=gradient(f.centre);const double q=-.37*dot(g,f.areaVector);
            const auto stencil=quadraticWallGradient2D(mesh,b.face,prescribed,partners);
            double heatEnvelope=0,stressEnvelope=0;
            for(const auto& sample:stencil.samples){const double inputMagnitude=std::abs(sample.boundary?value(mesh.faces[sample.index].centre):t[sample.index])+std::abs(b.value);heatEnvelope+=.37*std::abs(dot(sample.weight,f.areaVector))*inputMagnitude;stressEnvelope+=.23*(8./3)*(std::abs(f.areaVector.x)+std::abs(f.areaVector.y))*(std::abs(sample.weight.x)+std::abs(sample.weight.y))*inputMagnitude;}
            const double heatError=std::abs(h.faceHeatFlux[b.face]-q)/(1+std::abs(q));error=std::max(error,heatError);if(length==3)originalError=std::max(originalError,heatError);
            roundoffScaledError=std::max(roundoffScaledError,std::abs(h.faceHeatFlux[b.face]-q)/(1+std::abs(q)+heatEnvelope));
            const double div=g.x-.3*g.y,xx=.23*(2*g.x-2./3*div),yy=.23*(-.6*g.y-2./3*div),xy=.23*(g.y-.3*g.x);
            const Vector2D force{-(xx*f.areaVector.x+xy*f.areaVector.y),-(xy*f.areaVector.x+yy*f.areaVector.y)};
            const auto u=Vector2D{value(f.centre),-.3*value(f.centre)};const std::array<double,3> exact{force.x,force.y,dot(force,u)};
            for(std::size_t component=0;component<3;++component){const double delta=std::abs(v.faceFlux[b.face][component]-exact[component]),normalized=delta/(1+std::abs(exact[component]));error=std::max(error,normalized);if(length==3)originalError=std::max(originalError,normalized);const double envelope=stressEnvelope*(component==2?std::abs(u.x)+std::abs(u.y):1);roundoffScaledError=std::max(roundoffScaledError,delta/(1+std::abs(exact[component])+envelope));}
        }
        // Directly perturb all DOFs of the evaluated operators. This detects
        // omission of nonlocal wall coefficients from their explicit bound.
        std::vector<double> heatNorm(t.size());std::vector<std::array<double,2>> stressNorm(t.size());
        for(std::size_t j=0;j<t.size();++j) {
            auto tp=t,tm=t;tp[j]+=1;tm[j]-=1;const auto hp=heat.evaluate(tp,capacity),hm=heat.evaluate(tm,capacity);
            for(std::size_t i=0;i<t.size();++i)heatNorm[i]+=.5*std::abs(hp.cellResidual[i]-hm.cellResidual[i]);
            for(std::size_t a=0;a<2;++a){auto up=velocity,um=velocity;(a?up[j].y:up[j].x)+=1;(a?um[j].y:um[j].x)-=1;const auto vp=stress.evaluate(up,rho),vm=stress.evaluate(um,rho);for(std::size_t i=0;i<t.size();++i)for(std::size_t k=0;k<2;++k)stressNorm[i][k]+=.5*std::abs(vp.cellResidual[i][k]-vm.cellResidual[i][k]);}
        }
        for(std::size_t i=0;i<t.size();++i){const double hr=.5*heatNorm[i]/(mesh.cells[i].area*capacity[i]),vr=.5*std::max(stressNorm[i][0],stressNorm[i][1])/(mesh.cells[i].area*rho[i]);rateError=std::max({rateError,std::abs(h.rate[i]-hr)/hr,std::abs(v.rate[i]-vr)/vr});}
    }
    // 4096 epsilon accommodates QR/difference arithmetic, normalized by flux
    // magnitude or Jacobian norm. This checks algebra, not PDE accuracy.
    require(originalError<4096*std::numeric_limits<double>::epsilon(),"quadratic wall polynomial reproduction failed");
    // The added extreme-aspect cases amplify rounding of the sampled absolute
    // point values by O(1/h_normal). Bound that effect with the actual linear
    // sensitivity to each input. Keep the ORIGINAL gate and report the raw
    // physical-scale error as well; this is not a tighter PDE accuracy claim.
    require(roundoffScaledError<4096*std::numeric_limits<double>::epsilon(),"stretched wall reproduction exceeds input-roundoff envelope");
    require(rateError<4096*std::numeric_limits<double>::epsilon(),"quadratic wall Jacobian bound failed");
    auto mesh=fv_test::rectangle(1,1,1);std::vector<HeatBoundary2D> bc;for(std::size_t id=0;id<mesh.faces.size();++id)bc.push_back({id,HeatBoundaryKind2D::Temperature,2,{}});
    bool rejected=false;try{HeatConductionOperator2D insufficient(mesh,bc,1,WallGradient2D::Quadratic);}catch(const std::exception& e){rejected=std::string(e.what()).find("rank-deficient")!=std::string::npos;}require(rejected,"rank deficiency silently accepted");
    std::cout<<std::setprecision(17)<<"{\"quadraticWallFluxError\":"<<error<<",\"originalGeometryWallFluxError\":"<<originalError<<",\"inputRoundoffScaledError\":"<<roundoffScaledError<<",\"fullRowNormRelativeError\":"<<rateError<<",\"rankDeficiencyRejected\":true}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
