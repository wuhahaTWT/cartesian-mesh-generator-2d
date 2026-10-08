#include "FlowSolverDetail2D.hpp"
#include <algorithm>
#include <map>
#include <limits>
namespace cartmesh2d::fv::solver_detail {
Vec rhieChowFlux(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,const FlowResult2D& r,
    const Vec& ra,const Vec& oldU,const Vec& oldV,const std::vector<Vector2D>& gu,
    const std::vector<Vector2D>& gv,const std::vector<Vector2D>& gup,const std::vector<Vector2D>& gvp,
    const std::vector<Vector2D>& gp,const std::vector<Vector2D>& forceGradient,const Vec& oldFluxDefect,
    bool previous,double timeStep,Vec& df) {
    const auto nf=m.faces.size();
        Vec predicted(nf);
        for(std::size_t id=0;id<nf;++id){const auto&f=m.faces[id];const auto i=f.owner;
            const double rf=interpolate(f,ra);df[id]=rf*f.transmissibility;
            if(f.neighbour){const auto j=*f.neighbour;const double w=f.neighbourWeight;
                const Point2D point{m.cells[i].centre.x*(1-w)+m.cells[j].centre.x*w,m.cells[i].centre.y*(1-w)+m.cells[j].centre.y*w};
                const auto skew=f.centre-point;
                const double uf=interpolate(f,r.u)+dot(interpolateGradient(f,gup),skew),vf=interpolate(f,r.v)+dot(interpolateGradient(f,gvp),skew);
                const Vector2D rag{(1-w)*ra[i]*forceGradient[i].x+w*ra[j]*forceGradient[j].x,(1-w)*ra[i]*forceGradient[i].y+w*ra[j]*forceGradient[j].y};
                predicted[id]=uf*f.areaVector.x+vf*f.areaVector.y+dot(rag,f.areaVector)-rf*(f.transmissibility*(r.p[j]-r.p[i])+dot(interpolateGradient(f,gp),f.correction));
                if (previous) {
                    const double oldUf=interpolate(f,oldU)+dot(interpolateGradient(f,gu),skew);
                    const double oldVf=interpolate(f,oldV)+dot(interpolateGradient(f,gv),skew);
                    // Backward Euler transports the accepted old face flux.
                    // Correct both its cell-interpolation defect and the inner
                    // under-relaxation defect; otherwise the dt -> 0 response
                    // depends on the arbitrary inner relaxation factor.
                    predicted[id]+=rf/timeStep*oldFluxDefect[id]
                        +(1-c.velocityRelaxation)*(r.flux[id]-oldUf*f.areaVector.x-oldVf*f.areaVector.y);
                } else {
                    // The face equation must retain the same implicit
                    // relaxation as the cell momentum equation. Otherwise
                    // the converged Rhie--Chow flux depends on alphaU.
                    const double oldUf=interpolate(f,oldU)+dot(interpolateGradient(f,gu),skew);
                    const double oldVf=interpolate(f,oldV)+dot(interpolateGradient(f,gv),skew);
                    predicted[id]+=(1-c.velocityRelaxation)*
                        (r.flux[id]-oldUf*f.areaVector.x-oldVf*f.areaVector.y);
                }
            }else if(b.role[id]==Role::Outlet || b.role[id]==Role::Opening || b.role[id]==Role::Farfield){
                const double pressureDifference = c.scenario == "custom"
                    ? f.transmissibility*(b.p[id]-r.p[i]) : -f.transmissibility*r.p[i];
                predicted[id]=r.u[i]*f.areaVector.x+r.v[i]*f.areaVector.y+ra[i]*dot(forceGradient[i],f.areaVector)-rf*(pressureDifference+dot(gp[i],f.correction));
                if (previous) predicted[id]+=rf/timeStep*oldFluxDefect[id]
                    +(1-c.velocityRelaxation)*(r.flux[id]-oldU[i]*f.areaVector.x-oldV[i]*f.areaVector.y);
                else predicted[id]+=(1-c.velocityRelaxation)*
                    (r.flux[id]-oldU[i]*f.areaVector.x-oldV[i]*f.areaVector.y);}
            else if(b.role[id]==Role::Inlet)predicted[id]=b.u[id]*f.areaVector.x+b.v[id]*f.areaVector.y;
            else predicted[id]=0;
        }
    return predicted;
}
}
