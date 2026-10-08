#include "FlowSolverDetail2D.hpp"
#include <algorithm>
#include <map>
#include <limits>
namespace cartmesh2d::fv::solver_detail {
void postprocessForces(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,FlowResult2D& r,
    const Vec& pressureBoundary,const std::vector<Vector2D>& gp,const std::vector<Vector2D>& gu,
    const std::vector<Vector2D>& gv,const std::vector<Vector2D>& stressCorrection){
    const auto nf=m.faces.size();
    const auto pf=detail::pressureFaceValues(m,r.p,gp,pressureBoundary,b.fixedP);
    const auto lu=c.convection==ConvectionScheme2D::LimitedLinearUpwind
        ? detail::faceReconstructionLimiter(m,r.u,gu,b.u,b.fixedU) : Vec{};
    const auto lv=c.convection==ConvectionScheme2D::LimitedLinearUpwind
        ? detail::faceReconstructionLimiter(m,r.v,gv,b.v,b.fixedV) : Vec{};
    const auto faceVelocity=c.convection==ConvectionScheme2D::FaceLimitedLinearUpwind
        ? detail::faceFrameVelocityValues(m,r.flux,r.u,r.v,gu,gv,b.u,b.v,b.fixedU,b.fixedV) : std::vector<Vector2D>{};
    r.faceMomentum.resize(nf);
    for(std::size_t id=0;id<nf;++id) {
        const auto& f=m.faces[id]; const auto i=f.owner;
        auto& fm=r.faceMomentum[id]; fm.pressure=finite(pf[id]);
        if (!f.neighbour && b.role[id]==Role::Outlet && r.flux[id]<0) {
            ++r.outletBackflowFaces;
            r.outletInflow=finite(r.outletInflow-r.flux[id]);
        }
        auto component=[&](const Vec& value,const std::vector<Vector2D>& g,
                           const Vec& bc,const std::vector<bool>& fixed,const Vec& limiter,bool y) {
            const bool normalInflow = !f.neighbour && !fixed[id] && r.flux[id]<0 &&
                (b.role[id]==Role::Farfield || b.role[id]==Role::Opening || (b.role[id]==Role::Outlet && c.outletBackflow==OutletBackflow2D::NormalInlet));
            const double faceValue=(!f.neighbour && fixed[id]) ? bc[id]
                : normalInflow ? value[i] : !faceVelocity.empty() ? (y?faceVelocity[id].y:faceVelocity[id].x)
                : detail::upwindFaceValue(m,id,r.flux[id],value,g,limiter);
            double diffusion=0;
            if(f.neighbour || fixed[id]) {
                const double other=f.neighbour ? value[*f.neighbour] : bc[id];
                diffusion=-faceNu(c,id)*(f.transmissibility*(other-value[i])+dot(interpolateGradient(f,g),f.correction));
            }
            return std::pair{finite(c.momentumInertia*r.flux[id]*faceValue),finite(diffusion)};
        };
        auto [ax,dx]=component(r.u,gu,b.u,b.fixedU,lu,false);
        auto [ay,dy]=component(r.v,gv,b.v,b.fixedV,lv,true);
        if (!stressCorrection.empty()) {dx+=stressCorrection[id].x;dy+=stressCorrection[id].y;}
        fm.advection={ax,ay}; fm.diffusion={dx,dy};
        fm.wall=!f.neighbour && (b.role[id]==Role::Wall || b.role[id]==Role::Lid);
        if (fm.wall) {
            r.wallForceX+=pf[id]*f.areaVector.x+dx;
            r.wallForceY+=pf[id]*f.areaVector.y+dy;
            r.wallViscousForceX+=dx; r.wallViscousForceY+=dy;
        }
        if(f.neighbour || f.patch!=BoundaryPatch2D::EmbeddedBoundary || b.role[id]!=Role::Wall) continue;
        const double px=pf[id]*f.areaVector.x,py=pf[id]*f.areaVector.y;
        r.pressureForceX+=px; r.pressureForceY+=py;
        r.discreteForceX+=px+dx; r.discreteForceY+=py+dy;
        // Preserve the old cell-gradient diagnostic for explicit comparison;
        // symmetric mode reports the actual shared-face stress above as force.
        r.reconstructedForceX+=px-faceNu(c,id)*(2*gu[i].x*f.areaVector.x+(gu[i].y+gv[i].x)*f.areaVector.y);
        r.reconstructedForceY+=py-faceNu(c,id)*((gu[i].y+gv[i].x)*f.areaVector.x+2*gv[i].y*f.areaVector.y);
    }
    r.forceX=c.viscousStress==ViscousStress2D::Symmetric?r.discreteForceX:r.reconstructedForceX;
    r.forceY=c.viscousStress==ViscousStress2D::Symmetric?r.discreteForceY:r.reconstructedForceY;
    finite(r.reconstructedForceX);finite(r.reconstructedForceY);
    finite(r.wallForceX);finite(r.wallForceY);
    finite(r.wallViscousForceX);finite(r.wallViscousForceY);
    finite(r.globalImbalance);finite(r.forceX);finite(r.forceY);
    finite(r.pressureForceX);finite(r.pressureForceY);finite(r.discreteForceX);finite(r.discreteForceY);
    if (c.scenario == "custom") {
        std::map<std::string, FlowWallLoad2D> loads;
        // Sort face IDs so input record ordering cannot change reductions.
        auto conditions = c.boundaryConditions;
        std::sort(conditions.begin(), conditions.end(), [](const auto& a, const auto& b) { return a.face < b.face; });
        for (const auto& condition : conditions) {
            const auto id = condition.face;
            const auto& momentum = r.faceMomentum[id];
            if (!momentum.wall) continue;
            const auto& face = m.faces[id];
            auto& load = loads[condition.name];
            load.name = condition.name;
            ++load.faces;
            load.length = finite(load.length + std::hypot(face.areaVector.x, face.areaVector.y));
            const auto pressure = face.areaVector * momentum.pressure;
            load.pressure.x = finite(load.pressure.x + pressure.x);
            load.pressure.y = finite(load.pressure.y + pressure.y);
            load.viscous.x = finite(load.viscous.x + momentum.diffusion.x);
            load.viscous.y = finite(load.viscous.y + momentum.diffusion.y);
            load.pressureTorque = finite(load.pressureTorque + face.centre.x*pressure.y - face.centre.y*pressure.x);
            load.viscousTorque = finite(load.viscousTorque + face.centre.x*momentum.diffusion.y - face.centre.y*momentum.diffusion.x);
        }
        for (auto& entry : loads) r.namedWallLoads.push_back(std::move(entry.second));
    }
}
}
