#pragma once

#include "cartmesh2d/fv/Incompressible2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cartmesh2d::fv {

// Compact C2 streamfunction; positive peakSpeed is counterclockwise. The
// complete support disk must lie strictly inside the final fluid domain.
// This is an initial condition, not a body force or a repeated time-step source.
struct FlowInitialVortex2D {
    Point2D centre;
    double radius = 0;
    double peakSpeed = 0;
};

inline void validateFlowInitialVortex2D(const FlowInitialVortex2D& vortex) {
    if (!std::isfinite(vortex.centre.x) || !std::isfinite(vortex.centre.y) ||
        !std::isfinite(vortex.radius) || vortex.radius<=0 || !std::isfinite(vortex.peakSpeed))
        throw std::runtime_error("Initial vortex requires finite centre/speed and positive radius");
}

struct FlowVortexSample2D { Vector2D velocity; double streamfunction=0; };
inline FlowVortexSample2D sampleInitialVortex2D(Point2D point, const FlowInitialVortex2D& vortex) {
    validateFlowInitialVortex2D(vortex);
    if (!std::isfinite(point.x) || !std::isfinite(point.y))
        throw std::runtime_error("Initial vortex sample point must be finite");
    const double x=(point.x-vortex.centre.x)/vortex.radius;
    const double y=(point.y-vortex.centre.y)/vortex.radius;
    const double radius=std::hypot(x,y);
    if (radius>=1) return {};
    const double s=1-radius*radius;
    const double scale=vortex.peakSpeed*(25*std::sqrt(5.)/16)*s*s;
    FlowVortexSample2D sample{{-scale*y,scale*x},
        vortex.peakSpeed*vortex.radius*(25*std::sqrt(5.)/96)*s*s*s};
    if (!std::isfinite(sample.velocity.x) || !std::isfinite(sample.velocity.y) || !std::isfinite(sample.streamfunction))
        throw std::runtime_error("Initial vortex numerical range exceeded");
    return sample;
}

inline FlowState2D withInitialVortex2D(const FvMesh2D& mesh, const FlowState2D& initial,
                                      const FlowInitialVortex2D& vortex) {
    validateFvMesh2D(mesh);
    validateFlowInitialVortex2D(vortex);
    if (initial.time!=0 || initial.u.size()!=mesh.cells.size() || initial.v.size()!=mesh.cells.size() ||
        initial.p.size()!=mesh.cells.size() || initial.flux.size()!=mesh.faces.size())
        throw std::runtime_error("Initial vortex requires a complete state at physical time zero");
    for (const auto* values:{&initial.u,&initial.v,&initial.p,&initial.flux})
        for (double value:*values) if (!std::isfinite(value))
            throw std::runtime_error("Initial vortex received nonfinite initial state");
    // Boundary faces are oriented with fluid to the left. Their winding sum
    // includes outer walls and holes, so centres inside solids are rejected.
    int winding=0;
    double clearance=std::numeric_limits<double>::infinity();
    const auto endpoints=[](const Face& face) {
        return std::pair{Point2D{face.centre.x+.5*face.areaVector.y,face.centre.y-.5*face.areaVector.x},
                         Point2D{face.centre.x-.5*face.areaVector.y,face.centre.y+.5*face.areaVector.x}};
    };
    for (const auto& face:mesh.faces) if (!face.neighbour) {
        const auto [a,b]=endpoints(face);
        const auto d=b-a, q=vortex.centre-a;
        const double length=std::hypot(d.x,d.y);
        const double along=std::clamp((q.x*(d.x/length)+q.y*(d.y/length))/length,0.,1.);
        clearance=std::min(clearance,std::hypot(q.x-along*d.x,q.y-along*d.y));
        const double cross=d.x*q.y-d.y*q.x;
        if (a.y<=vortex.centre.y && b.y>vortex.centre.y && cross>0) ++winding;
        if (a.y>vortex.centre.y && b.y<=vortex.centre.y && cross<0) --winding;
    }
    const double eps=16*TolerancePolicy{}.scale(std::max({1.,vortex.radius,std::abs(vortex.centre.x),std::abs(vortex.centre.y)}));
    if (winding!=1 || !(clearance>vortex.radius+eps))
        throw std::runtime_error("Initial vortex support disk must lie strictly inside fluid, clear of all walls and openings");
    auto result=initial;
    bool velocityResolved=false,fluxResolved=false;
    for (std::size_t i=0;i<mesh.cells.size();++i) {
        const auto sample=sampleInitialVortex2D(mesh.cells[i].centre,vortex);
        velocityResolved|=sample.velocity.x!=0 || sample.velocity.y!=0;
        result.u[i]+=sample.velocity.x;result.v[i]+=sample.velocity.y;
        if (!std::isfinite(result.u[i]) || !std::isfinite(result.v[i]))
            throw std::runtime_error("Initial vortex velocity overflow");
    }
    for (std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& face=mesh.faces[id];
        if (!face.neighbour) continue; // support is strictly away from boundaries
        const auto [a,b]=endpoints(face);
        const double flux=sampleInitialVortex2D(b,vortex).streamfunction-sampleInitialVortex2D(a,vortex).streamfunction;
        fluxResolved|=flux!=0;
        result.flux[id]+=flux;
        if (!std::isfinite(result.flux[id])) throw std::runtime_error("Initial vortex flux overflow");
    }
    if (vortex.peakSpeed!=0 && (!velocityResolved || !fluxResolved))
        throw std::runtime_error("Initial vortex is unresolved by cell centres or face vertices; refine the mesh or increase support radius");
    return result;
}
} // namespace cartmesh2d::fv
