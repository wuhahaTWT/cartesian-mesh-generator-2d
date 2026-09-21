#pragma once

#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>
#include <utility>

namespace cartmesh2d::fv::detail {

// Correct the normal derivative to the compact two-point/non-orthogonal
// derivative. Constant Dirichlet traces additionally have zero tangential
// derivative (stationary/moving straight walls, fixed slip-normal velocity).
inline Vector2D viscousFaceGradient(const FvMesh2D& m, std::size_t id,
    const std::vector<double>& value, const std::vector<Vector2D>& gradient,
    const std::vector<double>& boundary, const std::vector<bool>& fixed,
    const std::vector<bool>& constantTrace) {
    const auto& f=m.faces[id];
    const auto i=f.owner;
    const double area=std::hypot(f.areaVector.x,f.areaVector.y);
    const auto normal=f.areaVector*(1/area);
    auto g=gradient[i];
    if (f.neighbour) {
        const auto j=*f.neighbour;
        g={g.x*(1-f.neighbourWeight)+gradient[j].x*f.neighbourWeight,
           g.y*(1-f.neighbourWeight)+gradient[j].y*f.neighbourWeight};
    } else if (constantTrace[id]) {
        g={};
    }
    const double gn=g.x*normal.x+g.y*normal.y;
    if (!f.neighbour && !fixed[id]) return {g.x-normal.x*gn,g.y-normal.y*gn};
    const auto d=(f.neighbour?m.cells[*f.neighbour].centre:f.centre)-m.cells[i].centre;
    const double dn=d.x*normal.x+d.y*normal.y;
    if (!(dn>0)) throw std::runtime_error("Viscous face has non-positive normal distance");
    const double other=f.neighbour?value[*f.neighbour]:boundary[id];
    const double correction=(other-value[i]-g.x*d.x-g.y*d.y)/dn;
    return {g.x+normal.x*correction,g.y+normal.y*correction};
}

// Explicit addition to the existing -nu grad(U).S compact flux. It includes
// the transpose gradient and, at constant wall traces, removes the obsolete
// cell-gradient tangential correction. One shared value per face conserves
// both components across internal interfaces.
inline std::vector<Vector2D> symmetricViscousCorrection(const FvMesh2D& m,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<Vector2D>& gu, const std::vector<Vector2D>& gv,
    const std::vector<double>& bu, const std::vector<double>& bv,
    const std::vector<bool>& fu, const std::vector<bool>& fv,
    const std::vector<bool>& constantU, const std::vector<bool>& constantV, double nu,
    const std::vector<double>& faceViscosity = {}) {
    if (!faceViscosity.empty() && faceViscosity.size()!=m.faces.size())
        throw std::runtime_error("Viscous face coefficient dimensions mismatch");
    std::vector<Vector2D> result(m.faces.size());
    for (std::size_t id=0;id<m.faces.size();++id) {
        const double viscosity=faceViscosity.empty()?nu:faceViscosity[id];
        if (!std::isfinite(viscosity) || viscosity<=0)
            throw std::runtime_error("Viscous face coefficient must be finite positive");
        const auto& f=m.faces[id];
        const auto a=viscousFaceGradient(m,id,u,gu,bu,fu,constantU);
        const auto b=viscousFaceGradient(m,id,v,gv,bv,fv,constantV);
        result[id]={-viscosity*(a.x*f.areaVector.x+b.x*f.areaVector.y),
                    -viscosity*(a.y*f.areaVector.x+b.y*f.areaVector.y)};
        if (!f.neighbour) {
            if (constantU[id] && fu[id])
                result[id].x+=viscosity*(gu[f.owner].x*f.correction.x+gu[f.owner].y*f.correction.y);
            if (constantV[id] && fv[id])
                result[id].y+=viscosity*(gv[f.owner].x*f.correction.x+gv[f.owner].y*f.correction.y);
        }
    }
    return result;
}

// Reconstruction policy, not a mesh/solution acceptance tolerance. With
// inverse-distance-normalized rows this condition measures directional sample
// coverage. A chain of wall triangles may need more than two graph rings.
inline constexpr double pressureGradientConditionTarget2D = 16.;
inline constexpr std::size_t pressureGradientMaximumRings2D = 6;

// Fixed-mesh reconstruction geometry. Samples preserve the original row order,
// divisions and final 2x2 solve; only geometry/ring discovery is cached. Build a
// fresh stencil when the mesh or fixed-boundary mask changes. A solver owns its
// pressure stencil locally, so separate solves and rejected candidates cannot
// share mutable state. Boundary VALUES and cell fields remain live inputs.
struct FlowGradientStencil2D {
    struct Sample { std::size_t index; Vector2D direction; double length; bool boundary; bool zero; };
    struct Row { std::size_t begin,end; double xx,xy,yy,det; };
    std::vector<Sample> samples;
    std::vector<Row> rows;
    std::size_t boundarySize=0;
    std::vector<Vector2D> apply(const std::vector<double>& field,const std::vector<double>& boundary) const {
        if(field.size()!=rows.size() || boundary.size()!=boundarySize)
            throw std::runtime_error("Cached flow gradient field size differs from stencil");
        std::vector<Vector2D> result(rows.size());
        for(std::size_t i=0;i<rows.size();++i) {
            const auto& row=rows[i];double bx=0,by=0;
            for(std::size_t k=row.begin;k<row.end;++k) {
                const auto& sample=samples[k];
                const double value=sample.zero?0:((sample.boundary?boundary[sample.index]:field[sample.index])-field[i])/sample.length;
                bx+=sample.direction.x*value;by+=sample.direction.y*value;
            }
            result[i]={(row.yy*bx-row.xy*by)/row.det,(row.xx*by-row.xy*bx)/row.det};
            if(!std::isfinite(result[i].x) || !std::isfinite(result[i].y))
                throw std::runtime_error("Flow gradient numerical range exceeded");
        }
        return result;
    }
};

// Pressure at velocity boundaries is extrapolated from interior values.
// It is not a prescribed zero physical pressure gradient. Velocity slip/outflow
// retains the zero-normal row; pressure-correction face flux remains a separate BC.
inline std::vector<Vector2D> flowGradient(const FvMesh2D& m,
                               const std::vector<double>& u,
                               const std::vector<double>& bc,
                               const std::vector<bool>& fixed,
                               bool extrapolateUnknown = false,
                               FlowGradientStencil2D* capture = nullptr) {
    if(capture) {
        capture->rows.clear();capture->rows.reserve(m.cells.size());
        capture->samples.clear();capture->samples.reserve(2*m.faces.size());
        capture->boundarySize=m.faces.size();
    }
    std::vector<Vector2D> g(u.size());
    for (std::size_t i = 0; i < u.size(); ++i) {
        double xx = 0;
        double xy = 0;
        double yy = 0;
        double bx = 0;
        double by = 0;
        bool omittedBoundary=false;
        const auto rowBegin=capture?capture->samples.size():0;
        for (auto id : m.cells[i].faces) {
            const auto& f = m.faces[id];
            const auto j =
                f.owner == i ? f.neighbour : std::optional<std::size_t>(f.owner);
            Vector2D d;
            double value = 0;
            double sampleLength=1;
            if (j || fixed[id]) {
                d = (j ? m.cells[*j].centre : f.centre) - m.cells[i].centre;
                const double length = std::hypot(d.x, d.y);
                if (!(length > 0)) throw std::runtime_error("Flow gradient degenerate stencil");
                sampleLength=length;
                value = ((j ? u[*j] : bc[id]) - u[i]) / length;
                d = d * (1 / length);
            } else {
                if (extrapolateUnknown) {omittedBoundary=true;continue;}
                d = f.areaVector;
                const double length = std::hypot(d.x, d.y);
                d = d * (1 / length);
            }
            if(capture)capture->samples.push_back({j?*j:id,d,sampleLength,!j,!j&&!fixed[id]});
            xx += d.x * d.x;
            xy += d.x * d.y;
            yy += d.y * d.y;
            bx += d.x * value;
            by += d.y * value;
        }
        const auto fullRank = [&] {
            return xx * yy - xy * xy > 64 * std::numeric_limits<double>::epsilon() *
                                               (xx + yy) * (xx + yy);
        };
        if (extrapolateUnknown && (omittedBoundary || !fullRank())) {
            // Keep the complete second ring used by existing boundary
            // reconstruction, then expand only poorly conditioned stencils.
            // Add whole sorted rings for deterministic rotation-independent
            // selection. Use all available samples if the graph ends sooner;
            // the original numerical rank check still applies below.
            std::vector<std::size_t> visited{i}, frontier{i};
            for (std::size_t ring=1;ring<=pressureGradientMaximumRings2D;++ring) {
                std::vector<std::size_t> next;
                for (const auto j:frontier) for (const auto id:m.cells[j].faces) {
                    const auto& f=m.faces[id];if(!f.neighbour)continue;
                    const auto k=f.owner==j?*f.neighbour:f.owner;
                    if(!std::binary_search(visited.begin(),visited.end(),k))next.push_back(k);
                }
                std::sort(next.begin(),next.end());next.erase(std::unique(next.begin(),next.end()),next.end());
                if(next.empty())break;
                if(ring>1)for(const auto k:next){
                    auto d=m.cells[k].centre-m.cells[i].centre;const double length=std::hypot(d.x,d.y);
                    if (!(length>0)) throw std::runtime_error("Flow gradient degenerate extended stencil");
                    const double value=(u[k]-u[i])/length;d=d*(1/length);
                    if(capture)capture->samples.push_back({k,d,length,false,false});
                    xx+=d.x*d.x;xy+=d.x*d.y;yy+=d.y*d.y;bx+=d.x*value;by+=d.y*value;
                }
                visited.insert(visited.end(),next.begin(),next.end());std::sort(visited.begin(),visited.end());frontier=std::move(next);
                const double maximum=.5*(xx+yy+std::hypot(xx-yy,2*xy));
                if(ring>=2 && (xx*yy-xy*xy)>=maximum*maximum/pressureGradientConditionTarget2D)break;
            }
        }
        const double det = xx * yy - xy * xy;
        if (!fullRank())
            throw std::runtime_error("Flow gradient rank deficient; cannot reconstruct from the available stencil");
        if(capture)capture->rows.push_back({rowBegin,capture->samples.size(),xx,xy,yy,det});
        g[i] = {(yy * bx - xy * by) / det, (xx * by - xy * bx) / det};
        if (!std::isfinite(g[i].x) || !std::isfinite(g[i].y))
            throw std::runtime_error("Flow gradient numerical range exceeded");
    }
    return g;
}


inline FlowGradientStencil2D buildFlowGradientStencil2D(const FvMesh2D& mesh,
    const std::vector<bool>& fixed,bool extrapolateUnknown=false) {
    if(fixed.size()!=mesh.faces.size())throw std::runtime_error("Flow gradient boundary mask size differs from mesh");
    FlowGradientStencil2D result;
    (void)flowGradient(mesh,std::vector<double>(mesh.cells.size()),std::vector<double>(mesh.faces.size()),fixed,extrapolateUnknown,&result);
    return result;
}

// Velocity boundary types can change with outlet backflow. Compare the entire
// mask before each use and rebuild on any change; values never enter the cache.
// The mesh referenced by this per-solve object must remain immutable.
class ChangingFlowGradientStencil2D {
    const FvMesh2D& mesh_;
    std::vector<bool> fixed_;
    std::optional<FlowGradientStencil2D> stencil_;
public:
    explicit ChangingFlowGradientStencil2D(const FvMesh2D& mesh):mesh_(mesh) {}
    std::vector<Vector2D> apply(const std::vector<double>& field,const std::vector<double>& boundary,
                               const std::vector<bool>& fixed) {
        if(!stencil_ || fixed_!=fixed) {
            auto next=buildFlowGradientStencil2D(mesh_,fixed);
            stencil_=std::move(next);fixed_=fixed;
        }
        return stencil_->apply(field,boundary);
    }
};

// Internal operators: the caller validates the mesh and field extents once.
// A single reconstructed value is stored per shared face, never per incidence.
inline std::vector<double> pressureFaceValues(
    const FvMesh2D& mesh, const std::vector<double>& p,
    const std::vector<Vector2D>& gradient, const std::vector<double>& boundary,
    const std::vector<bool>& fixed) {
    std::vector<double> result(mesh.faces.size());
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& f = mesh.faces[id];
        const auto i = f.owner;
        if (f.neighbour) {
            const auto j = *f.neighbour;
            const double w = f.neighbourWeight;
            const Point2D point{(1-w)*mesh.cells[i].centre.x+w*mesh.cells[j].centre.x,
                                (1-w)*mesh.cells[i].centre.y+w*mesh.cells[j].centre.y};
            const Vector2D g{(1-w)*gradient[i].x+w*gradient[j].x,
                             (1-w)*gradient[i].y+w*gradient[j].y};
            const auto d = f.centre-point;
            result[id] = (1-w)*p[i]+w*p[j]+g.x*d.x+g.y*d.y;
        } else if (fixed[id]) {
            result[id] = boundary[id];
        } else {
            const auto d = f.centre-mesh.cells[i].centre;
            result[id] = p[i]+gradient[i].x*d.x+gradient[i].y*d.y;
        }
    }
    return result;
}

inline std::vector<Vector2D> conservativePressureGradient(
    const FvMesh2D& mesh, const std::vector<double>& facePressure) {
    std::vector<Vector2D> result(mesh.cells.size());
    for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto& f = mesh.faces[id];
        const Vector2D force{facePressure[id]*f.areaVector.x,
                             facePressure[id]*f.areaVector.y};
        result[f.owner].x += force.x;
        result[f.owner].y += force.y;
        if (f.neighbour) {
            result[*f.neighbour].x -= force.x;
            result[*f.neighbour].y -= force.y;
        }
    }
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i].x /= mesh.cells[i].area;
        result[i].y /= mesh.cells[i].area;
    }
    return result;
}

// Barth-Jespersen-style cell limiter, applied at every actual face centre.
// It bounds reconstructed face values by the cell's neighbour/Dirichlet stencil;
// it does not assert a maximum principle for the coupled Navier-Stokes solution.
inline std::vector<double> faceReconstructionLimiter(
    const FvMesh2D& mesh, const std::vector<double>& field,
    const std::vector<Vector2D>& gradient, const std::vector<double>& boundary,
    const std::vector<bool>& fixed) {
    std::vector<double> limiter(mesh.cells.size(), 1.);
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        double lo = field[i], hi = field[i];
        for (const auto id : mesh.cells[i].faces) {
            const auto& f = mesh.faces[id];
            if (f.neighbour) {
                const auto j = f.owner == i ? *f.neighbour : f.owner;
                lo = std::min(lo, field[j]); hi = std::max(hi, field[j]);
            } else if (fixed[id]) {
                lo = std::min(lo, boundary[id]); hi = std::max(hi, boundary[id]);
            }
        }
        for (const auto id : mesh.cells[i].faces) {
            const auto d = mesh.faces[id].centre-mesh.cells[i].centre;
            const double delta = gradient[i].x*d.x+gradient[i].y*d.y;
            if (delta > 0) limiter[i] = std::min(limiter[i], (hi-field[i])/delta);
            else if (delta < 0) limiter[i] = std::min(limiter[i], (lo-field[i])/delta);
        }
        limiter[i] = std::clamp(limiter[i], 0., 1.);
    }
    return limiter;
}

inline double upwindFaceValue(
    const FvMesh2D& mesh, std::size_t id, double flux,
    const std::vector<double>& field, const std::vector<Vector2D>& gradient,
    const std::vector<double>& limiter) {
    const auto& f = mesh.faces[id];
    const auto up = (f.neighbour && flux < 0) ? *f.neighbour : f.owner;
    if (limiter.empty()) return field[up];
    const auto d = f.centre-mesh.cells[up].centre;
    return field[up]+limiter[up]*(gradient[up].x*d.x+gradient[up].y*d.y);
}

// Limit momentum in the normal/tangent frame of each target face. Limiting
// global Cartesian components separately is nonlinear and depends on the
// coordinate frame. Here both scalar bounds and gradient increments are
// projected into a frame carried by the mesh. The projected reconstructions
// remain inside the upwind cell's neighbour/Dirichlet stencil bounds.
// This is a local Barth-Jespersen construction, not a maximum-principle claim
// for the coupled solution or a reuse of a shallow-water solver's algorithm.
inline std::vector<Vector2D> faceFrameVelocityValues(
    const FvMesh2D& mesh, const std::vector<double>& flux,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<Vector2D>& gu, const std::vector<Vector2D>& gv,
    const std::vector<double>& bu, const std::vector<double>& bv,
    const std::vector<bool>& fixedU, const std::vector<bool>& fixedV) {
    std::vector<Vector2D> result(mesh.faces.size());
    for (std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& target=mesh.faces[id];
        const auto up=target.neighbour && flux[id]<0 ? *target.neighbour : target.owner;
        const double length=std::hypot(target.areaVector.x,target.areaVector.y);
        const Vector2D normal=target.areaVector*(1/length), tangent{-normal.y,normal.x};
        const auto displacement=target.centre-mesh.cells[up].centre;
        const auto increment=[&](Vector2D d) {return Vector2D{gu[up].x*d.x+gu[up].y*d.y,
                                                             gv[up].x*d.x+gv[up].y*d.y};};
        const auto project=[](Vector2D a,Vector2D b) {return a.x*b.x+a.y*b.y;};
        const auto limited=[&](Vector2D axis) {
            double lo=0,hi=0;
            for (const auto fid:mesh.cells[up].faces) {
                const auto& f=mesh.faces[fid];
                Vector2D delta{};
                if (f.neighbour) {
                    const auto other=f.owner==up ? *f.neighbour : f.owner;
                    delta={u[other]-u[up],v[other]-v[up]};
                } else if (fixedU[fid] || fixedV[fid]) {
                    delta={fixedU[fid]?bu[fid]-u[up]:0, fixedV[fid]?bv[fid]-v[up]:0};
                } else continue;
                const double sample=project(delta,axis);
                lo=std::min(lo,sample);hi=std::max(hi,sample);
            }
            double phi=1;
            for (const auto fid:mesh.cells[up].faces) {
                const double change=project(increment(mesh.faces[fid].centre-mesh.cells[up].centre),axis);
                if (change>0) phi=std::min(phi,hi/change);
                else if (change<0) phi=std::min(phi,lo/change);
            }
            return std::clamp(phi,0.,1.)*project(increment(displacement),axis);
        };
        const double dn=limited(normal),dt=limited(tangent);
        result[id]={u[up]+normal.x*dn+tangent.x*dt,v[up]+normal.y*dn+tangent.y*dt};
    }
    return result;
}

} // namespace cartmesh2d::fv::detail
