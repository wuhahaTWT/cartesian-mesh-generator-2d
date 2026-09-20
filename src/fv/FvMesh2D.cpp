#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>

namespace cartmesh2d::fv {
namespace {
void require(bool test, std::string_view message) {
    // Successful per-face/cell checks must not allocate diagnostic strings.
    // A temporary string argument remains alive throughout this call.
    if (!test) throw std::runtime_error("FVM mesh: " + std::string(message));
}
bool finite(Point2D p) { return std::isfinite(p.x) && std::isfinite(p.y); }
}
FvMesh2D makeFvMesh2D(const TopologyMesh2D& topology, const TolerancePolicy& tol) {
    require(topology.valid() && !topology.cells.empty(), "invalid/empty topology");
    require(std::isfinite(tol.absolute) && std::isfinite(tol.relative) &&
            tol.absolute >= 0 && tol.relative >= 0, "invalid geometry tolerance");
    FvMesh2D mesh;
    mesh.cells.resize(topology.cells.size()); mesh.faces.resize(topology.edges.size());
    std::vector<int> ownerSeen(mesh.faces.size()), neighbourSeen(mesh.faces.size());
    std::map<std::pair<double,double>,std::size_t> coordinates;
    for (std::size_t i=0; i<topology.vertices.size(); ++i) {
        const auto& v=topology.vertices[i];
        require(v.id==i && finite(v.point), "non-contiguous vertex ID or nonfinite coordinate");
        require(coordinates.emplace(std::make_pair(v.point.x,v.point.y),i).second, "duplicate vertex coordinate");
    }
    std::set<std::pair<std::size_t,std::size_t>> endpointPairs;
    for (std::size_t i=0; i<topology.edges.size(); ++i) {
        const auto& e=topology.edges[i];
        require(e.id==i && e.v0<topology.vertices.size() && e.v1<topology.vertices.size() && e.v0!=e.v1,
                "invalid edge IDs/endpoints at face " + std::to_string(i));
        require(e.owner<mesh.cells.size() && (!e.neighbour || (*e.neighbour<mesh.cells.size() && *e.neighbour!=e.owner)), "invalid owner/neighbour");
        require(endpointPairs.emplace(std::minmax(e.v0,e.v1)).second, "duplicate edge");
        require(e.neighbour ? e.patch==BoundaryPatch2D::None :
                (e.patch==BoundaryPatch2D::DomainBoundary || e.patch==BoundaryPatch2D::EmbeddedBoundary), "unclassified/inconsistent patch");
        auto& f=mesh.faces[i]; f.owner=e.owner; f.neighbour=e.neighbour; f.patch=e.patch;
        const auto a=topology.vertices[e.v0].point, b=topology.vertices[e.v1].point;
        f.centre={a.x+(b.x-a.x)*0.5,a.y+(b.y-a.y)*0.5};
    }
    for (std::size_t i=0; i<topology.cells.size(); ++i) {
        const auto& cell=topology.cells[i]; auto& out=mesh.cells[i];
        require(cell.id==i && cell.vertices.size()>=3 && cell.edges.size()==cell.vertices.size(), "invalid cell loop");
        Polygon2D poly;
        for (auto v:cell.vertices) {
            require(v<topology.vertices.size(), "invalid cell vertex"); poly.vertices.push_back(topology.vertices[v].point);
        }
        require(BoundaryLoop(poly.vertices).diagnose(tol).valid(), "invalid polygon at cell " + std::to_string(i));
        const auto metric=evaluateSolverCellMetrics2D(poly,tol);
        require(metric.valid && poly.signedArea()>0 && metric.area>0 && finite(metric.centroid), "nonpositive or invalid cell geometry");
        const auto box=poly.bounds(); const double extent=std::max(box.max.x-box.min.x,box.max.y-box.min.y);
        require(std::isfinite(cell.geometryArea) && std::abs(cell.geometryArea-metric.area)<=tol.areaScale(extent), "stored area does not match final polygon");
        out.area=metric.area; out.centre=metric.centroid; out.faces=cell.edges;
        std::set<std::size_t> seen;
        Vector2D closure{}; double perimeter=0;
        for (std::size_t j=0; j<cell.edges.size(); ++j) {
            const auto id=cell.edges[j];
            require(id<mesh.faces.size() && seen.insert(id).second, "invalid/repeated cell edge");
            const auto& e=topology.edges[id]; const auto a=cell.vertices[j], b=cell.vertices[(j+1)%cell.vertices.size()];
            require((e.v0==a && e.v1==b)||(e.v0==b && e.v1==a), "edge does not match polygon loop");
            const auto d=topology.vertices[b].point-topology.vertices[a].point;
            const Vector2D s{d.y,-d.x};
            closure.x+=s.x; closure.y+=s.y; perimeter+=std::hypot(s.x,s.y);
            if (e.owner==i) { ++ownerSeen[id]; mesh.faces[id].areaVector=s; }
            else { require(e.neighbour && *e.neighbour==i,"cell references unrelated edge"); ++neighbourSeen[id]; }
        }
        const double closing=std::hypot(closure.x,closure.y);
        require(closing<=tol.scale(perimeter), "cell outward normals do not close");
        mesh.maxClosureError=std::max(mesh.maxClosureError,closing/perimeter);
    }
    for (std::size_t i=0; i<mesh.faces.size(); ++i) {
        auto& f=mesh.faces[i];
        require(ownerSeen[i]==1 && neighbourSeen[i]==(f.neighbour?1:0),"orphan/nonmanifold face incidence");
        if (f.neighbour) {
            const auto& c=topology.cells[*f.neighbour];
            const auto pos=static_cast<std::size_t>(std::find(c.edges.begin(),c.edges.end(),i)-c.edges.begin());
            const auto d=topology.vertices[c.vertices[(pos+1)%c.vertices.size()]].point-topology.vertices[c.vertices[pos]].point;
            require(std::hypot(f.areaVector.x+d.y,f.areaVector.y-d.x)<=tol.scale(std::hypot(d.x,d.y)), "shared face orientations are not opposite");
        }
        const auto o=mesh.cells[f.owner].centre;
        const auto d=(f.neighbour?mesh.cells[*f.neighbour].centre:f.centre)-o;
        const double sd=dot(f.areaVector,d), s2=squaredNorm(f.areaVector);
        require(std::isfinite(sd) && sd>0 && s2>0, "nonpositive normal centre distance at face " + std::to_string(i));
        f.transmissibility=s2/sd;
        f.correction={f.areaVector.x-f.transmissibility*d.x,f.areaVector.y-f.transmissibility*d.y};
        if (f.neighbour) {
            f.neighbourWeight=dot(f.areaVector,f.centre-o)/sd;
            require(f.neighbourWeight>0 && f.neighbourWeight<1,"face outside centre normal bracket");
        }
        require(std::isfinite(f.transmissibility) && std::isfinite(f.correction.x) && std::isfinite(f.correction.y),"nonfinite face coefficient");
    }
    const auto quality=evaluateSolverQuality2D(topology,{},tol);
    require(quality.valid(), "existing Solver quality gate failed ("+std::to_string(quality.issues.size())+" issues)");
    return mesh;
}

void validateFvMesh2D(const FvMesh2D& mesh) {
    require(!mesh.cells.empty() && !mesh.faces.empty(), "empty finite-volume mesh");
    std::vector<int> owners(mesh.faces.size()), neighbours(mesh.faces.size());
    for (const auto& f:mesh.faces) {
        require(f.owner<mesh.cells.size() && (!f.neighbour || (*f.neighbour<mesh.cells.size() && *f.neighbour!=f.owner)), "invalid cached owner/neighbour");
        require(finite(f.centre) && std::isfinite(f.areaVector.x) && std::isfinite(f.areaVector.y) &&
                std::hypot(f.areaVector.x,f.areaVector.y)>0 && std::isfinite(f.correction.x) && std::isfinite(f.correction.y) &&
                std::isfinite(f.transmissibility) && f.transmissibility>0, "invalid cached face geometry");
        require(f.neighbour ? (f.patch==BoundaryPatch2D::None && f.neighbourWeight>0 && f.neighbourWeight<1) :
                (f.patch==BoundaryPatch2D::DomainBoundary || f.patch==BoundaryPatch2D::EmbeddedBoundary), "invalid cached patch/weight");
    }
    for (std::size_t i=0;i<mesh.cells.size();++i) {
        const auto& c=mesh.cells[i];
        require(finite(c.centre) && std::isfinite(c.area) && c.area>0 && c.faces.size()>=3,"invalid cached cell geometry");
        for(auto id:c.faces) {
            require(id<mesh.faces.size(),"invalid cached face reference");
            const auto& f=mesh.faces[id];
            if(f.owner==i) ++owners[id];
            else {require(f.neighbour && *f.neighbour==i,"unrelated cached face");++neighbours[id];}
        }
    }
    for(std::size_t i=0;i<mesh.faces.size();++i)
        require(owners[i]==1 && neighbours[i]==(mesh.faces[i].neighbour?1:0),"invalid cached face incidence");
}

std::vector<Vector2D> reconstructGradient(const FvMesh2D& mesh,
        const std::vector<double>& u,const std::vector<double>& bc) {
    validateFvMesh2D(mesh);
    if (u.size()!=mesh.cells.size() || bc.size()!=mesh.faces.size()) throw std::invalid_argument("FVM gradient: wrong field size");
    std::vector<Vector2D> result(u.size());
    for (std::size_t i=0;i<u.size();++i) {
        double xx=0,xy=0,yy=0,bx=0,by=0;
        if (!std::isfinite(u[i])) throw std::runtime_error("FVM gradient: nonfinite field");
        for (auto id:mesh.cells[i].faces) {
            const auto& f=mesh.faces[id];
            const auto other=f.owner==i ? f.neighbour : std::optional<std::size_t>(f.owner);
            const auto d=(other?mesh.cells[*other].centre:f.centre)-mesh.cells[i].centre;
            const double delta=(other?u[*other]:bc[id])-u[i], r2=squaredNorm(d);
            if (!(r2>0) || !std::isfinite(delta)) throw std::runtime_error("FVM gradient: invalid sample");
            xx+=d.x*d.x/r2; xy+=d.x*d.y/r2; yy+=d.y*d.y/r2;
            bx+=d.x*delta/r2; by+=d.y*delta/r2;
        }
        const double det=xx*yy-xy*xy;
        if (!(det>64*std::numeric_limits<double>::epsilon()*(xx+yy)*(xx+yy)))
            throw std::runtime_error("FVM gradient: rank-deficient stencil at cell "+std::to_string(i));
        result[i]={(yy*bx-xy*by)/det,(xx*by-xy*bx)/det};
        if (!std::isfinite(result[i].x)||!std::isfinite(result[i].y)) throw std::runtime_error("FVM gradient: nonfinite reconstruction");
    }
    return result;
}
}
