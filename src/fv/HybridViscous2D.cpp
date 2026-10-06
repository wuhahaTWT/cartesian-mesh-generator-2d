#include "cartmesh2d/fv/HybridViscous2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(std::string("hybrid viscous: ")+message);}
double checked(double x){require(std::isfinite(x),"nonfinite arithmetic");return x;}
bool finite(Vector2D v){return std::isfinite(v.x)&&std::isfinite(v.y);}
double component(Vector2D v,std::size_t k){return k?v.y:v.x;}
}
struct HybridViscousOperator2D::Data {
    struct Face {
        std::size_t owner=0,start=0,count=0;
        std::optional<std::size_t> neighbour,partner;
        ViscousBoundaryKind2D kind=ViscousBoundaryKind2D::ZeroTraction;
        Vector2D value{},normal{};
        std::array<Vector2D,2> basis{{{1,0},{0,1}}};
    };
    std::vector<Face> faces;
    std::vector<HybridViscousCell2D> cells;
    std::vector<double> areas,cellBound,diagonal,off;
    std::unique_ptr<detail::SparsePattern2D> pattern;
    double viscosity=0;
    std::size_t unknowns=0;
};

HybridViscousOperator2D::HybridViscousOperator2D(const FvMesh2D& mesh,
    const std::vector<ViscousBoundary2D>& boundaries,double viscosity) {
    validateFvMesh2D(mesh);require(std::isfinite(viscosity)&&viscosity>0,"viscosity must be positive");
    auto data=std::make_shared<Data>();auto& d=*data;d.viscosity=viscosity;
    std::vector<const ViscousBoundary2D*> lookup(mesh.faces.size());
    for(const auto& b:boundaries) {
        require(b.face<lookup.size()&&!mesh.faces[b.face].neighbour&&!lookup[b.face],"invalid or duplicate boundary");
        require(b.kind==ViscousBoundaryKind2D::Velocity||b.kind==ViscousBoundaryKind2D::Slip||b.kind==ViscousBoundaryKind2D::ZeroTraction||b.kind==ViscousBoundaryKind2D::Periodic,"unknown boundary kind");
        require(finite(b.velocity),"nonfinite boundary velocity");
        require(b.kind==ViscousBoundaryKind2D::Velocity||(b.velocity.x==0&&b.velocity.y==0),"inactive boundary velocity must be zero");
        require((b.kind==ViscousBoundaryKind2D::Periodic)==b.partner.has_value(),"periodic pairing mismatch");lookup[b.face]=&b;
    }
    d.faces.resize(mesh.faces.size());
    std::vector<std::size_t> parent(mesh.cells.size());std::iota(parent.begin(),parent.end(),0);
    const auto root=[&](std::size_t c){while(parent[c]!=c){parent[c]=parent[parent[c]];c=parent[c];}return c;};
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id];auto& face=d.faces[id];const auto* bc=lookup[id];
        require(f.neighbour||bc,"missing boundary");face.owner=f.owner;face.neighbour=f.neighbour;
        const double length=std::hypot(f.areaVector.x,f.areaVector.y);
        face.normal={f.areaVector.x/length,f.areaVector.y/length};
        if(bc){face.kind=bc->kind;face.value=bc->velocity;face.partner=bc->partner;}
        if(face.partner) {
            const auto p=*face.partner;
            require(p<lookup.size()&&p!=id&&lookup[p]&&lookup[p]->partner==id,"nonreciprocal periodic pair");
            const auto& other=mesh.faces[p];
            require(std::hypot(f.areaVector.x+other.areaVector.x,f.areaVector.y+other.areaVector.y)<=1e-10*length,"periodic normals differ");
            parent[root(f.owner)]=root(other.owner);
        }else if(f.neighbour)parent[root(f.owner)]=root(*f.neighbour);
        if(face.partner&&id>*face.partner){face.start=d.faces[*face.partner].start;face.count=2;continue;}
        face.start=d.unknowns;
        if(f.neighbour||face.partner||face.kind==ViscousBoundaryKind2D::ZeroTraction)face.count=2;
        else if(face.kind==ViscousBoundaryKind2D::Slip){face.count=1;face.basis[0]={-face.normal.y,face.normal.x};}
        d.unknowns+=face.count;
    }
    std::vector<bool> anchored(mesh.cells.size());
    for(const auto& face:d.faces)if(!face.neighbour&&!face.partner&&face.kind==ViscousBoundaryKind2D::Velocity)anchored[root(face.owner)]=true;
    for(std::size_t c=0;c<mesh.cells.size();++c)require(anchored[root(c)],"each connected component needs a prescribed-velocity face to remove rigid trace modes");

    std::vector<std::pair<std::size_t,std::size_t>> connections;
    for(std::size_t c=0;c<mesh.cells.size();++c) {
        d.cells.emplace_back(mesh,c,viscosity);d.areas.push_back(mesh.cells[c].area);
        const auto& faces=d.cells.back().faces();const auto& a=d.cells.back().tractionMatrix();const auto n=2*faces.size();
        std::array<double,4> cellBlock{};
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)cellBlock[2*(i%2)+j%2]+=a[i*n+j];
        // Static condensation minimizes the nonnegative local energy, hence
        // 0 <= K_cell <= blockdiag(E' A E). This block row bound controls the
        // mass-scaled momentum spectrum; it is NOT a pointwise positivity bound.
        d.cellBound.push_back(checked(std::max(std::abs(cellBlock[0])+std::abs(cellBlock[1]),std::abs(cellBlock[2])+std::abs(cellBlock[3]))));
        for(auto f:faces)for(auto g:faces)for(std::size_t i=0;i<d.faces[f].count;++i)for(std::size_t j=0;j<d.faces[g].count;++j) {
            const auto row=d.faces[f].start+i,column=d.faces[g].start+j;
            if(row<column)connections.emplace_back(row,column);
        }
    }
    if(d.unknowns) {
        d.pattern=std::make_unique<detail::SparsePattern2D>(d.unknowns,connections);
        detail::SparseSystem2D matrix(*d.pattern);
        for(const auto& cell:d.cells) {
            const auto& faces=cell.faces();const auto& a=cell.tractionMatrix();const auto n=2*faces.size();
            for(std::size_t i=0;i<faces.size();++i)for(std::size_t j=0;j<faces.size();++j) {
                const auto& f=d.faces[faces[i]];const auto& g=d.faces[faces[j]];
                for(std::size_t p=0;p<f.count;++p)for(std::size_t q=0;q<g.count;++q) {
                    const auto row=f.start+p,column=g.start+q;if(row<column)continue;
                    double value=0;
                    for(std::size_t k=0;k<2;++k)for(std::size_t l=0;l<2;++l)value+=component(f.basis[p],k)*a[(2*i+k)*n+2*j+l]*component(g.basis[q],l)/viscosity;
                    if(row==column)matrix.diag[row]+=value;
                    else {matrix.add(row,column,value);matrix.add(column,row,value);}
                }
            }
        }
        for(double value:matrix.diag)require(std::isfinite(value)&&value>0,"nonpositive trace diagonal");
        d.diagonal=std::move(matrix.diag);d.off=std::move(matrix.off);
    }
    data_=std::move(data);
}

std::size_t HybridViscousOperator2D::traceUnknowns() const{return data_->unknowns;}
void HybridViscousOperator2D::visitTraceMatrix(const std::function<void(std::size_t,std::size_t,double)>& visit) const {
    const auto& d=*data_;
    for(std::size_t i=0;i<d.unknowns;++i){visit(i,i,d.viscosity*d.diagonal[i]);for(auto k=d.pattern->rows[i];k<d.pattern->rows[i+1];++k)visit(i,d.pattern->columns[k],d.viscosity*d.off[k]);}
}
HybridViscousResult2D HybridViscousOperator2D::evaluate(const std::vector<Vector2D>& u,
    const std::vector<double>& rho,double tolerance) const {
    const auto& d=*data_;require(u.size()==d.cells.size()&&rho.size()==u.size(),"velocity or density dimensions differ");
    require(std::isfinite(tolerance)&&tolerance>0&&tolerance<=1e-2,"invalid relative tolerance");
    for(std::size_t c=0;c<u.size();++c)require(finite(u[c])&&std::isfinite(rho[c])&&rho[c]>0,"invalid velocity or density");
    HybridViscousResult2D out;out.trace.resize(d.faces.size());out.faceFlux.resize(d.faces.size());
    out.cellResidual.resize(u.size());out.cellDissipation.resize(u.size());out.cellMechanicalHeating.resize(u.size());out.rate.resize(u.size());out.traceUnknowns=d.unknowns;
    const auto reference=u.front();
    for(std::size_t id=0;id<d.faces.size();++id) {
        const auto& f=d.faces[id];
        if(f.count==2)out.trace[id]=reference;
        else if(f.count==1){const auto t=f.basis[0];const double tangential=dot(reference,t);out.trace[id]={t.x*tangential,t.y*tangential};}
        else out.trace[id]=f.value;
    }
    double velocityScale=0;
    for(std::size_t c=0;c<u.size();++c)for(auto id:d.cells[c].faces())velocityScale=std::max({velocityScale,std::abs(out.trace[id].x-u[c].x),std::abs(out.trace[id].y-u[c].y)});
    if(d.unknowns&&velocityScale>0) {
        detail::SparseSystem2D system(*d.pattern);system.diag=d.diagonal;system.off=d.off;
        for(std::size_t c=0;c<u.size();++c) {
            const auto& faces=d.cells[c].faces();const auto& a=d.cells[c].tractionMatrix();const auto n=2*faces.size();
            for(std::size_t i=0;i<faces.size();++i) {
                Vector2D traction{};
                for(std::size_t j=0;j<faces.size();++j)for(std::size_t k=0;k<2;++k) {
                    const double delta=(component(out.trace[faces[j]],k)-component(u[c],k))/velocityScale;
                    traction.x+=a[(2*i)*n+2*j+k]/d.viscosity*delta;
                    traction.y+=a[(2*i+1)*n+2*j+k]/d.viscosity*delta;
                }
                const auto& f=d.faces[faces[i]];for(std::size_t k=0;k<f.count;++k)system.rhs[f.start+k]-=dot(f.basis[k],traction);
            }
        }
        detail::LinearVector2D x(d.unknowns);detail::LinearWorkspace2D work(d.unknowns);
        out.iterations=system.solvePressure(x,work,detail::LinearPressureMethod2D::Jacobi,tolerance);
        detail::LinearVector2D residual(d.unknowns);
        for(std::size_t i=0;i<d.unknowns;++i)residual[i]=system.compensatedResidualRow(i,x);
        const double physicalScale=checked(d.viscosity*velocityScale);
        out.traceResidualNorm=checked(physicalScale*detail::linearNorm(residual));
        out.traceResidualTarget=checked(physicalScale*(1e-13+tolerance*detail::linearNorm(system.rhs)));
        for(std::size_t id=0;id<d.faces.size();++id) {
            const auto& f=d.faces[id];
            for(std::size_t k=0;k<f.count;++k){const double value=velocityScale*x[f.start+k];out.trace[id].x+=f.basis[k].x*value;out.trace[id].y+=f.basis[k].y*value;}
            require(finite(out.trace[id]),"nonfinite solved trace");
        }
    }
    std::vector<std::array<double,3>> other(d.faces.size());
    for(std::size_t c=0;c<u.size();++c) {
        const auto& cell=d.cells[c];std::vector<Vector2D> trace;for(auto f:cell.faces())trace.push_back(out.trace[f]);
        const auto local=cell.evaluate(u[c],trace);out.cellDissipation[c]=local.dissipation;out.dissipation+=local.dissipation;
        out.rate[c]=checked(d.cellBound[c]/(rho[c]*d.areas[c]));
        for(std::size_t i=0;i<cell.faces().size();++i) {
            const auto id=cell.faces()[i];if(d.faces[id].owner==c)out.faceFlux[id]=local.outwardFlux[i];else other[id]=local.outwardFlux[i];
        }
    }
    for(std::size_t id=0;id<d.faces.size();++id) {
        const auto& f=d.faces[id];auto& flux=out.faceFlux[id];
        if(f.partner&&id>*f.partner)continue;
        if(f.neighbour||f.partner) {
            const auto q=f.partner?out.faceFlux[*f.partner]:other[id];
            out.maximumTractionJump=std::max(out.maximumTractionJump,std::hypot(flux[0]+q[0],flux[1]+q[1]));
            if(f.partner)for(std::size_t k=0;k<3;++k)out.faceFlux[*f.partner][k]=-flux[k];
        }else if(f.kind==ViscousBoundaryKind2D::ZeroTraction) {
            out.maximumBoundaryConstraintResidual=std::max(out.maximumBoundaryConstraintResidual,std::hypot(flux[0],flux[1]));flux={};
        }else if(f.kind==ViscousBoundaryKind2D::Slip) {
            const auto t=f.basis[0];out.maximumBoundaryConstraintResidual=std::max(out.maximumBoundaryConstraintResidual,std::abs(t.x*flux[0]+t.y*flux[1]));
            const double normal=f.normal.x*flux[0]+f.normal.y*flux[1];flux={normal*f.normal.x,normal*f.normal.y,0};
        }
    }
    for(std::size_t id=0;id<d.faces.size();++id) {
        const auto& f=d.faces[id];for(std::size_t k=0;k<3;++k){out.cellResidual[f.owner][k]+=out.faceFlux[id][k];if(f.neighbour)out.cellResidual[*f.neighbour][k]-=out.faceFlux[id][k];}
    }
    for(std::size_t c=0;c<u.size();++c) {
        const auto& r=out.cellResidual[c];const double heat=checked(-r[2]+u[c].x*r[0]+u[c].y*r[1]);
        out.cellMechanicalHeating[c]=heat;out.mechanicalHeating+=heat;
        out.maximumCellWorkDefect=std::max(out.maximumCellWorkDefect,std::abs(heat-out.cellDissipation[c]));
    }
    checked(out.dissipation);checked(out.mechanicalHeating);return out;
}
} // namespace cartmesh2d::fv
