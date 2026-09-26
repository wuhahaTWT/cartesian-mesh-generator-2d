#include "cartmesh2d/fv/ViscousStress2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(std::string("viscous stress: ")+message);}
double checked(double x){require(std::isfinite(x),"non-finite arithmetic");return x;}
// Sparse linear forms are used only during preparation of the complete 2x2
// momentum Jacobian, including cross derivatives and nonorthogonal corrections.
using Form=std::map<std::size_t,double>;
void add(Form& a,const Form& b,double s){for(const auto& [j,v]:b)a[j]+=s*v;}
using Gradient=std::array<Vector2D,2>; // rows grad(u), grad(v)
}
ViscousStressOperator2D::ViscousStressOperator2D(const FvMesh2D& mesh,
    const std::vector<ViscousBoundary2D>& boundaries,double viscosity,WallGradient2D wallGradient):viscosity_(viscosity) {
    validateFvMesh2D(mesh);require(std::isfinite(viscosity)&&viscosity>0,"dynamic viscosity must be positive");
    require(wallGradient==WallGradient2D::Linear||wallGradient==WallGradient2D::Quadratic,"invalid wall gradient scheme");
    std::vector<const ViscousBoundary2D*> lookup(mesh.faces.size());
    for(const auto& b:boundaries) {
        require(b.face<lookup.size()&&!mesh.faces[b.face].neighbour&&!lookup[b.face],"invalid/duplicate boundary");
        require(b.kind==ViscousBoundaryKind2D::Velocity||b.kind==ViscousBoundaryKind2D::Slip||b.kind==ViscousBoundaryKind2D::ZeroTraction||b.kind==ViscousBoundaryKind2D::Periodic,"unknown boundary kind");
        require(std::isfinite(b.velocity.x)&&std::isfinite(b.velocity.y),"non-finite boundary velocity");
        require(b.kind==ViscousBoundaryKind2D::Velocity||(b.velocity.x==0&&b.velocity.y==0),"inactive velocity must be zero");
        require((b.kind==ViscousBoundaryKind2D::Periodic)==b.partner.has_value(),"periodic pairing mismatch");
        lookup[b.face]=&b;
    }
    faces_.resize(mesh.faces.size());gradients_.resize(mesh.cells.size());
    for(const auto& c:mesh.cells)areas_.push_back(c.area);
    for(std::size_t id=0;id<faces_.size();++id) {
        const auto& f=mesh.faces[id];auto& a=faces_[id];const auto* b=lookup[id];
        require(f.neighbour||b,"missing viscous boundary");
        a.owner=f.owner;a.neighbour=f.neighbour;a.area=f.areaVector;
        const double length=std::hypot(a.area.x,a.area.y);a.normal={a.area.x/length,a.area.y/length};
        a.ownerOffset={f.centre.x-mesh.cells[f.owner].centre.x,f.centre.y-mesh.cells[f.owner].centre.y};
        a.d=a.ownerOffset;
        if(b){a.kind=b->kind;a.value=b->velocity;a.partner=b->partner;}
        if(a.partner) {
            const auto p=*a.partner;
            require(p<lookup.size()&&p!=id&&lookup[p]&&lookup[p]->partner==id,"nonreciprocal periodic pair");
            const auto& other=mesh.faces[p];
            require(std::hypot(a.area.x+other.areaVector.x,a.area.y+other.areaVector.y)<=1e-10*length,"periodic normals differ");
            a.neighbour=other.owner;
            a.neighbourOffset={other.centre.x-mesh.cells[other.owner].centre.x,other.centre.y-mesh.cells[other.owner].centre.y};
        }else if(a.neighbour)a.neighbourOffset={f.centre.x-mesh.cells[*a.neighbour].centre.x,f.centre.y-mesh.cells[*a.neighbour].centre.y};
        if(a.neighbour)a.d={a.ownerOffset.x-a.neighbourOffset.x,a.ownerOffset.y-a.neighbourOffset.y};
        a.normalDistance=checked(dot(a.normal,a.d));require(a.normalDistance>0,"nonpositive normal distance");
        a.weight=a.neighbour?dot(a.normal,a.ownerOffset)/a.normalDistance:0;
        require(!a.neighbour||(a.weight>0&&a.weight<1),"face outside centre bracket");
    }
    wallGradients_.resize(faces_.size());
    if(wallGradient==WallGradient2D::Quadratic) {
        std::vector<bool> prescribed(faces_.size());std::vector<std::optional<std::size_t>> partners(faces_.size());
        for(std::size_t id=0;id<faces_.size();++id){prescribed[id]=!faces_[id].neighbour&&faces_[id].kind==ViscousBoundaryKind2D::Velocity;partners[id]=faces_[id].partner;}
        for(std::size_t id=0;id<faces_.size();++id)if(prescribed[id]){wallGradients_[id]=quadraticWallGradient2D(mesh,id,prescribed,partners);++quadraticWalls_;}
    }
    std::vector<std::array<Form,4>> jac(mesh.cells.size());
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        double xx=0,xy=0,yy=0;const auto c=mesh.cells[i].centre;
        for(auto id:mesh.cells[i].faces) {
            const auto& f=mesh.faces[id];const auto& a=faces_[id];Sample s;s.value=a.value;s.kind=a.kind;s.normal=a.normal;
            Vector2D d{f.centre.x-c.x,f.centre.y-c.y};
            if(f.neighbour){s.cell=i==f.owner?*f.neighbour:f.owner;d={mesh.cells[*s.cell].centre.x-c.x,mesh.cells[*s.cell].centre.y-c.y};}
            else if(a.partner){s.cell=a.neighbour;d=a.d;}
            else if(a.kind!=ViscousBoundaryKind2D::Velocity){const double dn=dot(d,a.normal);d={dn*a.normal.x,dn*a.normal.y};}
            const double norm2=checked(dot(d,d));require(norm2>0,"zero stencil distance");
            xx+=d.x*d.x/norm2;xy+=d.x*d.y/norm2;yy+=d.y*d.y/norm2;
            s.weight={d.x/norm2,d.y/norm2};gradients_[i].push_back(s);
        }
        const double det=checked(xx*yy-xy*xy),trace=xx+yy;
        require(det>64*std::numeric_limits<double>::epsilon()*trace*trace,"rank-deficient velocity gradient stencil");
        for(auto& s:gradients_[i]) {
            s.weight={(yy*s.weight.x-xy*s.weight.y)/det,(xx*s.weight.y-xy*s.weight.x)/det};
            const double n[2]={s.normal.x,s.normal.y},w[2]={s.weight.x,s.weight.y};
            for(std::size_t k=0;k<2;++k)for(std::size_t a=0;a<2;++a) {
                auto& form=jac[i][2*k+a];
                if(s.cell){form[2*(*s.cell)+k]+=w[a];form[2*i+k]-=w[a];}
                else if(s.kind==ViscousBoundaryKind2D::Velocity)form[2*i+k]-=w[a];
                else if(s.kind==ViscousBoundaryKind2D::Slip)for(std::size_t l=0;l<2;++l)form[2*i+l]-=w[a]*n[k]*n[l];
            }
        }
    }
    std::vector<Form> rows(2*mesh.cells.size());
    for(std::size_t id=0;id<faces_.size();++id) {
        const auto& f=faces_[id];if((f.partner&&id>*f.partner)||(!f.neighbour&&f.kind==ViscousBoundaryKind2D::ZeroTraction))continue;
        std::array<Form,4> g;const double n[2]={f.normal.x,f.normal.y};
        for(std::size_t k=0;k<2;++k) {
            for(std::size_t a=0;a<2;++a){add(g[2*k+a],jac[f.owner][2*k+a],1-f.weight);if(f.neighbour)add(g[2*k+a],jac[*f.neighbour][2*k+a],f.weight);}
            if(f.neighbour||f.kind==ViscousBoundaryKind2D::Velocity) {
                Form delta{{2*f.owner+k,-1}};if(f.neighbour)delta[2*(*f.neighbour)+k]+=1;
                add(delta,g[2*k],-f.d.x);add(delta,g[2*k+1],-f.d.y);
                for(std::size_t a=0;a<2;++a)add(g[2*k+a],delta,n[a]/f.normalDistance);
            }
        }
        if(wallGradients_[id]) {
            g={};
            for(const auto& sample:wallGradients_[id]->samples)if(!sample.boundary)for(std::size_t k=0;k<2;++k){g[2*k][2*sample.index+k]+=sample.weight.x;g[2*k+1][2*sample.index+k]+=sample.weight.y;}
        }
        Form tx,ty;add(tx,g[0],4./3*f.area.x);add(tx,g[3],-2./3*f.area.x);add(tx,g[1],f.area.y);add(tx,g[2],f.area.y);
        add(ty,g[1],f.area.x);add(ty,g[2],f.area.x);add(ty,g[3],4./3*f.area.y);add(ty,g[0],-2./3*f.area.y);
        if(!f.neighbour&&f.kind==ViscousBoundaryKind2D::Slip){Form normal;add(normal,tx,n[0]);add(normal,ty,n[1]);tx.clear();ty.clear();add(tx,normal,n[0]);add(ty,normal,n[1]);}
        for(std::size_t k=0;k<2;++k){const auto& t=k==0?tx:ty;add(rows[2*f.owner+k],t,-viscosity_);if(f.neighbour)add(rows[2*(*f.neighbour)+k],t,viscosity_);}
    }
    rowNorm_.resize(areas_.size());
    for(std::size_t i=0;i<areas_.size();++i)for(std::size_t k=0;k<2;++k){double sum=0;for(const auto& [j,v]:rows[2*i+k])sum+=std::abs(v);rowNorm_[i]=std::max(rowNorm_[i],checked(.5*sum));}
}
ViscousStressResult2D ViscousStressOperator2D::evaluate(const std::vector<Vector2D>& u,const std::vector<double>& rho) const {
    require(u.size()==areas_.size()&&rho.size()==u.size(),"velocity/density size mismatch");
    for(std::size_t i=0;i<u.size();++i)require(std::isfinite(u[i].x)&&std::isfinite(u[i].y)&&std::isfinite(rho[i])&&rho[i]>0,"invalid velocity or density");
    ViscousStressResult2D out;out.faceFlux.resize(faces_.size());out.cellResidual.resize(u.size());out.rate.resize(u.size());
    std::vector<Gradient> gradients(u.size());
    for(std::size_t i=0;i<u.size();++i) {
        for(const auto& s:gradients_[i]) {
            Vector2D delta{};
            if(s.cell)delta={u[*s.cell].x-u[i].x,u[*s.cell].y-u[i].y};
            else if(s.kind==ViscousBoundaryKind2D::Velocity)delta={s.value.x-u[i].x,s.value.y-u[i].y};
            else if(s.kind==ViscousBoundaryKind2D::Slip){const double un=dot(u[i],s.normal);delta={-un*s.normal.x,-un*s.normal.y};}
            gradients[i][0].x+=delta.x*s.weight.x;gradients[i][0].y+=delta.x*s.weight.y;
            gradients[i][1].x+=delta.y*s.weight.x;gradients[i][1].y+=delta.y*s.weight.y;
        }
        for(const auto& g:gradients[i]){checked(g.x);checked(g.y);}out.rate[i]=checked(rowNorm_[i]/(areas_[i]*rho[i]));
    }
    for(std::size_t id=0;id<faces_.size();++id) {
        const auto& f=faces_[id];if((f.partner&&id>*f.partner)||(!f.neighbour&&f.kind==ViscousBoundaryKind2D::ZeroTraction))continue;
        auto g=gradients[f.owner];Vector2D uf=f.value;
        if(f.neighbour) {
            for(std::size_t k=0;k<2;++k){g[k].x=(1-f.weight)*g[k].x+f.weight*gradients[*f.neighbour][k].x;g[k].y=(1-f.weight)*g[k].y+f.weight*gradients[*f.neighbour][k].y;}
            uf={(1-f.weight)*(u[f.owner].x+dot(gradients[f.owner][0],f.ownerOffset))+f.weight*(u[*f.neighbour].x+dot(gradients[*f.neighbour][0],f.neighbourOffset)),
                (1-f.weight)*(u[f.owner].y+dot(gradients[f.owner][1],f.ownerOffset))+f.weight*(u[*f.neighbour].y+dot(gradients[*f.neighbour][1],f.neighbourOffset))};
        }
        if(f.neighbour||f.kind==ViscousBoundaryKind2D::Velocity) {
            const auto other=f.neighbour?u[*f.neighbour]:f.value;
            const double delta[2]={other.x-u[f.owner].x,other.y-u[f.owner].y};
            for(std::size_t k=0;k<2;++k){const double correction=(delta[k]-dot(g[k],f.d))/f.normalDistance;g[k].x+=correction*f.normal.x;g[k].y+=correction*f.normal.y;}
        }
        if(wallGradients_[id]) {
            g={};
            for(const auto& sample:wallGradients_[id]->samples) {
                const auto value=sample.boundary?faces_[sample.index].value:u[sample.index];const double delta[2]={value.x-f.value.x,value.y-f.value.y};
                for(std::size_t k=0;k<2;++k){g[k].x+=sample.weight.x*delta[k];g[k].y+=sample.weight.y*delta[k];}
            }
        }
        // Planar Newtonian gas, Stokes hypothesis. The 2/3 coefficient follows
        // the molecular constitutive law, not the dimension of the mesh.
        const double div=g[0].x+g[1].y,xx=viscosity_*(2*g[0].x-2./3*div),yy=viscosity_*(2*g[1].y-2./3*div),xy=viscosity_*(g[0].y+g[1].x);
        Vector2D traction{xx*f.area.x+xy*f.area.y,xy*f.area.x+yy*f.area.y};
        if(!f.neighbour&&f.kind==ViscousBoundaryKind2D::Slip){const double normal=dot(traction,f.normal);traction={normal*f.normal.x,normal*f.normal.y};uf={};}
        const std::array<double,3> flux{checked(-traction.x),checked(-traction.y),checked(-dot(uf,traction))};
        out.faceFlux[id]=flux;
        for(std::size_t k=0;k<3;++k){out.cellResidual[f.owner][k]+=flux[k];if(f.neighbour)out.cellResidual[*f.neighbour][k]-=flux[k];if(f.partner)out.faceFlux[*f.partner][k]=-flux[k];}
    }
    for(const auto& r:out.cellResidual)for(double v:r)checked(v);
    return out;
}
} // namespace cartmesh2d::fv
