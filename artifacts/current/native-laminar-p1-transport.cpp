// RT1-conforming test/load and conservative upwind transport on native fans.
// Uses same coupled hybrid stress, pressure and exact original face topology.
// On each fan triangle: -int R(u).(beta.grad R(v)); internal radial fluxes
// use beta_n R(u)_up [R(v)], original faces beta_n u_up (R(v)-v_F).
// The natural outlet adds beta_n u_F.v_F. beta=R(previous u) has the same
// normal moments and P1 divergence as the common pressure B block.
// Hence c(beta;u,u)=1/2 int div(beta)|R(u)|^2 + upwind jump dissipation
// minus 1/2 prescribed-face beta_n|u_F|^2 plus 1/2 outlet beta_n|u_F|^2.
// This is static pressure, NOT a Bernoulli variable. Positive outlet flow is
// not assumed by the implementation; backflow can inject energy and is reported.
// R is affine-exact but not a general quadratic velocity reconstruction:
// lost quadratic-channel exactness on irregular coarse cells is retained.
// Boundary mode "open" retains the earlier natural total-stress outlet.
// Boundary mode "pressure" uses the product decomposition: fixed static
// pressure via the weak pseudo-traction, retaining the transpose-gradient
// symmetric correction. This does not separately impose strong p and d_n u
// traces or prove equivalence to the collocated discrete stencil.
#define CARTMESH_P1_OSEEN_NO_MAIN
#include "native-laminar-p1-oseen.cpp"
#include "cartmesh2d/fv/detail/CompatibleFlowTransport2D.hpp"

// beta.n is P1 on every RT1 edge. Split at its exact interior zero so
// upwind/absolute-value integrals are polynomial on each interval, rather
// than sampling a nonsmooth flux with an unsplit Gauss rule. No fitted cutoff.
std::vector<std::pair<double,double>> upwindRule(int order,double left,double right){
    std::vector<double> ends{0,1};
    if((left<0&&right>0)||(left>0&&right<0))ends.insert(ends.begin()+1,left/(left-right));
    std::vector<std::pair<double,double>> out;auto rule=gauss(order);
    for(std::size_t i=1;i<ends.size();++i)for(auto [z,w]:rule)out.push_back({ends[i-1]+(ends[i]-ends[i-1])*z,(ends[i]-ends[i-1])*w});
    return out;
}
Vector2D difference(Vector2D a,Vector2D b){return {a.x-b.x,a.y-b.y};}
struct Transport : Oseen {
    Transport(const Fixture& f,int t,const std::string& problem,double nu,bool nonlinear,int outletModel,int order,const Vec& previous):Oseen(f,t,problem,nu,nonlinear,outletModel,order,previous){
        int m=e.a.m;for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j){e.matrix(i,j)-=convection(i,j);convection(i,j)=0;}
        std::vector<bool> open(f.mesh.cells[t].faces.size());
        for(std::size_t id=0;id<open.size();++id)open[id]=boundaryKind(f.mesh.faces[f.mesh.cells[t].faces[id]],outletModel,problem)==1;
        convection=transportMatrix(f.mesh,t,e.a,lift,beta,open,order);
        for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j)e.matrix(i,j)+=convection(i,j);
        e.condense();
    }
    Vector2D transported(const Vec& v,int k,Point2D p)const{Vector2D u{};for(std::size_t j=0;j<2*e.a.m;++j){auto r=lift.value(k,j,p);u.x+=v[j]*r.x;u.y+=v[j]*r.y;}return u;}
    double energyIdentity(const Vec& v)const{
        int m=e.a.m;double result=0;
        for(std::size_t k=0;k<lift.tri.size();++k){const auto& tr=lift.tri[k];for(auto q:tr.quadrature(order)){auto phi=e.a.basis.phi(q.p);double div=0;for(int j=0;j<3;++j)for(int l=0;l<m;++l)div+=phi[j]*(e.a.gx(j,l)*beta[l]+e.a.gy(j,l)*beta[m+l]);auto u=transported(v,int(k),q.p);result+=.5*q.w*div*dot(u,u);}
            int prev=(int(k)+int(lift.tri.size())-1)%int(lift.tri.size());auto d=tr.a-tr.c;Vector2D S{d.y,-d.x};
            for(auto [z,w]:upwindRule(order,dot(betaAt(int(k),tr.c),S),dot(betaAt(int(k),tr.a),S))){Point2D p{tr.c.x+z*d.x,tr.c.y+z*d.y};double flux=dot(betaAt(int(k),p),S);auto jump=difference(transported(v,int(k),p),transported(v,prev,p));result+=.5*w*std::abs(flux)*dot(jump,jump);}
            int l=lift.faceLocal[k];const auto& face=f.mesh.faces[f.mesh.cells[t].faces[l]];double sign=face.owner==std::size_t(t)?1:-1;
            auto fluxEnd=[&](double ss){return sign*((beta[3+2*l]+ss*beta[4+2*l])*face.areaVector.x+(beta[m+3+2*l]+ss*beta[m+4+2*l])*face.areaVector.y);};
            for(auto [z,w]:upwindRule(order,fluxEnd(-.5),fluxEnd(.5))){double ss=z-.5;auto sf=face.areaVector;Point2D p{face.centre.x-ss*sf.y,face.centre.y+ss*sf.x};double flux=sign*((beta[3+2*l]+ss*beta[4+2*l])*sf.x+(beta[m+3+2*l]+ss*beta[m+4+2*l])*sf.y);Vector2D uf{v[3+2*l]+ss*v[4+2*l],v[m+3+2*l]+ss*v[m+4+2*l]};auto jump=difference(transported(v,int(k),p),uf);result+=.5*w*(std::abs(flux)*dot(jump,jump)+(boundaryKind(face,outletModel,problem)==1?1:-1)*flux*dot(uf,uf));}
        }return result;
    }
};
#ifndef CARTMESH_P1_TRANSPORT_NO_MAIN
int main(int argc,char** argv){return runOseen<Transport>(argc,argv);}
#endif
