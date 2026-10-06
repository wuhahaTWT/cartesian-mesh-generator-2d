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
// This research natural stress outlet differs from the collocated product's
// pressure Dirichlet / zero-normal-gradient velocity boundary stencil.
#define CARTMESH_P1_OSEEN_NO_MAIN
#include "native-laminar-p1-oseen.cpp"

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
    Transport(const Fixture& f,int t,const std::string& problem,double nu,bool nonlinear,bool open,int order,const Vec& previous):Oseen(f,t,problem,nu,nonlinear,open,order,previous){
        int m=e.a.m;for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j){e.matrix(i,j)-=convection(i,j);convection(i,j)=0;}
        for(std::size_t k=0;k<lift.tri.size();++k){const auto& tr=lift.tri[k];
            for(auto q:tr.quadrature(order)){auto z=tr.ref(q.p);auto b=betaAt(int(k),q.p);Vec r0(2*m),r1(2*m),d0(2*m),d1(2*m);auto shape=tr.shape(q.p);
                std::array<Vector2D,8> dx{{{0,0},{1,0},{0,0},{0,0},{0,1},{0,0},{2*z.x,z.y},{z.y,0}}},dy{{{0,0},{0,0},{1,0},{0,0},{0,0},{0,1},{0,z.x},{z.x,2*z.y}}};
                double bx=(tr.j1.y*b.x-tr.j1.x*b.y)/(tr.h*tr.det),by=(-tr.j0.y*b.x+tr.j0.x*b.y)/(tr.h*tr.det);
                for(int l=0;l<8;++l){Vector2D d{bx*dx[l].x+by*dy[l].x,bx*dx[l].y+by*dy[l].y};d={(tr.j0.x*d.x+tr.j1.x*d.y)/tr.det,(tr.j0.y*d.x+tr.j1.y*d.y)/tr.det};
                    for(int j=0;j<2*m;++j){double c=lift.coefficients(8*k+l,j);r0[j]+=c*shape[l].x;r1[j]+=c*shape[l].y;d0[j]+=c*d.x;d1[j]+=c*d.y;}}
                for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j)convection(i,j)-=q.w*(d0[i]*r0[j]+d1[i]*r1[j]);
            }
            // Each internal radial edge once, normal outward from triangle k.
            int prev=(int(k)+int(lift.tri.size())-1)%int(lift.tri.size());auto d=tr.a-tr.c;Vector2D S{d.y,-d.x};
            for(auto [z,w]:upwindRule(order,dot(betaAt(int(k),tr.c),S),dot(betaAt(int(k),tr.a),S))){Point2D p{tr.c.x+z*d.x,tr.c.y+z*d.y};double flux=dot(betaAt(int(k),p),S);std::vector<Vector2D> left(2*m),right(2*m);
                for(int j=0;j<2*m;++j){left[j]=lift.value(int(k),j,p);right[j]=lift.value(prev,j,p);}
                for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j)convection(i,j)+=w*flux*dot(difference(left[i],right[i]),flux>=0?left[j]:right[j]);
            }
            // Original polygon face remains a shared two-component P1 unknown.
            int l=lift.faceLocal[k];const auto& face=f.mesh.faces[f.mesh.cells[t].faces[l]];double sign=face.owner==std::size_t(t)?1:-1;
            auto fluxEnd=[&](double ss){return sign*((beta[3+2*l]+ss*beta[4+2*l])*face.areaVector.x+(beta[m+3+2*l]+ss*beta[m+4+2*l])*face.areaVector.y);};
            for(auto [z,w]:upwindRule(order,fluxEnd(-.5),fluxEnd(.5))){double ss=z-.5;auto Sface=face.areaVector;Point2D p{face.centre.x-ss*Sface.y,face.centre.y+ss*Sface.x};double flux=sign*((beta[3+2*l]+ss*beta[4+2*l])*Sface.x+(beta[m+3+2*l]+ss*beta[m+4+2*l])*Sface.y);std::vector<Vector2D> inside(2*m),facev(2*m);
                for(int j=0;j<2*m;++j)inside[j]=lift.value(int(k),j,p);
                facev[3+2*l].x=1;facev[4+2*l].x=ss;facev[m+3+2*l].y=1;facev[m+4+2*l].y=ss;
                for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j){convection(i,j)+=w*flux*dot(difference(inside[i],facev[i]),flux>=0?inside[j]:facev[j]);if(boundaryKind(face,open)==1)convection(i,j)+=w*flux*dot(facev[i],facev[j]);}
            }
        }
        for(int i=0;i<2*m;++i)for(int j=0;j<2*m;++j)e.matrix(i,j)+=convection(i,j);
        Mat ii(8,8);for(int i=0;i<8;++i)for(int j=0;j<8;++j)ii(i,j)=e.matrix(e.inside[i],e.inside[j]);DenseLU lu(ii.v,8);
        for(std::size_t j=0;j<e.outside.size();++j){Vec r(8);for(int i=0;i<8;++i)r[i]=e.matrix(e.inside[i],e.outside[j]);r=lu.solve(r);for(int i=0;i<8;++i)e.eliminated(i,int(j))=r[i];}Vec r(8);for(int i=0;i<8;++i)r[i]=e.rhs[e.inside[i]];e.loadInternal=lu.solve(r);
    }
    Vector2D transported(const Vec& v,int k,Point2D p)const{Vector2D u{};for(int j=0;j<2*e.a.m;++j){auto r=lift.value(k,j,p);u.x+=v[j]*r.x;u.y+=v[j]*r.y;}return u;}
    double energyIdentity(const Vec& v)const{
        int m=e.a.m;double result=0;
        for(std::size_t k=0;k<lift.tri.size();++k){const auto& tr=lift.tri[k];for(auto q:tr.quadrature(order)){auto phi=e.a.basis.phi(q.p);double div=0;for(int j=0;j<3;++j)for(int l=0;l<m;++l)div+=phi[j]*(e.a.gx(j,l)*beta[l]+e.a.gy(j,l)*beta[m+l]);auto u=transported(v,int(k),q.p);result+=.5*q.w*div*dot(u,u);}
            int prev=(int(k)+int(lift.tri.size())-1)%int(lift.tri.size());auto d=tr.a-tr.c;Vector2D S{d.y,-d.x};
            for(auto [z,w]:upwindRule(order,dot(betaAt(int(k),tr.c),S),dot(betaAt(int(k),tr.a),S))){Point2D p{tr.c.x+z*d.x,tr.c.y+z*d.y};double flux=dot(betaAt(int(k),p),S);auto jump=difference(transported(v,int(k),p),transported(v,prev,p));result+=.5*w*std::abs(flux)*dot(jump,jump);}
            int l=lift.faceLocal[k];const auto& face=f.mesh.faces[f.mesh.cells[t].faces[l]];double sign=face.owner==std::size_t(t)?1:-1;
            auto fluxEnd=[&](double ss){return sign*((beta[3+2*l]+ss*beta[4+2*l])*face.areaVector.x+(beta[m+3+2*l]+ss*beta[m+4+2*l])*face.areaVector.y);};
            for(auto [z,w]:upwindRule(order,fluxEnd(-.5),fluxEnd(.5))){double ss=z-.5;auto sf=face.areaVector;Point2D p{face.centre.x-ss*sf.y,face.centre.y+ss*sf.x};double flux=sign*((beta[3+2*l]+ss*beta[4+2*l])*sf.x+(beta[m+3+2*l]+ss*beta[m+4+2*l])*sf.y);Vector2D uf{v[3+2*l]+ss*v[4+2*l],v[m+3+2*l]+ss*v[m+4+2*l]};auto jump=difference(transported(v,int(k),p),uf);result+=.5*w*(std::abs(flux)*dot(jump,jump)+(boundaryKind(face,open)==1?1:-1)*flux*dot(uf,uf));}
        }return result;
    }
};
int main(int argc,char** argv){return runOseen<Transport>(argc,argv);}
