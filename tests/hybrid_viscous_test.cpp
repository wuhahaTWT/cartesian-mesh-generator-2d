#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/HybridViscous2D.hpp"
#include <iomanip>
#include <iostream>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}require(rejected,"invalid hybrid viscous input accepted");}
double affineError=0,workError=0,reciprocityError=0,balanceError=0,tractionError=0,rawWorkDefect=0,rawTractionJump=0;
std::size_t maximumIterations=0;
std::vector<ViscousBoundary2D> prescribed(const FvMesh2D& mesh,const std::function<Vector2D(Point2D)>& value) {
    std::vector<ViscousBoundary2D> b;
    for(std::size_t f=0;f<mesh.faces.size();++f)if(!mesh.faces[f].neighbour)b.push_back({f,ViscousBoundaryKind2D::Velocity,value(mesh.faces[f].centre),{}});
    return b;
}
void sharedChecks(const FvMesh2D& mesh,const std::vector<Vector2D>& u,const HybridViscousResult2D& out,
    const std::vector<ViscousBoundary2D>& b,double mu) {
    std::array<double,3> sum{},boundary{},scale{};double speed=0;
    for(auto v:u)speed=std::max(speed,std::hypot(v.x,v.y));for(auto v:out.trace)speed=std::max(speed,std::hypot(v.x,v.y));
    // Normalize with the actual constitutive input scales, including when
    // rigid rotation makes the exact output vanish. Dividing SI defects by
    // '1 + near-zero output' mixes units and fails under a change of units.
    const double forceScale=mu*speed,workScaleInput=mu*speed*speed;
    const auto ratio=[](double error,double envelope){require(envelope>0||error==0,"zero error scale");return envelope>0?error/envelope:0;};
    std::vector<std::array<double,3>> actual(u.size());std::vector<const ViscousBoundary2D*> lookup(mesh.faces.size());for(const auto& bc:b)lookup[bc.face]=&bc;
    for(std::size_t f=0;f<mesh.faces.size();++f) {
        const auto& face=mesh.faces[f];const auto& q=out.faceFlux[f];const auto* bc=lookup[f];
        if(bc&&bc->partner)for(std::size_t k=0;k<3;++k)require(q[k]==-out.faceFlux[*bc->partner][k],"periodic flux not exactly opposite");
        if(bc&&bc->kind==ViscousBoundaryKind2D::ZeroTraction)for(double v:q)require(v==0,"open boundary exports traction/work");
        if(bc&&bc->kind==ViscousBoundaryKind2D::Slip)require(q[2]==0,"static slip wall exports work");
        else workError=std::max(workError,ratio(std::abs(q[2]-out.trace[f].x*q[0]-out.trace[f].y*q[1]),workScaleInput+std::abs(q[2])+std::abs(out.trace[f].x*q[0])+std::abs(out.trace[f].y*q[1])));
        for(std::size_t k=0;k<3;++k) {
            scale[k]+=std::abs(q[k]);actual[face.owner][k]+=q[k];if(face.neighbour)actual[*face.neighbour][k]-=q[k];
            if(!face.neighbour&&!(bc&&bc->partner))boundary[k]+=q[k];
        }
    }
    double mechanical=0,dissipation=0,workScale=0;
    for(std::size_t c=0;c<u.size();++c) {
        for(std::size_t k=0;k<3;++k){require(actual[c][k]==out.cellResidual[c][k],"residual does not use exported shared flux");sum[k]+=actual[c][k];}
        const double h=-actual[c][2]+u[c].x*actual[c][0]+u[c].y*actual[c][1];
        require(h==out.cellMechanicalHeating[c],"heating not based on shared momentum and work");
        require(out.cellDissipation[c]>=0&&out.rate[c]>0,"invalid local dissipation or spectral bound");
        mechanical+=h;dissipation+=out.cellDissipation[c];workScale+=std::abs(actual[c][2])+std::abs(u[c].x*actual[c][0])+std::abs(u[c].y*actual[c][1])+out.cellDissipation[c];
    }
    for(std::size_t k=0;k<3;++k)balanceError=std::max(balanceError,ratio(std::abs(sum[k]-boundary[k]),(k==2?workScaleInput:forceScale)+scale[k]));
    rawWorkDefect=std::max(rawWorkDefect,std::abs(mechanical-dissipation));rawTractionJump=std::max(rawTractionJump,std::max(out.maximumTractionJump,out.maximumBoundaryConstraintResidual));
    workError=std::max(workError,ratio(std::abs(mechanical-dissipation),workScaleInput+workScale));
    tractionError=std::max(tractionError,ratio(std::max(out.maximumTractionJump,out.maximumBoundaryConstraintResidual),forceScale+scale[0]+scale[1]));
    maximumIterations=std::max(maximumIterations,out.iterations);
}
void affine(double mu,double lengthScale,double speedScale,double angle) {
    auto mesh=fv_test::rectangle(5,4,2,true);
    const auto rotate=[&](Vector2D p){return Vector2D{std::cos(angle)*p.x-std::sin(angle)*p.y,std::sin(angle)*p.x+std::cos(angle)*p.y};};
    for(auto& c:mesh.cells){const auto p=rotate({c.centre.x,c.centre.y});c.centre={lengthScale*p.x,lengthScale*p.y};c.area*=lengthScale*lengthScale;}
    for(auto& f:mesh.faces){const auto p=rotate({f.centre.x,f.centre.y}),s=rotate(f.areaVector),v=rotate(f.correction);f.centre={lengthScale*p.x,lengthScale*p.y};f.areaVector={lengthScale*s.x,lengthScale*s.y};f.correction={lengthScale*v.x,lengthScale*v.y};}
    for(const auto a:{std::array<double,4>{0,0,0,0},{0,-2,2,0},{.4,-.7,1.2,-.9}}) {
        const auto value=[&](Point2D p){return Vector2D{speedScale*(.3+(a[0]*p.x+a[1]*p.y)/lengthScale),speedScale*(-.6+(a[2]*p.x+a[3]*p.y)/lengthScale)};};
        const auto b=prescribed(mesh,value);HybridViscousOperator2D op(mesh,b,mu);
        std::vector<Vector2D> u;for(const auto& c:mesh.cells)u.push_back(value(c.centre));
        const auto out=op.evaluate(u,std::vector<double>(u.size(),1),2e-14);sharedChecks(mesh,u,out,b,mu);
        const double xx=mu*speedScale/lengthScale*(2*a[0]-2./3*(a[0]+a[3])),yy=mu*speedScale/lengthScale*(2*a[3]-2./3*(a[0]+a[3])),xy=mu*speedScale/lengthScale*(a[1]+a[2]);
        for(std::size_t f=0;f<mesh.faces.size();++f) {
            const auto s=mesh.faces[f].areaVector,v=value(mesh.faces[f].centre);
            const std::array<double,3> expected{-xx*s.x-xy*s.y,-xy*s.x-yy*s.y,-v.x*(xx*s.x+xy*s.y)-v.y*(xy*s.x+yy*s.y)};
            for(std::size_t k=0;k<3;++k)affineError=std::max(affineError,std::abs(out.faceFlux[f][k]-expected[k])/(mu*speedScale*(k==2?speedScale:1)+std::abs(expected[k])));
        }
        if(a==std::array<double,4>{0,0,0,0})require(out.dissipation==0&&out.iterations==0,"uniform translation not preserved exactly");
    }
}
void arbitraryAndReciprocity() {
    auto mesh=fv_test::rectangle(5,4,2,true);const auto b=prescribed(mesh,[](Point2D){return Vector2D{};});HybridViscousOperator2D op(mesh,b,.37);
    std::vector<Vector2D> u,v;std::vector<double> rho(mesh.cells.size(),1.2);
    for(std::size_t c=0;c<mesh.cells.size();++c){u.push_back({std::sin(double(c)),std::cos(double(3*c))});v.push_back({std::cos(double(2*c)),std::sin(double(5*c))});}
    const auto a=op.evaluate(u,rho,2e-14),z=op.evaluate(v,rho,2e-14);sharedChecks(mesh,u,a,b,.37);sharedChecks(mesh,v,z,b,.37);
    double uv=0,vu=0,envelope=0,bound=0,quadratic=0;
    for(std::size_t c=0;c<u.size();++c){const double x=u[c].x*z.cellResidual[c][0]+u[c].y*z.cellResidual[c][1],y=v[c].x*a.cellResidual[c][0]+v[c].y*a.cellResidual[c][1];uv+=x;vu+=y;envelope+=std::abs(x)+std::abs(y);quadratic+=u[c].x*a.cellResidual[c][0]+u[c].y*a.cellResidual[c][1];bound+=rho[c]*mesh.cells[c].area*a.rate[c]*dot(u[c],u[c]);}
    reciprocityError=std::abs(uv-vu)/(1+envelope);require(quadratic>=0&&quadratic<=bound,"condensed momentum operator violates its energy bound");
    const auto n=op.traceUnknowns();std::vector<double> matrix(n*n);op.visitTraceMatrix([&](auto i,auto j,double x){matrix[i*n+j]=x;});
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)require(matrix[i*n+j]==matrix[j*n+i],"trace matrix not exactly symmetric");
    // Small native dense Cholesky of the ACTUAL assembled trace matrix. This
    // checks anchoring/assembly on this mesh, not arbitrary-grid qualification.
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<=i;++j){double x=matrix[i*n+j];for(std::size_t k=0;k<j;++k)x-=matrix[i*n+k]*matrix[j*n+k];if(i==j){require(x>0,"anchored trace matrix not positive definite");matrix[i*n+j]=std::sqrt(x);}else matrix[i*n+j]=x/matrix[j*n+j];}
}
void boundariesAndKernel() {
    const auto mesh=fv_test::rectangle(4,4,2);std::vector<ViscousBoundary2D> b;
    for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour) {
        const auto p=mesh.faces[id].centre;
        if(p.x==0||p.x==2){std::optional<std::size_t> partner;for(std::size_t j=0;j<mesh.faces.size();++j)if(!mesh.faces[j].neighbour&&mesh.faces[j].centre.y==p.y&&mesh.faces[j].centre.x==2-p.x)partner=j;require(partner.has_value(),"missing test periodic partner");b.push_back({id,ViscousBoundaryKind2D::Periodic,{},partner});}
        else b.push_back({id,ViscousBoundaryKind2D::Velocity,{.3+.7*p.y,0},{}});
    }
    std::vector<Vector2D> u;for(const auto& c:mesh.cells)u.push_back({.3+.7*c.centre.y,0});std::vector<double> rho(u.size(),1);
    HybridViscousOperator2D periodic(mesh,b,.37);sharedChecks(mesh,u,periodic.evaluate(u,rho,2e-14),b,.37);
    auto mixed=b;for(auto& bc:mixed)if(bc.partner){bc.partner.reset();bc.kind=ViscousBoundaryKind2D::ZeroTraction;}else if(mesh.faces[bc.face].centre.y==0){bc.kind=ViscousBoundaryKind2D::Slip;bc.velocity={};}
    HybridViscousOperator2D open(mesh,mixed,.37);sharedChecks(mesh,u,open.evaluate(u,rho,2e-14),mixed,.37);
    auto noAnchor=b;for(auto& bc:noAnchor)if(!bc.partner){bc.kind=ViscousBoundaryKind2D::Slip;bc.velocity={};}
    rejects([&]{HybridViscousOperator2D unsupported(mesh,noAnchor,.37);});
    // Exact failure mechanism: alternating local rigid rotations share their
    // face velocity on an orthogonal grid while every cell velocity is zero.
    // These traces have no local strain/stabilization energy, and slip walls
    // do not constrain them. Keep this counterexample when improving the form.
    std::vector<Vector2D> trace(mesh.faces.size());std::vector<bool> assigned(mesh.faces.size());double dissipation=0;
    for(std::size_t c=0;c<mesh.cells.size();++c) {
        const double omega=((c%4+c/4)%2)?-1:1;std::vector<Vector2D> values;HybridViscousCell2D local(mesh,c,.37);
        for(auto f:local.faces()){const auto p=mesh.faces[f].centre,q=mesh.cells[c].centre;const Vector2D t{-omega*(p.y-q.y),omega*(p.x-q.x)};if(assigned[f])require(trace[f].x==t.x&&trace[f].y==t.y,"checkerboard traces not shared");trace[f]=t;assigned[f]=true;values.push_back(t);}dissipation+=local.evaluate({},values).dissipation;
    }
    require(dissipation==0,"rigid trace counterexample changed");
    auto bad=b;bad.pop_back();rejects([&]{HybridViscousOperator2D incomplete(mesh,bad,.37);});
    rho[0]=0;rejects([&]{(void)periodic.evaluate(u,rho);});
    const auto one=fv_test::rectangle(1,1,1);const auto fixed=prescribed(one,[](Point2D p){return Vector2D{p.y,0};});HybridViscousOperator2D single(one,fixed,.37);require(single.traceUnknowns()==0,"single prescribed cell has free trace");sharedChecks(one,{{.5,0}},single.evaluate({{.5,0}},{1}),fixed,.37);
}
}
int main(){try {
    for(double mu:{1e-5,.37,1e5})for(double length:{1e-3,1.,1e3})for(double speed:{1e-4,1.,1e4})affine(mu,length,speed,.731);
    arbitraryAndReciprocity();boundariesAndKernel();
    std::cout<<std::setprecision(17)<<"{\"affineRelativeError\":"<<affineError<<",\"sharedMechanicalWorkError\":"<<workError<<",\"momentumReciprocityError\":"<<reciprocityError<<",\"conservativeBalanceError\":"<<balanceError<<",\"tractionEquilibriumError\":"<<tractionError<<",\"maximumPhysicalWorkDefect\":"<<rawWorkDefect<<",\"maximumPhysicalTractionJump\":"<<rawTractionJump<<",\"maximumIterations\":"<<maximumIterations<<"}\n";
    // Dimensionless algebra checks: 65536 eps allows the accumulated global
    // trace-solve error at requested 2e-14 relative tolerance. Physical accuracy,
    // time stability, curved-grid cost and full Euler coupling remain separate.
    const double budget=65536*std::numeric_limits<double>::epsilon();
    require(affineError<budget,"affine global traction/work error");require(workError<budget,"shared mechanical heating differs from dissipation");require(reciprocityError<budget,"condensed momentum operator not reciprocal");require(balanceError<budget,"shared flux not conservative");require(tractionError<budget,"trace constraints not satisfied");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
