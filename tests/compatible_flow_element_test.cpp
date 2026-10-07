#include "cartmesh2d/fv/detail/CompatibleFlowElement2D.hpp"
#include "fixtures/PolygonMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::fv::detail::compatible;
namespace {
double fieldError=0, pressureError=0, internalError=0, gradientLoadError=0;
std::size_t cellsChecked=0;
void require(bool ok,const char* what){if(!ok)throw std::runtime_error(what);}
// Unit-size fixtures: U=L=nu=1, kinematic p uses U^2. This is a roundoff
// regression bound for small local sums/solves, not a physical accuracy gate.
// 4096 eps per local DOF leaves room for the local saddle factorization.
void near(double actual,double expected,std::size_t dofs,const char* what){
    const double bound=4096*std::numeric_limits<double>::epsilon()*double(dofs)*(1+std::abs(expected));
    require(std::isfinite(actual)&&std::abs(actual-expected)<=bound,what);
}
Vector2D velocity(Point2D p){return {p.x*p.x-p.y*p.y,-2*p.x*p.y};}
double pressure(Point2D p){return 1+2*p.x-3*p.y;}
double phi(Point2D p){return p.x*p.x*p.x+p.x*p.y*p.y;}
double diameter(const FvMesh2D& mesh,std::size_t cell){
    std::vector<Point2D> v;
    for(auto id:mesh.cells[cell].faces){const auto& f=mesh.faces[id];
        v.push_back({f.centre.x+f.areaVector.y/2,f.centre.y-f.areaVector.x/2});
        v.push_back({f.centre.x-f.areaVector.y/2,f.centre.y+f.areaVector.x/2});}
    double result=0;for(auto a:v)for(auto b:v)result=std::max(result,std::hypot(a.x-b.x,a.y-b.y));return result;
}
Vec project(const FvMesh2D& mesh,std::size_t cell,const P1Local& a,bool hydrostatic){
    Vec state(2*a.m+3);std::array<Vec,3> moment{Vec(3),Vec(3),Vec(3)};
    for(const auto& q:a.q){const auto b=a.basis.phi(q.p);const auto u=hydrostatic?Vector2D{}:velocity(q.p);
        for(std::size_t j=0;j<3;++j){moment[0][j]+=q.w*b[j]*u.x;moment[1][j]+=q.w*b[j]*u.y;moment[2][j]+=q.w*b[j]*(hydrostatic?phi(q.p):pressure(q.p));}}
    const DenseLU mass(a.mass.v,3);
    for(std::size_t c=0;c<3;++c){const auto coeff=mass.solve(moment[c]);for(std::size_t j=0;j<3;++j)state[(c<2?c*a.m:2*a.m)+j]=coeff[j];}
    if(!hydrostatic)for(std::size_t l=0;l<mesh.cells[cell].faces.size();++l){const auto& f=mesh.faces[mesh.cells[cell].faces[l]];
        for(auto [z,w]:gauss(6)){const double s=z-.5;const auto u=velocity({f.centre.x-s*f.areaVector.y,f.centre.y+s*f.areaVector.x});
            state[3+2*l]+=w*u.x;state[4+2*l]+=12*w*s*u.x;state[a.m+3+2*l]+=w*u.y;state[a.m+4+2*l]+=12*w*s*u.y;}}
    return state;
}
void check(const FvMesh2D& mesh){
    for(std::size_t cell=0;cell<mesh.cells.size();++cell){
        P1System e(P1Local(mesh,cell,diameter(mesh,cell),6));const auto& a=e.a;
        Lift lift(mesh,cell,a,6);const auto state=project(mesh,cell,a,false);
        const auto n=e.rhs.size();double area=0;for(const auto& q:a.q)area+=q.w;
        near(area,mesh.cells[cell].area,n,"original polygon area changed");
        e.rhs=liftedBodyForce(a,lift,6,[](Point2D){return Vector2D{2,-3};});
        e.condense();Vec retained;for(auto id:e.outside)retained.push_back(state[id]);
        const auto recovered=e.recover(retained);
        for(auto id:e.inside){double residual=-e.rhs[id];for(std::size_t j=0;j<n;++j)residual+=e.matrix(id,j)*state[j];
            internalError=std::max(internalError,std::abs(residual));near(residual,0,n,"quadratic Stokes internal equilibrium");
            near(recovered[id],state[id],n,"condensed analytic state differs");}
        for(const auto& q:a.q){const auto theta=a.basis.theta(q.p);const auto b=a.basis.phi(q.p);const auto exact=velocity(q.p);Vector2D got{};double p=0;
            for(std::size_t j=0;j<a.m;++j)for(std::size_t k=0;k<6;++k){got.x+=theta[k]*a.potential(k,j)*recovered[j];got.y+=theta[k]*a.potential(k,j)*recovered[a.m+j];}
            for(std::size_t k=0;k<3;++k)p+=b[k]*recovered[2*a.m+k];
            fieldError=std::max({fieldError,std::abs(got.x-exact.x),std::abs(got.y-exact.y)});pressureError=std::max(pressureError,std::abs(p-pressure(q.p)));
            near(got.x,exact.x,n,"quadratic x velocity reconstruction");near(got.y,exact.y,n,"quadratic y velocity reconstruction");near(p,pressure(q.p),n,"affine pressure reconstruction");}
        near(lift.traceResidual,0,n,"RT1 radial normal trace continuity");near(lift.divResidual,0,n,"RT1 prescribed divergence");
        const auto hydro=project(mesh,cell,a,true);
        const auto load=liftedBodyForce(a,lift,6,[](Point2D p){return Vector2D{3*p.x*p.x+p.y*p.y,2*p.x*p.y};});
        // grad(phi) must be balanced by P1-projected pressure for every
        // interior velocity test; this cannot be replaced by a small speed.
        for(std::size_t c=0;c<2;++c)for(std::size_t k=0;k<3;++k){const auto row=c*a.m+k;double r=-load[row];
            for(std::size_t j=0;j<n;++j)r+=e.matrix(row,j)*hydro[j];gradientLoadError=std::max(gradientLoadError,std::abs(r));near(r,0,n,"gradient body force lost pressure compatibility");}
        bool rejected=false;try{e.recover(Vec(e.outside.size()-1));}catch(const std::invalid_argument&){rejected=true;}require(rejected,"invalid retained size accepted");
        rejected=false;try{liftedBodyForce(a,lift,6,[](Point2D){return Vector2D{std::numeric_limits<double>::quiet_NaN(),0};});}catch(const std::runtime_error&){rejected=true;}require(rejected,"nonfinite force accepted");
        // A failed reassembly must not make the previous local factor look
        // valid for a new candidate; the owning solver keeps its accepted field.
        for(auto j:e.inside)e.matrix(e.inside.front(),j)=0;
        rejected=false;try{e.condense();}catch(const std::runtime_error&){rejected=true;}require(rejected,"singular interior block accepted");
        rejected=false;try{e.recover(retained);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"failed reassembly exposed stale recovery");
        ++cellsChecked;
    }
}
}
int main()try{
    check(makeFvMesh2D(test::fromPolygons({{{{0,0},{1,0},{1.7,1},{.7,1}}}})));
    check(makeFvMesh2D(test::fromPolygons({{{{0,.127},{1,.437},{1,1},{0,1}}}})));
    check(makeFvMesh2D(test::fromPolygons({{{{0,0},{1,0},{1,1},{1,2},{0,2}}},{{{1,0},{2,0},{2,1},{1,1}}},{{{1,1},{2,1},{2,2},{1,2}}}})));
    check(makeFvMesh2D(test::fromPolygons({{{{0,0},{1,0},{0,1}}},{{{1,0},{1,1},{0,1}}}})));
    bool rejected=false;try{gauss(2);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"underintegrated polynomial rule accepted");
    rejected=false;try{DenseLU(Vec(3),2);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"invalid local matrix size accepted");
    std::cout<<std::setprecision(17)<<"{\"cells\":"<<cellsChecked<<",\"velocityMax\":"<<fieldError<<",\"pressureMax\":"<<pressureError<<",\"internalEquilibriumMax\":"<<internalError<<",\"gradientLoadBalanceMax\":"<<gradientLoadError<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
