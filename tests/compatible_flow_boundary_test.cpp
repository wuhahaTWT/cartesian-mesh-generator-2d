#include "cartmesh2d/fv/CompatibleFlowBoundary2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/fv/detail/CompatibleFlowElement2D.hpp"
#include "fixtures/PolygonMesh2D.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::fv::detail::compatible;
namespace {
void require(bool ok,const char* m){if(!ok)throw std::runtime_error(m);}
struct Case {
    double angle;bool symmetry,outlet;double nu=.2;FvMesh2D mesh;
    Vector2D rotate(Vector2D p)const{const auto c=std::cos(angle),s=std::sin(angle);return {c*p.x-s*p.y,s*p.x+c*p.y};}
    Point2D coordinates(Point2D p)const{const auto c=std::cos(angle),s=std::sin(angle);return {c*p.x+s*p.y,-s*p.x+c*p.y};}
    Vector2D baseU(Point2D p)const{return symmetry?Vector2D{1-p.y*p.y,0}:Vector2D{.2+.3*p.x+.7*p.y,-.1-.4*p.x-.3*p.y};}
    Vector2D velocity(Point2D p)const{return rotate(baseU(coordinates(p)));}
    double pressure(Point2D p)const{p=coordinates(p);return symmetry?2*nu*(1-p.x):1+2*p.x-3*p.y;}
    Vector2D traction(Point2D p,Vector2D normal)const{
        const auto q=coordinates(p);const auto n=coordinates({normal.x,normal.y});const auto pr=pressure(p);
        return rotate(symmetry?Vector2D{-pr*n.x-2*nu*q.y*n.y,-2*nu*q.y*n.x-pr*n.y}:Vector2D{(.6*nu-pr)*n.x+.3*nu*n.y,.3*nu*n.x+(-.6*nu-pr)*n.y});
    }
    Case(double a,bool symmetric,bool open):angle(a),symmetry(symmetric),outlet(open){
        constexpr int n=3;std::vector<Polygon2D> polygons;auto node=[&](int i,int j){const auto p=rotate({double(i)/n,double(j)/n});return Point2D{p.x,p.y};};
        for(int j=0;j<n;++j)for(int i=0;i<n;++i)polygons.push_back({{node(i,j),node(i+1,j),node(i+1,j+1),node(i,j+1)}});
        mesh=makeFvMesh2D(cartmesh2d::test::fromPolygons(polygons));
    }
    CompatibleFlowControls2D controls()const{
        CompatibleFlowControls2D c;c.viscosity=nu;c.equation=symmetry?CompatibleEquation2D::Stokes:CompatibleEquation2D::NavierStokes;
        c.globalization=symmetry?CompatibleGlobalization2D::Backtracking:CompatibleGlobalization2D::PseudoTime;
        if(outlet)c.pressureInverse=CompatiblePressureInverse2D::DiagonalSchur;
        if(!symmetry)c.acceleration=[this](Point2D p){const auto u=baseU(coordinates(p));return rotate({2+.3*u.x+.7*u.y,-3-.4*u.x-.3*u.y});};
        for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour){const auto& f=mesh.faces[id];const auto x=coordinates(f.centre);const double length=std::hypot(f.areaVector.x,f.areaVector.y);const Vector2D n{f.areaVector.x/length,f.areaVector.y/length};CompatibleBoundary2D b;b.face=id;
            if((!symmetry&&(std::abs(x.y)<1e-12||std::abs(x.y-1)<1e-12))||(symmetry&&std::abs(x.y)<1e-12)){
                b.kind=CompatibleBoundaryKind2D::NormalVelocity;
                if(!symmetry){b.normalVelocity=[this,n](Point2D p){return dot(velocity(p),n);};b.tangentialTraction=[this,n](Point2D p){return dot(traction(p,n),Vector2D{-n.y,n.x});};}
            }else if(outlet&&std::abs(x.x-1)<1e-12){b.kind=CompatibleBoundaryKind2D::Traction;b.value=[this,n](Point2D p){return traction(p,n);};}
            else b.value=[this](Point2D p){return velocity(p);};c.boundaries.push_back(std::move(b));}
        return c;
    }
};
double diameter(const FvMesh2D& mesh,std::size_t t){std::vector<Point2D> v;for(auto id:mesh.cells[t].faces){const auto& f=mesh.faces[id];v.push_back({f.centre.x+f.areaVector.y/2,f.centre.y-f.areaVector.x/2});v.push_back({f.centre.x-f.areaVector.y/2,f.centre.y+f.areaVector.x/2});}double h=0;for(auto p:v)for(auto q:v)h=std::max(h,std::hypot(p.x-q.x,p.y-q.y));return h;}
void backflow() {
    const Case fixture(0,false,true);const auto& mesh=fixture.mesh;
    CompatibleFlowControls2D c;c.viscosity=.2;c.equation=CompatibleEquation2D::Stokes;c.globalization=CompatibleGlobalization2D::Backtracking;
    for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour){
        const auto& face=mesh.faces[id];CompatibleBoundary2D b;b.face=id;
        if(face.centre.x<1e-12)b.value=[](Point2D p){return Vector2D{p.y-.1,0};};
        else if(face.centre.x>1-1e-12){b.kind=CompatibleBoundaryKind2D::Traction;b.value=[](Point2D){return Vector2D{0,.2};};}
        else {b.kind=CompatibleBoundaryKind2D::NormalVelocity;b.tangentialTraction=[](Point2D){return -.2;};}
        c.boundaries.push_back(std::move(b));
    }
    const auto allowed=solveCompatibleIncompressible2D(mesh,c);require(allowed.converged(),"Backflow reference shear failed");
    bool negativeEndpoint=false;
    for(auto& b:c.boundaries)if(b.kind==CompatibleBoundaryKind2D::Traction){const auto& u=allowed.lastAccepted->faces[b.face];require(u[0]>0,"Reference must have positive outlet face means");negativeEndpoint=negativeEndpoint||u[0]-.5*std::abs(u[1])<-.05;b.rejectBackflow=true;}
    require(negativeEndpoint,"Reference lacks reversed P1 endpoint");
    const auto rejected=solveCompatibleIncompressible2D(mesh,c);
    require(rejected.stop==CompatibleFlowStop2D::BoundaryFailure&&!rejected.lastAccepted&&rejected.lastRejected,"Outlet accepted backflow or lost rejected trace");
    require(!rejected.iterations.back().accepted,"Rejected backflow labelled accepted");
    std::cout<<"{\"backflowPositiveMeanNegativeEndpoint\":true,\"trials\":"<<rejected.iterations.back().trials<<"}\n";
}
void adapter() {
    // Real ring of eight cells around an omitted solid square, with both
    // physical wall and far-field boundary components.
    std::vector<Polygon2D> polygons;
    for(int j=0;j<3;++j)for(int i=0;i<3;++i)if(i!=1||j!=1)
        polygons.push_back({{{double(i)/3,double(j)/3},{double(i+1)/3,double(j)/3},{double(i+1)/3,double(j+1)/3},{double(i)/3,double(j+1)/3}}});
    auto topology=cartmesh2d::test::fromPolygons(polygons);
    for(auto& e:topology.edges)if(!e.neighbour){const auto a=topology.vertices[e.v0].point,b=topology.vertices[e.v1].point;const double x=(a.x+b.x)/2,y=(a.y+b.y)/2;
        if(x>0&&x<1&&y>0&&y<1)e.patch=BoundaryPatch2D::EmbeddedBoundary;}
    const auto mesh=makeFvMesh2D(topology);FlowControls2D c;c.scenario="external";c.speed=1;
    const auto explicitBoundary=explicitFlowBoundaryPreset2D(mesh,c);
    const auto compatible=compatibleFlowBoundaries2D(mesh,explicitBoundary,c.speed);
    std::size_t symmetry=0,outlet=0;
    for(const auto& b:compatible){symmetry+=b.kind==CompatibleBoundaryKind2D::NormalVelocity;outlet+=b.rejectBackflow;}
    require(symmetry==6&&outlet==3,"External preset lost symmetry/outlet policy");
    auto unsupported=explicitBoundary;
    for(auto& b:unsupported)if(b.kind==FlowBoundaryKind2D::PressureOutlet)b.kind=FlowBoundaryKind2D::PressureOpening;
    bool rejected=false;try{(void)compatibleFlowBoundaries2D(mesh,unsupported,c.speed);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"Dynamic pressure-opening policy silently substituted");
}
void run(double angle,bool symmetry,bool outlet){const Case c(angle,symmetry,outlet);const auto controls=c.controls();const auto result=solveCompatibleIncompressible2D(c.mesh,controls);
    if(!result.converged())throw std::runtime_error(result.reason);const auto& state=*result.lastAccepted;double ue=0,pe=0,normalError=0;const auto gauge=outlet?0:c.pressure(c.mesh.cells.back().centre);
    for(std::size_t t=0;t<c.mesh.cells.size();++t){const P1Local a(c.mesh,t,diameter(c.mesh,t),6);Vec v(2*a.m);for(std::size_t k=0;k<2;++k)for(std::size_t j=0;j<3;++j)v[k*a.m+j]=state.cells[t][3*k+j];
        for(std::size_t f=0;f<c.mesh.cells[t].faces.size();++f)for(std::size_t k=0;k<2;++k)for(std::size_t j=0;j<2;++j)v[k*a.m+3+2*f+j]=state.faces[c.mesh.cells[t].faces[f]][2*k+j];
        for(const auto& q:a.q){const auto phi=a.basis.phi(q.p);const auto theta=a.basis.theta(q.p);Vector2D got{};double p=0;
            for(std::size_t j=0;j<a.m;++j)for(std::size_t k=0;k<6;++k){got.x+=theta[k]*a.potential(k,j)*v[j];got.y+=theta[k]*a.potential(k,j)*v[a.m+j];}for(std::size_t j=0;j<3;++j)p+=phi[j]*state.cells[t][6+j];const auto expected=c.velocity(q.p);
            ue=std::max(ue,std::hypot(got.x-expected.x,got.y-expected.y));pe=std::max(pe,std::abs(p-c.pressure(q.p)+gauge));}
    }
    for(const auto& b:controls.boundaries)if(b.kind==CompatibleBoundaryKind2D::NormalVelocity){const auto& f=c.mesh.faces[b.face];const auto& u=state.faces[b.face];const double length=std::hypot(f.areaVector.x,f.areaVector.y);for(double s:{-.5,.5}){
        const Point2D point{f.centre.x-s*f.areaVector.y,f.centre.y+s*f.areaVector.x};const double actual=((u[0]+s*u[1])*f.areaVector.x+(u[2]+s*u[3])*f.areaVector.y)/length;
        normalError=std::max(normalError,std::abs(actual-(b.normalVelocity?b.normalVelocity(point):0)));}}
    // Same normalized polynomial regression allowance as the existing API
    // tests; not a new geometry or physical-accuracy acceptance threshold.
    require(ue<1e-7&&pe<1e-7&&normalError<1e-12,"Rotated mixed-boundary analytic field/normal trace failed");
    const auto& last=*result.iterations.back().metrics;
    std::cout<<std::setprecision(17)<<"{\"angle\":"<<angle<<",\"symmetry\":"<<(symmetry?"true":"false")<<",\"outlet\":"<<(outlet?"true":"false")<<",\"iterations\":"<<result.iterations.size()<<",\"velocityMax\":"<<ue<<",\"pressureMax\":"<<pe<<",\"normalVelocityMax\":"<<normalError<<",\"cellMomentum\":"<<last.cellMomentum<<",\"faceMomentum\":"<<last.faceMomentum<<",\"divergence\":"<<last.divergence<<"}\n";
}
}
int main()try{backflow();adapter();run(0,false,true);run(.61,false,true);run(.61,false,false);run(0,true,true);run(.61,true,true);
    Case c(.61,false,true);auto options=c.controls();for(auto& b:options.boundaries)if(b.kind==CompatibleBoundaryKind2D::NormalVelocity){b.value=[](Point2D){return Vector2D{};};break;}bool rejected=false;try{(void)solveCompatibleIncompressible2D(c.mesh,options);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"Conflicting mixed boundary input silently ignored");return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
