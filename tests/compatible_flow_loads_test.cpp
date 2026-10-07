#include "cartmesh2d/fv/CompatibleIncompressible2D.hpp"
#include "cartmesh2d/fv/detail/CompatibleFlowElement2D.hpp"
#include "fixtures/PolygonMesh2D.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void check(bool yes,const char* message){if(!yes)throw std::runtime_error(message);}
struct Fixture {
    double angle,U,L;int mode;FvMesh2D mesh;
    Vector2D rotate(Vector2D p)const{const auto c=std::cos(angle),s=std::sin(angle);return {c*p.x-s*p.y,s*p.x+c*p.y};}
    Vector2D unrotate(Vector2D p)const{const auto c=std::cos(angle),s=std::sin(angle);return {c*p.x+s*p.y,-s*p.x+c*p.y};}
    Point2D local(Point2D p)const{const auto q=unrotate({p.x/L,p.y/L});return {q.x,q.y};}
    Vector2D u(Point2D p)const{return mode==3?Vector2D{}:Vector2D{.2+.3*p.x+.7*p.y,-.1-.4*p.x-.3*p.y};}
    double p(Point2D q)const{return 1+2*q.x-3*q.y;}
    Vector2D velocity(Point2D point)const{const auto a=rotate(u(local(point)));return {U*a.x,U*a.y};}
    Vector2D traction(Point2D point,Vector2D normal,bool pseudo=false,double gauge=0)const {
        const auto n=unrotate(normal);const double pressure=p(local(point))-gauge,nu=mode==3?0:.2;
        const auto a=rotate(pseudo?Vector2D{(.3*nu-pressure)*n.x+.7*nu*n.y,-.4*nu*n.x+(-.3*nu-pressure)*n.y}:
                                        Vector2D{(.6*nu-pressure)*n.x+.3*nu*n.y,.3*nu*n.x+(-.6*nu-pressure)*n.y});
        return {U*U*a.x,U*U*a.y};
    }
    Fixture(double a,double speed,double length,int kind):angle(a),U(speed),L(length),mode(kind) {
        std::vector<Polygon2D> polygons;const auto node=[&](int i,int j){const auto q=rotate({L*(i+.25*j)/3,L*j/3});return Point2D{q.x,q.y};};
        for(int j=0;j<3;++j)for(int i=0;i<3;++i)polygons.push_back({{node(i,j),node(i+1,j),node(i+1,j+1),node(i,j+1)}});
        mesh=makeFvMesh2D(cartmesh2d::test::fromPolygons(polygons));
    }
    CompatibleFlowControls2D controls()const {
        CompatibleFlowControls2D c;c.referenceLength=L;c.referenceVelocity=U;c.viscosity=.2*U*L;
        if(mode==3){c.equation=CompatibleEquation2D::Stokes;c.globalization=CompatibleGlobalization2D::Backtracking;}
        if(mode==1||mode==2)c.pressureInverse=CompatiblePressureInverse2D::DiagonalSchur;
        c.acceleration=[this](Point2D point){const auto v=u(local(point));const auto f=rotate({2+.3*v.x+.7*v.y,-3-.4*v.x-.3*v.y});return Vector2D{U*U/L*f.x,U*U/L*f.y};};
        for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour){const auto& face=mesh.faces[id];const auto q=local(face.centre);const double length=std::hypot(face.areaVector.x,face.areaVector.y);const Vector2D n{face.areaVector.x/length,face.areaVector.y/length};CompatibleBoundary2D b;b.face=id;
            if((mode==1||mode==2)&&q.x-.25*q.y>1-1e-12){b.kind=mode==1?CompatibleBoundaryKind2D::Traction:CompatibleBoundaryKind2D::PseudoTraction;b.value=[this,n](Point2D x){return traction(x,n,mode==2);};}
            else if(mode==1&&(q.y<1e-12||q.y>1-1e-12)){b.kind=CompatibleBoundaryKind2D::NormalVelocity;b.normalVelocity=[this,n](Point2D x){return dot(velocity(x),n);};b.tangentialTraction=[this,n](Point2D x){return dot(traction(x,n),Vector2D{-n.y,n.x});};}
            else b.value=[this](Point2D x){return velocity(x);};c.boundaries.push_back(std::move(b));
        }
        return c;
    }
};
void run(double angle,double U,double L,int mode) {
    const Fixture f(angle,U,L,mode);const auto c=f.controls();const auto result=solveCompatibleIncompressible2D(f.mesh,c);
    check(result.converged()&&result.checkpoint,"Analytic load control did not converge");const Point2D origin{.13*L,-.27*L};
    const auto report=evaluateCompatibleFlowLoads2D(f.mesh,c,*result.checkpoint,origin);
    check(report.acceptedIterations==result.checkpoint->acceptedIterations(),"Load source lost accepted iteration identity");
    check(report.absolutePressureReference==(mode==1||mode==2),"Load pressure reference mislabelled");
    const double gauge=report.absolutePressureReference?0:f.p(f.local(f.mesh.cells.back().centre));
    double error=0,pressureError=0,torqueError=0;
    for(const auto& load:report.boundaries){const auto& face=f.mesh.faces[load.face];const auto S=face.areaVector;const double length=std::hypot(S.x,S.y);const Vector2D normal{S.x/length,S.y/length};
        std::array<Vector2D,2> expected{},pressure{};double torque=0;
        for(auto [z,w]:detail::compatible::gauss(6)){const double s=z-.5;const Point2D x{face.centre.x-s*S.y,face.centre.y+s*S.x};const auto traction=f.traction(x,normal,false,gauge);const double p=U*U*(f.p(f.local(x))-gauge);
            torque+=w*length*((x.x-origin.x)*traction.y-(x.y-origin.y)*traction.x);
            for(int k=0;k<2;++k){const double weight=w*length*(k?s:1);expected[k].x+=weight*traction.x;expected[k].y+=weight*traction.y;pressure[k].x-=weight*p*normal.x;pressure[k].y-=weight*p*normal.y;}}
        for(int k=0;k<2;++k){error=std::max(error,std::hypot(expected[k].x-load.tractionMoments[k].x,expected[k].y-load.tractionMoments[k].y)/(U*U*length));pressureError=std::max(pressureError,std::hypot(pressure[k].x-load.pressureMoments[k].x,pressure[k].y-load.pressureMoments[k].y)/(U*U*length));}
        torqueError=std::max(torqueError,std::abs(torque-load.torqueOnFluid)/(U*U*L*L));
    }
    // Existing normalized 1e-7 polynomial regression allowance: larger than
    // solver's 1e-9 stop, not a new spatial/physical accuracy threshold.
    const auto balance=std::hypot(report.momentumImbalance.x,report.momentumImbalance.y)/(U*U*L);
    check(error<1e-7&&pressureError<1e-7&&torqueError<1e-7&&balance<1e-7,"Analytic traction, pressure, torque or momentum balance failed");
    const Point2D shift{.31*L,-.17*L};const auto shifted=evaluateCompatibleFlowLoads2D(f.mesh,c,*result.checkpoint,{origin.x+shift.x,origin.y+shift.y});
    check(std::abs(shifted.boundaryTorqueOnFluid-report.boundaryTorqueOnFluid+shift.x*report.boundaryTraction.y-shift.y*report.boundaryTraction.x)/(U*U*L*L)<1e-12,"Torque origin transformation failed");
    if(mode==1){auto wrong=c;wrong.viscosity*=2;bool rejected=false;try{(void)evaluateCompatibleFlowLoads2D(f.mesh,wrong,*result.checkpoint);}catch(const std::invalid_argument&){rejected=true;}check(rejected,"Wrong physical context accepted for loads");
        wrong=c;wrong.stopRequested=[]{return true;};rejected=false;try{(void)evaluateCompatibleFlowLoads2D(f.mesh,wrong,*result.checkpoint);}catch(const std::runtime_error&){rejected=true;}check(rejected,"Cancelled load evaluation returned a report");
        auto limited=c;limited.maximumIterations=1;const auto partial=solveCompatibleIncompressible2D(f.mesh,limited);check(!partial.converged()&&partial.checkpoint,"Missing incomplete accepted control");
        const auto last=evaluateCompatibleFlowLoads2D(f.mesh,limited,*partial.checkpoint);check(last.acceptedIterations==1,"Incomplete checkpoint load source lost");
        std::ostringstream saved;writeCompatibleFlowCheckpoint2D(saved,*partial.checkpoint);std::istringstream stream(saved.str());const auto restored=readCompatibleFlowCheckpoint2D(stream,f.mesh);
        const auto readback=evaluateCompatibleFlowLoads2D(f.mesh,limited,restored);check(readback.boundaries.size()==last.boundaries.size(),"Reloaded load topology differs");
        for(std::size_t i=0;i<last.boundaries.size();++i)for(int k=0;k<2;++k)check(last.boundaries[i].tractionMoments[k].x==readback.boundaries[i].tractionMoments[k].x&&last.boundaries[i].tractionMoments[k].y==readback.boundaries[i].tractionMoments[k].y,"Checkpoint load moments changed bits");
    }
    std::cout<<std::setprecision(17)<<"{\"mode\":"<<mode<<",\"angle\":"<<angle<<",\"U\":"<<U<<",\"L\":"<<L<<",\"normalizedTractionMomentError\":"<<error<<",\"normalizedPressureMomentError\":"<<pressureError<<",\"normalizedTorqueError\":"<<torqueError<<",\"normalizedMomentumImbalance\":"<<balance<<"}\n";
}
}
int main()try {for(int mode:{0,1,2,3})run(.61,1,1,mode);run(.37,3.25,2.3,1);run(0,1,1,2);return 0;}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
