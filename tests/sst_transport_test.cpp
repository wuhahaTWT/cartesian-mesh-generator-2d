#include "cartmesh2d/fv/Sst2003m2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void check(bool ok,const char* message) { if (!ok) throw std::runtime_error(message); }
void near(double a,double b,double tolerance,const char* message) { check(std::abs(a-b)<=tolerance,message); }
template<class F> void rejects(F&& f,const std::string& reason) {
    try { f(); } catch(const std::exception& e) {
        check(std::string(e.what()).find(reason)!=std::string::npos,"unexpected rejection reason"); return;
    }
    throw std::runtime_error("invalid SST input accepted");
}
FvMesh2D oneCell() {
    TopologyMesh2D t;
    t.vertices={{0,{0,0}},{1,{1,0}},{2,{1,1}},{3,{0,1}}};
    for(std::size_t i=0;i<4;++i) t.edges.push_back({i,i,(i+1)%4,0,{},BoundaryPatch2D::DomainBoundary});
    TopologyCell2D c; c.id=0; c.vertices={0,1,2,3};c.edges={0,1,2,3};c.geometryArea=1;
    t.cells.push_back(c);return makeFvMesh2D(t);
}
void coefficients() {
    Sst2003mPoint2D p{.02,4,1e-5,1e-4,0,{0,0},{0,0}};
    auto c=evaluateSst2003m2D(p);
    near(c.f1,1,0,"inner F1");near(c.f2,1,0,"inner F2");
    near(c.beta,.075,0,"inner beta");near(c.gamma,5./9.,0,"2003 gamma");
    near(c.turbulentViscosity,.005,1e-17,"unlimited turbulent viscosity");
    near(c.lossRateK,.36,1e-15,"k decay rate");near(c.lossRateOmega,.3,1e-15,"omega decay rate");
    p.strainMagnitude=100;
    c=evaluateSst2003m2D(p);
    near(c.turbulentViscosity,.000062,1e-18,"strain limited viscosity");
    near(c.productionK,.072,1e-15,"factor ten production limiter");
    near(c.productionOmega,(5./9.)*.072/.000062,1e-10,"omega ALSO uses production limiter");
    // Pure rigid rotation has zero strain despite nonzero vorticity: S=0 must
    // not trigger the old SST vorticity-based limiter.
    p.strainMagnitude=0;p.wallDistance=100;p.gradientK={1,0};p.gradientOmega={.1,0};
    c=evaluateSst2003m2D(p);
    check(c.f1<1e-8,"outer blend");near(c.beta,.0828,1e-10,"outer beta");
    near(c.crossDiffusion,2*(1-c.f1)*.856*.1/4,1e-15,"positive cross diffusion");
    near(c.sourceOmega,c.crossDiffusion,1e-15,"positive cross remains explicit");
    near(c.lossRateOmega,c.beta*4,1e-15,"positive cross not double counted");
    p.gradientOmega={-.1,0};c=evaluateSst2003m2D(p);
    check(c.crossDiffusion<0,"negative cross retained");
    near(c.sourceOmega,0,0,"negative cross not explicit");
    near(c.sourceOmega-c.lossRateOmega*p.omega,
         c.productionOmega-c.beta*p.omega*p.omega+c.crossDiffusion,1e-15,"split preserves nonlinear omega source");
    // Intermediate blend, not just saturated 0/1 branches, from independent
    // evaluation of published functions with d=1, k=.02, omega=4, grad dot=.1.
    p.wallDistance=1;p.gradientOmega={.1,0};c=evaluateSst2003m2D(p);
    const double a=std::min(std::max(std::sqrt(.02)/(.09*4),500*1e-5/4),4*.856*.02/(2*.856*.1/4));
    near(c.f1,std::tanh(std::pow(a,4)),1e-15,"intermediate F1");
    near(c.f2,std::tanh(std::pow(2*std::sqrt(.02)/(.09*4),2)),1e-15,"intermediate F2");
    p.k=0;p.strainMagnitude=3;c=evaluateSst2003m2D(p);
    near(c.turbulentViscosity,0,0,"zero k no artificial floor");
    near(c.productionK,0,0,"zero k production");check(std::isfinite(c.productionOmega),"zero k continuous omega production");
    p.k=-1;rejects([&]{(void)evaluateSst2003m2D(p);},"nonnegative");p.k=.02;
    p.omega=0;rejects([&]{(void)evaluateSst2003m2D(p);},"positive");p.omega=4;
    p.wallDistance=0;rejects([&]{(void)evaluateSst2003m2D(p);},"positive");p.wallDistance=1;
    p.gradientK.x=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{(void)evaluateSst2003m2D(p);},"range");
}
FrozenSst2003mProblem2D problem(const FvMesh2D& m) {
    FrozenSst2003mProblem2D p;const auto n=m.cells.size(),nf=m.faces.size();
    p.k.assign(n,.02);p.omega.assign(n,4);p.wallDistance.assign(n,1e-4);
    p.strainMagnitude.assign(n,0);p.gradientK.resize(n);p.gradientOmega.resize(n);
    p.volumeFlux.assign(nf,0);
    p.boundaryK.assign(nf,{ScalarBoundaryKind2D::DiffusiveFlux,0,{}});p.boundaryOmega=p.boundaryK;
    return p;
}
void invalidTransport() {
    const auto m=oneCell();auto p=problem(m);const auto oldK=p.k,oldW=p.omega;
    rejects([&]{(void)solveFrozenSst2003mTransport2D(m,p,{},oldK,{},.1);},"both previous");
    p.wallDistance.clear();rejects([&]{(void)solveFrozenSst2003mTransport2D(m,p);},"dimensions");p=problem(m);
    p.boundaryOmega[0]={ScalarBoundaryKind2D::Value,0,{}};
    rejects([&]{(void)solveFrozenSst2003mTransport2D(m,p);},"boundary value");p=problem(m);
    ScalarTransportControls2D c;c.maxCorrections=1;
    rejects([&]{(void)solveFrozenSst2003mTransport2D(m,p,c,oldK,oldW,.1);},"did not converge");
    p.boundaryK[0]={ScalarBoundaryKind2D::DiffusiveFlux,1e4,{}};
    rejects([&]{(void)solveFrozenSst2003mTransport2D(m,p,{},oldK,oldW,.1);},"negative k after");
}
void frozenSourceSplit() {
    const auto m=oneCell();
    for(double sign:{-1.,1.}) {
        auto p=problem(m);p.wallDistance[0]=100;p.strainMagnitude[0]=3;
        p.gradientK[0]={1,0};p.gradientOmega[0]={sign*.1,0};
        const auto oldK=p.k,oldW=p.omega;const double dt=.2;
        const auto r=solveFrozenSst2003mTransport2D(m,p,{},oldK,oldW,dt);
        const auto c=r.coefficients[0];
        near(r.k.values[0],(oldK[0]+dt*c.sourceK)/(1+dt*c.lossRateK),1e-9,"frozen k source/sink solve");
        near(r.omega.values[0],(oldW[0]+dt*c.sourceOmega)/(1+dt*c.lossRateOmega),1e-9,"frozen omega source/sink solve");
        near(r.k.sinkIntegrals[0],c.lossRateK*r.k.values[0],1e-12,"reported k loss");
        near(r.omega.sinkIntegrals[0],c.lossRateOmega*r.omega.values[0],1e-12,"reported omega loss");
    }
}
void reconstructedBoundaryDerivatives() {
    const auto m=oneCell();auto p=problem(m);
    // Affine k=.02+.004*x; omega=4+.4*x, values on vertical sides,
    // zero normal derivative on horizontal sides.
    p.k[0]=.022;p.omega[0]=4.2;
    std::vector<SstVelocityBoundary2D> bc(m.faces.size());
    for(std::size_t id=0;id<m.faces.size();++id) {
        const auto& f=m.faces[id];const bool fixed=std::abs(f.areaVector.x)>0;
        if(fixed) {
            p.boundaryK[id]={ScalarBoundaryKind2D::Value,.02+.004*f.centre.x,{}};
            p.boundaryOmega[id]={ScalarBoundaryKind2D::Value,4+.4*f.centre.x,{}};
        }
        bc[id]={{f.centre.x,2},fixed,false};
    }
    const auto g=reconstructSst2003mGradients2D(m,p,{{.5,2}},bc);
    near(g.k[0].x,.004,1e-15,"mixed boundary k derivative");near(g.k[0].y,0,0,"zero flux k derivative");
    near(g.omega[0].x,.4,1e-14,"mixed omega derivative");
    near(g.strainMagnitude[0],std::sqrt(2.),1e-14,"normal strain includes diagonal tensor terms");
    rejects([&]{(void)reconstructSst2003mGradients2D(m,p,{},bc);},"dimensions");
    p.k[0]=-1;rejects([&]{(void)reconstructSst2003mGradients2D(m,p,{{.5,2}},bc);},"nonnegative");
}
// Homogeneous no-shear decay is an ODE limit of the TWO PDEs. Wall distance is
// prescribed to saturate F1=1, not claimed to represent a physical wall domain.
// Converge the frozen iterations against the ORIGINAL nonlinear cell residual.
void decay(const FvMesh2D& mesh,double dt,int steps,const std::string& prefix) {
    auto p=problem(mesh);ScalarTransportControls2D control;
    control.relativeTolerance=1e-12;control.absoluteTolerance=1e-14;control.cellTolerance=1e-12;
    std::ofstream history,fields,faces;
    if(!prefix.empty()) {
        history.open(prefix+".history.csv");fields.open(prefix+".cells.csv");faces.open(prefix+".faces.csv");
        check(bool(history)&&bool(fields)&&bool(faces),"output open failed");
        history<<std::setprecision(17)<<"step,time,dt,outerIterations,k,omega,maxNonlinearRateResidual\n";
        fields<<std::setprecision(17)<<"cell,x,y,area,k,omega,previousK,previousOmega\n";
        faces<<std::setprecision(17)<<"face,kFlux,omegaFlux\n";
    }
    double discreteK=.02,discreteW=4;FrozenSst2003mResult2D r;
    for(int step=1;step<=steps;++step) {
        const auto oldK=p.k,oldW=p.omega;double residual=0;int outer=0;
        for(outer=1;outer<=80;++outer) {
            r=solveFrozenSst2003mTransport2D(mesh,p,control,oldK,oldW,dt);
            p.k=r.k.values;p.omega=r.omega.values;
            std::vector<double> rk(mesh.cells.size()),rw(rk.size());
            for(std::size_t i=0;i<rk.size();++i) {
                rk[i]=mesh.cells[i].area*((p.k[i]-oldK[i])/dt+.09*p.omega[i]*p.k[i]);
                rw[i]=mesh.cells[i].area*((p.omega[i]-oldW[i])/dt+.075*p.omega[i]*p.omega[i]);
            }
            for(std::size_t id=0;id<mesh.faces.size();++id) {
                const auto& f=mesh.faces[id];const double qk=r.k.advectiveFlux[id]+r.k.diffusiveFlux[id];
                const double qw=r.omega.advectiveFlux[id]+r.omega.diffusiveFlux[id];
                rk[f.owner]+=qk;rw[f.owner]+=qw;
                if(f.neighbour){rk[*f.neighbour]-=qk;rw[*f.neighbour]-=qw;}
            }
            residual=0;
            for(std::size_t i=0;i<rk.size();++i)
                residual=std::max({residual,std::abs(rk[i])/mesh.cells[i].area,std::abs(rw[i])/mesh.cells[i].area});
            if(residual<1e-9)break;
        }
        check(outer<=80,"nonlinear SST decay did not converge");
        discreteW=2*discreteW/(1+std::sqrt(1+4*dt*.075*discreteW));
        discreteK/=1+dt*.09*discreteW;
        for(std::size_t i=0;i<mesh.cells.size();++i) {
            near(p.omega[i],discreteW,1e-8,"omega backward Euler analytic root");
            near(p.k[i],discreteK,1e-9,"k backward Euler analytic root");
        }
        if(history.is_open())history<<step<<','<<step*dt<<','<<dt<<','<<outer<<','<<p.k[0]<<','<<p.omega[0]<<','<<residual<<'\n';
        if(step==steps&&fields.is_open()) {
            for(std::size_t i=0;i<mesh.cells.size();++i) {
                const auto& c=mesh.cells[i];fields<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','
                    <<p.k[i]<<','<<p.omega[i]<<','<<oldK[i]<<','<<oldW[i]<<'\n';
            }
            for(std::size_t id=0;id<mesh.faces.size();++id)
                faces<<id<<','<<r.k.advectiveFlux[id]+r.k.diffusiveFlux[id]<<','
                    <<r.omega.advectiveFlux[id]+r.omega.diffusiveFlux[id]<<'\n';
        }
    }
    std::cout<<"cells="<<mesh.cells.size()<<" dt="<<dt<<" steps="<<steps
        <<" k="<<p.k[0]<<" omega="<<p.omega[0]<<'\n';
}
}
int main(int argc,char** argv) {
    try {
        coefficients();invalidTransport();frozenSourceSplit();reconstructedBoundaryDerivatives();
        if(argc==1) {decay(oneCell(),.2,4,"");return 0;}
        check(argc==5,"usage: sst_transport_tests [mesh.cm2d outputPrefix dt steps]");
        const auto read=readCm2dTopology(argv[1]);check(read.valid(),"invalid real mesh");
        std::size_t used=0;const std::string dtArg=argv[3],stepsArg=argv[4];
        const double dt=std::stod(dtArg,&used);check(used==dtArg.size()&&std::isfinite(dt)&&dt>0,"invalid dt");
        const int steps=std::stoi(stepsArg,&used);check(used==stepsArg.size()&&steps>0&&steps<=100,"invalid steps");
        decay(makeFvMesh2D(read.topology),dt,steps,argv[2]);return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
