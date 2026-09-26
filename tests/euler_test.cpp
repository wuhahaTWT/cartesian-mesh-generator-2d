#include "cartmesh2d/fv/Euler2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
EulerStepControls2D activeControls;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
template<class F> void rejects(F fn){bool rejected=false;try{fn();}catch(const std::exception&){rejected=true;}require(rejected,"invalid Euler input accepted");}
FvMesh2D rectangle(int nx=24, int ny=8, double length=4, bool warped=false) {
    TopologyMesh2D t;
    for (int j=0;j<=ny;++j) for (int i=0;i<=nx;++i) {
        double x=length*i/nx,y=double(j)/ny;
        if(warped && i>0 && i<nx && j>0 && j<ny) {
            const double pi=std::acos(-1.),amplitude=.2*std::min(length/nx,1./ny);
            x+=amplitude*std::sin(2*pi*length*i/nx/length)*std::sin(pi*double(j)/ny);
            y+=amplitude*std::sin(pi*length*i/nx/length)*std::sin(2*pi*double(j)/ny);
        }
        t.vertices.push_back({t.vertices.size(),{x,y}});
    }
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> edges;
    for (int j=0;j<ny;++j) for (int i=0;i<nx;++i) {
        TopologyCell2D cell;
        cell.id=t.cells.size();cell.geometryArea=length/nx/ny;
        const auto a=static_cast<std::size_t>(j*(nx+1)+i);
        cell.vertices={a,a+1,a+static_cast<std::size_t>(nx)+2,a+static_cast<std::size_t>(nx)+1};
        if(warped) {
            cell.geometryArea=0;
            for(std::size_t k=0;k<4;++k) {
                const auto p=t.vertices[cell.vertices[k]].point,q=t.vertices[cell.vertices[(k+1)%4]].point;
                cell.geometryArea+=(p.x*q.y-p.y*q.x)/2;
            }
        }
        for (std::size_t k=0;k<4;++k) {
            const auto x=cell.vertices[k],y=cell.vertices[(k+1)%4];
            const auto [it,inserted]=edges.emplace(std::minmax(x,y),t.edges.size());
            if (inserted) t.edges.push_back({it->second,x,y,cell.id,{},BoundaryPatch2D::DomainBoundary});
            else {t.edges[it->second].neighbour=cell.id;t.edges[it->second].patch=BoundaryPatch2D::None;}
            cell.edges.push_back(it->second);
        }
        t.cells.push_back(cell);
    }
    return makeFvMesh2D(t);
}

std::vector<EulerBoundary2D> boundaries(const FvMesh2D& mesh,EulerBoundaryKind2D kind,EulerPrimitive2D reference={}) {
    std::vector<EulerBoundary2D> result;
    for(std::size_t i=0;i<mesh.faces.size();++i)if(!mesh.faces[i].neighbour)
        result.push_back({i,kind,reference,{},"boundary"});
    return result;
}
EulerState2D constant(const FvMesh2D& mesh,EulerPrimitive2D q) {
    return {0,0,std::vector<EulerConservative2D>(mesh.cells.size(),eulerConservative2D(q))};
}
void same(const EulerState2D& a,const EulerState2D& b,double tolerance=1e-12) {
    require(a.cells.size()==b.cells.size(),"state size differs");
    for(std::size_t i=0;i<a.cells.size();++i)for(std::size_t k=0;k<4;++k)
        require(std::abs(a.cells[i][k]-b.cells[i][k])<tolerance,"conserved state differs");
}
void restAndFreeStream() {
    const auto mesh=rectangle(12,6,2,true);
    for(auto kind:{EulerBoundaryKind2D::SlipWall,EulerBoundaryKind2D::Farfield,EulerBoundaryKind2D::Transmissive}) {
        const EulerPrimitive2D q{1.3,kind==EulerBoundaryKind2D::SlipWall?0:.7,kind==EulerBoundaryKind2D::SlipWall?0:-.2,2};
        const auto bc=boundaries(mesh,kind,q);const auto initial=constant(mesh,q);auto state=initial;
        for(int step=0;step<40;++step) {
            const auto next=advanceEuler2D(mesh,bc,{},state,activeControls);
            require(next.acousticCourant<=.400000000000001&&next.maximumCellBalanceError<1e-14,"CFL / conservation failure");
            state=next.state;
        }
        same(initial,state);
    }
    for(double speed:{-5.,5.}) {
        const EulerPrimitive2D q{1,speed,.3,1};const auto state=constant(mesh,q);
        same(state,advanceEuler2D(mesh,boundaries(mesh,EulerBoundaryKind2D::Farfield,q),{},state,activeControls).state);
    }
}
std::vector<EulerBoundary2D> periodic(const FvMesh2D& mesh) {
    auto bc=boundaries(mesh,EulerBoundaryKind2D::Periodic);
    for(auto& b:bc) {
        const auto& f=mesh.faces[b.face];const bool x=std::abs(f.areaVector.x)>std::abs(f.areaVector.y);b.name=x?"periodic-x":"periodic-y";
        for(const auto& other:bc) {
            const auto& g=mesh.faces[other.face];
            if(other.face!=b.face&&std::hypot(f.areaVector.x+g.areaVector.x,f.areaVector.y+g.areaVector.y)<1e-12&&
               std::abs(x?f.centre.y-g.centre.y:f.centre.x-g.centre.x)<1e-12){b.partner=other.face;break;}
        }
    }
    return bc;
}
void periodicConservation() {
    const auto mesh=rectangle(20,10,2,true);const auto bc=periodic(mesh);auto state=constant(mesh,{1,.3,-.2,1});
    const double pi=std::acos(-1.);
    for(std::size_t i=0;i<state.cells.size();++i) {
        const auto c=mesh.cells[i].centre;
        state.cells[i]=eulerConservative2D({1+.2*std::sin(pi*c.x)*std::cos(2*pi*c.y),.3+.1*std::cos(pi*c.x),-.2,1});
    }
    EulerConservative2D before{},after{};
    for(std::size_t i=0;i<state.cells.size();++i)for(std::size_t k=0;k<4;++k)before[k]+=mesh.cells[i].area*state.cells[i][k];
    for(int i=0;i<80;++i) {
        const auto next=advanceEuler2D(mesh,bc,{},state,activeControls);
        for(const auto& b:bc)for(std::size_t k=0;k<4;++k)require(next.faceFlux[b.face][k]==-next.faceFlux[*b.partner][k],"periodic action-reaction differs");
        state=next.state;
    }
    for(std::size_t i=0;i<state.cells.size();++i)for(std::size_t k=0;k<4;++k)after[k]+=mesh.cells[i].area*state.cells[i][k];
    for(std::size_t k=0;k<4;++k)require(std::abs(after[k]-before[k])<2e-14,"periodic total conservation failed");
    auto bad=bc;bad[0].partner=bad[0].face;rejects([&]{validateEulerBoundaries2D(mesh,bad);});
    bad=bc;bad[0].name="unpaired";rejects([&]{validateEulerBoundaries2D(mesh,bad);});
}
FvMesh2D rotated(FvMesh2D mesh,double angle) {
    const auto vector=[&](Vector2D p){return Vector2D{p.x*std::cos(angle)-p.y*std::sin(angle),p.x*std::sin(angle)+p.y*std::cos(angle)};};
    for(auto& c:mesh.cells){const auto r=vector({c.centre.x,c.centre.y});c.centre={r.x,r.y};}
    for(auto& f:mesh.faces){const auto r=vector({f.centre.x,f.centre.y});f.centre={r.x,r.y};f.areaVector=vector(f.areaVector);f.correction=vector(f.correction);}
    return mesh;
}
void wallAndUnitScaling() {
    const auto mesh=rotated(rectangle(10,6,2,true),.37);
    const auto bc=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    auto state=constant(mesh,{1,0,0,1}),scaled=state;
    const double rhoScale=1.225,pressureScale=101325,velocityScale=std::sqrt(pressureScale/rhoScale);
    const std::array<double,4> units{rhoScale,rhoScale*velocityScale,rhoScale*velocityScale,pressureScale};
    for(std::size_t i=0;i<state.cells.size();++i) {
        const auto c=mesh.cells[i].centre;
        state.cells[i]=eulerConservative2D({1+.1*std::sin(c.x),.4+.05*std::sin(c.y),-.2,1+.1*std::cos(c.y)});
        for(std::size_t k=0;k<4;++k)scaled.cells[i][k]=state.cells[i][k]*units[k];
    }
    double maxError=0;
    for(unsigned step=0;step<40;++step) {
        auto controls=activeControls;controls.maximumStep=.002;
        const auto next=advanceEuler2D(mesh,bc,{},state,controls);
        controls.maximumStep/=velocityScale;
        const auto si=advanceEuler2D(mesh,bc,{},scaled,controls);
        for(const auto* result:{&next,&si}) {
            require(result->boundaryFlux[0]==0&&result->boundaryFlux[3]==0,"closed wall mass/energy leak");
            for(const auto& b:bc) {
                const auto& f=result->faceFlux[b.face];const auto n=mesh.faces[b.face].areaVector;
                require(f[0]==0&&f[3]==0,"wall flux symmetry not exact");
                require(std::abs(f[1]*n.y-f[2]*n.x)<=8*std::numeric_limits<double>::epsilon()*std::hypot(n.x,n.y)*std::hypot(f[1],f[2]),"wall has tangential traction");
            }
        }
        for(std::size_t i=0;i<state.cells.size();++i)for(std::size_t k=0;k<4;++k)
            maxError=std::max(maxError,std::abs(si.state.cells[i][k]/units[k]-next.state.cells[i][k])/(1+std::abs(next.state.cells[i][k])));
        state=next.state;scaled=si.state;
    }
    // Roundoff-only equivalence after 40 accepted steps, normalized back to the
    // original units; this is not a physical solution-accuracy requirement.
    require(maxError<256*std::numeric_limits<double>::epsilon()*state.steps,"Euler depends on chosen physical units");
    std::cout<<"rotated wall SI scaling: max normalized field difference "<<maxError<<'\n';
}
void shockRotationAndCheckpoint() {
    const auto mesh=rectangle(100,4,1);const double angle=.37;const auto rotatedMesh=rotated(mesh,angle);
    auto bc=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    for(auto& b:bc)if(std::abs(mesh.faces[b.face].areaVector.x)>0){b.kind=EulerBoundaryKind2D::Transmissive;b.name="ends";}
    EulerState2D state=constant(mesh,{}),rotatedState;
    for(std::size_t i=0;i<state.cells.size();++i)state.cells[i]=eulerConservative2D(mesh.cells[i].centre.x<.5?EulerPrimitive2D{1,0,0,1}:EulerPrimitive2D{.125,0,0,.1});
    rotatedState=state;const auto initial=state;double localRotationError=0;
    while(state.time<.2) {
        auto controls=activeControls;controls.maximumStep=.2-state.time;
        const auto next=advanceEuler2D(mesh,bc,{},state,controls);
        if(activeControls.order==2&&activeControls.fluxScheme==EulerFluxScheme2D::Hllc) {
            auto exactRotated=state;
            for(auto& q:exactRotated.cells){const auto x=q[1],y=q[2];q[1]=x*std::cos(angle)-y*std::sin(angle);q[2]=x*std::sin(angle)+y*std::cos(angle);}
            auto cap=controls;cap.maximumStep=.9*next.step;
            const auto referenceStep=advanceEuler2D(mesh,bc,{},state,cap),rotatedStep=advanceEuler2D(rotatedMesh,bc,{},exactRotated,cap);
            double local=0;
            for(std::size_t i=0;i<state.cells.size();++i) {
                const auto& q=referenceStep.state.cells[i];const auto& r=rotatedStep.state.cells[i];
                const EulerConservative2D expected{q[0],q[1]*std::cos(angle)-q[2]*std::sin(angle),q[1]*std::sin(angle)+q[2]*std::cos(angle),q[3]};
                for(std::size_t k=0;k<4;++k)local=std::max(local,std::abs(expected[k]-r[k]));
            }
            if(local>localRotationError){localRotationError=local;if(local>1e-10)std::cerr<<"local rotation: step "<<state.steps<<" error "<<local<<" dt "<<referenceStep.step<<' '<<rotatedStep.step<<'\n';}
        }
        controls.maximumStep=next.step;
        rotatedState=advanceEuler2D(rotatedMesh,bc,{},rotatedState,controls).state;state=next.state;
        require(next.minimumDensity>0&&next.minimumPressure>0,"Sod positivity failed");
        if(state.steps==50) {
            std::ostringstream stream;writeEulerCheckpoint2D(stream,mesh,bc,{},state,"Sod");
            std::istringstream input(stream.str());const auto restored=readEulerCheckpoint2D(input,mesh,bc,{},"Sod");
            same(state,restored,1e-15);require(state.time==restored.time&&state.steps==restored.steps,"restart counter differs");
            same(advanceEuler2D(mesh,bc,{},state,controls).state,advanceEuler2D(mesh,bc,{},restored,controls).state,1e-15);
            rejects([&]{std::istringstream in(stream.str());(void)readEulerCheckpoint2D(in,mesh,bc,{1.3,287.05},"Sod");});
            rejects([&]{std::istringstream in(stream.str()+"extra");(void)readEulerCheckpoint2D(in,mesh,bc,{},"Sod");});
            rejects([&]{std::istringstream in(stream.str());(void)readEulerCheckpoint2D(in,rotatedMesh,bc,{},"Sod");});
        }
    }
    double rotationError=0;
    for(std::size_t i=0;i<state.cells.size();++i) {
        const auto a=eulerPrimitive2D(state.cells[i]),b=eulerPrimitive2D(rotatedState.cells[i]);
        rotationError=std::max(rotationError,std::max({std::abs(a.density-b.density),std::abs(a.pressure-b.pressure),
            std::abs(a.u*std::cos(angle)-a.v*std::sin(angle)-b.u),std::abs(a.u*std::sin(angle)+a.v*std::cos(angle)-b.v)}));
    }
    std::cout<<"rotation diagnostic "<<std::setprecision(17)<<rotationError<<", max single-step "<<localRotationError<<", times "<<state.time<<' '<<rotatedState.time<<std::endl;
    // Keep the original first-order tolerance. The added nonlinear, two-stage
    // path uses a roundoff accumulation budget proportional to accepted steps,
    // in these nondimensional Sod variables; this is not a physical error gate.
    const double rotationBudget=activeControls.order==1?1e-12:256*std::numeric_limits<double>::epsilon()*static_cast<double>(state.steps);
    require(rotationError<rotationBudget,"Sod rotation covariance failed");
    require(localRotationError<1e-12,"single-step rotation covariance failed");
    std::cout<<"Sod 400 cells: "<<state.steps<<" steps, rotated state max difference "<<rotationError<<'\n';
    auto controls=EulerStepControls2D{};controls.minimumStep=.1;controls.maximumStep=1;
    rejects([&]{(void)advanceEuler2D(mesh,bc,{},initial,controls);});
    require(initial.time==0&&initial.steps==0,"failed step changed input");
}
void strongWaves() {
    const auto mesh=rectangle(80,2,1);
    const auto bc=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    for(bool expansion:{false,true}) {
        auto state=constant(mesh,{});
        for(std::size_t i=0;i<state.cells.size();++i) {
            const bool left=mesh.cells[i].centre.x<.5;
            state.cells[i]=eulerConservative2D(expansion?EulerPrimitive2D{1,left?-2.:2.,0,.4}:
                EulerPrimitive2D{1,0,0,left?1000.:.01});
        }
        double energy=0,mass=0;
        for(std::size_t i=0;i<state.cells.size();++i){energy+=mesh.cells[i].area*state.cells[i][3];mass+=mesh.cells[i].area*state.cells[i][0];}
        for(int j=0;j<200;++j) {
            const auto next=advanceEuler2D(mesh,bc,{},state,activeControls);state=next.state;
            require(next.minimumDensity>0&&next.minimumPressure>0,"strong wave positivity failed");
            require(std::abs(next.afterIntegral[3]-energy)<1e-12*energy&&std::abs(next.afterIntegral[0]-mass)<1e-13,"closed blast/expansion lost energy or mass");
        }
    }
}
void faceFluxIdentities() {
    const auto state=eulerConservative2D({1.2,.4,-.3,2});
    for(auto scheme:{EulerFluxScheme2D::Rusanov,EulerFluxScheme2D::Hllc}) {
        const auto f=eulerFaceFlux2D(state,state,{.6,.8},{},scheme);const double vn=.4*.6-.3*.8;
        const EulerConservative2D exact{state[0]*vn,state[1]*vn+1.2,state[2]*vn+1.6,(state[3]+2)*vn};
        for(std::size_t k=0;k<4;++k)require(std::abs(f.integratedFlux[k]-exact[k])<3e-15,"consistent physical flux failed");
        for(int i=0;i<40;++i) {
            const auto l=eulerConservative2D({.1+.07*i,-3+.17*i,1-.09*i,.2+.13*i});
            const auto r=eulerConservative2D({3-.03*i,2-.19*i,-.7+.04*i,4-.09*i});
            const auto a=eulerFaceFlux2D(l,r,{.3,.4},{},scheme),b=eulerFaceFlux2D(r,l,{-.3,-.4},{},scheme);
            for(std::size_t k=0;k<4;++k)require(std::abs(a.integratedFlux[k]+b.integratedFlux[k])<1e-12*(1+std::abs(a.integratedFlux[k])),"flux owner reversal failed");
            require(a.hllcFallback==b.hllcFallback&&std::abs(a.waveSpeed-b.waveSpeed)<=8*std::numeric_limits<double>::epsilon()*a.waveSpeed,"wave/fallback owner reversal failed");
        }
    }
    for(double velocity:{-5.,0.,5.}) {
        const auto l=eulerConservative2D({1,velocity,2,1}),r=eulerConservative2D({3,velocity,-1,1});
        const auto f=eulerFaceFlux2D(l,r,{1,0},{},EulerFluxScheme2D::Hllc);const auto& q=velocity>=0?l:r;
        const EulerConservative2D exact{q[0]*velocity,q[1]*velocity+1,q[2]*velocity,(q[3]+1)*velocity};
        for(std::size_t k=0;k<4;++k)require(std::abs(f.integratedFlux[k]-exact[k])<2e-13,"HLLC contact/shear resolution failed");
        require(!f.hllcFallback,"admissible contact unexpectedly fell back");
    }
    const auto expansion=eulerFaceFlux2D(eulerConservative2D({1,-2,0,.4}),eulerConservative2D({1,2,0,.4}),{1,0},{},EulerFluxScheme2D::Hllc);
    require(expansion.hllcFallback,"non-positive HLLC star pressure did not report fallback");
    rejects([&]{(void)eulerFaceFlux2D(state,state,{0,0},{},EulerFluxScheme2D::Hllc);});
    rejects([&]{(void)eulerFaceFlux2D(state,state,{1,0},{},static_cast<EulerFluxScheme2D>(99));});
}
void stationaryContact() {
    const auto mesh=rectangle(40,6,2);const auto bc=periodic(mesh);auto state=constant(mesh,{});
    for(std::size_t i=0;i<state.cells.size();++i)state.cells[i]=eulerConservative2D({mesh.cells[i].centre.x<1?1.:3.,0,0,1});
    const auto initial=state;auto controls=activeControls;controls.fluxScheme=EulerFluxScheme2D::Hllc;
    for(int j=0;j<50;++j){const auto step=advanceEuler2D(mesh,bc,{},state,controls);require(step.hllcFallbackEvaluations==0,"contact fallback");state=step.state;}
    same(initial,state,2e-14);
}
void perturbedNormalShock() {
    const double gamma=1.4,amplitude=1e-5;
    for(double mach:{3.,10.}) {
        const auto mesh=rectangle(48,12,2);auto bc=periodic(mesh);
        const EulerPrimitive2D left{1,mach*std::sqrt(gamma),0,1};
        const double ratio=(gamma+1)*mach*mach/((gamma-1)*mach*mach+2);
        const EulerPrimitive2D right{ratio,left.u/ratio,0,1+2*gamma/(gamma+1)*(mach*mach-1)};
        for(auto& b:bc)if(b.name=="periodic-x") {
            const bool inlet=mesh.faces[b.face].areaVector.x<0;
            b.kind=EulerBoundaryKind2D::Farfield;b.partner.reset();b.reference=inlet?left:right;b.name=inlet?"inlet":"outlet";
        }
        auto state=constant(mesh,left);
        for(std::size_t i=0;i<state.cells.size();++i) {
            const auto c=mesh.cells[i].centre;auto q=c.x<1?left:right;
            if(std::abs(c.x-1)<.1)q.density*=1+amplitude*((i/48)%2==0?1:-1);
            state.cells[i]=eulerConservative2D(q);
        }
        EulerStepControls2D controls;controls.fluxScheme=EulerFluxScheme2D::Hllc;controls.order=2;
        double maximumVariation=0;
        const double end=.3/left.u; // fixed convected distance at either Mach number
        while(state.time<end) {
            controls.maximumStep=end-state.time;const auto next=advanceEuler2D(mesh,bc,{},state,controls);state=next.state;
            require(next.minimumDensity>0&&next.minimumPressure>0&&next.minimumContactRestoration<.9,"shock guard/positivity not active");
            for(int i=0;i<48;++i) {
                double low=INFINITY,high=-INFINITY;
                for(int j=0;j<12;++j){const double rho=state.cells[static_cast<std::size_t>(j*48+i)][0];low=std::min(low,rho);high=std::max(high,rho);}
                maximumVariation=std::max(maximumVariation,(high-low)/ratio);
            }
        }
        std::cout<<"perturbed normal shock: Mach="<<mach<<", maximum transverse density span/postshock density="<<maximumVariation<<'\n';
        // Input density sawtooth spans 2e-5 of the postshock density. Bound its
        // growth by one decade over the stated window; this is not a universal
        // shock-stability claim and does not tune the product shock sensor.
        require(maximumVariation<20*amplitude,"shock transverse perturbation grew by more than one decade");
    }
}
void smoothAccuracy() {
    const double pi=std::acos(-1.),end=.3;std::array<double,2> last{};
    for(unsigned order:{1u,2u}) {
        double prev=0;
        for(int nx:{20,40,80}) {
            const int ny=nx/2;const auto mesh=rectangle(nx,ny,2);const auto bc=periodic(mesh);
            auto state=constant(mesh,{});const double sinc=std::sin(pi/nx)/(pi/nx)*std::sin(pi/ny)/(pi/ny);
            const auto density=[&](Point2D c,double t){return 1+.2*sinc*std::sin(pi*(c.x-.3*t))*std::cos(2*pi*(c.y+.2*t));};
            for(std::size_t i=0;i<state.cells.size();++i)state.cells[i]=eulerConservative2D({density(mesh.cells[i].centre,0),.3,-.2,1});
            EulerStepControls2D controls;controls.fluxScheme=EulerFluxScheme2D::Hllc;controls.order=order;
            while(state.time<end){controls.maximumStep=end-state.time;state=advanceEuler2D(mesh,bc,{},state,controls).state;}
            double error=0;
            for(std::size_t i=0;i<state.cells.size();++i)error+=mesh.cells[i].area*std::abs(state.cells[i][0]-density(mesh.cells[i].centre,end));
            error/=2;
            std::cout<<"smooth entropy wave: order="<<order<<", nx="<<nx<<", area-weighted density L1="<<error;
            if(prev>0)std::cout<<", observed order="<<std::log2(prev/error);
            std::cout<<'\n';
            if(nx==80) {
                // A smooth periodic transport case tests the implemented order.
                // A BJ limiter reduces accuracy near extrema; 1.4 is a bounded
                // development regression floor, not a general CFD accuracy gate.
                require(prev/error>(order==2?std::pow(2.,1.4):1.5),"smooth grid refinement lost expected order");
                last[order-1]=error;
            }
            prev=error;
        }
    }
    require(last[1]<.3*last[0],"second-order smooth transport did not improve resolution");
}
int main() {
    try {
        const EulerPrimitive2D p{1.225,320,-12,101325};const auto round=eulerPrimitive2D(eulerConservative2D(p));
        require(std::abs(round.pressure-p.pressure)<1e-10&&round.density==p.density,"ideal gas round trip failed");
        rejects([]{(void)eulerConservative2D({-1,0,0,1});});
        rejects([]{(void)eulerPrimitive2D({1,10,0,1});});
        rejects([]{validateIdealGas2D({1,287});});
        faceFluxIdentities();
        for(auto scheme:{EulerFluxScheme2D::Rusanov,EulerFluxScheme2D::Hllc})for(unsigned order:{1u,2u}) {
            activeControls.fluxScheme=scheme;activeControls.order=order;
            std::cout<<"scheme="<<(scheme==EulerFluxScheme2D::Hllc?"hllc":"rusanov")<<", order="<<order<<std::endl;
            restAndFreeStream();periodicConservation();wallAndUnitScaling();shockRotationAndCheckpoint();strongWaves();stationaryContact();
        }
        perturbedNormalShock();smoothAccuracy();
        std::cout<<"Euler conservation, acoustic CFL, boundaries, rotation and checkpoint passed\n";return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
