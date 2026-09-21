#include "cartmesh2d/fv/Euler2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
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
            const auto next=advanceEuler2D(mesh,bc,{},state,{});
            require(next.acousticCourant<=.400000000000001&&next.maximumCellBalanceError<1e-14,"CFL / conservation failure");
            state=next.state;
        }
        same(initial,state);
    }
    for(double speed:{-5.,5.}) {
        const EulerPrimitive2D q{1,speed,.3,1};const auto state=constant(mesh,q);
        same(state,advanceEuler2D(mesh,boundaries(mesh,EulerBoundaryKind2D::Farfield,q),{},state,{}).state);
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
        const auto next=advanceEuler2D(mesh,bc,{},state,{});
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
void shockRotationAndCheckpoint() {
    const auto mesh=rectangle(100,4,1);const double angle=.37;const auto rotatedMesh=rotated(mesh,angle);
    auto bc=boundaries(mesh,EulerBoundaryKind2D::SlipWall);
    for(auto& b:bc)if(std::abs(mesh.faces[b.face].areaVector.x)>0){b.kind=EulerBoundaryKind2D::Transmissive;b.name="ends";}
    EulerState2D state=constant(mesh,{}),rotatedState;
    for(std::size_t i=0;i<state.cells.size();++i)state.cells[i]=eulerConservative2D(mesh.cells[i].centre.x<.5?EulerPrimitive2D{1,0,0,1}:EulerPrimitive2D{.125,0,0,.1});
    rotatedState=state;const auto initial=state;
    while(state.time<.2) {
        EulerStepControls2D controls;controls.maximumStep=.2-state.time;
        const auto next=advanceEuler2D(mesh,bc,{},state,controls);
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
    require(rotationError<1e-12,"Sod rotation covariance failed");
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
            const auto next=advanceEuler2D(mesh,bc,{},state,{});state=next.state;
            require(next.minimumDensity>0&&next.minimumPressure>0,"strong wave positivity failed");
            require(std::abs(next.afterIntegral[3]-energy)<1e-12*energy&&std::abs(next.afterIntegral[0]-mass)<1e-13,"closed blast/expansion lost energy or mass");
        }
    }
}
int main() {
    try {
        const EulerPrimitive2D p{1.225,320,-12,101325};const auto round=eulerPrimitive2D(eulerConservative2D(p));
        require(std::abs(round.pressure-p.pressure)<1e-10&&round.density==p.density,"ideal gas round trip failed");
        rejects([]{(void)eulerConservative2D({-1,0,0,1});});
        rejects([]{(void)eulerPrimitive2D({1,10,0,1});});
        rejects([]{validateIdealGas2D({1,287});});
        restAndFreeStream();periodicConservation();shockRotationAndCheckpoint();strongWaves();
        std::cout<<"Euler conservation, acoustic CFL, boundaries, rotation and checkpoint passed\n";return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
