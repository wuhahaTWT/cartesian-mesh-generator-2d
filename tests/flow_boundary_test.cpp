#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/fv/FlowCheckpoint2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }

FvMesh2D rectangle(int nx=24, int ny=8, double length=4) {
    TopologyMesh2D t;
    for (int j=0;j<=ny;++j) for (int i=0;i<=nx;++i)
        t.vertices.push_back({t.vertices.size(),{length*i/nx, double(j)/ny}});
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> edges;
    for (int j=0;j<ny;++j) for (int i=0;i<nx;++i) {
        TopologyCell2D cell;
        cell.id=t.cells.size();cell.geometryArea=length/nx/ny;
        const auto a=static_cast<std::size_t>(j*(nx+1)+i);
        cell.vertices={a,a+1,a+static_cast<std::size_t>(nx)+2,a+static_cast<std::size_t>(nx)+1};
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

FlowControls2D conditions(const FvMesh2D& mesh, bool cavity=false) {
    FlowControls2D c;c.scenario="custom";c.nu=.2;c.tolerance=1e-10;c.maxIterations=3000;
    c.pressurePreconditioner=PressurePreconditioner2D::Aggregation;
    for (std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& face=mesh.faces[id];if(face.neighbour)continue;
        FlowBoundaryCondition2D b{id,FlowBoundaryKind2D::Wall,{},0,"wall"};
        if(cavity) {
            if(face.areaVector.y>0) b={id,FlowBoundaryKind2D::MovingWall,{1,0},0,"lid"};
        } else if(face.areaVector.x<0) b={id,FlowBoundaryKind2D::VelocityInlet,{1,0},0,"inlet"};
        else if(face.areaVector.x>0) b={id,FlowBoundaryKind2D::PressureOutlet,{},0,"outlet"};
        c.boundaryConditions.push_back(b);
    }
    return c;
}

Vector2D rotate(Vector2D p, double angle) {
    return {p.x*std::cos(angle)-p.y*std::sin(angle),p.x*std::sin(angle)+p.y*std::cos(angle)};
}
Point2D rotatePoint(Point2D p,double angle) {const auto r=rotate({p.x,p.y},angle);return {r.x,r.y};}

FvMesh2D rotated(FvMesh2D m,double angle) {
    for(auto& cell:m.cells)cell.centre=rotatePoint(cell.centre,angle);
    for(auto& f:m.faces) {
        f.centre=rotatePoint(f.centre,angle);f.areaVector=rotate(f.areaVector,angle);
        f.correction=rotate(f.correction,angle);
    }
    validateFvMesh2D(m);return m;
}

void compare(const FlowResult2D& a,const FlowResult2D& b,double angle=0,double offset=0) {
    require(a.converged&&b.converged,"solution did not converge");
    double velocityError=0,pressureError=0,fluxError=0;
    for(std::size_t i=0;i<a.u.size();++i) {
        const auto velocity=rotate({a.u[i],a.v[i]},angle);
        velocityError=std::max(velocityError,std::hypot(b.u[i]-velocity.x,b.v[i]-velocity.y));
        pressureError=std::max(pressureError,std::abs(b.p[i]-a.p[i]-offset));
    }
    for(std::size_t i=0;i<a.flux.size();++i)fluxError=std::max(fluxError,std::abs(b.flux[i]-a.flux[i]));
    std::cout<<"max errors u="<<velocityError<<" p="<<pressureError<<" flux="<<fluxError<<'\n';
    require(velocityError<2e-8&&pressureError<2e-8&&fluxError<2e-9,"rotation/pressure-gauge equivalence failed");
}

template<class F> void rejects(F function) {
    bool rejected=false;
    try {function();} catch(const std::runtime_error&) {rejected=true;}
    require(rejected,"invalid custom boundary was accepted");
}

int main() {
    try {
        const auto mesh=rectangle();const auto control=conditions(mesh);
        const auto baseline=solveIncompressible2D(mesh,control);
        std::ostringstream boundaryFile;writeFlowBoundaryConditions2D(boundaryFile,mesh,control);
        std::istringstream boundaryInput(boundaryFile.str());auto fromFile=control;
        fromFile.boundaryConditions=readFlowBoundaryConditions2D(boundaryInput,mesh,control);
        compare(baseline,solveIncompressible2D(mesh,fromFile));
        rejects([&]{std::istringstream changed(boundaryFile.str());readFlowBoundaryConditions2D(changed,rotated(mesh,.3),control);});
        rejects([&]{std::istringstream changed(boundaryFile.str()+"EXTRA");readFlowBoundaryConditions2D(changed,mesh,control);});
        auto preset=control;preset.scenario="duct";preset.boundaryConditions.clear();
        compare(baseline,solveIncompressible2D(mesh,preset));
        auto templateControl=control;templateControl.boundaryConditions=explicitFlowBoundaryPreset2D(mesh,preset);
        compare(baseline,solveIncompressible2D(mesh,templateControl));
        auto shifted=control;
        for(auto& b:shifted.boundaryConditions)if(b.kind==FlowBoundaryKind2D::PressureOutlet)b.pressure=7.25;
        compare(baseline,solveIncompressible2D(mesh,shifted),0,7.25);
        const auto initial=initialIncompressibleState2D(mesh,control);
        const auto step=advanceIncompressible2D(mesh,control,initial,.03);
        compare(step,advanceIncompressible2D(mesh,shifted,initialIncompressibleState2D(mesh,shifted),.03),0,7.25);
        for (double angle:{std::acos(-1.)/2,.63}) {
            const auto turnedMesh=rotated(mesh,angle);auto turned=control;
            for(auto& b:turned.boundaryConditions)b.velocity=rotate(b.velocity,angle);
            compare(baseline,solveIncompressible2D(turnedMesh,turned),angle);
            compare(step,advanceIncompressible2D(turnedMesh,turned,initialIncompressibleState2D(turnedMesh,turned),.03),angle);
        }
        const auto square=rectangle(10,10,1);auto closed=conditions(square,true);
        auto lid=closed;lid.scenario="cavity";lid.boundaryConditions.clear();
        compare(solveIncompressible2D(square,lid),solveIncompressible2D(square,closed));
        auto bad=control;bad.boundaryConditions.pop_back();rejects([&]{validateFlowBoundaryConditions2D(mesh,bad);});
        bad=control;bad.boundaryConditions.push_back(bad.boundaryConditions.front());rejects([&]{validateFlowBoundaryConditions2D(mesh,bad);});
        bad=control;bad.boundaryConditions.front().face=mesh.faces.size();rejects([&]{validateFlowBoundaryConditions2D(mesh,bad);});
        bad=control;bad.boundaryConditions.front().name="";rejects([&]{validateFlowBoundaryConditions2D(mesh,bad);});
        bad=control;bad.boundaryConditions.front().velocity.x=std::numeric_limits<double>::quiet_NaN();rejects([&]{validateFlowBoundaryConditions2D(mesh,bad);});
        bad=control;bad.scenario="duct";rejects([&]{solveIncompressible2D(mesh,bad);});
        bad=control;bad.outletBackflow=OutletBackflow2D::NormalInlet;rejects([&]{validateFlowBoundaryConditions2D(mesh,bad);});
        bad=control;
        for(auto& b:bad.boundaryConditions)if(b.kind==FlowBoundaryKind2D::Wall) {
            b.kind=FlowBoundaryKind2D::MovingWall;b.velocity=mesh.faces[b.face].areaVector;break;
        }
        rejects([&]{validateFlowBoundaryConditions2D(mesh,bad);});
        const FlowState2D accepted{step.time,step.u,step.v,step.p,step.flux};
        std::ostringstream output;writeFlowCheckpoint2D(output,mesh,control,accepted);
        require(output.str().starts_with("CARTMESH2D_FLOW_CHECKPOINT 4\n"),"explicit checkpoint version missing");
        auto reordered=control;std::reverse(reordered.boundaryConditions.begin(),reordered.boundaryConditions.end());
        std::istringstream saved(output.str());const auto restored=readFlowCheckpoint2D(saved,mesh,reordered);
        require(restored.time==accepted.time && restored.u==accepted.u && restored.v==accepted.v &&
                restored.p==accepted.p && restored.flux==accepted.flux,"explicit checkpoint changed accepted state");
        const auto continuous=advanceIncompressible2D(mesh,control,accepted,.03);
        compare(continuous,advanceIncompressible2D(mesh,reordered,restored,.03));
        rejects([&]{std::istringstream changed(output.str());readFlowCheckpoint2D(changed,mesh,shifted);});
        bad=control;bad.boundaryConditions.front().name="different-wall";
        rejects([&]{std::istringstream changed(output.str());readFlowCheckpoint2D(changed,mesh,bad);});
        auto viscous=control;viscous.faceViscosity.assign(mesh.faces.size(),.2);
        std::ostringstream variable;writeFlowCheckpoint2D(variable,mesh,viscous,accepted);
        std::istringstream variableInput(variable.str());readFlowCheckpoint2D(variableInput,mesh,viscous);
        rejects([&]{std::istringstream missingViscosity(variable.str());readFlowCheckpoint2D(missingViscosity,mesh,control);});
        std::cout<<"explicit boundary rotation, pressure gauge, closed lid and validation passed\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
