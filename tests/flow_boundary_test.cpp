#include <cstdlib>
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

void steadyRelaxationIndependence() {
    // Non-orthogonal cells expose the face/cell interpolation defect that
    // vanishes for the linear pressure of an orthogonal Poiseuille mesh.
    const auto mesh=rectangle(18,6,3,true);auto control=conditions(mesh);
    control.nu=.1;control.tolerance=1e-11;control.maxIterations=8000;
    control.convection=ConvectionScheme2D::FaceLimitedLinearUpwind;
    for(auto& b:control.boundaryConditions) {
        const auto& face=mesh.faces[b.face];
        if(b.kind==FlowBoundaryKind2D::VelocityInlet || b.kind==FlowBoundaryKind2D::PressureOutlet)
            b={b.face,FlowBoundaryKind2D::PressureOpening,{},face.areaVector.x<0?.6:0,b.name};
        else if(face.areaVector.y>0)b={b.face,FlowBoundaryKind2D::Symmetry,{},0,"symmetry"};
    }
    const auto baseline=solveIncompressible2D(mesh,control);
#ifdef __APPLE__
    auto direct=control;direct.pressurePreconditioner=PressurePreconditioner2D::SystemCholesky;direct.profile=true;
    const auto directResult=solveIncompressible2D(mesh,direct);
    compare(baseline,directResult);
    require(directResult.performance.pressureCholeskyBuilds>0,"Warped channel did not use Cholesky");
#endif
    auto shortRun=control;shortRun.maxIterations=30;
    const auto partial=solveIncompressible2D(mesh,shortRun);
    require(!partial.converged,"continuation fixture unexpectedly converged");
    FlowInitialGuess2D guess{partial.u,partial.v,partial.p,partial.flux};
    shortRun.maxIterations=1;
    const auto resumed=solveIncompressibleFromGuess2D(mesh,shortRun,guess);
    shortRun.maxIterations=31;
    const auto uninterrupted=solveIncompressible2D(mesh,shortRun);
    double difference=0;
    for(const auto fields:{std::make_pair(&resumed.u,&uninterrupted.u),std::make_pair(&resumed.v,&uninterrupted.v),std::make_pair(&resumed.p,&uninterrupted.p),std::make_pair(&resumed.flux,&uninterrupted.flux)})
        for(std::size_t i=0;i<fields.first->size();++i)difference=std::max(difference,std::abs((*fields.first)[i]-(*fields.second)[i]));
    require(difference<2e-10,"cell and face initial iterate lost a SIMPLE step");
    require(!resumed.converged,"finite partial iterate bypassed convergence gates");
    auto invalid=guess;invalid.flux.pop_back();
    rejects([&]{(void)solveIncompressibleFromGuess2D(mesh,shortRun,invalid);});
    invalid=guess;invalid.flux[0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{(void)solveIncompressibleFromGuess2D(mesh,shortRun,invalid);});
    invalid=guess;
    for(const auto& b:control.boundaryConditions)if(b.kind==FlowBoundaryKind2D::Symmetry) {invalid.flux[b.face]=.1;break;}
    rejects([&]{(void)solveIncompressibleFromGuess2D(mesh,shortRun,invalid);});
    for(double alpha:{.8,.9}) {
        auto changed=control;changed.velocityRelaxation=alpha;
        compare(baseline,solveIncompressible2D(mesh,changed));
    }
    auto accelerated=control;accelerated.steadyAcceleration=SteadyAcceleration2D::Anderson;
    const auto faster=solveIncompressible2D(mesh,accelerated);
    compare(baseline,faster);
    require(faster.performance.accelerationAccepted>0,"Anderson never accepted a candidate");
    require(faster.history.size()<baseline.history.size(),"Anderson did not accelerate the warped channel");
    compare(faster,solveIncompressible2D(rotated(mesh,std::acos(-1.)/2),accelerated),std::acos(-1.)/2);
    auto shifted=accelerated;
    for(auto& b:shifted.boundaryConditions)if(b.kind==FlowBoundaryKind2D::PressureOpening)b.pressure+=7.25;
    compare(faster,solveIncompressible2D(mesh,shifted),0,7.25);
    rejects([&]{(void)advanceIncompressible2D(mesh,accelerated,initialIncompressibleState2D(mesh,control),.02);});
    const auto cavity=rectangle(10,10,1);auto closed=conditions(cavity,true);
    const auto original=solveIncompressible2D(cavity,closed);
#ifdef __APPLE__
    auto closedDirect=closed;closedDirect.pressurePreconditioner=PressurePreconditioner2D::SystemCholesky;
    compare(original,solveIncompressible2D(cavity,closedDirect));
#endif
    closed.steadyAcceleration=SteadyAcceleration2D::Anderson;
    compare(original,solveIncompressible2D(cavity,closed));
}

void pressureOpenings() {
    double previousError=0;
    for (int n:{8,16}) {
        const auto mesh=rectangle(4*n,n);auto control=conditions(mesh);
        control.nu=.1;
        control.convection=ConvectionScheme2D::FaceLimitedLinearUpwind;
        for(auto& b:control.boundaryConditions) {
            if(b.kind==FlowBoundaryKind2D::VelocityInlet || b.kind==FlowBoundaryKind2D::PressureOutlet) {
                b.kind=FlowBoundaryKind2D::PressureOpening;b.velocity={};
                b.pressure=mesh.faces[b.face].areaVector.x<0 ? 4.8 : 0;
            }
        }
        const auto result=solveIncompressible2D(mesh,control);
        require(result.converged,"pressure-driven flow did not converge");
        double error=0,pressureError=0,flux=0;
        for(std::size_t i=0;i<mesh.cells.size();++i) {
            const auto point=mesh.cells[i].centre;
            const double exact=6*point.y*(1-point.y);
            error+=mesh.cells[i].area*std::pow(result.u[i]-exact,2);
            pressureError=std::max(pressureError,std::abs(result.p[i]-4.8*(1-point.x/4)));
        }
        error=std::sqrt(error/4);
        for(std::size_t i=0;i<mesh.faces.size();++i)
            if(!mesh.faces[i].neighbour && mesh.faces[i].areaVector.x>0)flux+=result.flux[i];
        std::cout<<"pressure opening n="<<n<<" U L2="<<error<<" p max="<<pressureError<<" flow="<<flux<<'\n';
        require(error<.05 && pressureError<1e-7 && std::abs(flux-1)<.05,"pressure-driven Poiseuille reference failed");
        if(previousError>0)require(previousError/error>3.8,"pressure-driven spatial order below expected second order");
        previousError=error;
        if(n!=8)continue;
        auto shifted=control;
        for(auto& b:shifted.boundaryConditions)if(b.kind==FlowBoundaryKind2D::PressureOpening)b.pressure+=7.25;
        compare(result,solveIncompressible2D(mesh,shifted),0,7.25);
        compare(result,solveIncompressible2D(rotated(mesh,std::acos(-1.)/2),control),std::acos(-1.)/2);
        auto equal=control;
        for(auto& b:equal.boundaryConditions)if(b.kind==FlowBoundaryKind2D::PressureOpening)b.pressure=7.25;
        const auto rest=solveIncompressible2D(mesh,equal);
        require(rest.converged,"equal-pressure state did not converge");
        for(std::size_t i=0;i<rest.u.size();++i)
            require(std::hypot(rest.u[i],rest.v[i])<1e-10 && std::abs(rest.p[i]-7.25)<1e-9,"equal-pressure opening creates flow");
        auto reverse=control;
        for(auto& b:reverse.boundaryConditions)if(b.kind==FlowBoundaryKind2D::PressureOpening)b.pressure=4.8-b.pressure;
        const auto reversed=solveIncompressible2D(mesh,reverse);
        require(reversed.converged,"reverse-pressure flow failed");
        for(std::size_t i=0;i<result.u.size();++i)
            require(std::abs(reversed.u[i]+result.u[i])<2e-8,"reverse-pressure velocity is not symmetric");
        auto rejectBackflow=reverse;
        for(auto& b:rejectBackflow.boundaryConditions)
            if(b.kind==FlowBoundaryKind2D::PressureOpening && mesh.faces[b.face].areaVector.x>0)
                b.kind=FlowBoundaryKind2D::PressureOutlet;
        rejects([&]{(void)solveIncompressible2D(mesh,rejectBackflow);});
        const auto first=advanceIncompressible2D(mesh,control,initialIncompressibleState2D(mesh,control),.03);
        require(first.converged,"pressure startup failed");
        FlowState2D state{first.time,first.u,first.v,first.p,first.flux};
        std::ostringstream file;writeFlowCheckpoint2D(file,mesh,control,state);
        std::istringstream stream(file.str());const auto restored=readFlowCheckpoint2D(stream,mesh,control);
        const auto whole=advanceIncompressible2D(mesh,control,state,.03);
        const auto resumed=advanceIncompressible2D(mesh,control,restored,.03);
        require(whole.converged && resumed.converged && whole.u==resumed.u && whole.v==resumed.v &&
                whole.p==resumed.p && whole.flux==resumed.flux,"pressure restart changed state");
        rejects([&]{validateFlowBoundaryConditions2D(rotated(mesh,.3),control);});
        auto bad=control;
        for(auto& b:bad.boundaryConditions)if(b.kind==FlowBoundaryKind2D::PressureOpening){b.velocity.x=1;break;}
        rejects([&]{validateFlowBoundaryConditions2D(mesh,bad);});
    }
}

void symmetryBoundaries() {
    double previous=0;
    for(int n:{8,16}) {
        const auto mesh=rectangle(4*n,n);auto c=conditions(mesh);
        c.nu=.1;c.maxIterations=10000;c.convection=ConvectionScheme2D::FaceLimitedLinearUpwind;
        for(auto& b:c.boundaryConditions) {
            const auto& f=mesh.faces[b.face];
            if(b.kind==FlowBoundaryKind2D::VelocityInlet || b.kind==FlowBoundaryKind2D::PressureOutlet)
                b={b.face,FlowBoundaryKind2D::PressureOpening,{},f.areaVector.x<0?1.2:0,b.name};
            else if(f.areaVector.y>0)b={b.face,FlowBoundaryKind2D::Symmetry,{},0,"symmetry"};
        }
        const auto result=solveIncompressible2D(mesh,c);require(result.converged,"half-channel did not converge");
        double error=0;
        for(std::size_t i=0;i<mesh.cells.size();++i) {
            const auto y=mesh.cells[i].centre.y;
            error+=mesh.cells[i].area*std::pow(result.u[i]-1.5*y*(2-y),2);
        }
        error=std::sqrt(error/4);std::cout<<"half-channel n="<<n<<" velocity L2="<<error<<'\n';
        require(error<.01,"half-channel analytic velocity failed");
        if(previous>0)require(previous/error>3.8,"half-channel spatial order failed");previous=error;
        for(const auto& b:c.boundaryConditions)if(b.kind==FlowBoundaryKind2D::Symmetry)
            require(result.flux[b.face]==0,"symmetry allows penetration");
        if(n!=8)continue;
        compare(result,solveIncompressible2D(rotated(mesh,std::acos(-1.)/2),c),std::acos(-1.)/2);
        rejects([&]{validateFlowBoundaryConditions2D(rotated(mesh,.3),c);});
        auto bad=c;
        for(auto& b:bad.boundaryConditions)if(b.kind==FlowBoundaryKind2D::Symmetry){b.pressure=1;break;}
        rejects([&]{validateFlowBoundaryConditions2D(mesh,bad);});
        // With both walls frictionless, the pressure gradient gives exactly
        // uniform acceleration dU/dt=Delta(p/rho)/L, not a steady balance.
        auto plug=c;
        for(auto& b:plug.boundaryConditions)if(b.kind==FlowBoundaryKind2D::Wall)
            b={b.face,FlowBoundaryKind2D::Symmetry,{},0,"bottom-symmetry"};
        auto state=initialIncompressibleState2D(mesh,plug);
        for(int step=0;step<3;++step){
            const auto next=advanceIncompressible2D(mesh,plug,state,.02);
            require(next.converged,"frictionless acceleration step failed");
            for(std::size_t i=0;i<next.u.size();++i)
                require(std::abs(next.u[i]-.3*next.time)<2e-8 && std::abs(next.v[i])<2e-8,"free-slip acceleration reference failed");
            state={next.time,next.u,next.v,next.p,next.flux};
        }
        std::ostringstream checkpoint;writeFlowCheckpoint2D(checkpoint,mesh,plug,state);
        std::istringstream input(checkpoint.str());const auto loaded=readFlowCheckpoint2D(input,mesh,plug);
        require(loaded.u==state.u && loaded.flux==state.flux,"symmetry checkpoint changed state");
    }
}

int main() {
#ifdef __APPLE__
    setenv("VECLIB_MAXIMUM_THREADS", "1", 1);
#endif
    try {
        steadyRelaxationIndependence();
        pressureOpenings();
        symmetryBoundaries();
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
