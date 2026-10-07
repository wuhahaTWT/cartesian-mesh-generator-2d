// Reproducible API comparison using existing native manufactured/mesh/state
// adapters. Only this research executable knows case names or writes files.
#define CARTMESH_NEWTON_NO_MAIN
#include "native-laminar-newton.cpp"
#include "cartmesh2d/fv/CompatibleIncompressible2D.hpp"
int main(int argc,char** argv)try {
    if(argc!=11&&argc!=12)throw std::runtime_error("usage: compatible-solver mesh n problem nu stokes|ns closed|traction|pseudo-traction|pressure output-prefix pseudo-step iterations mass|schur [linear-restarts]");
    const auto start=std::chrono::steady_clock::now();const auto fixture=readFixture(argv[1],std::stoi(argv[2]));const auto& mesh=fixture.mesh;
    const std::string problem=argv[3],boundary=argv[6],prefix=argv[7];const double nu=std::stod(argv[4]),step=std::stod(argv[8]);const bool nonlinear=std::string(argv[5])=="ns";
    if((!nonlinear&&std::string(argv[5])!="stokes")||(boundary!="closed"&&boundary!="traction"&&boundary!="pseudo-traction"&&boundary!="pressure")||(problem=="cylinder"&&boundary!="pressure"))throw std::runtime_error("unsupported research adapter; API requires explicit boundary data");
    if(std::filesystem::exists(prefix+".json")||std::filesystem::exists(prefix+".accepted.state"))throw std::runtime_error("existing API evidence output");
    CompatibleFlowControls2D c;c.viscosity=nu;c.maximumIterations=std::stoull(argv[9]);if(argc==12)c.maximumLinearRestarts=std::stoull(argv[11]);c.equation=nonlinear?CompatibleEquation2D::NavierStokes:CompatibleEquation2D::Stokes;
    c.globalization=step>0?CompatibleGlobalization2D::PseudoTime:CompatibleGlobalization2D::Backtracking;if(step>0)c.initialPseudoStep=step;
    c.pressureInverse=std::string(argv[10])=="schur"?CompatiblePressureInverse2D::DiagonalSchur:CompatiblePressureInverse2D::ViscousMass;
    c.acceleration=[=](Point2D p){return forceAt(p,problem,nu,nonlinear);};
    for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour){const auto& face=mesh.faces[id];CompatibleBoundary2D b;b.face=id;
        if(problem=="cylinder"){const double round=128*std::numeric_limits<double>::epsilon()*std::hypot(face.areaVector.x,face.areaVector.y);
            if(face.patch==BoundaryPatch2D::EmbeddedBoundary)b.value=[](Point2D){return Vector2D{};};
            else if(face.patch!=BoundaryPatch2D::DomainBoundary)throw std::runtime_error("cylinder requires explicit domain/embedded patches");
            else if(std::abs(face.areaVector.y)<=round&&face.areaVector.x>0)b.kind=CompatibleBoundaryKind2D::PseudoTraction;
            else if(std::abs(face.areaVector.y)<=round)b.value=[](Point2D){return Vector2D{1,0};};
            else if(std::abs(face.areaVector.x)<=round)b.kind=CompatibleBoundaryKind2D::NormalVelocity;
            else throw std::runtime_error("cylinder far boundary must be axis aligned");
        }
        else if(boundary!="closed"&&face.areaVector.x>0&&std::abs(face.areaVector.y)<128*std::numeric_limits<double>::epsilon()*std::hypot(face.areaVector.x,face.areaVector.y)){
            b.kind=boundary=="traction"?CompatibleBoundaryKind2D::Traction:CompatibleBoundaryKind2D::PseudoTraction;auto n=face.areaVector;const auto length=std::hypot(n.x,n.y);n.x/=length;n.y/=length;
            b.value=[=](Point2D p){const auto e=exactAt(p,problem);Vector2D value{-e.p*n.x,-e.p*n.y};
                if(boundary=="traction"){const auto gu=e.gradient[0],gv=e.gradient[1];value.x+=nu*(2*gu.x*n.x+(gu.y+gv.x)*n.y);value.y+=nu*((gu.y+gv.x)*n.x+2*gv.y*n.y);}return value;};
        }else b.value=[=](Point2D p){return problem=="noslip"||problem=="noslip-sheared"?Vector2D{}:exactAt(p,problem).u;};c.boundaries.push_back(std::move(b));}
    const auto r=solveCompatibleIncompressible2D(mesh,c);const auto save=[&](const auto& state,const std::string& name){if(!state)return;Vec out;for(const auto& cell:state->cells)out.insert(out.end(),cell.begin(),cell.end());for(const auto& face:state->faces)out.insert(out.end(),face.begin(),face.end());writeState(prefix+name+".state",fixture,out);};
    save(r.seed,".seed");save(r.lastAccepted,".accepted");save(r.lastRejected,".rejected");
    std::ostringstream s;s<<std::setprecision(17)<<"{\"completed\":"<<(r.converged()?"true":"false")<<",\"stop\":"<<int(r.stop)<<",\"reason\":\""<<r.reason<<"\",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"cells\":"<<mesh.cells.size()<<",\"faces\":"<<mesh.faces.size()<<",\"iterations\":[";
    for(std::size_t i=0;i<r.iterations.size();++i){const auto& x=r.iterations[i];const auto m=x.metrics.value_or(CompatibleFlowMetrics2D{});if(i)s<<',';s<<"{\"iteration\":"<<x.iteration<<",\"restarts\":"<<x.linearRestarts<<",\"products\":"<<x.matrixProducts<<",\"linearRelativeResidual\":"<<x.linearRelativeResidual<<",\"originalEquationsEvaluated\":"<<(x.metrics?"true":"false")<<",\"accepted\":"<<(x.accepted?"true":"false")<<",\"cellMomentum\":"<<m.cellMomentum<<",\"faceMomentum\":"<<m.faceMomentum<<",\"divergence\":"<<m.divergence<<",\"residualNorm\":"<<m.residualNorm<<",\"stateChange\":"<<m.stateChange<<",\"trials\":"<<x.trials<<",\"alpha\":"<<x.alpha<<",\"pseudoStep\":"<<x.pseudoStep<<'}';}s<<"]}\n";
    std::ofstream log(prefix+".json");log<<s.str();log.close();if(!log)throw std::runtime_error("API result write failed");std::cout<<s.str();return r.converged()?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
