#include "cartmesh2d/fv/Euler2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
volatile std::sig_atomic_t stopped=0;
void stop(int){stopped=1;}
void require(bool b,const std::string& message){if(!b)throw std::runtime_error(message);}
double number(const std::string& value) {
    std::size_t end=0;const double result=std::stod(value,&end);
    require(end==value.size()&&std::isfinite(result),"expected a finite number: "+value);return result;
}
std::size_t count(const std::string& value,std::size_t upper) {
    const double n=number(value);require(n>=1&&n<=static_cast<double>(upper)&&std::floor(n)==n,"invalid integer limit");
    return static_cast<std::size_t>(n);
}
std::string quote(const std::string& s) {
    std::ostringstream out;out<<'"';
    for(char raw:s) {
        const auto c=static_cast<unsigned char>(raw);
        if(c=='"'||c=='\\')out<<'\\'<<c;
        else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<static_cast<int>(c)<<std::dec;
        else out<<c;
    }
    return out.str()+'"';
}
std::ofstream output(const std::string& path) {
    require(!std::filesystem::is_symlink(path),"refusing symlink output: "+path);
    const auto parent=std::filesystem::path(path).parent_path();if(!parent.empty())std::filesystem::create_directories(parent);
    std::ofstream out(path);out.exceptions(std::ios::badbit|std::ios::failbit);out<<std::setprecision(17);return out;
}
std::string kindName(EulerBoundaryKind2D kind) {
    switch(kind) {
    case EulerBoundaryKind2D::SlipWall:return "slip-wall";
    case EulerBoundaryKind2D::Transmissive:return "transmissive";
    case EulerBoundaryKind2D::Farfield:return "farfield";
    case EulerBoundaryKind2D::Periodic:return "periodic";
    }
    throw std::runtime_error("invalid Euler boundary kind");
}
void boundaryFile(std::ostream& out,const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& bc) {
    out<<"CM2D_EULER_BOUNDARY 1\nMESH "<<mesh.cells.size()<<' '<<mesh.faces.size()<<'\n';
    for(const auto& b:bc) {
        const auto& f=mesh.faces[b.face];
        out<<b.face<<' '<<f.owner<<' '<<f.centre.x<<' '<<f.centre.y<<' '<<f.areaVector.x<<' '<<f.areaVector.y<<' '<<kindName(b.kind)<<' '<<b.reference.density<<' '<<b.reference.u<<' '<<b.reference.v<<' '<<b.reference.pressure<<' ';
        if(b.partner)out<<*b.partner;else out<<'-';
        out<<' '<<std::quoted(b.name)<<'\n';
    }
    out<<"END\n";
}
std::vector<EulerBoundary2D> readBoundaries(const std::string& path,const FvMesh2D& mesh) {
    std::ifstream in(path);std::string line;require(static_cast<bool>(std::getline(in,line))&&line=="CM2D_EULER_BOUNDARY 1","invalid Euler boundary header");
    const auto faces=mesh.faces.size();
    require(static_cast<bool>(std::getline(in,line)),"missing boundary mesh binding");
    std::ostringstream meshLine;meshLine<<"MESH "<<mesh.cells.size()<<' '<<faces;
    require(line==meshLine.str(),"boundary mesh size mismatch");
    std::vector<EulerBoundary2D> result;bool ended=false;
    while(std::getline(in,line)) {
        if(line=="END"){ended=true;break;}
        std::istringstream row(line);std::string id,kind,partner,trailing;EulerBoundary2D b;std::size_t owner=0;double x=0,y=0,sx=0,sy=0;
        require(static_cast<bool>(row>>id>>owner>>x>>y>>sx>>sy>>kind>>b.reference.density>>b.reference.u>>b.reference.v>>b.reference.pressure>>partner>>std::quoted(b.name))&&!(row>>trailing),"invalid Euler boundary row");
        const auto faceNumber=number(id);require(faceNumber>=0&&faceNumber<static_cast<double>(faces)&&faceNumber==std::floor(faceNumber),"invalid boundary face");b.face=static_cast<std::size_t>(faceNumber);
        const auto& face=mesh.faces[b.face];
        require(owner==face.owner&&x==face.centre.x&&y==face.centre.y&&sx==face.areaVector.x&&sy==face.areaVector.y,"boundary mesh geometry mismatch");
        if(kind=="slip-wall")b.kind=EulerBoundaryKind2D::SlipWall;
        else if(kind=="transmissive")b.kind=EulerBoundaryKind2D::Transmissive;
        else if(kind=="farfield")b.kind=EulerBoundaryKind2D::Farfield;
        else if(kind=="periodic")b.kind=EulerBoundaryKind2D::Periodic;
        else throw std::runtime_error("unknown Euler boundary kind: "+kind);
        if(partner!="-") {const auto n=number(partner);require(n>=0&&n<static_cast<double>(faces)&&n==std::floor(n),"invalid periodic partner");b.partner=static_cast<std::size_t>(n);}
        result.push_back(b);require(result.size()<=faces,"too many boundary rows");
    }
    std::string extra;require(ended&&!(in>>extra),"Euler boundary terminator missing or trailing data");return result;
}
}
int main(int argc,char** argv) {
 try {
    std::string meshPath,prefix,problem="sod",boundaryPath,exportBoundary,restart;
    IdealGas2D gas;EulerPrimitive2D reference{1,0,0,1};EulerStepControls2D controls;
    double endTime=.2,split=.5,beta=5,maximumSeconds=180;std::size_t maximumSteps=100000,checkpointEvery=25;
    bool specifiedReference=false;
    for(int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        if(arg=="--help") {
            std::cout<<"Native 2D ideal-gas Euler; inviscid, first-order Rusanov / forward Euler.\n"
                "--mesh FINAL.solver.cm2d --output PREFIX --case sod|uniform|external|vortex|custom\n"
                "--end-time .2 --max-step 1 --min-step 1e-14 --cfl .4 (0 < CFL <= .45)\n"
                "--gamma 1.4 --gas-r 287.05 --density 1 --u 0 --v 0 --pressure 1\n"
                "Sod: left (rho,p), right (.125*rho,.1*p), u=v=0, split x fraction --split .5.\n"
                "Vortex: fixed ambient rho=p=1, u=v=1, --vortex-strength 5; periodic rectangle width/height >= 20.\n"
                "External: embedded slip walls, characteristic outer farfield. Uniform: all farfield.\n"
                "Custom: --boundary FILE, full explicit face coverage; uniform initial state.\n"
                "--export-boundaries FILE exports the selected preset without solving.\n"
                "--restart PREFIX.checkpoint --checkpoint-every 25 --max-steps 100000 --max-seconds 180\n"
                "SI density kg/m3, absolute pressure Pa, velocity m/s, R J/(kg K); unit-depth integrals.\n"
                "Outputs .json .fields.json .cells.csv .faces.csv .history.csv .vtk .checkpoint .boundaries\n";return 0;
        }
        require(i+1<argc,"missing value for "+arg);const std::string value=argv[++i];
        if(arg=="--mesh")meshPath=value;else if(arg=="--output")prefix=value;else if(arg=="--case")problem=value;
        else if(arg=="--boundary")boundaryPath=value;else if(arg=="--export-boundaries")exportBoundary=value;
        else if(arg=="--restart")restart=value;
        else if(arg=="--end-time")endTime=number(value);else if(arg=="--max-step")controls.maximumStep=number(value);
        else if(arg=="--min-step")controls.minimumStep=number(value);else if(arg=="--cfl")controls.acousticCourant=number(value);
        else if(arg=="--gamma")gas.gamma=number(value);else if(arg=="--gas-r")gas.gasConstant=number(value);
        else if(arg=="--density"){reference.density=number(value);specifiedReference=true;}
        else if(arg=="--u"){reference.u=number(value);specifiedReference=true;}
        else if(arg=="--v"){reference.v=number(value);specifiedReference=true;}
        else if(arg=="--pressure"){reference.pressure=number(value);specifiedReference=true;}
        else if(arg=="--split")split=number(value);else if(arg=="--vortex-strength")beta=number(value);
        else if(arg=="--max-seconds")maximumSeconds=number(value);
        else if(arg=="--max-steps")maximumSteps=count(value,10000000);
        else if(arg=="--checkpoint-every")checkpointEvery=count(value,10000);
        else throw std::runtime_error("unknown argument "+arg);
    }
    require(!meshPath.empty()&&(!prefix.empty()||!exportBoundary.empty()),"--mesh and --output (or --export-boundaries) required");
    require(meshPath.ends_with(".solver.cm2d")&&!meshPath.ends_with(".failed.solver.cm2d"),"requires final *.solver.cm2d");
    require(problem=="sod"||problem=="uniform"||problem=="external"||problem=="vortex"||problem=="custom","unknown Euler case");
    require((problem=="custom")==!boundaryPath.empty(),"only custom Euler case requires --boundary");
    require(endTime>0&&maximumSeconds>0&&controls.minimumStep>0&&controls.maximumStep>=controls.minimumStep&&
        controls.acousticCourant>0&&controls.acousticCourant<=.45,"invalid Euler time controls");
    validateIdealGas2D(gas);(void)eulerConservative2D(reference,gas);
    require(problem=="sod"||(split==.5),"--split is only supported for Sod");
    require(problem=="vortex"||(beta==5),"--vortex-strength is only supported for vortex");
    if(problem=="sod")require(split>0&&split<1&&reference.u==0&&reference.v==0,"Sod requires split in (0,1) and zero initial velocity");
    if(problem=="vortex")require(!specifiedReference&&beta>0&&beta<10,"vortex has fixed ambient state; strength must be in (0,10)");
    auto read=readCm2dTopology(meshPath);require(read.valid(),read.error);const auto mesh=makeFvMesh2D(read.topology);
    double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY;
    for(const auto& v:read.topology.vertices){xmin=std::min(xmin,v.point.x);xmax=std::max(xmax,v.point.x);ymin=std::min(ymin,v.point.y);ymax=std::max(ymax,v.point.y);}
    const double scale=std::max(xmax-xmin,ymax-ymin);const double geometryTolerance=1e-10*scale;
    if(problem=="vortex")require(xmax-xmin>=20&&ymax-ymin>=20,"vortex periodic domain must be at least 20 by 20 to suppress artificial tail jumps");
    std::vector<EulerBoundary2D> bc;
    if(problem=="custom")bc=readBoundaries(boundaryPath,mesh);
    else for(std::size_t i=0;i<mesh.faces.size();++i) {
        const auto& f=mesh.faces[i];if(f.neighbour)continue;
        EulerBoundary2D b{i,EulerBoundaryKind2D::Farfield,reference,{},"farfield"};
        if(problem=="external"&&f.patch==BoundaryPatch2D::EmbeddedBoundary){b.kind=EulerBoundaryKind2D::SlipWall;b.name="body";}
        else if(problem=="sod"||problem=="vortex") {
            const double length=std::hypot(f.areaVector.x,f.areaVector.y);
            const bool x=std::abs(f.areaVector.x)>std::abs(f.areaVector.y);
            require(std::abs(x?f.areaVector.y:f.areaVector.x)<1e-12*length,"Sod/vortex preset requires an axis-aligned rectangle");
            require(x?(std::abs(f.centre.x-xmin)<geometryTolerance||std::abs(f.centre.x-xmax)<geometryTolerance):
                      (std::abs(f.centre.y-ymin)<geometryTolerance||std::abs(f.centre.y-ymax)<geometryTolerance),"Sod/vortex preset does not support holes");
            b.kind=problem=="vortex"?EulerBoundaryKind2D::Periodic:(x?EulerBoundaryKind2D::Transmissive:EulerBoundaryKind2D::SlipWall);
            b.name=problem=="vortex"?(x?"periodic-x":"periodic-y"):(x?"ends":"walls");
        } else require(problem=="external"||f.patch!=BoundaryPatch2D::EmbeddedBoundary,"use external or custom for embedded walls");
        bc.push_back(b);
    }
    if(problem=="vortex")for(auto& b:bc) {
        const auto& f=mesh.faces[b.face];const bool x=b.name=="periodic-x";
        for(const auto& other:bc) {
            const auto& g=mesh.faces[other.face];
            if(other.face!=b.face&&other.name==b.name&&std::hypot(f.areaVector.x+g.areaVector.x,f.areaVector.y+g.areaVector.y)<1e-10*std::hypot(f.areaVector.x,f.areaVector.y)&&
               std::abs(x?f.centre.y-g.centre.y:f.centre.x-g.centre.x)<geometryTolerance) {
                require(!b.partner,"ambiguous periodic pairing");b.partner=other.face;
            }
        }
    }
    validateEulerBoundaries2D(mesh,bc,gas);
    if(!exportBoundary.empty()) {
        require(restart.empty(),"boundary export cannot restart");
        require(std::filesystem::weakly_canonical(exportBoundary)!=std::filesystem::weakly_canonical(meshPath)&&
            (boundaryPath.empty()||std::filesystem::weakly_canonical(exportBoundary)!=std::filesystem::weakly_canonical(boundaryPath)),"boundary export collides with input");
        auto out=output(exportBoundary);boundaryFile(out,mesh,bc);return 0;
    }
    EulerState2D state;state.cells.reserve(mesh.cells.size());
    for(const auto& c:mesh.cells) {
        auto q=reference;
        if(problem=="sod"&&c.centre.x>=xmin+split*(xmax-xmin)){q.density*=.125;q.pressure*=.1;}
        if(problem=="vortex") {
            const double x=c.centre.x-.5*(xmin+xmax),y=c.centre.y-.5*(ymin+ymax),r2=x*x+y*y;
            const double factor=beta/(2*std::numbers::pi)*std::exp((1-r2)/2);
            const double temperature=1-(gas.gamma-1)*beta*beta/(8*gas.gamma*std::numbers::pi*std::numbers::pi)*std::exp(1-r2);
            require(temperature>0,"vortex strength gives non-positive temperature");
            q={std::pow(temperature,1/(gas.gamma-1)),1-y*factor,1+x*factor,0};q.pressure=std::pow(q.density,gas.gamma);
        }
        state.cells.push_back(eulerConservative2D(q,gas));
    }
    if(!restart.empty()) {std::ifstream in(restart);require(static_cast<bool>(in),"cannot read Euler restart");state=readEulerCheckpoint2D(in,mesh,bc,gas,problem);}
    require(endTime>state.time,"end time must exceed accepted restart time");
    const std::vector<std::string> suffixes={".json",".fields.json",".cells.csv",".faces.csv",".history.csv",".vtk",".checkpoint",".checkpoint.tmp",".boundaries"};
    for(const auto& suffix:suffixes) {
        const auto path=std::filesystem::weakly_canonical(prefix+suffix);
        require(path!=std::filesystem::weakly_canonical(meshPath)&&
            (boundaryPath.empty()||path!=std::filesystem::weakly_canonical(boundaryPath)),"output collides with input");
        require(!std::filesystem::is_symlink(prefix+suffix),"refusing symlink output");
    }
    const auto save=[&]{auto out=output(prefix+".checkpoint.tmp");writeEulerCheckpoint2D(out,mesh,bc,gas,state,problem);out.close();std::filesystem::rename(prefix+".checkpoint.tmp",prefix+".checkpoint");};
    save();auto boundaryOutput=output(prefix+".boundaries");boundaryFile(boundaryOutput,mesh,bc);boundaryOutput.close();
    auto history=output(prefix+".history.csv");history<<"step,time,dt,acousticCourant,minimumDensity,minimumPressure,cellBalanceError,rejectedCandidates,mass,momentumX,momentumY,totalEnergy,boundaryMass,boundaryMomentumX,boundaryMomentumY,boundaryEnergy,balanceMass,balanceMomentumX,balanceMomentumY,balanceEnergy\n";
    const auto started=std::chrono::steady_clock::now();const double initialTime=state.time;const auto initialSteps=state.steps;
    std::signal(SIGINT,stop);std::signal(SIGTERM,stop);std::optional<EulerStepResult2D> last;
    std::size_t rejected=0;std::string status="target_reached",failure;
    try {
        while(state.time<endTime) {
            require(!stopped,"calculation cancelled; accepted state retained");
            require(state.steps-initialSteps<maximumSteps,"accepted-step budget exhausted");
            require(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<maximumSeconds,"wall-time budget exhausted");
            auto step=controls;step.maximumStep=std::min(step.maximumStep,endTime-state.time);
            last=advanceEuler2D(mesh,bc,gas,state,step);state=last->state;rejected+=last->rejectedCandidates;
            history<<state.steps<<','<<state.time<<','<<last->step<<','<<last->acousticCourant<<','<<last->minimumDensity<<','<<last->minimumPressure<<','<<last->maximumCellBalanceError<<','<<last->rejectedCandidates;
            for(const auto& array:{last->afterIntegral,last->boundaryFlux,last->balanceError})for(double v:array)history<<','<<v;
            history<<'\n';
            if((state.steps-initialSteps)%checkpointEvery==0){save();history.flush();}
        }
    }catch(const std::exception& e){status=stopped?"cancelled":"failed";failure=e.what();}
    save();history.close();
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    auto cells=output(prefix+".cells.csv"),fields=output(prefix+".fields.json");
    cells<<"cell,x,y,area,rho,u,v,p,rhoU,rhoV,rhoE,temperature,mach,previousRho,previousRhoU,previousRhoV,previousRhoE\n";
    fields<<"{\"format\":\"cartmesh2d-euler-v1\",\"time\":"<<state.time<<",\"cells\":[\n";
    for(std::size_t i=0;i<state.cells.size();++i) {
        const auto& c=mesh.cells[i];const auto& q=state.cells[i];const auto p=eulerPrimitive2D(q,gas);
        const double temperature=p.pressure/(p.density*gas.gasConstant),mach=std::hypot(p.u,p.v)/eulerSoundSpeed2D(p,gas);
        require(std::isfinite(temperature)&&temperature>0&&std::isfinite(mach),"unrepresentable exported temperature or Mach number");
        cells<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','<<p.density<<','<<p.u<<','<<p.v<<','<<p.pressure<<','<<q[1]<<','<<q[2]<<','<<q[3]<<','<<temperature<<','<<mach;
        for(double v:last?last->previousCells[i]:q)cells<<','<<v;
        cells<<'\n';if(i)fields<<",\n";
        fields<<"{\"id\":"<<i<<",\"rho\":"<<p.density<<",\"u\":"<<p.u<<",\"v\":"<<p.v<<",\"p\":"<<p.pressure
            <<",\"rhoE\":"<<q[3]<<",\"temperature\":"<<temperature<<",\"mach\":"<<mach<<'}';
    }
    fields<<"\n]}\n";cells.close();fields.close();
    auto faces=output(prefix+".faces.csv");faces<<"face,owner,neighbour,x,y,sx,sy,kind,partner,waveSpeed,mass,momentumX,momentumY,energy\n";
    std::vector<const EulerBoundary2D*> lookup(mesh.faces.size());for(const auto& b:bc)lookup[b.face]=&b;
    for(std::size_t i=0;i<mesh.faces.size();++i) {
        const auto& f=mesh.faces[i];faces<<i<<','<<f.owner<<',';if(f.neighbour)faces<<*f.neighbour;else faces<<-1;
        faces<<','<<f.centre.x<<','<<f.centre.y<<','<<f.areaVector.x<<','<<f.areaVector.y<<','<<(lookup[i]?kindName(lookup[i]->kind):"internal")<<',';
        if(lookup[i]&&lookup[i]->partner)faces<<*lookup[i]->partner;else faces<<-1;
        faces<<','<<(last?last->faceWaveSpeed[i]:0);
        for(double value:last?last->faceFlux[i]:EulerConservative2D{})faces<<','<<value;
        faces<<'\n';
    }
    faces.close();std::string vtkError;require(writeLegacyVtk2D(read.topology,prefix+".vtk",&vtkError),vtkError);
    std::ofstream vtk(prefix+".vtk",std::ios::app);vtk.exceptions(std::ios::badbit|std::ios::failbit);vtk<<std::setprecision(17);
    for(const std::string field:{"density","pressure","temperature","mach"}) {
        vtk<<"SCALARS "<<field<<" double 1\nLOOKUP_TABLE default\n";
        for(const auto& q:state.cells){const auto p=eulerPrimitive2D(q,gas);vtk<<(field=="density"?p.density:field=="pressure"?p.pressure:field=="temperature"?p.pressure/(p.density*gas.gasConstant):std::hypot(p.u,p.v)/eulerSoundSpeed2D(p,gas))<<'\n';}
    }
    vtk<<"VECTORS velocity double\n";for(const auto& q:state.cells)vtk<<q[1]/q[0]<<' '<<q[2]/q[0]<<" 0\n";vtk.close();
    auto summary=output(prefix+".json");
    summary<<"{\n\"solver\":\"native 2D ideal-gas Euler\",\"method\":\"first-order Rusanov / forward Euler\",\"case\":"<<quote(problem)
        <<",\"status\":"<<quote(status)<<",\"failure\":"<<quote(failure)<<",\"targetReached\":"<<(status=="target_reached"?"true":"false")
        <<",\"cells\":"<<mesh.cells.size()<<",\"faces\":"<<mesh.faces.size()<<",\"gamma\":"<<gas.gamma<<",\"gasConstant\":"<<gas.gasConstant
        <<",\"time\":"<<state.time<<",\"requestedEndTime\":"<<endTime<<",\"initialTime\":"<<initialTime<<",\"steps\":"<<state.steps
        <<",\"acceptedSteps\":"<<state.steps-initialSteps<<",\"rejectedCandidates\":"<<rejected<<",\"lastStep\":"<<(last?last->step:0)
        <<",\"cflLimit\":"<<controls.acousticCourant<<",\"maximumStep\":"<<controls.maximumStep<<",\"minimumStep\":"<<controls.minimumStep
        <<",\"maximumSteps\":"<<maximumSteps<<",\"maximumSeconds\":"<<maximumSeconds<<",\"checkpointEvery\":"<<checkpointEvery<<",\"elapsedSeconds\":"<<elapsed
        <<",\"nativeTopologyRevalidated\":true,\"solverQualityPassed\":true,\"externalCheckMesh\":\"not run\""
        <<",\"units\":{\"rho\":\"kg/m3\",\"p\":\"Pa absolute\",\"rhoE\":\"J/m3\",\"temperature\":\"K\",\"flux\":\"outward-owner per unit depth\"}"
        <<",\"initialization\":\"cell-centre sampling; restart replaces complete conserved state\",\"scope\":\"inviscid ideal gas; no viscous stress, heat conduction or turbulence; target time is not steady convergence\"\n}\n";
    summary.close();std::cout<<status<<": "<<mesh.cells.size()<<" cells, t="<<std::setprecision(17)<<state.time<<", "<<state.steps-initialSteps<<" accepted steps, "<<elapsed<<" s\n";
    if(!failure.empty())std::cerr<<failure<<'\n';return status=="target_reached"?0:2;
 }catch(const std::exception& e){std::cerr<<"cartmesh2d_euler_cli: "<<e.what()<<'\n';return 1;}
}
