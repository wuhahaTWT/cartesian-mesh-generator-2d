#include "cartmesh2d/fv/ScalarTransport2D.hpp"
#include "cartmesh2d/fv/FlowCheckpoint2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
#include <stdexcept>
using namespace cartmesh2d;
namespace {
void require(bool b,const std::string& s) {if(!b) throw std::runtime_error(s);}
double number(const std::string& s) {
    std::size_t used=0; const auto v=std::stod(s,&used);
    require(used==s.size()&&std::isfinite(v),"expected finite number: "+s); return v;
}
std::string quote(const std::string& s) {
    std::ostringstream out; out<<'"';
    for (char raw:s) {
        const auto ch=static_cast<unsigned char>(raw);
        if(ch=='"'||ch=='\\') out<<'\\'<<ch;
        else if(ch<32) out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<static_cast<int>(ch)<<std::dec;
        else out<<ch;
    }
    return out.str()+'"';
}
std::ofstream output(const std::string& prefix,const char* suffix) {
    std::ofstream out(prefix+suffix); out.exceptions(std::ios::badbit|std::ios::failbit);
    out<<std::setprecision(17); return out;
}
fv::FlowState2D carrier(const std::string& path,const fv::FvMesh2D& mesh) {
    std::ifstream in(path); require(bool(in),"cannot open flow checkpoint");
    // Recover only physical configuration, then run the existing strict full
    // geometry/incidence/state parser from the beginning. No CSV-order guessing.
    std::string line; fv::FlowControls2D c;
    std::getline(in,line); require(line=="CARTMESH2D_FLOW_CHECKPOINT 1","invalid flow checkpoint header");
    std::getline(in,line); require(line=="DISCRETIZATION Euler-RC-v2","unsupported flow checkpoint discretization");
    std::getline(in,line); std::istringstream config(line); std::string token,scheme,stress;
    require(bool(config>>token>>std::quoted(c.scenario)>>c.nu>>c.speed>>scheme>>stress>>c.manufacturedPressureSlope)&&token=="CONFIG",
        "invalid flow checkpoint configuration");
    require(scheme=="upwind"||scheme=="limited-linear","invalid carrier convection");
    require(stress=="symmetric"||stress=="laplacian","invalid carrier stress");
    c.convection=scheme=="upwind"?fv::ConvectionScheme2D::Upwind:fv::ConvectionScheme2D::LimitedLinearUpwind;
    c.viscousStress=stress=="symmetric"?fv::ViscousStress2D::Symmetric:fv::ViscousStress2D::Laplacian;
    in.clear(); in.seekg(0); return fv::readFlowCheckpoint2D(in,mesh,c);
}
std::vector<fv::ScalarBoundary2D> boundaries(const std::string& path,const fv::FvMesh2D& mesh) {
    std::ifstream in(path); require(bool(in),"cannot open scalar boundary CSV");
    std::string line; std::getline(in,line);
    if(!line.empty()&&line.back()=='\r') line.pop_back();
    require(line=="face,type,value,inflowValue","boundary CSV header must be face,type,value,inflowValue");
    std::vector<fv::ScalarBoundary2D> bc(mesh.faces.size()); std::vector<bool> seen(mesh.faces.size(),false);
    while(std::getline(in,line)) {
        if(!line.empty()&&line.back()=='\r') line.pop_back();
        require(!line.empty(),"empty boundary CSV row");
        std::vector<std::string> fields; std::size_t pos=0;
        for(;;) {const auto end=line.find(',',pos); fields.push_back(line.substr(pos,end-pos)); if(end==std::string::npos)break;pos=end+1;}
        require(fields.size()==4,"boundary CSV row requires four fields");
        const double index=number(fields[0]);
        require(index>=0&&index<static_cast<double>(mesh.faces.size())&&index==std::floor(index),"invalid boundary face ID");
        const auto id=static_cast<std::size_t>(index);
        require(!mesh.faces[id].neighbour&&!seen[id],"duplicate or internal boundary face"); seen[id]=true;
        require(fields[1]=="value"||fields[1]=="flux","boundary type must be value or flux");
        bc[id].kind=fields[1]=="value"?fv::ScalarBoundaryKind2D::Value:fv::ScalarBoundaryKind2D::DiffusiveFlux;
        bc[id].value=number(fields[2]);
        if(!fields[3].empty()) bc[id].inflowValue=number(fields[3]);
    }
    for(std::size_t id=0;id<mesh.faces.size();++id)
        require(mesh.faces[id].neighbour.has_value()||seen[id],"missing scalar boundary face");
    return bc;
}
}
int main(int argc,char**argv) {
    std::string prefix; bool started=false;
    try {
        std::string meshPath,flowPath,bcPath,verification;
        double diffusivity=.01,source=0,initial=0,dt=0,speed=1; std::size_t steps=1;
        fv::ScalarTransportControls2D controls;
        for(int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if(arg=="--help") {
                std::cout<<"Native constant-property passive scalar / temperature transport.\n"
                    "--mesh FINAL.solver.cm2d --flow-checkpoint FLOW.checkpoint --boundary BC.csv --output PREFIX\n"
                    "BC.csv: face,type,value,inflowValue; every boundary face, type=value|flux.\n"
                    "flux is outward -D grad(s).n per length; negative carrier flux requires inflow value.\n"
                    "--diffusivity .01 --source 0 --initial 0 --convection upwind|limited-linear\n"
                    "--dt DT --steps N: backward Euler on FROZEN carrier flux (not coupled evolving flow).\n"
                    "Without dt: steady. Thermal D=k/(rho cp), source=Q/(rho cp), boundary flux=q/(rho cp).\n"
                    "Verification only: --verification sine|decay (unit-square decay), --speed 1 for sine.\n"
                    "Outputs .json .vtk .cells.csv .faces.csv .history.csv. No scalar restart/UI yet.\n";
                return 0;
            }
            require(i+1<argc,"missing option value"); const std::string value=argv[++i];
            if(arg=="--mesh")meshPath=value; else if(arg=="--output")prefix=value;
            else if(arg=="--flow-checkpoint")flowPath=value; else if(arg=="--boundary")bcPath=value;
            else if(arg=="--verification")verification=value;
            else if(arg=="--diffusivity")diffusivity=number(value);
            else if(arg=="--source")source=number(value); else if(arg=="--initial")initial=number(value);
            else if(arg=="--dt")dt=number(value); else if(arg=="--speed")speed=number(value);
            else if(arg=="--steps"||arg=="--max-corrections") {
                const auto n=number(value); require(n>=1&&n<=100000&&n==std::floor(n),"invalid iteration count");
                if(arg=="--steps")steps=static_cast<std::size_t>(n);else controls.maxCorrections=static_cast<std::size_t>(n);
            } else if(arg=="--convection") {
                require(value=="upwind"||value=="limited-linear","unknown convection scheme");
                controls.convection=value=="upwind"?fv::ConvectionScheme2D::Upwind:fv::ConvectionScheme2D::LimitedLinearUpwind;
            } else throw std::runtime_error("unknown option: "+arg);
        }
        require(!meshPath.empty()&&!prefix.empty(),"--mesh and --output required");
        require(meshPath.ends_with(".solver.cm2d")&&!meshPath.ends_with(".failed.solver.cm2d"),"requires final solver mesh");
        require(dt>=0&&diffusivity>0&&(dt>0||steps==1),"invalid time settings or diffusivity");
        require(verification.empty()||verification=="sine"||verification=="decay","unknown verification");
        require(verification.empty()?(!flowPath.empty()&&!bcPath.empty()):(flowPath.empty()&&bcPath.empty()),"choose explicit carrier/boundary files OR verification");
        require(verification!="sine"||dt==0,"sine verification is steady");
        require(verification!="decay"||dt>0,"decay verification needs dt");
        auto read=readCm2dTopology(meshPath); require(read.valid(),read.error);
        const auto mesh=fv::makeFvMesh2D(read.topology);
        fv::ScalarTransportProblem2D p; p.diffusivity=diffusivity;
        std::vector<fv::ScalarBoundary2D> bc;
        double carrierTime=0,time=0;
        const double pi=std::numbers::pi;
        const auto exact=[&](Point2D x) {return std::sin(pi*x.x)*std::sin(pi*x.y)*std::exp(verification=="decay"?-2*pi*pi*diffusivity*time:0.);};
        if(verification.empty()) {
            const auto state=carrier(flowPath,mesh); p.volumeFlux=state.flux; carrierTime=state.time;
            bc=boundaries(bcPath,mesh);
            p.boundary=[&](std::size_t id,const fv::Face&){return bc[id];};
            p.source=[&](Point2D){return source;};
        } else {
            if(verification=="decay") {
                double area=0;
                for(const auto& cell:mesh.cells) {
                    require(cell.centre.x>0&&cell.centre.x<1&&cell.centre.y>0&&cell.centre.y<1,"decay requires interior unit-square cells");
                    area+=cell.area;
                }
                require(std::abs(area-1)<1e-10,"decay requires a full unit-square domain");
                for(const auto& f:mesh.faces) if(!f.neighbour)
                    require(std::abs(f.centre.x)<1e-12||std::abs(f.centre.x-1)<1e-12||std::abs(f.centre.y)<1e-12||std::abs(f.centre.y-1)<1e-12,
                        "decay verification requires the unit square boundary");
            }
            p.volumeFlux.resize(mesh.faces.size());
            for(std::size_t id=0;id<mesh.faces.size();++id) p.volumeFlux[id]=verification=="decay"?0:speed*mesh.faces[id].areaVector.x;
            p.boundary=[&](std::size_t,const fv::Face& f){return fv::ScalarBoundary2D{fv::ScalarBoundaryKind2D::Value,exact(f.centre),{}};};
            p.source=[&](Point2D x){return verification=="decay"?0:2*diffusivity*pi*pi*exact(x)+speed*pi*std::cos(pi*x.x)*std::sin(pi*x.y);};
        }
        const auto parent=std::filesystem::path(prefix).parent_path();if(!parent.empty())std::filesystem::create_directories(parent);
        {auto pending=output(prefix,".json");pending<<"{\"status\":\"running\",\"converged\":false}\n";} started=true;
        auto history=output(prefix,".history.csv");history<<"step,time,iteration,linearIterations,residualNorm,relativeResidual,maxCellImbalance,maxDiagonalScaledImbalance\n";
        std::vector<double> previous;
        if(dt>0) {previous.assign(mesh.cells.size(),initial);if(verification=="decay") for(std::size_t i=0;i<previous.size();++i)previous[i]=exact(mesh.cells[i].centre);}
        fv::ScalarTransportResult2D result;
        for(std::size_t step=1;step<=steps;++step) {
            time=dt*static_cast<double>(step); require(std::isfinite(time),"physical time overflow");
            result=fv::solveScalarTransport2D(mesh,p,controls,previous,dt);
            for(const auto& h:result.history)history<<step<<','<<time<<','<<h.iteration<<','<<h.linearIterations<<','<<h.residualNorm<<','<<h.relativeResidual<<','<<h.maxCellImbalance<<','<<h.maxDiagonalScaledImbalance<<'\n';
            if(!result.converged)break;
            if(step<steps)previous=result.values;
        }
        auto cells=output(prefix,".cells.csv");cells<<"cell,x,y,area,value,previous,sourceIntegral,temporalIntegral,exact\n";
        double error2=0,area=0;
        for(std::size_t i=0;i<mesh.cells.size();++i) {
            const auto& cell=mesh.cells[i];
            cells<<i<<','<<cell.centre.x<<','<<cell.centre.y<<','<<cell.area<<','<<result.values[i]<<',';
            if(dt>0)cells<<previous[i];else cells<<"nan";
            cells<<','<<result.sourceIntegrals[i]<<','<<result.temporalIntegrals[i]<<',';
            if(!verification.empty()) {const double ex=exact(cell.centre);cells<<ex;error2+=cell.area*std::pow(result.values[i]-ex,2);area+=cell.area;}else cells<<"nan";
            cells<<'\n';
        }
        auto faces=output(prefix,".faces.csv");faces<<"face,owner,neighbour,volumeFlux,advectiveFlux,diffusiveFlux\n";
        for(std::size_t id=0;id<mesh.faces.size();++id) {
            const auto& f=mesh.faces[id]; faces<<id<<','<<f.owner<<',';
            if(f.neighbour)faces<<*f.neighbour;else faces<<-1;
            faces<<','<<p.volumeFlux[id]<<','<<result.advectiveFlux[id]<<','<<result.diffusiveFlux[id]<<'\n';
        }
        std::string error;require(writeLegacyVtk2D(read.topology,prefix+".vtk",&error),error);
        std::ofstream vtk(prefix+".vtk",std::ios::app);vtk.exceptions(std::ios::badbit|std::ios::failbit);vtk<<std::setprecision(17)<<"SCALARS scalar double 1\nLOOKUP_TABLE default\n";
        for(double v:result.values)vtk<<v<<'\n';
        const auto& h=result.history.back();auto json=output(prefix,".json");
        json<<"{\n\"format\":\"cartmesh2d-scalar-transport-v1\",\n\"status\":"<<quote(result.converged?"converged":"correction-limit")
            <<",\n\"converged\":"<<(result.converged?"true":"false")<<",\n\"mesh\":"<<quote(meshPath)
            <<",\n\"carrierCheckpoint\":"<<quote(flowPath)<<",\n\"carrierTime\":"<<carrierTime
            <<",\n\"boundaryFile\":"<<quote(bcPath)<<",\n\"verification\":"<<quote(verification)
            <<",\n\"scope\":\"constant-property passive scalar; frozen carrier flux; no buoyancy or thermal feedback\""
            <<",\n\"cells\":"<<mesh.cells.size()<<",\n\"faces\":"<<mesh.faces.size()<<",\n\"time\":"<<time<<",\n\"timeStep\":"<<dt
            <<",\n\"verificationSpeed\":"<<speed<<",\n\"constantSource\":"<<source<<",\n\"initialValue\":"<<initial
            <<",\n\"diffusivity\":"<<diffusivity<<",\n\"convection\":"<<quote(controls.convection==fv::ConvectionScheme2D::Upwind?"upwind":"limited-linear")
            <<",\n\"relativeTolerance\":"<<controls.relativeTolerance<<",\n\"absoluteTolerance\":"<<controls.absoluteTolerance<<",\n\"cellTolerance\":"<<controls.cellTolerance
            <<",\n\"residualNorm\":"<<h.residualNorm<<",\n\"maxCellImbalance\":"<<h.maxCellImbalance<<",\n\"maxDiagonalScaledImbalance\":"<<h.maxDiagonalScaledImbalance
            <<",\n\"boundaryFlux\":"<<result.boundaryFlux<<",\n\"sourceIntegral\":"<<result.sourceIntegral<<",\n\"temporalIntegral\":"<<result.temporalIntegral
            <<",\n\"globalBalance\":"<<result.globalBalance<<",\n\"maxCarrierImbalance\":"<<result.maxCarrierImbalance
            <<",\n\"minValue\":"<<result.minValue<<",\n\"maxValue\":"<<result.maxValue<<",\n\"maxCourant\":"<<result.maxCourant<<",\n\"l2Error\":";
        if(!verification.empty())json<<std::sqrt(error2/area);else json<<"null";
        json<<",\n\"nativeTopologyRevalidated\":true,\n\"solverQualityPassed\":true,\n\"externalCheckMesh\":\"not run\"\n}\n";
        cells.close();faces.close();vtk.close();history.close();json.close();
        std::cout<<(result.converged?"Converged":"Correction limit")<<": "<<mesh.cells.size()<<" cells, time="<<time<<", scalar range=["<<result.minValue<<','<<result.maxValue<<"], balance="<<result.globalBalance<<'\n';
        return result.converged?0:2;
    } catch(const std::exception& e) {
        if(started)try {auto out=output(prefix,".json");out<<"{\"status\":\"failed\",\"converged\":false,\"error\":"<<quote(e.what())<<"}\n";}catch(...){}
        std::cerr<<"cartmesh2d_transport_cli: "<<e.what()<<'\n';return 1;
    }
}
