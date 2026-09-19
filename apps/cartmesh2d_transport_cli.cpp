#include "cartmesh2d/fv/ScalarTransport2D.hpp"
#include "cartmesh2d/fv/ThermalCheckpoint2D.hpp"
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
    std::getline(in,line);
    require(line=="CARTMESH2D_FLOW_CHECKPOINT 1"||line=="CARTMESH2D_FLOW_CHECKPOINT 2","invalid flow checkpoint header");
    const bool v2=line.ends_with(" 2");
    std::getline(in,line); require(line=="DISCRETIZATION Euler-RC-v2","unsupported flow checkpoint discretization");
    std::getline(in,line); std::istringstream config(line); std::string token,scheme,stress;
    std::string backflow;
    require(bool(config>>token>>std::quoted(c.scenario)>>c.nu>>c.speed>>scheme>>stress>>c.manufacturedPressureSlope)
        && (!v2 || bool(config>>backflow)) && token=="CONFIG",
        "invalid flow checkpoint configuration");
    require(scheme=="upwind"||scheme=="limited-linear","invalid carrier convection");
    require(stress=="symmetric"||stress=="laplacian","invalid carrier stress");
    c.convection=scheme=="upwind"?fv::ConvectionScheme2D::Upwind:fv::ConvectionScheme2D::LimitedLinearUpwind;
    c.viscousStress=stress=="symmetric"?fv::ViscousStress2D::Symmetric:fv::ViscousStress2D::Laplacian;
    if (v2) {
        require(backflow=="reject"||backflow=="normal-inlet","invalid carrier outlet backflow model");
        c.outletBackflow=backflow=="normal-inlet"?fv::OutletBackflow2D::NormalInlet:fv::OutletBackflow2D::Reject;
    }
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
    std::string prefix; bool started=false; double acceptedTime=0;
    try {
        std::string meshPath,flowPath,bcPath,verification,evolve,restart;
        fv::FlowControls2D flowControls; flowControls.tolerance=1e-8;
        bool explicitVelocityRelaxation=false;
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
                    "Evolving flow: --evolve-flow external|channel|cavity --boundary BC.csv --dt DT --steps N\n"
                    "  --flow-nu .01 --flow-speed 1 --flow-tolerance 1e-8 --flow-max-iterations 1500\n"
                    "  --flow-velocity-relaxation 0.6: evolving carrier only; (0,1], larger may be unstable.\n"
                    "  --flow-convection upwind|limited-linear --pressure-preconditioner ic0|aggregation\n"
                    "  --outlet-backflow reject|normal-inlet (default reject)\n"
                    "  --restart PREFIX.thermal.checkpoint: resume both fields, same physical setup.\n"
                    "  --verification thermal-vortex: analytic evolving vortex/scalar decay on unit square.\n"
                    "Outputs .json .vtk .cells.csv .faces.csv .history.csv; evolving mode adds joint checkpoint.\n";
                return 0;
            }
            require(i+1<argc,"missing option value"); const std::string value=argv[++i];
            if(arg=="--mesh")meshPath=value; else if(arg=="--output")prefix=value;
            else if(arg=="--flow-checkpoint")flowPath=value; else if(arg=="--boundary")bcPath=value;
            else if(arg=="--verification")verification=value;
            else if(arg=="--evolve-flow")evolve=value;
            else if(arg=="--restart")restart=value;
            else if(arg=="--flow-nu")flowControls.nu=number(value);
            else if(arg=="--flow-speed")flowControls.speed=number(value);
            else if(arg=="--flow-tolerance")flowControls.tolerance=number(value);
            else if(arg=="--flow-velocity-relaxation") {
                flowControls.velocityRelaxation=number(value);
                require(flowControls.velocityRelaxation>0&&flowControls.velocityRelaxation<=1,
                        "flow-velocity-relaxation must be in (0,1]");
                explicitVelocityRelaxation=true;
            }
            else if(arg=="--flow-max-iterations") {
                const double n=number(value);require(n>=1&&n<=100000&&n==std::floor(n),"invalid flow iteration count");
                flowControls.maxIterations=static_cast<std::size_t>(n);
            } else if(arg=="--flow-convection") {
                require(value=="upwind"||value=="limited-linear","unknown flow convection");
                flowControls.convection=value=="upwind"?fv::ConvectionScheme2D::Upwind:fv::ConvectionScheme2D::LimitedLinearUpwind;
            } else if(arg=="--outlet-backflow") {
                require(value=="reject"||value=="normal-inlet","unknown outlet backflow model");
                flowControls.outletBackflow=value=="normal-inlet"?fv::OutletBackflow2D::NormalInlet:fv::OutletBackflow2D::Reject;
            } else if(arg=="--pressure-preconditioner") {
                require(value=="ic0"||value=="aggregation","unknown pressure preconditioner");
                flowControls.pressurePreconditioner=value=="ic0"?fv::PressurePreconditioner2D::IncompleteCholesky0:fv::PressurePreconditioner2D::Aggregation;
            }
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
        require(verification.empty()||verification=="sine"||verification=="decay"||verification=="thermal-vortex","unknown verification");
        if(verification=="thermal-vortex") {
            require(evolve.empty()&&flowPath.empty()&&bcPath.empty(),"thermal-vortex defines its own carrier and boundaries");
            evolve="taylor-green";
        }
        const bool evolving=!evolve.empty();
        require(!explicitVelocityRelaxation||evolving,"flow-velocity-relaxation option requires evolving flow");
        if(evolving) {
            require(dt>0&&flowPath.empty()&&(verification.empty()||verification=="thermal-vortex"),"evolving flow requires dt and cannot use a frozen carrier or other verification");
            require(evolve=="channel"||evolve=="cavity"||evolve=="external"||evolve=="taylor-green","unsupported evolving flow case");
            require(evolve!="taylor-green"||verification=="thermal-vortex","use thermal-vortex verification for Taylor-Green");
            flowControls.scenario=evolve;
        }
        require(restart.empty()||evolving,"joint restart requires evolving flow");
        require(verification.empty()?((!flowPath.empty()||evolving)&&!bcPath.empty()):(flowPath.empty()&&bcPath.empty()),"choose explicit carrier/boundary files OR verification");
        require(verification!="sine"||dt==0,"sine verification is steady");
        require(verification!="decay"||dt>0,"decay verification needs dt");
        auto read=readCm2dTopology(meshPath); require(read.valid(),read.error);
        const auto mesh=fv::makeFvMesh2D(read.topology);
        fv::ScalarTransportProblem2D p; p.diffusivity=diffusivity;
        std::vector<fv::ScalarBoundary2D> bc;
        double carrierTime=0,time=0;
        fv::ThermalSetup2D thermalSetup; fv::ThermalFlowState2D state;
        const double pi=std::numbers::pi;
        const auto exact=[&](Point2D x) {return std::sin(pi*x.x)*std::sin(pi*x.y)*std::exp(verification!="sine"?-2*pi*pi*diffusivity*time:0.);};
        if(verification.empty()) {
            if(!evolving) {const auto frozen=carrier(flowPath,mesh); p.volumeFlux=frozen.flux; carrierTime=frozen.time;}
            bc=boundaries(bcPath,mesh);
            p.boundary=[&](std::size_t id,const fv::Face&){return bc[id];};
            p.source=[&](Point2D){return source;};
        } else {
            if(verification=="decay"||verification=="thermal-vortex") {
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
            for(std::size_t id=0;id<mesh.faces.size();++id) p.volumeFlux[id]=verification!="sine"?0:speed*mesh.faces[id].areaVector.x;
            p.boundary=[&](std::size_t,const fv::Face& f){return fv::ScalarBoundary2D{fv::ScalarBoundaryKind2D::Value,exact(f.centre),{}};};
            p.source=[&](Point2D x){return verification!="sine"?0:2*diffusivity*pi*pi*exact(x)+speed*pi*std::cos(pi*x.x)*std::sin(pi*x.y);};
        }
        if(evolving) {
            thermalSetup.diffusivity=diffusivity;
            thermalSetup.sourceDensity.assign(mesh.cells.size(),verification.empty()?source:0.);
            thermalSetup.boundary=verification.empty()?bc:std::vector<fv::ScalarBoundary2D>(mesh.faces.size());
            fv::validateThermalSetup2D(mesh,thermalSetup);
            if(restart.empty()) {
                state.flow=fv::initialIncompressibleState2D(mesh,flowControls);
                state.scalar.assign(mesh.cells.size(),initial);
                if(verification=="thermal-vortex")for(std::size_t i=0;i<mesh.cells.size();++i)state.scalar[i]=exact(mesh.cells[i].centre);
            } else {
                std::ifstream in(restart);require(bool(in),"cannot open joint restart");
                state=fv::readThermalCheckpoint2D(in,mesh,flowControls,thermalSetup,controls);
            }
            acceptedTime=state.flow.time;
        }
        const auto parent=std::filesystem::path(prefix).parent_path();if(!parent.empty())std::filesystem::create_directories(parent);
        {auto pending=output(prefix,".json");pending<<"{\"status\":\"running\",\"converged\":false}\n";} started=true;
        auto history=output(prefix,".history.csv");history<<"step,time,iteration,linearIterations,residualNorm,relativeResidual,maxCellImbalance,maxDiagonalScaledImbalance\n";
        std::vector<double> previous;
        if(dt>0) {previous.assign(mesh.cells.size(),initial);if(verification=="decay") for(std::size_t i=0;i<previous.size();++i)previous[i]=exact(mesh.cells[i].centre);}
        fv::ScalarTransportResult2D result;
        std::ofstream thermalHistory;
        const auto saveAccepted=[&](const fv::ThermalFlowState2D& accepted) {
            auto out=output(prefix,".thermal.checkpoint.tmp");
            fv::writeThermalCheckpoint2D(out,mesh,flowControls,thermalSetup,controls,accepted);
            out.close();std::filesystem::rename(prefix+".thermal.checkpoint.tmp",prefix+".thermal.checkpoint");
        };
        if(evolving) {
            saveAccepted(state);
            thermalHistory=output(prefix,".thermal-history.csv");
            thermalHistory<<"step,time,accepted,flowIterations,flowMomentumResidual,flowContinuity,scalarIterations,scalarResidual,heatContent,scalarGlobalBalance,maxCourant\n";
        }
        for(std::size_t step=1;step<=steps;++step) {
            time=evolving?state.flow.time+dt:dt*static_cast<double>(step);
            require(std::isfinite(time),"physical time overflow");
            if(evolving) {
                previous=state.scalar;
                const auto attempt=fv::advanceThermalFlow2D(mesh,flowControls,thermalSetup,controls,state,dt);
                result=attempt.scalar;
                double heat=0;
                for(std::size_t i=0;i<result.values.size();++i)heat+=mesh.cells[i].area*result.values[i];
                if(attempt.accepted) {
                    saveAccepted(*attempt.accepted);state=*attempt.accepted;acceptedTime=state.flow.time;
                }
                const auto& fh=attempt.flow.history.back();
                thermalHistory<<step<<','<<time<<','<<(attempt.accepted?1:0)<<','<<attempt.flow.history.size()<<','<<fh.momentumResidual<<','<<fh.continuity
                    <<','<<result.history.size()<<',';
                if(result.history.empty())thermalHistory<<"nan";else thermalHistory<<result.history.back().residualNorm;
                thermalHistory<<','<<heat<<','<<result.globalBalance<<','<<attempt.flow.maxCourant<<'\n';thermalHistory.flush();
                if(!attempt.accepted) {
                    auto failed=output(prefix,".json");
                    failed<<"{\"status\":\"step-not-accepted\",\"converged\":false,\"failedStage\":"
                        <<quote(attempt.flow.converged?"scalar":"flow")<<",\"acceptedTime\":"<<acceptedTime<<",\"attemptedTime\":"<<time<<"}\n";
                    return 2;
                }
                std::cout<<std::setprecision(17)<<"{\"type\":\"thermal-time-step\",\"accepted\":1,\"step\":"<<step
                    <<",\"time\":"<<time<<",\"flowIterations\":"<<attempt.flow.history.size()
                    <<",\"momentumResidual\":"<<fh.momentumResidual<<",\"continuity\":"<<fh.continuity
                    <<",\"scalarIterations\":"<<result.history.size()<<",\"scalarResidual\":"<<result.history.back().residualNorm
                    <<",\"heatContent\":"<<heat<<",\"globalBalance\":"<<result.globalBalance
                    <<",\"maxCourant\":"<<attempt.flow.maxCourant<<"}\n"<<std::flush;
                p.volumeFlux=state.flow.flux;carrierTime=state.flow.time;
            } else result=fv::solveScalarTransport2D(mesh,p,controls,previous,dt);
            for(const auto& h:result.history)history<<step<<','<<time<<','<<h.iteration<<','<<h.linearIterations<<','<<h.residualNorm<<','<<h.relativeResidual<<','<<h.maxCellImbalance<<','<<h.maxDiagonalScaledImbalance<<'\n';
            history.flush();
            if(!result.converged)break;
            if(!evolving&&step<steps)previous=result.values;
        }
        if(evolving) {
            // Diagnostic carrier export is not the authoritative joint restart.
            flowPath=prefix+".carrier.checkpoint";
            auto out=output(prefix,".carrier.checkpoint");fv::writeFlowCheckpoint2D(out,mesh,flowControls,state.flow);out.close();
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
            <<",\n\"scope\":"<<quote(evolving?"synchronized new-time carrier and passive scalar; no buoyancy or thermal feedback":"constant-property passive scalar; frozen carrier flux; no buoyancy or thermal feedback")
            <<",\n\"evolvingFlow\":"<<(evolving?"true":"false")<<",\n\"acceptedTime\":"<<acceptedTime
            <<",\n\"restart\":"<<quote(restart)<<",\n\"flowCase\":"<<quote(evolve)
            <<",\n\"cells\":"<<mesh.cells.size()<<",\n\"faces\":"<<mesh.faces.size()<<",\n\"time\":"<<time<<",\n\"timeStep\":"<<dt
            <<",\n\"verificationSpeed\":"<<speed<<",\n\"constantSource\":"<<source<<",\n\"initialValue\":"<<initial
            <<",\n\"flowNu\":"<<flowControls.nu<<",\n\"flowSpeed\":"<<flowControls.speed
            <<",\n\"flowVelocityRelaxation\":"<<flowControls.velocityRelaxation
            <<",\n\"flowConvection\":"<<quote(flowControls.convection==fv::ConvectionScheme2D::Upwind?"upwind":"limited-linear")
            <<",\n\"outletBackflow\":"<<quote(flowControls.outletBackflow==fv::OutletBackflow2D::NormalInlet?"normal-inlet":"reject")
            <<",\n\"steps\":"<<steps
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
        if(started)try {auto out=output(prefix,".json");out<<"{\"status\":\"failed\",\"converged\":false,\"acceptedTime\":"<<acceptedTime<<",\"error\":"<<quote(e.what())<<"}\n";}catch(...){}
        std::cerr<<"cartmesh2d_transport_cli: "<<e.what()<<'\n';return 1;
    }
}
