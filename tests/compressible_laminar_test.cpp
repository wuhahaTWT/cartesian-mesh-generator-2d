#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
struct Evidence {double temperature=0,velocity=0,pressure=0,heat=0,work=0,rate=0,balance=0;std::size_t steps=0,rejected=0;};
Evidence couette(int n,const std::filesystem::path& output) {
    const auto mesh=fv_test::rectangle(2,n,1);const IdealGas2D gas{1.4,1};
    const EulerTransport2D transport{.1,.1};const double speed=.5,wallTemperature=1,endTime=60;
    std::vector<EulerBoundary2D> bc;
    for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour) {
        const auto& f=mesh.faces[id];EulerBoundary2D b{id,EulerBoundaryKind2D::NoSlipWall,{}, {},"lower"};
        if(f.areaVector.x!=0) {
            b.kind=EulerBoundaryKind2D::Periodic;b.name="periodic-x";
            for(std::size_t j=0;j<mesh.faces.size();++j) {
                const auto& g=mesh.faces[j];
                if(!g.neighbour&&g.areaVector.x==-f.areaVector.x&&g.centre.y==f.centre.y)b.partner=j;
            }
        }else {
            b.thermalKind=HeatBoundaryKind2D::Temperature;b.thermalValue=wallTemperature;
            if(f.areaVector.y>0){b.name="upper";b.wallVelocity={speed,0};}
        }
        bc.push_back(b);
    }
    EulerState2D state{0,0,std::vector<EulerConservative2D>(mesh.cells.size(),eulerConservative2D({1,0,0,1},gas))};
    EulerStepControls2D controls;controls.fluxScheme=EulerFluxScheme2D::Hllc;controls.order=2;
    controls.wallGradient=WallGradient2D::Quadratic;controls.endTime=endTime;
    const EulerStepper2D solver(mesh,bc,gas,transport,controls.wallGradient);
    Evidence e;EulerStepResult2D last;double integratedBoundaryEnergy=0;
    // Cold-start test uses the actual coupled conservative solver. Analytic
    // steady Couette is compared only after physical evolution, never installed
    // as the accepted state. No duplicate Python equation audit is involved.
    while(state.time<endTime) {
        last=solver.advance(state,controls);integratedBoundaryEnergy+=last.step*last.boundaryFlux[3];
        e.rejected+=last.rejectedCandidates;state=last.state;
        require(state.steps<200000,"Couette step budget exhausted");
    }
    const double a=transport.dynamicViscosity*speed*speed/(2*transport.thermalConductivity);
    const double peakRise=a/4,mid=wallTemperature+peakRise;
    const double inverseTemperatureIntegral=2/std::sqrt(a*mid)*std::atanh(.5*std::sqrt(a/mid));
    const double exactPressure=gas.gasConstant/inverseTemperatureIntegral;
    const double wallPower=transport.dynamicViscosity*speed*speed;
    double totalEnergy=0;
    for(std::size_t i=0;i<state.cells.size();++i) {
        const auto q=eulerPrimitive2D(state.cells[i],gas);const double y=mesh.cells[i].centre.y;
        e.temperature=std::max(e.temperature,std::abs(q.pressure/q.density-(wallTemperature+a*y*(1-y)))/peakRise);
        e.velocity=std::max(e.velocity,std::hypot(q.u-speed*y,q.v)/speed);
        e.pressure=std::max(e.pressure,std::abs(q.pressure/exactPressure-1));
        totalEnergy+=mesh.cells[i].area*state.cells[i][3];
        for(std::size_t k=0;k<4;++k)e.rate=std::max(e.rate,std::abs(state.cells[i][k]-last.previousCells[i][k])/last.step);
    }
    e.heat=std::abs(last.boundaryHeat/wallPower-1);e.work=std::abs(last.boundaryViscousWork/wallPower+1);
    e.balance=std::abs(totalEnergy-1/(gas.gamma-1)+integratedBoundaryEnergy)/(totalEnergy+1/(gas.gamma-1)+std::abs(integratedBoundaryEnergy));
    e.steps=state.steps;
    // Retain the previous Couette physical targets: 1e-6 of the analytic
    // temperature rise / wall speed / wall power; pressure is judged by
    // refinement separately. The 8192 eps budget checks accumulated arithmetic.
    require(e.temperature<1e-6&&e.velocity<1e-6,"cold Couette did not reach the analytic velocity/temperature target");
    require(e.heat<1e-6&&e.work<1e-6,"cold Couette wall heat/work accuracy failed");
    require(e.rate<1e-7,"Couette still evolving at the requested physical endpoint");
    require(e.balance<8192*std::numeric_limits<double>::epsilon(),"Couette energy budget failed");
    if(!output.empty()) {
        std::ofstream field(output/("couette-"+std::to_string(n)+".csv"));field.exceptions(std::ios::badbit|std::ios::failbit);
        field<<std::setprecision(17)<<"x,y,rho,u,v,p,T,T_exact,u_exact\n";
        for(std::size_t i=0;i<state.cells.size();++i){const auto q=eulerPrimitive2D(state.cells[i],gas);const auto c=mesh.cells[i].centre;field<<c.x<<','<<c.y<<','<<q.density<<','<<q.u<<','<<q.v<<','<<q.pressure<<','<<q.pressure/q.density<<','<<wallTemperature+a*c.y*(1-c.y)<<','<<speed*c.y<<'\n';}
        std::ofstream checkpoint(output/("couette-"+std::to_string(n)+".checkpoint"));
        writeEulerCheckpoint2D(checkpoint,mesh,bc,gas,state,"couette-cold",transport);
    }
    // Rebuild the integrator from the physical checkpoint and verify the next
    // accepted step, including the new clock, is identical.
    std::stringstream saved;writeEulerCheckpoint2D(saved,mesh,bc,gas,state,"couette-cold",transport);
    const auto restored=readEulerCheckpoint2D(saved,mesh,bc,gas,"couette-cold",transport);
    auto continued=controls;continued.endTime=endTime+.001;
    const auto direct=solver.advance(state,continued),restart=EulerStepper2D(mesh,bc,gas,transport,controls.wallGradient).advance(restored,continued);
    require(direct.state.cells==restart.state.cells&&direct.state.time==restart.state.time&&direct.state.steps==restart.state.steps,"Couette restart changed the accepted trajectory");
    return e;
}
}
int main(int argc,char** argv){try {
    require(argc==1||argc==3,"usage: compressible_laminar_tests [--output NEW_DIRECTORY]");
    std::filesystem::path output;
    if(argc==3){require(std::string(argv[1])=="--output","unknown argument");output=argv[2];require(!std::filesystem::exists(output),"refusing to overwrite an existing evidence directory");std::filesystem::create_directories(output);}
    std::array<Evidence,3> evidence{};std::ostringstream report;report<<std::setprecision(17)<<"{\"case\":\"cold-couette\",\"endTime\":60,\"wallMach\":"<<.5/std::sqrt(1.4)<<",\"reynolds\":5,\"mu\":0.1,\"k\":0.1,\"wallSpeed\":0.5,\"wallTemperature\":1,\"gasR\":1,\"gamma\":1.4,\"levels\":[";
    for(std::size_t level=0;level<3;++level) {
        const int n=8*(1<<level);const auto e=couette(n,output);evidence[level]=e;
        if(level)report<<',';
        report<<"{\"ny\":"<<n<<",\"steps\":"<<e.steps<<",\"rejected\":"<<e.rejected<<",\"temperatureRiseRelativeLinf\":"<<e.temperature<<",\"wallSpeedRelativeLinf\":"<<e.velocity<<",\"pressureRelativeLinf\":"<<e.pressure<<",\"wallHeatRelativeError\":"<<e.heat<<",\"wallWorkRelativeError\":"<<e.work<<",\"finalConservativeRateLinf\":"<<e.rate<<",\"energyBudgetRelativeError\":"<<e.balance<<'}';
        std::cerr<<"cold Couette ny="<<n<<" accepted="<<e.steps<<" pressure error="<<e.pressure<<'\n';
    }
    const double order=std::log2(evidence[1].pressure/evidence[2].pressure);
    require(order>1.8&&evidence[2].pressure<.0001,"Couette pressure refinement target failed");
    report<<"],\"pressureOrder\":"<<order<<",\"checkpointNextStepIdentical\":true}\n";
    std::cout<<report.str();if(!output.empty()){std::ofstream file(output/"couette.json");file.exceptions(std::ios::badbit|std::ios::failbit);file<<report.str();}
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
