// Native branch-certificate probe for ordinary internal, manufactured and closed flows.
// Reuses the product mesh fixture and product solver; it introduces no
// independent discretized equations or acceptance thresholds.
#define main inherited_boundary_test_main
#include "../../tests/flow_boundary_test.cpp"
#undef main

#include <chrono>
#include <fstream>
#include <iomanip>

namespace {
using namespace cartmesh2d;
using namespace cartmesh2d::fv;

struct FieldDifference {
    double velocityRms = 0;
    double velocityMaximum = 0;
    double pressureRms = 0;
    double pressureMaximum = 0;
    double pressureGaugeOffset = 0;
};

FieldDifference compareFields(const FvMesh2D& mesh,const FlowResult2D& a,const FlowResult2D& b) {
    FieldDifference d;double area=0;
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        area+=mesh.cells[i].area;
        d.pressureGaugeOffset+=mesh.cells[i].area*(a.p[i]-b.p[i]);
    }
    d.pressureGaugeOffset/=area;
    double velocity2=0,pressure2=0;
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        const double du=a.u[i]-b.u[i],dv=a.v[i]-b.v[i];
        const double velocity=std::hypot(du,dv);
        const double pressure=std::abs(a.p[i]-b.p[i]-d.pressureGaugeOffset);
        velocity2+=mesh.cells[i].area*velocity*velocity;
        pressure2+=mesh.cells[i].area*pressure*pressure;
        d.velocityMaximum=std::max(d.velocityMaximum,velocity);
        d.pressureMaximum=std::max(d.pressureMaximum,pressure);
    }
    d.velocityRms=std::sqrt(velocity2/area);
    d.pressureRms=std::sqrt(pressure2/area);
    return d;
}

void writeField(const std::string& path,const FvMesh2D& mesh,const FlowResult2D& r) {
    std::ofstream out(path);out<<std::setprecision(17)<<"id,x,y,area,u,v,p\n";
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        const auto& cell=mesh.cells[i];
        out<<i<<','<<cell.centre.x<<','<<cell.centre.y<<','<<cell.area<<','
           <<r.u[i]<<','<<r.v[i]<<','<<r.p[i]<<'\n';
    }
}

void writeRun(std::ostream& out,const char* name,const FlowControls2D& c,
    const FlowResult2D& r,double seconds) {
    const auto& last=r.history.back();
    out<<std::quoted(name)<<":{\"nu\":"<<c.nu<<",\"tolerance\":"<<c.tolerance
       <<",\"converged\":"<<(r.converged?"true":"false")
       <<",\"evaluations\":"<<last.iteration<<",\"seconds\":"<<seconds
       <<",\"continuity\":"<<last.continuity
       <<",\"momentumResidual\":"<<last.momentumResidual
       <<",\"velocityChange\":"<<last.velocityChange
       <<",\"pressureChange\":"<<last.pressureChange
       <<",\"maximumSpeedRatio\":"<<r.fieldAmplitude.maximumSpeedRatio
       <<",\"pressureRangeRatio\":"<<r.fieldAmplitude.pressureRangeRatio<<'}';
}
}

int main(int argc,char** argv) {
    if(argc!=5) {
        std::cerr<<"case n target-tolerance output-prefix\n";
        return 1;
    }
    const std::string kind=argv[1],prefix=argv[4];
    if(kind!="channel"&&kind!="manufactured"&&kind!="cavity")
        throw std::runtime_error("case must be channel, manufactured or cavity");
    const int n=std::stoi(argv[2]);
    const auto mesh=rectangle(kind=="channel"?4*n:n,n,kind=="channel"?4:1,
        kind=="manufactured");
    FlowControls2D target;target.scenario=kind;target.nu=kind=="cavity"?.01:.1;
    target.speed=1;target.tolerance=std::stod(argv[3]);target.maxIterations=4000;target.profile=true;
    auto guide=target;guide.nu=10*target.nu;guide.tolerance=std::sqrt(target.tolerance);
    std::ofstream record(prefix+".json");record<<std::setprecision(17);
    try {
        const auto timed=[&](const auto& solve) {
            const auto start=std::chrono::steady_clock::now();
            auto result=solve();
            return std::pair{std::move(result),std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start).count()};
        };
        auto direct=timed([&]{return solveIncompressible2D(mesh,target);});
        auto guided=timed([&]{return solveIncompressible2D(mesh,guide);});
        FlowInitialGuess2D initial{guided.first.u,guided.first.v,guided.first.p,guided.first.flux};
        auto certified=timed([&]{return solveIncompressibleFromGuess2D(mesh,target,initial);});
        const auto difference=compareFields(mesh,direct.first,certified.first);
        writeField(prefix+"-direct.cells.csv",mesh,direct.first);
        writeField(prefix+"-guide.cells.csv",mesh,guided.first);
        writeField(prefix+"-certified.cells.csv",mesh,certified.first);
        record<<"{\"schema\":1,\"case\":"<<std::quoted(kind)<<",\"n\":"<<n
              <<",\"cells\":"<<mesh.cells.size()<<",\"defaultMethod\":\"newton-krylov\""
              <<",\"defaultConvection\":\"face-limited-linear\",\"runs\":{";
        writeRun(record,"direct",target,direct.first,direct.second);record<<',';
        writeRun(record,"guide",guide,guided.first,guided.second);record<<',';
        writeRun(record,"certifiedTarget",target,certified.first,certified.second);
        record<<"},\"directVsCertifiedTarget\":{\"velocityRms_m_s\":"<<difference.velocityRms
              <<",\"velocityMaximum_m_s\":"<<difference.velocityMaximum
              <<",\"pressureRms_m2_s2\":"<<difference.pressureRms
              <<",\"pressureMaximum_m2_s2\":"<<difference.pressureMaximum
              <<",\"pressureGaugeOffset_m2_s2\":"<<difference.pressureGaugeOffset
              <<"},\"evaluationCostRatio\":"
              <<double(guided.first.history.back().iteration+certified.first.history.back().iteration)
                    /direct.first.history.back().iteration<<"}\n";
        return direct.first.converged&&guided.first.converged&&certified.first.converged?0:2;
    } catch(const std::exception& e) {
        record<<"{\"schema\":1,\"case\":"<<std::quoted(kind)
              <<",\"error\":"<<std::quoted(e.what())<<"}\n";
        return 1;
    }
}
