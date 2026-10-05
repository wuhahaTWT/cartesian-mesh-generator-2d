// Native solver accuracy/cost probe; reuses the existing native mesh fixture.
// No independent discretized equations or acceptance thresholds are introduced.
#define main inherited_boundary_test_main
#include "../../tests/flow_boundary_test.cpp"
#undef main
#include "cartmesh2d/fv/ManufacturedFlow2D.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
int main(int argc,char** argv) {
    if(argc!=7){std::cerr<<"case n tolerance scheme method output-prefix\n";return 1;}
    const std::string kind=argv[1],scheme=argv[4],method=argv[5],prefix=argv[6];
    const int n=std::stoi(argv[2]);const double length=kind=="channel"?4:1;
    const auto mesh=rectangle(kind=="channel"?4*n:n,n,length,kind=="manufactured");
    FlowControls2D c;c.scenario=kind;c.nu=kind=="cavity"?.01:.1;c.speed=1;
    c.tolerance=std::stod(argv[3]);c.maxIterations=3000;c.profile=true;
    if(scheme=="face-limited-linear")c.convection=ConvectionScheme2D::FaceLimitedLinearUpwind;
    else if(scheme!="upwind")throw std::runtime_error("unknown scheme");
    if(method=="none")c.steadyAcceleration=SteadyAcceleration2D::None;
    else if(method!="default")throw std::runtime_error("unknown method");
    const auto start=std::chrono::steady_clock::now();
    std::ofstream record(prefix+".json");record<<std::setprecision(17);
    try {
        const auto r=solveIncompressible2D(mesh,c);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream field(prefix+".cells.csv");field<<std::setprecision(17)<<"id,x,y,area,u,v,p\n";
        double area=0,velocity2=0,pressure2=0,velocityMax=0,pressureMax=0,pOffset=0;
        if(kind=="manufactured")for(const auto& cell:mesh.cells){const auto exact=manufacturedFlow2D(cell.centre,c.speed,c.nu);pOffset+=cell.area*(r.p[&cell-mesh.cells.data()]-exact.pressure);area+=cell.area;}
        if(area>0)pOffset/=area;area=0;
        for(const auto& cell:mesh.cells){const auto i=&cell-mesh.cells.data();field<<i<<','<<cell.centre.x<<','<<cell.centre.y<<','<<cell.area<<','<<r.u[i]<<','<<r.v[i]<<','<<r.p[i]<<'\n';
            if(kind=="channel"||kind=="manufactured"){
                const auto exact=kind=="manufactured"?manufacturedFlow2D(cell.centre,c.speed,c.nu):ManufacturedFlowSample2D{{4*cell.centre.y*(1-cell.centre.y),0},8*c.nu*(length-cell.centre.x),{}};
                const double ve=std::hypot(r.u[i]-exact.velocity.x,r.v[i]-exact.velocity.y),pe=std::abs(r.p[i]-exact.pressure-pOffset);
                velocity2+=cell.area*ve*ve;pressure2+=cell.area*pe*pe;velocityMax=std::max(velocityMax,ve);pressureMax=std::max(pressureMax,pe);
            }area+=cell.area;
        }
        const auto& last=r.history.back();
        record<<"{\"case\":"<<std::quoted(kind)<<",\"n\":"<<n<<",\"cells\":"<<mesh.cells.size()<<",\"nu\":"<<c.nu<<",\"tolerance\":"<<c.tolerance<<",\"scheme\":"<<std::quoted(scheme)<<",\"method\":"<<std::quoted(method)
            <<",\"converged\":"<<(r.converged?"true":"false")<<",\"evaluations\":"<<last.iteration<<",\"seconds\":"<<seconds<<",\"continuity\":"<<last.continuity<<",\"momentumResidual\":"<<last.momentumResidual<<",\"velocityChange\":"<<last.velocityChange<<",\"pressureChange\":"<<last.pressureChange
            ;
        if(kind!="cavity") record<<",\"velocityRms_m_s\":"<<std::sqrt(velocity2/area)<<",\"velocityMax_m_s\":"<<velocityMax<<",\"pressureRms_m2_s2\":"<<std::sqrt(pressure2/area)<<",\"pressureMax_m2_s2\":"<<pressureMax<<",\"pressureGaugeOffset_m2_s2\":"<<pOffset;
        record<<"}\n";
        return r.converged?0:2;
    } catch(const std::exception& e){record<<"{\"case\":"<<std::quoted(kind)<<",\"n\":"<<n<<",\"error\":"<<std::quoted(e.what())<<"}\n";return 1;}
}
