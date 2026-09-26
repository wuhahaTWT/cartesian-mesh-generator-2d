#include "cartmesh2d/immersed/CartesianFlow2D.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>

using namespace cartmesh2d;
namespace im=cartmesh2d::immersed;
namespace {
volatile std::sig_atomic_t interrupted=0;
void cancel(int) { interrupted=1; }
std::string jsonString(const std::string& value) {
    std::ostringstream s;s<<'"';
    for(char raw:value) {
        const auto c=static_cast<unsigned char>(raw);
        if(c=='"'||c=='\\')s<<'\\'<<c;
        else if(c<32)s<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<static_cast<unsigned>(c)<<std::dec;
        else s<<c;
    }
    s<<'"';return s.str();
}
double number(const std::string& s) { std::size_t n=0;const double x=std::stod(s,&n);if(n!=s.size()||!std::isfinite(x))throw std::runtime_error("invalid number: "+s);return x; }
std::size_t integer(const std::string& s) {
    if(s.empty()||s.find_first_not_of("0123456789")!=std::string::npos)throw std::runtime_error("invalid unsigned integer: "+s);
    std::size_t n=0;const auto x=std::stoull(s,&n);if(n!=s.size()||x>std::numeric_limits<std::size_t>::max())throw std::runtime_error("integer overflow");return static_cast<std::size_t>(x);
}
std::vector<BoundaryLoop> readSolids(const std::filesystem::path& path) {
    std::ifstream in(path);if(!in)throw std::runtime_error("cannot read XY boundary");
    std::vector<BoundaryLoop> loops;std::vector<Point2D> points;std::string line;
    const auto flush=[&](){if(!points.empty()){loops.emplace_back(points);points.clear();}};
    while(std::getline(in,line)) {
        const auto start=line.find_first_not_of(" \r\t");
        if(start==std::string::npos){flush();continue;}
        if(line[start]=='#') {
            // This entry has uniform stationary-wall physics. Named non-wall
            // product patches cannot silently acquire different boundary roles.
            if(line.find("cartmesh2d-loop",start)!=std::string::npos)
                throw std::runtime_error("immersed XY accepts plain solid loops only; named boundary metadata is not supported");
            continue;
        }
        std::istringstream row(line);Point2D p;std::string extra;
        if(!(row>>p.x>>p.y)||!std::isfinite(p.x)||!std::isfinite(p.y)||(row>>extra && extra[0]!='#'))throw std::runtime_error("invalid XY coordinate row");
        points.push_back(p);
    }
    if(!in.eof())throw std::runtime_error("failed reading XY boundary");flush();
    if(loops.empty())throw std::runtime_error("empty XY boundary");return loops;
}
std::ofstream output(const std::filesystem::path& p) {
    std::ofstream out(p);out.exceptions(std::ios::badbit|std::ios::failbit);out<<std::setprecision(17);return out;
}
void exportResult(const std::filesystem::path& dir,const im::Grid& g,const im::Result& r,const std::string& caseName) {
    const auto start=std::chrono::steady_clock::now();const auto& c=g.controls;const auto& s=r.state;const auto& m=r.metrics;
    std::array<std::size_t,3> counts{};for(auto v:g.classification)++counts[v];
    auto cells=output(dir/"cells.csv"),us=output(dir/"u.csv"),vs=output(dir/"v.csv"),walls=output(dir/"walls.csv");
    cells<<"i,j,x,y,u,v,pressure_fluctuation,divergence,classification,solid_mask\n";
    us<<"i,j,x,y,u,solid_mask,wall_force\n";vs<<"i,j,x,y,v,solid_mask,wall_force\n";
    walls<<"x,y,nx,ny,u,v,segment_weight\n";
    std::vector<double> centreU(s.u.size()),centreV(s.u.size()),div(s.u.size());
    for(std::size_t j=0;j<c.ny;++j)for(std::size_t i=0;i<c.nx;++i) {
        const auto k=g.index(i,j),ip=(i+1)%c.nx;
        const double x=(static_cast<double>(i)+.5)*g.dx(),y=(static_cast<double>(j)+.5)*g.dy();
        centreU[k]=.5*(s.u[k]+s.u[g.index(ip,j)]);centreV[k]=.5*(s.v[k]+s.v[g.index(i,j+1)]);
        div[k]=(s.u[g.index(ip,j)]-s.u[k])/g.dx()+(s.v[g.index(i,j+1)]-s.v[k])/g.dy();
        cells<<i<<','<<j<<','<<x<<','<<y<<','<<centreU[k]<<','<<centreV[k]<<','<<s.p[k]<<','<<div[k]<<','<<g.classification[k]<<','<<g.maskCell[k]<<'\n';
        us<<i<<','<<j<<','<<static_cast<double>(i)*g.dx()<<','<<y<<','<<s.u[k]<<','<<g.maskU[k]<<','<<s.wallForceU[k]<<'\n';
    }
    for(std::size_t j=0;j<=c.ny;++j)for(std::size_t i=0;i<c.nx;++i) {
        const auto k=g.index(i,j);
        vs<<i<<','<<j<<','<<(static_cast<double>(i)+.5)*g.dx()<<','<<static_cast<double>(j)*g.dy()<<','<<s.v[k]<<','<<g.maskV[k]<<','<<s.wallForceV[k]<<'\n';
    }
    for(const auto& w:im::sampleWalls(g,s))walls<<w.point.x<<','<<w.point.y<<','<<w.normal.x<<','<<w.normal.y<<','<<w.velocity.x<<','<<w.velocity.y<<','<<w.length<<'\n';
    cells.close();us.close();vs.close();walls.close();
    auto markers=output(dir/"wall-markers.csv");markers<<"x,y,segment_weight,weight,multiplier_u,multiplier_v\n";
    for(std::size_t k=0;k<g.wallStencils.size();++k) {
        const auto& w=g.wallStencils[k];markers<<w.point.x<<','<<w.point.y<<','<<w.length<<','<<w.weight<<','<<s.wallMultipliers[2*k]<<','<<s.wallMultipliers[2*k+1]<<'\n';
    }
    markers.close();
    auto vtk=output(dir/"field.vtk");
    vtk<<"# vtk DataFile Version 3.0\nCartesian immersed analysis lattice including auxiliary solid cells\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS "<<c.nx+1<<' '<<c.ny+1<<" 1\nX_COORDINATES "<<c.nx+1<<" double\n";
    for(std::size_t i=0;i<=c.nx;++i)vtk<<static_cast<double>(i)*g.dx()<<' ';
    vtk<<"\nY_COORDINATES "<<c.ny+1<<" double\n";
    for(std::size_t j=0;j<=c.ny;++j)vtk<<static_cast<double>(j)*g.dy()<<' ';
    vtk<<"\nZ_COORDINATES 1 double\n0\nCELL_DATA "<<c.nx*c.ny<<"\nVECTORS velocity double\n";
    for(std::size_t k=0;k<centreU.size();++k)vtk<<centreU[k]<<' '<<centreV[k]<<" 0\n";
    const auto scalar=[&](const char* name,const auto& values) { vtk<<"SCALARS "<<name<<" double 1\nLOOKUP_TABLE default\n";for(auto a:values)vtk<<a<<'\n'; };
    scalar("pressure_fluctuation",s.p);scalar("divergence",div);scalar("classification",g.classification);scalar("solid_mask",g.maskCell);vtk.close();
    auto geometry=output(dir/"boundary.xy");for(const auto& loop:g.solids) {for(auto p:loop.vertices())geometry<<p.x<<' '<<p.y<<'\n';geometry<<'\n';}geometry.close();
    const double exportSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    auto report=output(dir/"summary.json");
    report<<"{\n\"format\":\"cartmesh2d-immersed-flow-v1\",\n\"solver\":\"periodic-mac-immersed\",\n\"case\":"<<jsonString(caseName)<<",\n\"experimental\":true,\n\"product_solver_ready\":false,\n"
          <<"\"boundary_conditions\":{\"x\":\"periodic velocity and pressure fluctuation\",\"y\":\"stationary no-slip\",\"solid\":"<<jsonString(c.surfacePenalty?"coupled surface penalty plus Brinkman":"finite Brinkman penalty")<<"},\n"
          <<"\"converged\":"<<(r.converged?"true":"false")<<",\n\"stop_reason\":"<<jsonString(r.stopReason)<<",\n\"error\":"<<jsonString(r.error)<<",\n\"state_status\":"<<jsonString(s.steps?"last-accepted-iterate":"initial-only")<<",\n"
          <<"\"steps\":"<<s.steps<<",\n\"pseudo_time\":"<<s.pseudoTime<<",\n\"grid\":{\"nx\":"<<c.nx<<",\"ny\":"<<c.ny<<",\"length\":"<<c.length<<",\"height\":"<<c.height<<",\"dx\":"<<g.dx()<<",\"dy\":"<<g.dy()<<",\"cells\":"<<c.nx*c.ny<<",\"classification_counts\":["<<counts[0]<<','<<counts[1]<<','<<counts[2]<<"],\"solid_area_geometry\":"<<g.solidArea<<",\"fluid_area_geometry\":"<<c.length*c.height-g.solidArea<<"},\n"
          <<"\"controls\":{\"viscosity\":"<<c.viscosity<<",\"drive\":"<<c.drive<<",\"penalty_time\":"<<c.penaltyTime<<",\"mask_half_width_cells\":"<<c.maskHalfWidthCells<<",\"max_steps\":"<<c.maxSteps<<",\"max_time_step\":"<<c.maxTimeStep<<",\"steady_tolerance\":"<<c.steadyTolerance<<",\"continuity_tolerance\":"<<c.continuityTolerance<<",\"linear_tolerance\":"<<c.linearTolerance<<",\"linear_solver\":"<<jsonString(c.pressureSolver==im::PressureSolver::SystemCholesky?"cholesky":(c.pressureSolver==im::PressureSolver::Jacobi?"jacobi":"ic0"))<<",\"wall_method\":"<<jsonString(c.surfacePenalty?"surface-penalty":"brinkman")<<",\"wall_quadrature\":\"piecewise-gauss3\",\"wall_penalty_time\":"<<c.wallPenaltyTime<<"},\n"
          <<"\"normalization\":{\"reference_velocity\":"<<c.drive*c.height*c.height/(12*c.viscosity)<<",\"momentum\":\"max abs(nu Laplacian(u)-div(uu)-grad(p)-chi*u/eta+wall_force+drive) / drive\",\"continuity\":\"max abs(div(u))*H/Uref, including pressure reference cell\",\"field_change\":\"max abs(delta velocity)/Uref\",\"wall_constraint\":\"max abs(Ju - wall_penalty_time*multiplier/sqrt(weight))/Uref; algebraic wall law, not physical no-slip error\",\"pressure_linear\":\"true L2 residual of pressure or coupled pressure/wall system <= 1e-13 + linear_tolerance*L2(rhs)\",\"pressure\":\"kinematic periodic fluctuation; total p = fluctuation - drive*x\"},\n"
          <<"\"metrics\":{\"momentum\":"<<m.momentum<<",\"continuity\":"<<m.continuity<<",\"field_change\":"<<m.fieldChange<<",\"pressure_linear_residual\":"<<m.pressureLinearResidual<<",\"pressure_linear_rhs_norm\":"<<m.pressureLinearRhsNorm<<",\"last_dt\":"<<m.dt<<",\"mean_velocity\":"<<m.meanVelocity<<",\"max_speed\":"<<m.maxSpeed<<",\"flux_spread\":"<<m.fluxSpread<<",\"wall_speed_max\":"<<m.wallSpeed<<",\"wall_normal_speed_max\":"<<m.wallNormalSpeed<<",\"wall_tangential_speed_max\":"<<m.wallTangentialSpeed<<",\"deep_solid_speed_max\":"<<m.deepSolidSpeed<<",\"penalty_drag_per_density\":"<<m.penaltyDrag<<",\"penalty_lift_per_density\":"<<m.penaltyLift<<",\"surface_drag_per_density\":"<<m.surfaceDrag<<",\"surface_lift_per_density\":"<<m.surfaceLift<<",\"surface_power_per_density\":"<<m.surfacePower<<",\"wall_constraint_residual\":"<<m.wallConstraintResidual<<",\"marker_speed_max\":"<<m.markerSpeed<<",\"wall_normal_flux_net\":"<<m.wallNormalFluxNet<<",\"wall_normal_flux_abs\":"<<m.wallNormalFluxAbs<<",\"channel_wall_drag_per_density\":"<<m.channelWallDrag<<",\"force_balance\":"<<m.forceBalance<<",\"channel_relative_l2\":";
    if(g.solids.empty())report<<m.channelRelativeL2;else report<<"null";
    report<<"},\n\"timing\":{\"grid_seconds\":"<<g.gridSeconds<<",\"boundary_seconds\":"<<g.boundarySeconds<<",\"solve_seconds\":"<<r.solveSeconds<<",\"pressure_seconds\":"<<r.pressureSeconds<<",\"export_seconds\":"<<exportSeconds<<",\"pressure_iterations\":"<<r.pressureIterations<<"},\n"
          <<"\"qualification\":{\"physical_accuracy\":\"not-qualified\",\"mesh_independence\":\"not-run\",\"external_checkMesh\":\"not-applicable: full analysis lattice, no conformal fluid polyMesh\",\"transient_accuracy\":\"not-qualified; pseudo-time used to reach steady state\"}\n}\n";
    report.close();
}
}
int main(int argc,char** argv) {
    try {
        im::Controls c;std::string caseName="channel",boundary,outputDir;
        for(int i=1;i<argc;++i) {
            const std::string key=argv[i];
            if(key=="--help") { std::cout<<"Experimental native 2D uniform Cartesian MAC immersed steady flow\n"
                <<"--output NEW_DIRECTORY --case channel|cylinder|custom [--boundary SOLID.xy]\n"
                <<"--nx N --ny N --length L --height H --nu KINEMATIC_VISCOSITY --drive ACCELERATION\n"
                <<"--penalty-time SECONDS --mask-half-width CELLS --max-steps N --max-time-step SECONDS\n"
                <<"--steady-tolerance R --continuity-tolerance R --linear-tolerance R --linear-solver auto|ic0|jacobi|cholesky\n"
                <<"--wall-method brinkman|surface-penalty --wall-penalty-time SECONDS (opt-in, finite slip)\n"
                <<"Auto uses IC0 for Brinkman and Jacobi for surface penalty; coupled IC0 is rejected.\n"
                <<"x is periodic, y walls are no-slip. Pressure output is the periodic kinematic fluctuation.\n"
                <<"Solid auxiliary fields are retained and labelled; this is not a CM2D fluid mesh.\n"
                <<"Exit 0: steady-converged, 2: budget exhausted/candidate failed, 130: cancelled, 1: input/output error.\n";return 0; }
            if(i+1>=argc)throw std::runtime_error("missing value for "+key);
            const std::string value=argv[++i];
            if(key=="--output")outputDir=value;else if(key=="--case")caseName=value;else if(key=="--boundary")boundary=value;
            else if(key=="--nx")c.nx=integer(value);else if(key=="--ny")c.ny=integer(value);
            else if(key=="--max-steps")c.maxSteps=integer(value);else if(key=="--length")c.length=number(value);
            else if(key=="--height")c.height=number(value);else if(key=="--nu")c.viscosity=number(value);
            else if(key=="--drive")c.drive=number(value);else if(key=="--penalty-time")c.penaltyTime=number(value);
            else if(key=="--mask-half-width")c.maskHalfWidthCells=number(value);else if(key=="--max-time-step")c.maxTimeStep=number(value);
            else if(key=="--steady-tolerance")c.steadyTolerance=number(value);else if(key=="--continuity-tolerance")c.continuityTolerance=number(value);
            else if(key=="--linear-tolerance")c.linearTolerance=number(value);else if(key=="--linear-solver") {
                if(value=="auto")c.pressureSolver=im::PressureSolver::Automatic;
                else if(value=="ic0")c.pressureSolver=im::PressureSolver::IC0;
                else if(value=="jacobi")c.pressureSolver=im::PressureSolver::Jacobi;
                else if(value=="cholesky")c.pressureSolver=im::PressureSolver::SystemCholesky;
                else throw std::runtime_error("unknown linear solver");
            } else if(key=="--wall-method") {
                if(value!="brinkman"&&value!="surface-penalty")throw std::runtime_error("unknown wall method");c.surfacePenalty=value=="surface-penalty";
            } else if(key=="--wall-penalty-time")c.wallPenaltyTime=number(value);
            else throw std::runtime_error("unknown argument: "+key);
        }
        if(outputDir.empty())throw std::runtime_error("--output NEW_DIRECTORY is required");
        if(std::filesystem::exists(outputDir))throw std::runtime_error("output directory already exists; choose a new directory");
        std::vector<BoundaryLoop> solids;
        if(caseName=="custom") { if(boundary.empty())throw std::runtime_error("custom case requires --boundary");solids=readSolids(boundary); }
        else {
            if(!boundary.empty())throw std::runtime_error("--boundary requires --case custom");
            if(caseName=="cylinder") {
                std::vector<Point2D> vertices;
                for(int k=0;k<128;++k) {const double a=2*std::numbers::pi*k/128;vertices.push_back({c.length/3+.15*c.height*std::cos(a),.5*c.height+.15*c.height*std::sin(a)});}
                solids.emplace_back(vertices);
            } else if(caseName!="channel")throw std::runtime_error("unknown case");
        }
#ifdef __APPLE__
        // Match the shared factorization contract before its first Accelerate
        // call. This is local to this CLI process; product defaults are untouched.
        if(c.pressureSolver==im::PressureSolver::SystemCholesky && setenv("VECLIB_MAXIMUM_THREADS","1",1)!=0)
            throw std::runtime_error("cannot initialize deterministic system Cholesky threads");
#endif
        const auto grid=im::makeGrid(c,solids);std::filesystem::create_directories(outputDir);
        auto history=output(std::filesystem::path(outputDir)/"history.csv");
        history<<"step,pseudo_time,dt,momentum,continuity,field_change,pressure_linear_residual,pressure_linear_rhs_norm,pressure_iterations,wall_constraint_residual,marker_speed\n";
        std::signal(SIGINT,cancel);std::signal(SIGTERM,cancel);
        const auto result=im::solve(grid,[&](const im::State& s,const im::Metrics& m) {
            history<<s.steps<<','<<s.pseudoTime<<','<<m.dt<<','<<m.momentum<<','<<m.continuity<<','<<m.fieldChange<<','<<m.pressureLinearResidual<<','<<m.pressureLinearRhsNorm<<','<<m.pressureIterations<<','<<m.wallConstraintResidual<<','<<m.markerSpeed<<'\n';
            if(s.steps==1||s.steps%500==0)std::cout<<"step="<<s.steps<<" momentum="<<m.momentum<<" continuity="<<m.continuity<<std::endl;
            return !interrupted;
        });
        history.close();exportResult(outputDir,grid,result,caseName);
        std::cout<<"stop_reason="<<result.stopReason<<" steps="<<result.state.steps<<" momentum="<<result.metrics.momentum<<" continuity="<<result.metrics.continuity<<" wall_speed="<<result.metrics.wallSpeed<<" solve_seconds="<<result.solveSeconds<<'\n';
        if(!result.error.empty())std::cerr<<result.error<<'\n';
        return result.converged?0:(result.stopReason=="cancelled"?130:2);
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
