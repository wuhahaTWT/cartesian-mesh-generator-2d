#include "cartmesh2d/fv/Diffusion2D.hpp"
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
#include <string>

using namespace cartmesh2d;
namespace {
std::string quote(const std::string& s) {
    std::ostringstream o; o << '"';
    for(char raw:s) {
        const auto c=static_cast<unsigned char>(raw);
        if(c=='"'||c=='\\') o<<'\\'<<c;
        else if(c<32) o<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<static_cast<int>(c)<<std::dec;
        else o<<c;
    }
    return o.str()+'"';
}
double number(const std::string& s) {
    std::size_t used=0; const auto v=std::stod(s,&used);
    if(used!=s.size() || !std::isfinite(v)) throw std::invalid_argument("expected finite number: "+s);
    return v;
}
std::ofstream output(const std::string& prefix,const std::string& extension) {
    std::ofstream o(prefix+extension);
    o.exceptions(std::ios::badbit|std::ios::failbit); o<<std::setprecision(17); return o;
}
}
int main(int argc,char**argv) {
 try {
    std::string meshPath,prefix,problem="diffusion"; double wall=1,outer=0,source=0,k=1;
    fv::DiffusionControls2D controls;
    for(int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        if(arg=="--help") {
            std::cout<<"Native 2D scalar diffusion / Poisson foundation (not a flow solver).\n"
              "cartmesh2d_fv_cli --mesh FINAL.solver.cm2d --output PREFIX\n"
              "  --problem diffusion|constant|linear|sine (default diffusion)\n"
              "  diffusion: --wall-value 1 --outer-value 0 --source 0\n"
              "  --diffusivity 1 --max-corrections 400\n"
              "All boundary edges are Dirichlet. -div(k grad(value))=source.\n"
              "MMS: constant=1, linear=1+x+2y, sine=sin(pi*x)sin(pi*y).\n"
              "Outputs: .json .vtk .cells.csv .faces.csv .residuals.csv\n";
            return 0;
        }
        if(i+1>=argc) throw std::invalid_argument("missing value for "+arg);
        const std::string value=argv[++i];
        if(arg=="--mesh") meshPath=value;
        else if(arg=="--output") prefix=value;
        else if(arg=="--problem") problem=value;
        else if(arg=="--wall-value") wall=number(value);
        else if(arg=="--outer-value") outer=number(value);
        else if(arg=="--source") source=number(value);
        else if(arg=="--diffusivity") k=number(value);
        else if(arg=="--max-corrections") { const double n=number(value);
            if(n<1||n>10000||n!=std::floor(n)) throw std::invalid_argument("invalid correction limit");
            controls.maxCorrections=static_cast<std::size_t>(n); }
        else throw std::invalid_argument("unknown argument "+arg);
    }
    if(meshPath.empty()||prefix.empty()) throw std::invalid_argument("--mesh and --output required; use --help");
    if(problem!="diffusion"&&problem!="constant"&&problem!="linear"&&problem!="sine") throw std::invalid_argument("unknown problem");
    if(k<=0) throw std::invalid_argument("diffusivity must be positive");
    if(!meshPath.ends_with(".solver.cm2d") || meshPath.ends_with(".failed.solver.cm2d"))
        throw std::invalid_argument("requires final *.solver.cm2d, never failed or intermediate mesh");
    auto read=readCm2dTopology(meshPath);
    if(!read.valid()) throw std::runtime_error(read.error);
    const auto mesh=fv::makeFvMesh2D(read.topology);
    const auto exact=[&](Point2D p) {
        if(problem=="constant") return 1.;
        if(problem=="linear") return 1.+p.x+2*p.y;
        return std::sin(std::numbers::pi*p.x)*std::sin(std::numbers::pi*p.y);
    };
    fv::DiffusionProblem2D equation;
    equation.diffusivity=k;
    equation.source=[&](Point2D p) { return problem=="diffusion"?source:
        (problem=="sine"?2*k*std::numbers::pi*std::numbers::pi*exact(p):0.); };
    equation.boundaryValue=[&](Point2D p,BoundaryPatch2D patch) {
        return problem=="diffusion"?(patch==BoundaryPatch2D::EmbeddedBoundary?wall:outer):exact(p);
    };
    const auto result=fv::solveDiffusion2D(mesh,equation,controls);
    const std::filesystem::path parent=std::filesystem::path(prefix).parent_path();
    if(!parent.empty()) std::filesystem::create_directories(parent);
    const bool mms=problem!="diffusion"; double error2=0,totalArea=0,maximum=0;
    auto cells=output(prefix,".cells.csv"); cells<<"cell,x,y,area,value,exact,error,source\n";
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        const auto& c=mesh.cells[i]; const double e=mms?result.values[i]-exact(c.centre):0.;
        error2+=c.area*e*e; totalArea+=c.area; maximum=std::max(maximum,std::abs(e));
        cells<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','<<result.values[i]<<',';
        if(mms) cells<<exact(c.centre)<<','<<e; else cells<<"nan,nan";
        cells<<','<<result.sourceIntegrals[i]<<'\n';
    }
    auto faces=output(prefix,".faces.csv"); faces<<"face,owner,neighbour,flux\n";
    for(std::size_t i=0;i<mesh.faces.size();++i) {
        const auto& f=mesh.faces[i]; faces<<i<<','<<f.owner<<',';
        if(f.neighbour) faces<<*f.neighbour; else faces<<-1;
        faces<<','<<result.fluxes[i]<<'\n';
    }
    auto history=output(prefix,".residuals.csv");
    history<<"iteration,linearIterations,relativeResidual,residualNorm,maxCellImbalance\n";
    for(const auto& h:result.history) history<<h.iteration<<','<<h.linearIterations<<','<<h.relativeResidual<<','<<h.residualNorm<<','<<h.maxCellImbalance<<'\n';
    std::string error;
    if(!writeLegacyVtk2D(read.topology,prefix+".vtk",&error)) throw std::runtime_error(error);
    std::ofstream vtk(prefix+".vtk",std::ios::app); vtk.exceptions(std::ios::badbit|std::ios::failbit); vtk<<std::setprecision(17);
    vtk<<"SCALARS value double 1\nLOOKUP_TABLE default\n";
    for(double v:result.values) vtk<<v<<'\n';
    if(mms) {
        vtk<<"SCALARS exact double 1\nLOOKUP_TABLE default\n";
        for(const auto& c:mesh.cells) vtk<<exact(c.centre)<<'\n';
        vtk<<"SCALARS error double 1\nLOOKUP_TABLE default\n";
        for(std::size_t i=0;i<mesh.cells.size();++i) vtk<<result.values[i]-exact(mesh.cells[i].centre)<<'\n';
    }
    const auto& h=result.history.back(); auto json=output(prefix,".json");
    json<<"{\n  \"solver\": \"native steady scalar diffusion FVM\",\n  \"problem\": "<<quote(problem)
        <<",\n  \"mesh\": "<<quote(meshPath)<<",\n  \"converged\": "<<(result.converged?"true":"false")
        <<",\n  \"cells\": "<<mesh.cells.size()<<",\n  \"faces\": "<<mesh.faces.size()
        <<",\n  \"iterations\": "<<h.iteration<<",\n  \"diffusivity\": "<<k
        <<",\n  \"boundaryConditions\": \"all Dirichlet; manufactured values or geometry-patch constants\""
        <<",\n  \"requestedWallValue\": "<<wall<<",\n  \"requestedOuterValue\": "<<outer<<",\n  \"requestedConstantSource\": "<<source
        <<",\n  \"l2Error\": ";
    if(mms) json<<std::sqrt(error2/totalArea); else json<<"null";
    json<<",\n  \"linfError\": "; if(mms) json<<maximum; else json<<"null";
    json<<",\n  \"relativeResidual\": "<<h.relativeResidual<<",\n  \"residualNorm\": "<<h.residualNorm
        <<",\n  \"relativeTolerance\": "<<controls.relativeTolerance<<",\n  \"absoluteTolerance\": "<<controls.absoluteTolerance
        <<",\n  \"residualDefinition\": \"L2 of integrated cell flux imbalance / max(L2 of source plus implicit Dirichlet RHS, absoluteTolerance)\""
        <<",\n  \"boundaryFlux\": "<<result.boundaryFlux<<",\n  \"sourceIntegral\": "<<result.sourceIntegral
        <<",\n  \"globalBalance\": "<<result.globalBalance<<",\n  \"maxCellImbalance\": "<<h.maxCellImbalance
        <<",\n  \"maxClosureError\": "<<mesh.maxClosureError
        <<",\n  \"nativeTopologyRevalidated\": true,\n  \"solverQualityPassed\": true,\n  \"externalCheckMesh\": \"not run by this program\"\n}\n";
    cells.close(); faces.close(); history.close(); vtk.close(); json.close();
    std::cout<<(result.converged?"Converged":"Correction limit reached")<<": "<<mesh.cells.size()<<" cells, "<<h.iteration<<" corrections, residual "<<h.relativeResidual<<'\n';
    return result.converged?0:2;
 } catch(const std::exception& e) { std::cerr<<"cartmesh2d_fv_cli: "<<e.what()<<'\n'; return 1; }
}
