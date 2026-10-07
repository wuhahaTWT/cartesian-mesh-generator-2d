#include "cartmesh2d/fv/ScalarTransport2D.hpp"
#include "cartmesh2d/fv/WallGradient2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {
constexpr double diffusivity=.05, background=300., pulseAge=.1, endTime=.2;
struct Reference {
    bool pulse=false;
    double time=0;
    double value(Point2D p) const {
        const double radius2=p.x*p.x+p.y*p.y;
        return background+(pulse ? pulseAge/(pulseAge+time)*std::exp(-radius2/(4*diffusivity*(pulseAge+time)))
                                 : .5*std::log(radius2));
    }
    // Exact integral over the actual straight polygon face, not over a circle
    // with a guessed radius and not a midpoint approximation to the heat flux.
    double flux(const Face& f) const {
        if(!pulse) {
            const Point2D a{f.centre.x+f.areaVector.y/2,f.centre.y-f.areaVector.x/2};
            const Point2D b{f.centre.x-f.areaVector.y/2,f.centre.y+f.areaVector.x/2};
            return -diffusivity*std::atan2(a.x*b.y-a.y*b.x,a.x*b.x+a.y*b.y);
        }
        const long double length=std::hypot(f.areaVector.x,f.areaVector.y);
        const long double normal=(f.centre.x*f.areaVector.x+f.centre.y*f.areaVector.y)/length;
        const long double tangent=(-f.centre.x*f.areaVector.y+f.centre.y*f.areaVector.x)/length;
        const long double age=pulseAge+time, scale=2*std::sqrt(diffusivity*age);
        const long double integral=std::sqrt(std::numbers::pi_v<long double>*diffusivity*age)
            *(std::erf((tangent+length/2)/scale)-std::erf((tangent-length/2)/scale));
        return static_cast<double>(pulseAge/(2*age*age)*normal*std::exp(-normal*normal/(4*diffusivity*age))*integral);
    }
};
struct FluxError {
    double integral=0, absoluteError=0, maxDensityError=0;
    void add(double actual,double exact,double length) {
        integral+=actual;absoluteError+=std::abs(actual-exact);
        maxDensityError=std::max(maxDensityError,std::abs(actual-exact)/length);
    }
    void print(double reference,double absoluteReference) const {
        std::cout<<"{\"integral\":"<<integral<<",\"relativeIntegralError\":"<<(integral-reference)/std::abs(reference)
                 <<",\"relativeL1FluxError\":"<<absoluteError/absoluteReference
                 <<",\"maxFluxDensityError\":"<<maxDensityError<<'}';
    }
};
std::ofstream output(const std::string& prefix,const char* suffix) {
    std::ofstream stream(prefix+suffix);stream.exceptions(std::ios::failbit|std::ios::badbit);
    stream<<std::setprecision(17);return stream;
}
}

int main(int argc,char** argv) {
    try {
        if(argc<4 || argc>6)throw std::runtime_error("usage: thermal_spatial_probe FINAL.solver.cm2d PREFIX harmonic|pulse [steps] [unrestricted|bounded|bounded-spatial]");
        const std::string kind=argv[3],prefix=argv[2],mode=argc>5?argv[5]:"unrestricted";
        if(kind!="harmonic" && kind!="pulse")throw std::runtime_error("unknown reference case");
        if(mode!="unrestricted" && mode!="bounded" && mode!="bounded-spatial")throw std::runtime_error("unknown flux correction");
        const bool pulse=kind=="pulse";
        if(!pulse && (argc>4 || mode!="unrestricted"))throw std::runtime_error("harmonic probe is steady and unrestricted");
        std::size_t steps=100;
        if(argc>4) {std::size_t used=0;const std::string argument=argv[4];steps=std::stoull(argument,&used);
            if(used!=argument.size() || !steps || steps>100000)throw std::runtime_error("invalid step count");}
        const double dt=pulse?endTime/static_cast<double>(steps):0;
        auto read=readCm2dTopology(argv[1]);if(!read.valid())throw std::runtime_error(read.error);
        const auto mesh=makeFvMesh2D(read.topology);
        Reference exact{pulse,0};
        ScalarTransportProblem2D problem;problem.diffusivity=diffusivity;
        problem.volumeFlux.assign(mesh.faces.size(),0.);problem.source=[](Point2D){return 0.;};
        problem.boundary=[&](std::size_t,const Face& f){return ScalarBoundary2D{ScalarBoundaryKind2D::Value,exact.value(f.centre),{}};};
        ScalarTransportControls2D controls;controls.preconditioner=ScalarPreconditioner2D::ILU0;
        if(mode=="bounded")controls.fluxCorrection=ScalarFluxCorrection2D::Bounded;
        else if(mode=="bounded-spatial")controls.fluxCorrection=ScalarFluxCorrection2D::BoundedSpatial;
        ScalarTransportWorkspace2D workspace;
        auto history=output(prefix,".history.csv");
        history<<"step,time,dt,minValue,maxValue,globalBalance,maxDiagonalScaledImbalance,limitedFaces,minimumFluxCorrection\n";
        ScalarTransportResult2D result;
        std::vector<double> previous;for(const auto& cell:mesh.cells)previous.push_back(exact.value(cell.centre));
        double minimum=*std::min_element(previous.begin(),previous.end()),maximum=*std::max_element(previous.begin(),previous.end());
        for(std::size_t step=0;step<(pulse?steps:1);++step) {
            exact.time=pulse?static_cast<double>(step+1)*dt:0;
            result=solveScalarTransport2D(mesh,problem,controls,pulse?previous:std::vector<double>{},dt,&workspace);
            if(!result.converged)throw std::runtime_error("native reference solve did not converge at step "+std::to_string(step+1));
            previous=result.values;minimum=std::min(minimum,result.minValue);maximum=std::max(maximum,result.maxValue);
            history<<step+1<<','<<exact.time<<','<<dt<<','<<result.minValue<<','<<result.maxValue<<','<<result.globalBalance<<','
                   <<result.history.back().maxDiagonalScaledImbalance<<','<<result.limitedFaces<<','<<result.minimumFluxCorrection<<'\n';
        }
        std::vector<double> exactValues,exactPrevious;
        for(const auto& cell:mesh.cells) {exactValues.push_back(exact.value(cell.centre));Reference old{pulse,exact.time-dt};exactPrevious.push_back(old.value(cell.centre));}
        const auto onExact=evaluateScalarTransport2D(mesh,problem,exactValues,controls,pulse?exactPrevious:std::vector<double>{},dt);
        std::vector<bool> prescribed(mesh.faces.size());std::vector<std::optional<std::size_t>> partners(mesh.faces.size());
        for(std::size_t id=0;id<mesh.faces.size();++id)prescribed[id]=!mesh.faces[id].neighbour;
        double area=0,squaredError=0,maxError=0,heatError=0,wallReference=0,absoluteReference=0;
        std::size_t wallFaces=0;
        auto cells=output(prefix,".cells.csv"),faces=output(prefix,".faces.csv");
        cells<<"cell,x,y,area,exact,value,error\n";
        for(std::size_t id=0;id<mesh.cells.size();++id) {
            const auto& cell=mesh.cells[id];const double delta=result.values[id]-exactValues[id];
            area+=cell.area;squaredError+=cell.area*delta*delta;maxError=std::max(maxError,std::abs(delta));heatError+=cell.area*delta;
            cells<<id<<','<<cell.centre.x<<','<<cell.centre.y<<','<<cell.area<<','<<exactValues[id]<<','<<result.values[id]<<','<<delta<<'\n';
        }
        faces<<"face,x,y,length,exact,operatorOnExact,quadraticFitOnExact,solved\n";
        FluxError operatorError,fitError,solvedError;
        for(std::size_t id=0;id<mesh.faces.size();++id) {
            const auto& f=mesh.faces[id];if(f.neighbour || f.patch!=BoundaryPatch2D::EmbeddedBoundary)continue;
            const double length=std::hypot(f.areaVector.x,f.areaVector.y),target=exact.flux(f),wallValue=exact.value(f.centre);
            wallReference+=target;absoluteReference+=std::abs(target);++wallFaces;
            const auto stencil=quadraticWallGradient2D(mesh,id,prescribed,partners);Vector2D gradient{};
            for(const auto& sample:stencil.samples) {
                const double delta=(sample.boundary?exact.value(mesh.faces[sample.index].centre):exactValues[sample.index])-wallValue;
                gradient.x+=sample.weight.x*delta;gradient.y+=sample.weight.y*delta;
            }
            const double fitted=-diffusivity*dot(gradient,f.areaVector);
            operatorError.add(onExact.diffusiveFlux[id],target,length);fitError.add(fitted,target,length);
            solvedError.add(result.diffusiveFlux[id],target,length);
            faces<<id<<','<<f.centre.x<<','<<f.centre.y<<','<<length<<','<<target<<','<<onExact.diffusiveFlux[id]<<','<<fitted<<','<<result.diffusiveFlux[id]<<'\n';
        }
        if(!wallFaces || !std::isfinite(wallReference) || !absoluteReference)throw std::runtime_error("requires resolved embedded wall and nonzero reference heat flux");
        if(!pulse) {
            // This contour identity verifies an origin-containing hole. The
            // budget covers a floating-point sum, not a PDE accuracy threshold.
            const double expected=2*std::numbers::pi*diffusivity;
            const double roundoff=128*std::numeric_limits<double>::epsilon()*static_cast<double>(wallFaces)*expected;
            if(std::abs(wallReference-expected)>roundoff)throw std::runtime_error("harmonic probe requires a hole containing the origin");
        }
        std::cout<<std::setprecision(17)<<"{\"format\":\"cartmesh2d-thermal-spatial-probe-v1\",\"case\":\""<<kind<<"\",\"fluxCorrection\":\""<<mode
                 <<"\",\"cells\":"<<mesh.cells.size()<<",\"wallFaces\":"<<wallFaces<<",\"converged\":true,\"time\":"<<exact.time
                 <<",\"dt\":"<<dt<<",\"steps\":"<<(pulse?steps:0)<<",\"l2TemperatureError\":"<<std::sqrt(squaredError/area)
                 <<",\"maxTemperatureError\":"<<maxError<<",\"centroidHeatDifference\":"<<heatError
                 <<",\"minimumAcceptedTemperature\":"<<minimum<<",\"maximumAcceptedTemperature\":"<<maximum
                 <<",\"exactWallFlux\":"<<wallReference<<",\"globalBalance\":"<<result.globalBalance<<",\"operatorOnExact\":";
        operatorError.print(wallReference,absoluteReference);std::cout<<",\"quadraticFitOnExact\":";fitError.print(wallReference,absoluteReference);
        std::cout<<",\"solvedFlux\":";solvedError.print(wallReference,absoluteReference);std::cout<<"}\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 2;}
}
