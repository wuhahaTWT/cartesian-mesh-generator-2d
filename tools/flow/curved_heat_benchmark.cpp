#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/fv/HeatConduction2D.hpp"
#include "cartmesh2d/fv/detail/EulerNewton2D.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numbers>
#include <algorithm>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
struct Integrals {double temperature=0,inverse=0;};
Integrals analyticTriangle(Point2D a,Point2D b,Point2D c) {
    // Tensor Gauss-Legendre 8 on a Duffy triangle map, signed Jacobian.
    // This integrates the analytic reference only, never numerical fluxes.
    constexpr double x[]={-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498,
                          .1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
    constexpr double w[]={.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620,
                          .3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
    const double cross=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);Integrals out;
    for(unsigned i=0;i<8;++i)for(unsigned j=0;j<8;++j) {
        const double u=.5*(x[i]+1),v=.5*(x[j]+1);
        const Point2D p{a.x+u*(b.x-a.x)+u*v*(c.x-b.x),a.y+u*(b.y-a.y)+u*v*(c.y-b.y)};
        const double t=2+.2/std::log(2.)*std::log(std::hypot(p.x,p.y)/.5);
        const double weight=.25*w[i]*w[j]*u*cross;out.temperature+=weight*t;out.inverse+=weight/t;
    }return out;
}
Integrals refinedTriangle(Point2D a,Point2D b,Point2D c) {
    const Point2D ab{.5*(a.x+b.x),.5*(a.y+b.y)},ac{.5*(a.x+c.x),.5*(a.y+c.y)},bc{.5*(b.x+c.x),.5*(b.y+c.y)};
    Integrals result;for(const auto q:{analyticTriangle(a,ab,ac),analyticTriangle(ab,b,bc),analyticTriangle(ac,bc,c),analyticTriangle(ab,bc,ac)}){result.temperature+=q.temperature;result.inverse+=q.inverse;}
    return result;
}
}
// Native spatial experiment, not an independent reconstruction of the PDE.
// Solve the ACTUAL affine Fourier operator, starting at a nonanalytic constant.
// Coupled compressible relaxation is run separately by the existing Euler CLI.
int main(int argc,char** argv) {try {
    if(argc!=5)throw std::runtime_error("mesh prefix linear|quadratic trace|isothermal");
    const auto read=readCm2dTopology(argv[1]);if(!read.valid())throw std::runtime_error(read.error);
    const auto mesh=makeFvMesh2D(read.topology);const std::string prefix=argv[2],mode=argv[4];
    const auto scheme=std::string(argv[3])=="quadratic"?WallGradient2D::Quadratic:WallGradient2D::Linear;
    if(std::string(argv[3])!="linear"&&std::string(argv[3])!="quadratic")throw std::runtime_error("invalid wall scheme");
    if(mode!="trace"&&mode!="isothermal")throw std::runtime_error("invalid boundary mode");
    const double k=.37,beta=.2/std::log(2.);
    const auto exact=[&](Point2D p){return 2+beta*std::log(std::hypot(p.x,p.y)/.5);};
    std::vector<HeatBoundary2D> bc;double geometryTrace=0,perimeter=0,minimumArea=1e99;
    std::ofstream boundaries(prefix+".boundaries");boundaries<<std::setprecision(17)<<"CM2D_EULER_BOUNDARY 3\nMESH "<<mesh.cells.size()<<' '<<mesh.faces.size()<<'\n';
    for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour) {
        const auto& f=mesh.faces[id];const bool inner=dot(Vector2D{f.centre.x,f.centre.y},f.areaVector)<0;
        const double value=mode=="trace"?exact(f.centre):inner?2:2.2;
        bc.push_back({id,HeatBoundaryKind2D::Temperature,value,{}});
        const double length=std::hypot(f.areaVector.x,f.areaVector.y);
        geometryTrace+=length*std::abs(exact(f.centre)-(inner?2:2.2));perimeter+=length;
        boundaries<<id<<' '<<f.owner<<' '<<f.centre.x<<' '<<f.centre.y<<' '<<f.areaVector.x<<' '<<f.areaVector.y<<" no-slip-wall 1 0 0 2.1 - \""<<(inner?"inner":"outer")<<"\" temperature "<<value<<" 0 0\n";
    }boundaries<<"END\n";boundaries.close();
    const HeatConductionOperator2D op(mesh,bc,k,scheme);const std::size_t count=mesh.cells.size();
    std::vector<double> t(count,2.1),capacity(count,1),correction(count);
    const auto base=op.evaluate(t,capacity).cellResidual;
    std::vector<std::vector<std::pair<std::size_t,double>>> coefficients(count);std::vector<std::pair<std::size_t,std::size_t>> edges;
    for(std::size_t j=0;j<count;++j){t[j]+=1;const auto next=op.evaluate(t,capacity).cellResidual;t[j]-=1;for(std::size_t i=0;i<count;++i){const double a=next[i]-base[i];if(a!=0){coefficients[i].push_back({j,a});if(i!=j)edges.push_back({i,j});}}}
    std::ofstream matrix(prefix+".matrix.csv");matrix<<std::setprecision(17)<<"row,column,value,area\n";for(std::size_t i=0;i<count;++i)for(const auto& [j,a]:coefficients[i])matrix<<i<<','<<j<<','<<a<<','<<mesh.cells[i].area<<'\n';matrix.close();
    // A signed nonsymmetric operator must be measurable even when its
    // diagonal/ILU pivots are negative. GMRES diagnoses the actual operator;
    // solving such a steady system does NOT qualify it as stable diffusion.
    std::vector<double> diagonal(count);for(std::size_t i=0;i<count;++i)for(const auto& [j,a]:coefficients[i])diagonal[i]+=std::abs(a);
    std::size_t nonpositiveDiagonalRows=0;for(std::size_t i=0;i<count;++i){double d=0;for(const auto& [j,a]:coefficients[i])if(j==i)d=a;if(d<=0)++nonpositiveDiagonalRows;}
    const auto multiply=[&](const std::vector<double>& x){std::vector<double> y(count);for(std::size_t i=0;i<count;++i)for(const auto& [j,a]:coefficients[i])y[i]+=a*x[j]/diagonal[j];return y;};
    const auto solve=[&](const std::vector<double>& rhs,double tolerance){
        std::vector<double> x(count),r=rhs;const double target=tolerance*detail::newtonNorm(rhs);std::size_t work=0;
        for(unsigned iteration=0;iteration<14&&detail::newtonNorm(r)>target;++iteration){const auto delta=detail::eulerGmres(multiply,r,5000,work);for(std::size_t i=0;i<count;++i)x[i]+=delta[i];const auto ax=multiply(x);for(std::size_t i=0;i<count;++i)r[i]=rhs[i]-ax[i];}
        if(detail::newtonNorm(r)>std::max(target,2e-14*detail::newtonNorm(rhs)))throw std::runtime_error("benchmark GMRES tolerance exhausted");
        for(std::size_t i=0;i<count;++i)x[i]/=diagonal[i];return x;
    };
    std::vector<double> rhs=base;for(auto& x:rhs)x=-x;
    const auto loose=solve(rhs,1e-10);correction=solve(rhs,1e-12);
    double sensitivity=0;for(std::size_t i=0;i<count;++i){sensitivity=std::max(sensitivity,std::abs(loose[i]-correction[i]));t[i]+=correction[i];}
    for(int iteration=0;iteration<3;++iteration){rhs=op.evaluate(t,capacity).cellResidual;for(auto& x:rhs)x=-x;correction=solve(rhs,1e-10);for(std::size_t i=0;i<count;++i)t[i]+=correction[i];}
    const auto result=op.evaluate(t,capacity);
    std::ofstream cells(prefix+".heat.cells.csv");cells<<std::setprecision(17)<<"cell,x,y,area,temperature,analyticCentroidTemperature,analyticMeanTemperature,analyticHarmonicTemperature\n";
    double temperatureL1=0,volume=0,residual=0,nearError=0,nearArea=0;
    double meanError=0,harmonicError=0,representationError=0,quadratureSensitivity=0;
    for(std::size_t i=0;i<count;++i){const auto& c=mesh.cells[i];const double e=std::abs(t[i]-exact(c.centre));temperatureL1+=c.area*e;volume+=c.area;minimumArea=std::min(minimumArea,c.area);residual+=std::abs(result.cellResidual[i]);
        const bool wall=std::any_of(c.faces.begin(),c.faces.end(),[&](auto f){return !mesh.faces[f].neighbour;});if(wall){nearError+=c.area*e;nearArea+=c.area;}
        Integrals integral,refined;const auto& polygon=read.topology.cells[i].vertices;
        for(std::size_t v=0;v<polygon.size();++v) {
            const auto a=read.topology.vertices[polygon[v]].point,b=read.topology.vertices[polygon[(v+1)%polygon.size()]].point;
            const auto q=analyticTriangle(c.centre,a,b),r=refinedTriangle(c.centre,a,b);
            integral.temperature+=q.temperature;integral.inverse+=q.inverse;refined.temperature+=r.temperature;refined.inverse+=r.inverse;
        }
        const double mean=refined.temperature/c.area,harmonic=c.area/refined.inverse;
        meanError+=c.area*std::abs(t[i]-mean);harmonicError+=c.area*std::abs(t[i]-harmonic);
        representationError+=c.area*std::abs(harmonic-exact(c.centre));
        quadratureSensitivity=std::max(quadratureSensitivity,std::abs(harmonic-c.area/integral.inverse));
        cells<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','<<t[i]<<','<<exact(c.centre)<<','<<mean<<','<<harmonic<<'\n';
    }
    std::ofstream faces(prefix+".heat.faces.csv");faces<<std::setprecision(17)<<"face,x,y,sx,sy,temperatureBoundary,numericalHeat,analyticChordHeat,inner\n";
    double error=0,magnitude=0,innerHeat=0,outerHeat=0;
    for(const auto& boundary:bc){const auto id=boundary.face;const auto& f=mesh.faces[id];
        // Exact integral of grad(log r) dot outward S on the actual chord:
        // signed angle between its endpoints, not a face-centre approximation.
        const Point2D a{f.centre.x+.5*f.areaVector.y,f.centre.y-.5*f.areaVector.x},b{f.centre.x-.5*f.areaVector.y,f.centre.y+.5*f.areaVector.x};
        const double q=-k*beta*std::atan2(a.x*b.y-a.y*b.x,a.x*b.x+a.y*b.y);
        error+=std::abs(result.faceHeatFlux[id]-q);magnitude+=std::abs(q);const bool inner=dot(Vector2D{f.centre.x,f.centre.y},f.areaVector)<0;
        (inner?innerHeat:outerHeat)+=result.faceHeatFlux[id];faces<<id<<','<<f.centre.x<<','<<f.centre.y<<','<<f.areaVector.x<<','<<f.areaVector.y<<','<<boundary.value<<','<<result.faceHeatFlux[id]<<','<<q<<','<<inner<<'\n';
    }
    std::ofstream evidence(prefix+".json");evidence<<std::setprecision(17)<<"{\"scope\":\"native steady Fourier spatial isolation; not coupled compressible evolution\",\"cells\":"<<count<<",\"boundaryFaces\":"<<bc.size()<<",\"scheme\":\""<<argv[3]<<"\",\"boundaryMode\":\""<<mode<<"\",\"minimumArea\":"<<minimumArea<<",\"temperatureRelativeL1\":"<<temperatureL1/volume/.2<<",\"temperatureCellMeanRelativeL1\":"<<meanError/volume/.2<<",\"temperatureHarmonicMeanRelativeL1\":"<<harmonicError/volume/.2<<",\"harmonicCentroidRepresentationRelativeL1\":"<<representationError/volume/.2<<",\"analyticQuadratureSensitivityK\":"<<quadratureSensitivity<<",\"nearWallTemperatureRelativeL1\":"<<nearError/nearArea/.2<<",\"wallFluxRelativeL1\":"<<error/magnitude<<",\"innerHeat\":"<<innerHeat<<",\"outerHeat\":"<<outerHeat<<",\"analyticHeatMagnitude\":"<<2*std::numbers::pi*k*beta<<",\"geometryTraceRelativeL1\":"<<geometryTrace/perimeter/.2<<",\"algebraicResidualRelative\":"<<residual/magnitude<<",\"toleranceSensitivityK\":"<<sensitivity<<",\"nonpositiveDiffusionDiagonalRows\":"<<nonpositiveDiagonalRows<<",\"nonMonotoneRows\":"<<op.nonMonotoneRows()<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
