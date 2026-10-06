#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Contract the already-exported native face momentum fluxes with the difference
// of two accepted cell velocity fields on the identical final mesh. No second
// flow equation is assembled, no cell is omitted, and no pressure gauge shift
// is needed: a constant face pressure cancels by each polygon's closure.
int main(int argc,char** argv) {
    using namespace cartmesh2d;using namespace cartmesh2d::fv;
    if(argc!=5)throw std::runtime_error("usage: branch-energy mesh prefixA prefixB output.csv");
    const auto read=readCm2dTopology(argv[1]);if(!read.valid())throw std::runtime_error(read.error);
    const auto m=makeFvMesh2D(read.topology);const auto n=m.cells.size(),nf=m.faces.size();
    const auto loadCells=[&](const std::string& prefix) {
        std::ifstream in(prefix+".cells.csv");if(!in)throw std::runtime_error("missing cells");
        std::string line;std::getline(in,line);std::vector<std::array<double,3>> u(n);std::size_t count=0;
        while(std::getline(in,line)) {
            std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);
            std::size_t id;double x,y,a,p,speed;Vector2D v;
            if(!(s>>id>>x>>y>>a>>v.x>>v.y>>p>>speed)||id!=count||id>=n)throw std::runtime_error("invalid cells");
            const auto& c=m.cells[id];if(x!=c.centre.x||y!=c.centre.y||a!=c.area)throw std::runtime_error("different mesh");
            u[id]={v.x,v.y,p};++count;
        }
        if(count!=n)throw std::runtime_error("truncated cells");
        return u;
    };
    const auto loadForces=[&](const std::string& prefix, std::vector<double>& flux) {
        std::ifstream in(prefix+".faces.csv");if(!in)throw std::runtime_error("missing faces");
        std::string line;std::getline(in,line);std::vector<std::array<Vector2D,3>> force(n);std::size_t count=0;
        while(std::getline(in,line)) {
            std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);
            std::size_t id,owner;long long neighbour;double q,p,ax,ay,dx,dy;int wall;
            if(!(s>>id>>owner>>neighbour>>q>>p>>ax>>ay>>dx>>dy>>wall)||id!=count||id>=nf)throw std::runtime_error("invalid faces");
            flux.at(id)=q;
            const auto& f=m.faces[id];if(owner!=f.owner||neighbour!=(f.neighbour?static_cast<long long>(*f.neighbour):-1))throw std::runtime_error("different incidence");
            const std::array<Vector2D,3> value{{{ax,ay},{dx,dy},{p*f.areaVector.x,p*f.areaVector.y}}};
            for(std::size_t k=0;k<3;++k) {
                force[owner][k].x+=value[k].x;force[owner][k].y+=value[k].y;
                if(f.neighbour){force[*f.neighbour][k].x-=value[k].x;force[*f.neighbour][k].y-=value[k].y;}
            }
            ++count;
        }
        if(count!=nf)throw std::runtime_error("truncated faces");
        return force;
    };
    const auto a=loadCells(argv[2]),b=loadCells(argv[3]);
    std::vector<double> qa(nf),qb(nf);const auto fa=loadForces(argv[2],qa),fb=loadForces(argv[3],qb);
    std::vector<double> du(n),dv(n),dp(n),zero(nf),divCell(n),divFlux(n);
    std::vector<bool> fixedU(nf),fixedV(nf);
    double xmin=1e300,xmax=-1e300;
    for(const auto& f:m.faces)if(!f.neighbour){xmin=std::min(xmin,f.centre.x);xmax=std::max(xmax,f.centre.x);}
    for(std::size_t i=0;i<n;++i){du[i]=a[i][0]-b[i][0];dv[i]=a[i][1]-b[i][1];dp[i]=a[i][2]-b[i][2];}
    for(std::size_t id=0;id<nf;++id) {
        const auto& f=m.faces[id];if(f.neighbour)continue;
        const bool wall=f.patch==BoundaryPatch2D::EmbeddedBoundary,inlet=f.centre.x==xmin,outlet=f.centre.x==xmax;
        fixedU[id]=wall||inlet;fixedV[id]=!outlet;
    }
    const auto gu=detail::buildFlowGradientStencil2D(m,fixedU).apply(du,zero);
    const auto gv=detail::buildFlowGradientStencil2D(m,fixedV).apply(dv,zero);
    for(std::size_t id=0;id<nf;++id) {
        const auto& f=m.faces[id];const auto i=f.owner;double q=0;
        if(f.neighbour) {
            const auto j=*f.neighbour;const double w=f.neighbourWeight;
            const Point2D point{m.cells[i].centre.x*(1-w)+m.cells[j].centre.x*w,m.cells[i].centre.y*(1-w)+m.cells[j].centre.y*w};
            const auto skew=f.centre-point;
            const Vector2D ug{gu[i].x*(1-w)+gu[j].x*w,gu[i].y*(1-w)+gu[j].y*w};
            const Vector2D vg{gv[i].x*(1-w)+gv[j].x*w,gv[i].y*(1-w)+gv[j].y*w};
            q=(du[i]*(1-w)+du[j]*w+dot(ug,skew))*f.areaVector.x+(dv[i]*(1-w)+dv[j]*w+dot(vg,skew))*f.areaVector.y;
        } else if(f.centre.x==xmax) q=du[i]*f.areaVector.x+dv[i]*f.areaVector.y;
        divCell[i]+=q;divFlux[i]+=qa[id]-qb[id];
        if(f.neighbour){divCell[*f.neighbour]-=q;divFlux[*f.neighbour]-=qa[id]-qb[id];}
    }
    double interpolatedDivergencePower=0,conservativeDivergencePower=0,maxInterpolatedDiv=0,maxConservativeDiv=0;
    for(std::size_t i=0;i<n;++i){interpolatedDivergencePower+=dp[i]*divCell[i];conservativeDivergencePower+=dp[i]*divFlux[i];
        maxInterpolatedDiv=std::max(maxInterpolatedDiv,std::abs(divCell[i])/m.cells[i].area);maxConservativeDiv=std::max(maxConservativeDiv,std::abs(divFlux[i])/m.cells[i].area);}

    std::vector<bool> fixedP(nf);
    for(std::size_t id=0;id<nf;++id)fixedP[id]=!m.faces[id].neighbour && m.faces[id].centre.x==xmax;
    const auto gp=detail::buildFlowGradientStencil2D(m,fixedP,true).apply(dp,zero);
    std::array<double,4> pressurePieces{};
    double centralDivergencePower=0;
    for(std::size_t id=0;id<nf;++id) {
        const auto& f=m.faces[id];const auto i=f.owner;std::array<double,4> pf{};
        double contraction=-du[i]*f.areaVector.x-dv[i]*f.areaVector.y;
        if(f.neighbour) {
            const auto j=*f.neighbour;const double w=f.neighbourWeight;
            const Point2D point{m.cells[i].centre.x*(1-w)+m.cells[j].centre.x*w,m.cells[i].centre.y*(1-w)+m.cells[j].centre.y*w};
            const auto skew=f.centre-point;
            const Vector2D g{gp[i].x*(1-w)+gp[j].x*w,gp[i].y*(1-w)+gp[j].y*w};
            pf[0]=.5*(dp[i]+dp[j]);pf[1]=(w-.5)*(dp[j]-dp[i]);pf[2]=dot(g,skew);
            contraction+=du[j]*f.areaVector.x+dv[j]*f.areaVector.y;
            centralDivergencePower+=(dp[i]-dp[j])*(.5*(du[i]+du[j])*f.areaVector.x+.5*(dv[i]+dv[j])*f.areaVector.y);
        } else if(!fixedP[id]) {pf[0]=dp[i];pf[3]=dot(gp[i],f.centre-m.cells[i].centre);}
        else centralDivergencePower+=dp[i]*(du[i]*f.areaVector.x+dv[i]*f.areaVector.y);
        for(std::size_t k=0;k<4;++k)pressurePieces[k]+=pf[k]*contraction;
    }
    std::ofstream out(argv[4]);out<<std::setprecision(17)<<"cell,x,y,area,deltaKineticEnergy,advectionPower,viscousPower,pressurePower,powerClosure\n";
    std::array<double,3> total{};double energy=0,absoluteClosure=0,maxClosure=0;
    for(std::size_t i=0;i<n;++i) {
        const Vector2D delta{du[i],dv[i]};const auto& c=m.cells[i];
        const double e=.5*c.area*dot(delta,delta);energy+=e;std::array<double,3>w{};
        for(std::size_t k=0;k<3;++k){w[k]=-delta.x*(fa[i][k].x-fb[i][k].x)-delta.y*(fa[i][k].y-fb[i][k].y);total[k]+=w[k];}
        const double closure=w[0]+w[1]+w[2];absoluteClosure+=std::abs(closure);maxClosure=std::max(maxClosure,std::abs(closure));
        out<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','<<e<<','<<w[0]<<','<<w[1]<<','<<w[2]<<','<<closure<<'\n';
    }
    out.close();if(!out)throw std::runtime_error("output failed");
    std::cout<<std::setprecision(17)<<"{\"cells\":"<<n<<",\"differenceKineticEnergy_m4_s2\":"<<energy
      <<",\"advectionPower_m4_s3\":"<<total[0]<<",\"viscousPower_m4_s3\":"<<total[1]<<",\"pressurePower_m4_s3\":"<<total[2]
      <<",\"powerClosure_m4_s3\":"<<total[0]+total[1]+total[2]<<",\"sumAbsoluteCellPowerClosure_m4_s3\":"<<absoluteClosure
      <<",\"maximumCellPowerClosure_m4_s3\":"<<maxClosure
      <<",\"interpolatedDivergencePressurePower_m4_s3\":"<<interpolatedDivergencePower
      <<",\"conservativeDivergencePressurePower_m4_s3\":"<<conservativeDivergencePower
      <<",\"fluxDefectPressurePower_m4_s3\":"<<interpolatedDivergencePower-conservativeDivergencePower
      <<",\"pressureVelocityNonAdjointPower_m4_s3\":"<<total[2]-interpolatedDivergencePower
      <<",\"maximumInterpolatedDivergence_s1\":"<<maxInterpolatedDiv
      <<",\"maximumConservativeDivergence_s1\":"<<maxConservativeDiv
      <<",\"centralPressurePower_m4_s3\":"<<pressurePieces[0]
      <<",\"distanceWeightPressurePower_m4_s3\":"<<pressurePieces[1]
      <<",\"skewPressurePower_m4_s3\":"<<pressurePieces[2]
      <<",\"boundaryExtrapolationPressurePower_m4_s3\":"<<pressurePieces[3]
      <<",\"centralAdjointIdentityError_m4_s3\":"<<pressurePieces[0]-centralDivergencePower<<"}\n";
}
