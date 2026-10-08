// Native annulus diagnostic. Every update uses production FV and linear
// operators. Optional compact transpose blocks change only iteration splitting.
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include "FlowSolverDetail2D.hpp"
#include "FlowEquation2D.hpp"
#include "FlowPressure2D.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::fv::solver_detail;
Vector2D add(Vector2D a,Vector2D b){return {a.x+b.x,a.y+b.y};}

struct Step {
    FlowResult2D prediction, result;
    Vec ra, pc, predicted, correction, df;
    std::vector<Vector2D> gu, gv, gp, gauss, stress;
    double splittingResidualMaxError=0;
};

Step advance(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,
             const FlowInitialGuess2D& old,unsigned blockScope) {
    const auto n=m.cells.size(),nf=m.faces.size();
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    for(const auto& f:m.faces)if(f.neighbour)edges.emplace_back(f.owner,*f.neighbour);
    detail::SparsePattern2D pattern(n,edges);
    detail::LinearWorkspace2D work(n);
    System au(pattern),av(pattern),ap(pattern);
    Step s;s.result.u=old.u;s.result.v=old.v;s.result.p=old.p;s.result.flux=old.flux;
    FlowEquation2D equation(m,b.fixedP);
    const auto reconstruction=equation.reconstruct(s.result,c,b,b.p,old.flux);
    const auto& stencil=equation.pressureStencil();
    s.gu=reconstruction.gu;s.gv=reconstruction.gv;s.gp=reconstruction.gp;
    s.gauss=reconstruction.pressureForce;s.stress=reconstruction.stress;
    equation.assembleMomentum(au,av,s.result,c,b,reconstruction,old.flux,{});
    const double stop=.01*c.tolerance*c.speed*c.velocityRelaxation;
    for(std::size_t i=0;i<n;++i) {
        const double du=au.diag[i],dv=av.diag[i];
        au.diag[i]/=c.velocityRelaxation;av.diag[i]/=c.velocityRelaxation;
        au.rhs[i]+=(au.diag[i]-du)*old.u[i];av.rhs[i]+=(av.diag[i]-dv)*old.v[i];
    }
    s.ra.resize(n);for(std::size_t i=0;i<n;++i)s.ra[i]=m.cells[i].area/au.diag[i];
    if(blockScope) {
        std::vector<std::pair<std::size_t,std::size_t>> blockEdges;
        for(const auto& f:m.faces)if(f.neighbour)for(std::size_t a=0;a<2;++a)for(std::size_t k=0;k<2;++k)blockEdges.emplace_back(2*f.owner+a,2*(*f.neighbour)+k);
        for(std::size_t i=0;i<n;++i)blockEdges.emplace_back(2*i,2*i+1);
        detail::SparsePattern2D blockPattern(2*n,blockEdges);System a(blockPattern);
        detail::LinearWorkspace2D blockWork(2*n);Vec uv(2*n);
        for(std::size_t i=0;i<n;++i) {
            a.diag[2*i]=au.diag[i];a.diag[2*i+1]=av.diag[i];
            a.rhs[2*i]=au.rhs[i];a.rhs[2*i+1]=av.rhs[i];uv[2*i]=old.u[i];uv[2*i+1]=old.v[i];
            for(auto k=pattern.rows[i];k<pattern.rows[i+1];++k) {
                a.add(2*i,2*pattern.columns[k],au.off[k]);a.add(2*i+1,2*pattern.columns[k]+1,av.off[k]);
            }
        }
        for(const auto& f:m.faces)if(!f.neighbour||blockScope==2) {
            const auto i=f.owner;const double length=std::hypot(f.areaVector.x,f.areaVector.y);
            const auto normal=f.areaVector*(1/length);const double d=c.nu*f.transmissibility;
            const double xx=d*normal.x*normal.x,xy=d*normal.x*normal.y,yy=d*normal.y*normal.y;
            a.diag[2*i]+=xx;a.diag[2*i+1]+=yy;a.add(2*i,2*i+1,xy);a.add(2*i+1,2*i,xy);
            const double du=old.u[i]-(f.neighbour?old.u[*f.neighbour]:0),dv=old.v[i]-(f.neighbour?old.v[*f.neighbour]:0);
            a.rhs[2*i]+=xx*du+xy*dv;a.rhs[2*i+1]+=xy*du+yy*dv;
            if(f.neighbour) {
                const auto j=*f.neighbour;
                a.diag[2*j]+=xx;a.diag[2*j+1]+=yy;a.add(2*j,2*j+1,xy);a.add(2*j+1,2*j,xy);
                a.add(2*i,2*j,-xx);a.add(2*i,2*j+1,-xy);a.add(2*i+1,2*j,-xy);a.add(2*i+1,2*j+1,-yy);
                a.add(2*j,2*i,-xx);a.add(2*j,2*i+1,-xy);a.add(2*j+1,2*i,-xy);a.add(2*j+1,2*i+1,-yy);
                a.rhs[2*j]-=xx*du+xy*dv;a.rhs[2*j+1]-=xy*du+yy*dv;
            }
        }
        for(std::size_t i=0;i<n;++i) {
            s.splittingResidualMaxError=std::max(s.splittingResidualMaxError,std::abs(a.compensatedResidualRow(2*i,uv)-au.compensatedResidualRow(i,old.u)));
            s.splittingResidualMaxError=std::max(s.splittingResidualMaxError,std::abs(a.compensatedResidualRow(2*i+1,uv)-av.compensatedResidualRow(i,old.v)));
        }
        a.solve(uv,blockWork,stop);
        for(std::size_t i=0;i<n;++i){s.result.u[i]=uv[2*i];s.result.v[i]=uv[2*i+1];}
    } else {au.solve(s.result.u,work,stop);av.solve(s.result.v,work,stop);}
    s.prediction=s.result;
    const auto gup=equation.velocityGradient(s.result.u,b,false),gvp=equation.velocityGradient(s.result.v,b,true);
    s.df.resize(nf);
    s.predicted=rhieChowFlux(m,c,b,s.result,s.ra,old.u,old.v,s.gu,s.gv,gup,gvp,s.gp,s.gauss,{},false,0,s.df);
    s.pc.resize(n);const Vec zeros(nf);
    s.correction=solvePressureCorrection(m,c,b,s.result,ap,s.pc,s.ra,s.df,s.predicted,stencil,zeros,
        [&](const System& a,Vec& x,bool){return a.solvePressure(x,work,detail::LinearPressureMethod2D::IC0);});
    const auto gc=detail::conservativePressureGradient(m,detail::pressureFaceValues(m,s.pc,stencil.apply(s.pc,zeros),zeros,b.fixedP));
    for(std::size_t i=0;i<n;++i){s.result.u[i]-=s.ra[i]*gc[i].x;s.result.v[i]-=s.ra[i]*gc[i].y;s.result.p[i]+=c.pressureRelaxation*s.pc[i];}
    for(std::size_t id=0;id<nf;++id){const auto& f=m.faces[id];s.result.flux[id]=s.predicted[id];if(f.neighbour||b.fixedP[id])s.result.flux[id]+=s.df[id]*(s.pc[f.owner]-(f.neighbour?s.pc[*f.neighbour]:0))+s.correction[id];}
    return s;
}

void terms(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,
           const FlowInitialGuess2D& x,std::size_t cell,const std::string& label) {
    const auto gu=detail::flowGradient(m,x.u,b.u,b.fixedU),gv=detail::flowGradient(m,x.v,b.v,b.fixedV);
    const auto gp=detail::flowGradient(m,x.p,b.p,b.fixedP,true);
    const auto pf=detail::pressureFaceValues(m,x.p,gp,b.p,b.fixedP);
    const auto gg=detail::conservativePressureGradient(m,pf);
    const auto stress=detail::symmetricViscousCorrection(m,x.u,x.v,gu,gv,b.u,b.v,b.fixedU,b.fixedV,b.constantU,b.constantV,c.nu);
    std::array<Vector2D,5> sum{};
    std::cout<<"face_terms,"<<label<<",face,other,wall,transmissibility,normalDistance,skew,lap_x,lap_y,nonorth_x,nonorth_y,transpose_x,transpose_y,pressure_x,pressure_y,normal_x,normal_y\n";
    for(const auto id:m.cells[cell].faces) {
        const auto& f=m.faces[id];const double sign=f.owner==cell?1:-1;
        const auto other=f.neighbour?(f.owner==cell?*f.neighbour:f.owner):cell;
        const double d=c.nu*f.transmissibility,length=std::hypot(f.areaVector.x,f.areaVector.y);
        const Vector2D normal=f.areaVector*(sign/length);
        const double ou=f.neighbour?x.u[other]:b.u[id],ov=f.neighbour?x.v[other]:b.v[id];
        const Vector2D lap{d*(x.u[cell]-ou),d*(x.v[cell]-ov)};
        const Vector2D nonorth{-sign*c.nu*dot(interpolateGradient(f,gu),f.correction),-sign*c.nu*dot(interpolateGradient(f,gv),f.correction)};
        const Vector2D transpose=stress[id]*sign;
        const Vector2D pressure=f.areaVector*(sign*pf[id]);
        const std::size_t base=f.neighbour?0:2;
        sum[base]=add(add(sum[base],lap),nonorth);sum[base+1]=add(sum[base+1],transpose);sum[4]=add(sum[4],pressure);
        const auto delta=(f.neighbour?m.cells[*f.neighbour].centre:f.centre)-m.cells[f.owner].centre;
        const double distance=dot(delta,f.areaVector*(1/length));
        const auto ownCentre=m.cells[f.owner].centre,otherCentre=f.neighbour?m.cells[*f.neighbour].centre:ownCentre;
        const Point2D mid{ownCentre.x*(1-f.neighbourWeight)+otherCentre.x*f.neighbourWeight,ownCentre.y*(1-f.neighbourWeight)+otherCentre.y*f.neighbourWeight};
        const double skew=std::hypot(f.centre.x-mid.x,f.centre.y-mid.y);
        std::cout<<"face,"<<label<<','<<id<<','<<(f.neighbour?static_cast<long long>(other):-1)<<','<<!f.neighbour<<','<<f.transmissibility<<','<<distance<<','<<skew<<','<<lap.x<<','<<lap.y<<','<<nonorth.x<<','<<nonorth.y<<','<<transpose.x<<','<<transpose.y<<','<<pressure.x<<','<<pressure.y<<','<<normal.x<<','<<normal.y<<'\n';
    }
    const char* names[]={"interior_laplacian","interior_transpose","wall_laplacian","wall_transpose","gauss_pressure"};
    for(std::size_t i=0;i<sum.size();++i)std::cout<<"term,"<<label<<','<<names[i]<<','<<sum[i].x<<','<<sum[i].y<<'\n';
    std::cout<<"gradient,"<<label<<",LS,"<<gp[cell].x<<','<<gp[cell].y<<"\ngradient,"<<label<<",Gauss,"<<gg[cell].x<<','<<gg[cell].y<<'\n';
    Vector2D total{};for(const auto& v:sum)total=add(total,v);
    std::cout<<"residual,"<<label<<','<<total.x<<','<<total.y<<'\n';
}

FlowInitialGuess2D guess(const FlowResult2D& r){return {r.u,r.v,r.p,r.flux};}
double product(const FvMesh2D& m,const FlowInitialGuess2D& a,const FlowInitialGuess2D& b) {
    long double sum=0,area=0;for(std::size_t i=0;i<m.cells.size();++i){area+=m.cells[i].area;sum+=m.cells[i].area*(a.u[i]*b.u[i]+a.v[i]*b.v[i]+a.p[i]*b.p[i]);}return static_cast<double>(sum/area);
}

int main(int argc,char** argv) {
    if(argc<3||argc>8){std::cerr<<"Usage: momentum_terms_probe MESH BOUNDARIES [steps=120] [block=0:none,1:wall,2:all compact transpose] [cell=3710, zero based] [weight=1] [power|steady]\n";return 2;}
    try {
        const auto input=readCm2dTopology(argv[1]);ensure(input.valid(),input.error.c_str());const auto m=makeFvMesh2D(input.topology);
        FlowControls2D c;c.scenario="custom";c.nu=.1;c.speed=.5;c.tolerance=1e-8;c.maxIterations=1;c.momentumInertia=0;
        std::ifstream in(argv[2]);c.boundaryConditions=readFlowBoundaryConditions2D(in,m,c);const auto b=boundaries(m,c);
        ensure(b.closed,"Probe requires closed velocity-boundary domain");for(std::size_t id=0;id<m.faces.size();++id)ensure(m.faces[id].neighbour||(b.fixedU[id]&&b.fixedV[id]),"Probe requires prescribed wall velocity");
        const auto n=m.cells.size(),nf=m.faces.size(),steps=argc>3?std::stoul(argv[3]):120UL,cell=argc>5?std::stoul(argv[5]):3710UL;const unsigned block=argc>4?static_cast<unsigned>(std::stoul(argv[4])):0;const double weight=argc>6?std::stod(argv[6]):1.;const bool steady=argc>7&&std::string(argv[7])=="steady";ensure(cell<n,"Cell out of range");ensure(block<=2,"Unknown block scope");ensure(weight>0&&weight<=1,"Weight outside (0,1]");
        FlowInitialGuess2D zero;zero.u.resize(n);zero.v.resize(n);zero.p.resize(n);zero.flux.resize(nf);
        std::cout<<std::setprecision(17)<<"mesh,"<<n<<','<<nf<<",cell,"<<cell<<",area,"<<m.cells[cell].area<<",centre,"<<m.cells[cell].centre.x<<','<<m.cells[cell].centre.y<<'\n';
        const auto first=advance(m,c,b,zero,false);const auto native=solveIncompressibleFromGuess2D(m,c,zero);
        double delta=0;for(const auto& pair:{std::pair{&first.result.u,&native.u},std::pair{&first.result.v,&native.v},std::pair{&first.result.p,&native.p},std::pair{&first.result.flux,&native.flux}})for(std::size_t i=0;i<pair.first->size();++i)delta=std::max(delta,std::abs((*pair.first)[i]-(*pair.second)[i]));
        std::cout<<"production_step_max_error,"<<delta<<"\nfirst_predictor,"<<first.prediction.u[cell]<<','<<first.prediction.v[cell]<<"\nfirst_pc,"<<first.pc[cell]<<"\nfirst_accepted,"<<first.result.u[cell]<<','<<first.result.v[cell]<<','<<first.result.p[cell]<<'\n';
        terms(m,c,b,guess(first.result),cell,"first");
        if(steady) {
            auto x=zero;
            std::vector<std::pair<std::size_t,std::size_t>> edges;for(const auto& f:m.faces)if(f.neighbour)edges.emplace_back(f.owner,*f.neighbour);
            detail::SparsePattern2D pattern(n,edges);System au(pattern),av(pattern);
            std::cout<<"steady_header,iteration,momentum,continuity,velocityChange,pressureChange,worstCell\n";
            for(std::size_t it=1;it<=steps;++it) {
                const auto s=advance(m,c,b,x,block);auto next=guess(s.result);
                for(std::size_t i=0;i<n;++i){next.u[i]=(1-weight)*x.u[i]+weight*next.u[i];next.v[i]=(1-weight)*x.v[i]+weight*next.v[i];next.p[i]=(1-weight)*x.p[i]+weight*next.p[i];}
                for(std::size_t i=0;i<nf;++i)next.flux[i]=(1-weight)*x.flux[i]+weight*next.flux[i];
                const auto gu=detail::flowGradient(m,next.u,b.u,b.fixedU),gv=detail::flowGradient(m,next.v,b.v,b.fixedV),gp=detail::flowGradient(m,next.p,b.p,b.fixedP,true);
                const auto gg=detail::conservativePressureGradient(m,detail::pressureFaceValues(m,next.p,gp,b.p,b.fixedP));
                const auto stress=detail::symmetricViscousCorrection(m,next.u,next.v,gu,gv,b.u,b.v,b.fixedU,b.fixedV,b.constantU,b.constantV,c.nu);
                momentum(au,m,c,b,next.u,next.flux,gu,gg,{},stress,{},false,nullptr,0);momentum(av,m,c,b,next.v,next.flux,gv,gg,{},stress,{},true,nullptr,0);
                Vec div(n);for(std::size_t id=0;id<nf;++id){const auto& f=m.faces[id];div[f.owner]+=next.flux[id];if(f.neighbour)div[*f.neighbour]-=next.flux[id];}
                double mr=0,continuity=0,du=0,dp=0;std::size_t worst=0;
                for(std::size_t i=0;i<n;++i){const double residual=std::hypot(au.compensatedResidualRow(i,next.u),av.compensatedResidualRow(i,next.v))/((au.diag[i]+av.diag[i])*c.speed);if(residual>mr){mr=residual;worst=i;}continuity=std::max(continuity,std::abs(div[i])/(c.speed*std::sqrt(m.cells[i].area)));du=std::max(du,std::hypot(next.u[i]-x.u[i],next.v[i]-x.v[i])/c.speed);dp=std::max(dp,std::abs(next.p[i]-x.p[i])/(c.speed*c.speed));}
                const bool pass=mr<c.tolerance&&continuity<c.tolerance&&du<c.tolerance&&dp<c.tolerance;
                if(it==1||it%100==0||it==steps||pass)std::cout<<"steady,"<<it<<','<<mr<<','<<continuity<<','<<du<<','<<dp<<','<<worst<<'\n';
                x=std::move(next);
                if(pass||it==steps){terms(m,c,b,x,cell,"steady_last");std::cout<<"diagnostic_steady_threshold,"<<pass<<'\n';return pass?0:1;}
                if(mr>1e5){std::cout<<"diagnostic_steady_divergence,"<<it<<'\n';return 1;}
            }
        }
        auto homogeneous=b;std::fill(homogeneous.u.begin(),homogeneous.u.end(),0);std::fill(homogeneous.v.begin(),homogeneous.v.end(),0);std::fill(homogeneous.p.begin(),homogeneous.p.end(),0);
        auto x=zero;for(std::size_t i=0;i<n;++i){x.u[i]=std::sin(1.71*static_cast<double>(i));x.v[i]=std::cos(1.53*static_cast<double>(i));x.p[i]=std::sin(1.31*static_cast<double>(i));}x.p[0]=0;
        std::cout<<"power_header,iteration,normRatio,rayleigh,eigenResidual,worstCell\n";
        for(std::size_t it=1;it<=steps;++it) {
            const auto s=advance(m,c,homogeneous,x,block);auto next=guess(s.result);
            for(std::size_t i=0;i<n;++i){next.u[i]=(1-weight)*x.u[i]+weight*next.u[i];next.v[i]=(1-weight)*x.v[i]+weight*next.v[i];next.p[i]=(1-weight)*x.p[i]+weight*next.p[i];}
            for(std::size_t i=0;i<nf;++i)next.flux[i]=(1-weight)*x.flux[i]+weight*next.flux[i];
            const double norm=std::sqrt(product(m,next,next)),ray=product(m,x,next)/product(m,x,x);
            auto residual=next;for(std::size_t i=0;i<n;++i){residual.u[i]-=ray*x.u[i];residual.v[i]-=ray*x.v[i];residual.p[i]-=ray*x.p[i];}
            std::size_t worst=0;for(std::size_t i=1;i<n;++i)if(std::hypot(next.u[i],next.v[i])>std::hypot(next.u[worst],next.v[worst]))worst=i;
            if(it==1||it%10==0||it==steps)std::cout<<"power,"<<it<<','<<norm/std::sqrt(product(m,x,x))<<','<<ray<<','<<std::sqrt(product(m,residual,residual))/norm<<','<<worst<<'\n';
            if(it==steps) {
                std::cout<<"same_field_splitting_residual_max_error,"<<s.splittingResidualMaxError<<'\n';
                std::cout<<"mode_old,"<<x.u[cell]<<','<<x.v[cell]<<','<<x.p[cell]<<"\nmode_predictor,"<<s.prediction.u[cell]<<','<<s.prediction.v[cell]<<"\nmode_pc,"<<s.pc[cell]<<"\nmode_accepted,"<<next.u[cell]<<','<<next.v[cell]<<','<<next.p[cell]<<'\n';
                terms(m,c,homogeneous,x,cell,"mode_old");terms(m,c,homogeneous,next,cell,"mode_accepted");
                const auto gup=detail::flowGradient(m,s.prediction.u,homogeneous.u,homogeneous.fixedU),gvp=detail::flowGradient(m,s.prediction.v,homogeneous.v,homogeneous.fixedV);
                for(const auto id:m.cells[cell].faces){const auto& f=m.faces[id];if(!f.neighbour)continue;const auto i=f.owner,j=*f.neighbour;const double w=f.neighbourWeight;const Point2D point{m.cells[i].centre.x*(1-w)+m.cells[j].centre.x*w,m.cells[i].centre.y*(1-w)+m.cells[j].centre.y*w};const auto skew=f.centre-point;const double plain=interpolate(f,s.prediction.u)*f.areaVector.x+interpolate(f,s.prediction.v)*f.areaVector.y;const double velocitySkew=dot(interpolateGradient(f,gup),skew)*f.areaVector.x+dot(interpolateGradient(f,gvp),skew)*f.areaVector.y;const Vector2D rag{(1-w)*s.ra[i]*s.gauss[i].x+w*s.ra[j]*s.gauss[j].x,(1-w)*s.ra[i]*s.gauss[i].y+w*s.ra[j]*s.gauss[j].y};const double force=dot(rag,f.areaVector),pressure=-interpolate(f,s.ra)*(f.transmissibility*(x.p[j]-x.p[i])+dot(interpolateGradient(f,s.gp),f.correction));const double relaxation=s.predicted[id]-plain-velocitySkew-force-pressure;std::cout<<"rc,"<<id<<",plain,"<<plain<<",skew,"<<velocitySkew<<",gauss,"<<force<<",compact_pressure,"<<pressure<<",relaxation,"<<relaxation<<",predicted,"<<s.predicted[id]<<",accepted,"<<next.flux[id]<<'\n';}
            }
            for(auto* f:{&next.u,&next.v,&next.p,&next.flux})for(auto& v:*f)v/=norm;
            x=std::move(next);
        }
        return 0;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
