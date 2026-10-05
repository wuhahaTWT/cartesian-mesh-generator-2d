#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// Evaluate the product reconstruction on the accepted native mesh and fields.
// This is an operator diagnostic, not an independent discretized-equation audit.
int main(int argc,char** argv) {
    using namespace cartmesh2d;using namespace cartmesh2d::fv;
    if(argc!=3)throw std::runtime_error("usage: pressure-probe mesh.cm2d fields.cells.csv");
    const auto read=readCm2dTopology(argv[1]);if(!read.valid())throw std::runtime_error(read.error);
    const auto m=makeFvMesh2D(read.topology);
    std::vector<double> p(m.cells.size()),exact(p.size()),bc(m.faces.size());
    std::vector<bool> fixed(bc.size());
    std::ifstream in(argv[2]);if(!in)throw std::runtime_error("cannot read field");
    std::string line;std::getline(in,line);std::size_t count=0;
    while(std::getline(in,line)) {
        std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
        std::size_t i;double x,y,a,u,v;if(!(row>>i>>x>>y>>a>>u>>v)||i!=count||!(row>>p.at(i)))throw std::runtime_error("invalid field row");++count;
    }
    if(count!=p.size())throw std::runtime_error("field size mismatch");
    const auto reference=[](Point2D x){const double r2=x.x*x.x+x.y*x.y;return r2/18.-std::log(r2)/9.-1./(18.*r2);};
    double offset=0,area=0;for(std::size_t i=0;i<p.size();++i){exact[i]=reference(m.cells[i].centre);offset+=m.cells[i].area*(p[i]-exact[i]);area+=m.cells[i].area;}offset/=area;
    const auto stencil=detail::buildFlowGradientStencil2D(m,fixed,true);
    const auto ge=stencil.apply(exact,bc),gn=stencil.apply(p,bc);
    const auto fe=detail::pressureFaceValues(m,exact,ge,bc,fixed),fn=detail::pressureFaceValues(m,p,gn,bc,fixed);
    const auto de=detail::conservativePressureGradient(m,fe),dn=detail::conservativePressureGradient(m,fn);
    std::vector<double> faceExact(bc.size());for(std::size_t f=0;f<bc.size();++f)faceExact[f]=reference(m.faces[f].centre);
    const auto integrated=detail::conservativePressureGradient(m,faceExact);
    std::vector<double> ue(p.size()),ve(p.size()),bu(bc.size()),bv(bc.size());std::vector<bool> vf(bc.size()),constant(bc.size());
    for(std::size_t i=0;i<p.size();++i){auto x=m.cells[i].centre;const double q=(1./(x.x*x.x+x.y*x.y)-1)/3.;ue[i]=-x.y*q;ve[i]=x.x*q;}
    for(const auto& b:rotatingAnnulusBoundaryPreset2D(m,.5)){bu[b.face]=b.velocity.x;bv[b.face]=b.velocity.y;vf[b.face]=true;constant[b.face]=b.kind!=FlowBoundaryKind2D::SmoothMovingWall;}
    const auto vs=detail::buildFlowGradientStencil2D(m,vf);
    const auto viscousError=[&](const std::vector<double>& wallU,const std::vector<double>& wallV,const std::vector<bool>& trace){
    const auto gu=vs.apply(ue,wallU),gv=vs.apply(ve,wallV);
    std::vector<Vector2D> viscous(p.size());
    for(std::size_t f=0;f<bc.size();++f){
        const auto a=detail::viscousFaceGradient(m,f,ue,gu,wallU,vf,trace),b=detail::viscousFaceGradient(m,f,ve,gv,wallV,vf,trace);
        const auto& face=m.faces[f];const auto x=face.centre;const double k=1./(3*std::pow(x.x*x.x+x.y*x.y,2));
        const double ax=2*x.x*x.y*k,by=-ax,aybx=2*(x.y*x.y-x.x*x.x)*k;
        const Vector2D error{.1*(2*(a.x-ax)*face.areaVector.x+(a.y+b.x-aybx)*face.areaVector.y),.1*((a.y+b.x-aybx)*face.areaVector.x+2*(b.y-by)*face.areaVector.y)};
        viscous[face.owner].x+=error.x;viscous[face.owner].y+=error.y;
        if(face.neighbour){viscous[*face.neighbour].x-=error.x;viscous[*face.neighbour].y-=error.y;}
    }
    return viscous;};
    const auto viscous=viscousError(bu,bv,constant);
    // Counterfactual operator consistency only: analytic circular values on
    // polygon faces can have nonzero normal velocity and are NOT valid wall BCs.
    double normalLeak=0,traceDifference=0;
    for(std::size_t f=0;f<bc.size();++f){
        const auto& face=m.faces[f];const auto x=face.centre;const double q=(1./(x.x*x.x+x.y*x.y)-1)/3.;
        if(!face.neighbour){
            traceDifference=std::max(traceDifference,std::hypot(bu[f]+x.y*q,bv[f]-x.x*q));
            normalLeak=std::max(normalLeak,std::abs(-x.y*q*face.areaVector.x+x.x*q*face.areaVector.y)/std::hypot(face.areaVector.x,face.areaVector.y));
        }
        bu[f]=-x.y*q;bv[f]=x.x*q;
    }
    std::cerr<<std::setprecision(17)<<"{\"counterfactualMaxNormalVelocity_m_s\":"<<normalLeak<<",\"maximumWallTraceDifference_m_s\":"<<traceDifference<<"}\n";
    const auto smoothViscous=viscousError(bu,bv,std::vector<bool>(bc.size()));
    std::cout<<std::setprecision(17)<<"cell,x,y,area,p_error,samples,condition,wall_faces,exact_reconstruction_error,computed_gradient_error,ls_exact_error,velocity_condition,exact_viscous_error,counterfactual_smooth_trace_viscous_error\n";
    for(std::size_t i=0;i<p.size();++i){
        const auto& c=m.cells[i];const auto& row=stencil.rows[i];const double lambda=.5*(row.xx+row.yy+std::hypot(row.xx-row.yy,2*row.xy));
        std::size_t walls=0;for(auto f:c.faces)walls+=!m.faces[f].neighbour;
        const auto norm=[](Vector2D a,Vector2D b){return std::hypot(a.x-b.x,a.y-b.y);};
        const auto& vr=vs.rows[i];const double vl=.5*(vr.xx+vr.yy+std::hypot(vr.xx-vr.yy,2*vr.xy));
        std::cout<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','<<p[i]-exact[i]-offset<<','<<row.end-row.begin<<','<<lambda*lambda/row.det<<','<<walls<<','<<norm(de[i],integrated[i])<<','<<norm(dn[i],integrated[i])<<','<<norm(ge[i],integrated[i])<<','<<vl*vl/vr.det<<','<<std::hypot(viscous[i].x,viscous[i].y)/c.area<<','<<std::hypot(smoothViscous[i].x,smoothViscous[i].y)/c.area<<'\n';
    }
}
