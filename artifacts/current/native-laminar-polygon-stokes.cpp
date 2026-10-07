// Fixed actual polygon, smooth nonzero no-slip Stokes reference. Research only.
// Reuse the native operator/API and existing mesh/state adapters; no Python PDE.
#define CARTMESH_P1_OSEEN_NO_MAIN
#include "native-laminar-p1-oseen.cpp"
#include "cartmesh2d/fv/CompatibleIncompressible2D.hpp"

namespace polygon_stokes {
Vector2D operator+(Vector2D a,Vector2D b){return {a.x+b.x,a.y+b.y};}
Vector2D operator-(Vector2D a,Vector2D b){return {a.x-b.x,a.y-b.y};}
// Total-order-three Taylor coefficients, c_ij = d_x^i d_y^j f/(i!j!).
// Analytic polynomial differentiation; no finite-difference reference loads.
struct Jet {
    double c[4][4]{};
    explicit Jet(double value=0){c[0][0]=value;}
    static Jet affine(double value,double dx,double dy){Jet a(value);a.c[1][0]=dx;a.c[0][1]=dy;return a;}
};
Jet operator*(const Jet& a,const Jet& b){Jet r;
    for(int i=0;i<=3;++i)for(int j=0;j<=3-i;++j)
        for(int k=0;k<=i;++k)for(int l=0;l<=j;++l)r.c[i][j]+=a.c[k][l]*b.c[i-k][j-l];
    return r;
}
void checkJet(){
    // x^2*y^3 at (2,3) has explicit derivatives through total order three.
    const auto x=Jet::affine(2,1,0),y=Jet::affine(3,0,1),v=x*x*y*y*y;
    const std::array<std::array<double,4>,4> expected{{{{108,108,36,4}},{{108,108,36,0}},{{27,27,0,0}},{{0,0,0,0}}}};
    for(int i=0;i<=3;++i)for(int j=0;j<=3-i;++j)if(v.c[i][j]!=expected[i][j])throw std::runtime_error("analytic Taylor product control failed");
}
struct Reference {
    std::vector<Point2D> vertices;std::vector<std::array<double,3>> lines;
    double xmin=0,xmax=0,ymin=0,ymax=0,nu=1,amplitude=1,solidArea=0;
    explicit Reference(const std::string& xy,double padding,double viscosity):nu(viscosity){
        std::ifstream in(xy);Point2D v;while(in>>v.x>>v.y)vertices.push_back(v);
        if(!in.eof()||vertices.size()!=8||!std::isfinite(padding)||padding<=0)throw std::runtime_error("reference requires one eight-vertex polygon");
        if(!BoundaryLoop(vertices).diagnose().valid())throw std::runtime_error("invalid analytic polygon");
        xmin=xmax=vertices[0].x;ymin=ymax=vertices[0].y;
        double twiceArea=0;
        for(std::size_t i=0;i<vertices.size();++i){const auto a=vertices[i],b=vertices[(i+1)%vertices.size()];
            xmin=std::min(xmin,a.x);xmax=std::max(xmax,a.x);ymin=std::min(ymin,a.y);ymax=std::max(ymax,a.y);
            const double length=std::hypot(b.x-a.x,b.y-a.y);if(!(length>0))throw std::runtime_error("zero polygon edge");
            const double nx=(b.y-a.y)/length,ny=(a.x-b.x)/length;lines.push_back({nx,ny,-nx*a.x-ny*a.y});twiceArea+=a.x*b.y-b.x*a.y;
        }
        if(!(twiceArea>0))throw std::runtime_error("reference requires CCW polygon");
        solidArea=twiceArea/2;const double pad=padding*std::max(xmax-xmin,ymax-ymin);xmin-=pad;xmax+=pad;ymin-=pad;ymax+=pad;
        for(const auto& l:lines)for(auto p:vertices)if(l[0]*p.x+l[1]*p.y+l[2]>TolerancePolicy{}.scale(std::max(xmax-xmin,ymax-ymin)))throw std::runtime_error("reference requires convex polygon");
    }
    Exact at(Point2D p)const{
        // All line distances use the fixed 1 m reference length. Q is dimensionless.
        // psi = 1 m^2/s * [product_i l_i * Q]^2; amplitude is never fit to errors.
        Jet a(1);for(const auto& l:lines)a=a*Jet::affine(l[0]*p.x+l[1]*p.y+l[2],l[0],l[1]);
        const double hx=(xmax-xmin)/2,hy=(ymax-ymin)/2;
        a=a*Jet::affine((p.x-xmin)/hx,1/hx,0)*Jet::affine((xmax-p.x)/hx,-1/hx,0)
           *Jet::affine((p.y-ymin)/hy,0,1/hy)*Jet::affine((ymax-p.y)/hy,0,-1/hy);
        const auto psi=Jet(amplitude)*a*a;
        const double lapU=2*psi.c[2][1]+6*psi.c[0][3],lapV=-6*psi.c[3][0]-2*psi.c[1][2];
        return {{psi.c[0][1],-psi.c[1][0]},{2-nu*lapU,-3-nu*lapV},1+2*p.x-3*p.y,
                {{{psi.c[1][1],2*psi.c[0][2]},{-2*psi.c[2][0],-psi.c[1][1]}}}};
    }
    Vector2D traction(Point2D p,Vector2D n,double gauge)const{const auto e=at(p);const auto u=e.gradient[0],v=e.gradient[1];
        return {(2*nu*u.x-e.p+gauge)*n.x+nu*(u.y+v.x)*n.y,nu*(u.y+v.x)*n.x+(2*nu*v.y-e.p+gauge)*n.y};}
};
struct WallError {
    std::size_t faces=0;double length=0,pCell=0,pReaction=0,traction=0,shear=0,tractionMoment=0,projection=0,pCellMax=0,pReactionMax=0;
    Vector2D force{},exactForce{};
    void write(std::ostream& out)const{const auto rms=[&](double x){return std::sqrt(x/length);};
        out<<"{\"faces\":"<<faces<<",\"length\":"<<length<<",\"cellPressureRms\":"<<rms(pCell)<<",\"reactionPressureRms\":"<<rms(pReaction)
           <<",\"cellPressureEndpointMaximum\":"<<pCellMax<<",\"reactionPressureEndpointMaximum\":"<<pReactionMax<<",\"tractionRms\":"<<rms(traction)
           <<",\"shearRms\":"<<rms(shear)<<",\"tractionMomentRms\":"<<rms(tractionMoment)<<",\"bestP1TractionRms\":"<<rms(projection)
           <<",\"forceOnFluid\":["<<force.x<<','<<force.y<<"],\"exactForceOnFluid\":["<<exactForce.x<<','<<exactForce.y<<"]}";
    }
};
void checkedClose(std::ofstream& out){out.close();if(!out)throw std::runtime_error("reference evidence output failed");}
}
#ifndef CARTMESH_POLYGON_STOKES_NO_MAIN
int main(int argc,char** argv)try {
    using namespace polygon_stokes;
    if(argc!=7)throw std::runtime_error("usage: polygon-stokes mesh.solver.cm2d polygon.xy padding_fraction solve_order error_order fresh_prefix");
    checkJet();const auto start=std::chrono::steady_clock::now();const std::string prefix=argv[6];
    const int order=std::stoi(argv[4]),errorOrder=std::stoi(argv[5]);if(order<4||order>12||errorOrder<4||errorOrder>16)throw std::runtime_error("quadrature range");
    for(const auto* suffix:{".json",".cells.csv",".walls.csv",".iterations.csv",".accepted.state",".checkpoint"})
        if(std::filesystem::exists(prefix+suffix))throw std::runtime_error("existing reference evidence: "+prefix+suffix);
    const auto fixture=readFixture(argv[1],0);const auto& mesh=fixture.mesh;const Reference exact(argv[2],std::stod(argv[3]),1);
    CompatibleFlowControls2D c;c.viscosity=exact.nu;c.equation=CompatibleEquation2D::Stokes;c.globalization=CompatibleGlobalization2D::Backtracking;c.quadratureOrder=order;
    c.acceleration=[&](Point2D p){return exact.at(p).f;};
    double traceMaximum=0,normalViscousStressMaximum=0;std::array<double,4> domain{1e100,-1e100,1e100,-1e100};
    for(std::size_t id=0;id<mesh.faces.size();++id){const auto& face=mesh.faces[id];if(face.neighbour)continue;
        if(face.patch!=BoundaryPatch2D::EmbeddedBoundary&&face.patch!=BoundaryPatch2D::DomainBoundary)throw std::runtime_error("unclassified reference wall");
        if(face.patch==BoundaryPatch2D::EmbeddedBoundary){
            const Point2D a{face.centre.x+face.areaVector.y/2,face.centre.y-face.areaVector.x/2},b{face.centre.x-face.areaVector.y/2,face.centre.y+face.areaVector.x/2};
            bool onOriginalSegment=false;for(std::size_t j=0;j<exact.vertices.size();++j){const Segment2D segment{exact.vertices[j],exact.vertices[(j+1)%exact.vertices.size()]};
                onOriginalSegment=onOriginalSegment||(pointOnSegment(a,segment)&&pointOnSegment(b,segment));}
            if(!onOriginalSegment)throw std::runtime_error("embedded wall is not on one original polygon segment");
        }
        CompatibleBoundary2D b;b.face=id;c.boundaries.push_back(b);const double length=std::hypot(face.areaVector.x,face.areaVector.y);const auto n=face.areaVector*(1/length);
        auto samples=gauss(errorOrder);samples.push_back({0,0});samples.push_back({1,0});
        for(auto [z,w]:samples){(void)w;const double s=z-.5;const Point2D p{face.centre.x-s*face.areaVector.y,face.centre.y+s*face.areaVector.x};const auto e=exact.at(p);
            traceMaximum=std::max(traceMaximum,std::hypot(e.u.x,e.u.y));normalViscousStressMaximum=std::max(normalViscousStressMaximum,std::abs(2*exact.nu*(n.x*n.x*e.gradient[0].x+n.x*n.y*(e.gradient[0].y+e.gradient[1].x)+n.y*n.y*e.gradient[1].y)));
            if(face.patch==BoundaryPatch2D::DomainBoundary){domain[0]=std::min(domain[0],p.x);domain[1]=std::max(domain[1],p.x);domain[2]=std::min(domain[2],p.y);domain[3]=std::max(domain[3],p.y);}
        }
    }
    const std::array<double,4> referenceDomain{exact.xmin,exact.xmax,exact.ymin,exact.ymax};
    for(std::size_t j=0;j<4;++j)if(!TolerancePolicy{}.nearlyEqual(domain[j],referenceDomain[j],std::max(exact.xmax-exact.xmin,exact.ymax-exact.ymin)))throw std::runtime_error("mesh/reference outer domain mismatch");
    // Match the actual domain and retained fluid area using the existing geometric
    // tolerance policy. Report trace roundoff, without inventing a CFD error gate.
    double meshArea=0;for(const auto& cell:mesh.cells)meshArea+=cell.area;
    const double expectedArea=(exact.xmax-exact.xmin)*(exact.ymax-exact.ymin)-exact.solidArea;
    if(std::abs(meshArea-expectedArea)>TolerancePolicy{}.areaScale(std::max(exact.xmax-exact.xmin,exact.ymax-exact.ymin)))throw std::runtime_error("reference fluid area mismatch");
    std::ofstream iterations(prefix+".iterations.csv");iterations<<std::setprecision(17)<<"iteration,accepted,restarts,products,linear_residual,cell_momentum,face_momentum,divergence,state_change\n";
    const auto r=solveCompatibleIncompressible2D(mesh,c);const auto solveEnd=std::chrono::steady_clock::now();
    for(const auto& it:r.iterations){iterations<<it.iteration<<','<<it.accepted<<','<<it.linearRestarts<<','<<it.matrixProducts<<','<<it.linearRelativeResidual;
        if(it.metrics){const auto& m=*it.metrics;iterations<<','<<m.cellMomentum<<','<<m.faceMomentum<<','<<m.divergence<<','<<m.stateChange;}else iterations<<",,,,";iterations<<'\n';}checkedClose(iterations);
    auto save=[&](const auto& state,const std::string& name){if(!state)return;Vec values;for(const auto& cell:state->cells)values.insert(values.end(),cell.begin(),cell.end());for(const auto& face:state->faces)values.insert(values.end(),face.begin(),face.end());writeState(prefix+name+".state",fixture,values);};
    save(r.lastAccepted,".accepted");save(r.lastRejected,".rejected");if(r.checkpoint){std::ofstream cp(prefix+".checkpoint");writeCompatibleFlowCheckpoint2D(cp,*r.checkpoint);checkedClose(cp);}
    std::ofstream report(prefix+".json");report<<std::setprecision(17)<<"{\"converged\":"<<(r.converged()?"true":"false")<<",\"stop\":"<<int(r.stop)<<",\"reason\":\""<<r.reason
        <<"\",\"cells\":"<<mesh.cells.size()<<",\"faces\":"<<mesh.faces.size()<<",\"meshKey\":"<<meshKey(fixture)<<",\"solveQuadrature\":"<<order<<",\"errorQuadrature\":"<<errorOrder
        <<",\"viscosity\":"<<c.viscosity<<",\"amplitude\":"<<exact.amplitude<<",\"equationTolerance\":"<<c.equationTolerance<<",\"stateTolerance\":"<<c.stateTolerance
        <<",\"linearTolerance\":"<<c.linearTolerance<<",\"maximumLinearRestarts\":"<<c.maximumLinearRestarts<<",\"krylovDirections\":"<<c.krylovDirections
        <<",\"maximumIterations\":"<<c.maximumIterations<<",\"boundaryReferenceVelocityMaximum\":"<<traceMaximum<<",\"boundaryReferenceNormalViscousStressMaximum\":"<<normalViscousStressMaximum
        <<",\"solveAndReadSeconds\":"<<std::chrono::duration<double>(solveEnd-start).count();
    if(!r.converged()){report<<"}\n";checkedClose(report);std::cout<<"not converged; original failed state and iterations retained\n";return 2;}
    const auto& state=*r.lastAccepted;const auto loads=evaluateCompatibleFlowLoads2D(mesh,c,*r.checkpoint);
    const double gauge=exact.at(mesh.cells.back().centre).p;
    double area=0,uError=0,uP1Error=0,pError=0,pMeanError=0,pMaximum=0,gradientError=0,uNorm=0,pNorm=0,minArea=1e100,maxDiameter=0;
    std::ofstream cells(prefix+".cells.csv");cells<<std::setprecision(17)<<"cell,x,y,area,diameter,u0,ux,uy,v0,vx,vy,p0,px,py,velocity_p2_error_squared,velocity_p1_error_squared,pressure_error_squared,gradient_error_squared\n";
    for(std::size_t t=0;t<mesh.cells.size();++t){const auto& cell=mesh.cells[t];P1Local local(mesh,t,fixture.diameter[t],errorOrder);const auto m=local.m;Vec ux(m),uy(m);
        for(std::size_t j=0;j<3;++j){ux[j]=state.cells[t][j];uy[j]=state.cells[t][3+j];}
        for(std::size_t l=0;l<cell.faces.size();++l)for(std::size_t j=0;j<2;++j){ux[3+2*l+j]=state.faces[cell.faces[l]][j];uy[3+2*l+j]=state.faces[cell.faces[l]][2+j];}
        P6 rx{},ry{};for(std::size_t i=0;i<6;++i)for(std::size_t j=0;j<m;++j){rx[i]+=local.potential(i,j)*ux[j];ry[i]+=local.potential(i,j)*uy[j];}
        double eu=0,eu1=0,ep=0,eg=0;
        for(const auto& q:local.q){const auto e=exact.at(q.p);const auto theta=local.basis.theta(q.p);const auto phi=local.basis.phi(q.p);const auto grad=local.basis.grad(q.p);Vector2D u{},u1{},gu{},gv{};double p=0;
            for(std::size_t j=0;j<6;++j){u.x+=theta[j]*rx[j];u.y+=theta[j]*ry[j];gu.x+=grad[j].x*rx[j];gu.y+=grad[j].y*rx[j];gv.x+=grad[j].x*ry[j];gv.y+=grad[j].y*ry[j];}
            for(std::size_t j=0;j<3;++j){u1.x+=phi[j]*state.cells[t][j];u1.y+=phi[j]*state.cells[t][3+j];p+=phi[j]*state.cells[t][6+j];}
            const auto du=u-e.u,du1=u1-e.u,dgu=gu-e.gradient[0],dgv=gv-e.gradient[1];const double dp=p-(e.p-gauge);
            eu+=q.w*dot(du,du);eu1+=q.w*dot(du1,du1);ep+=q.w*dp*dp;pMeanError+=q.w*dp;eg+=q.w*(dot(dgu,dgu)+dot(dgv,dgv));uNorm+=q.w*dot(e.u,e.u);pNorm+=q.w*(e.p-gauge)*(e.p-gauge);
        }
        for(auto id:cell.faces){const auto& face=mesh.faces[id];for(double s:{-.5,.5}){const Point2D p{face.centre.x-s*face.areaVector.y,face.centre.y+s*face.areaVector.x};
            const auto phi=local.basis.phi(p);double pressure=0;for(std::size_t j=0;j<3;++j)pressure+=phi[j]*state.cells[t][6+j];
            pMaximum=std::max(pMaximum,std::abs(pressure-exact.at(p).p+gauge));}}
        area+=cell.area;minArea=std::min(minArea,cell.area);maxDiameter=std::max(maxDiameter,fixture.diameter[t]);uError+=eu;uP1Error+=eu1;pError+=ep;gradientError+=eg;
        cells<<t<<','<<cell.centre.x<<','<<cell.centre.y<<','<<cell.area<<','<<fixture.diameter[t];for(double x:state.cells[t])cells<<','<<x;cells<<','<<eu<<','<<eu1<<','<<ep<<','<<eg<<'\n';
    }checkedClose(cells);
    std::array<WallError,2> walls;std::ofstream wallCsv(prefix+".walls.csv");wallCsv<<std::setprecision(17)<<"face,owner,embedded,x,y,length,nx,ny,p_cell0,p_cells,p_reaction0,p_reactions,p_exact0,p_exacts,tx0,txs,ty0,tys,exact_tx0,exact_txs,exact_ty0,exact_tys\n";
    for(const auto& load:loads.boundaries){const auto& face=mesh.faces[load.face];const double length=std::hypot(face.areaVector.x,face.areaVector.y);const auto n=face.areaVector*(1/length);const bool embedded=face.patch==BoundaryPatch2D::EmbeddedBoundary;auto& w=walls[embedded?0:1];++w.faces;w.length+=length;
        std::array<Vector2D,2> numerical{load.tractionMoments[0]*(1/length),load.tractionMoments[1]*(12/length)},projected{};
        const auto& pc=state.cells[face.owner];const Basis basis{mesh.cells[face.owner].centre,fixture.diameter[face.owner],{}};const auto phi=basis.phi(face.centre);
        const std::array<double,2> cellP{phi[0]*pc[6]+phi[1]*pc[7]+phi[2]*pc[8],(-face.areaVector.y*pc[7]+face.areaVector.x*pc[8])/basis.h};
        const std::array<double,2> reactionP{-dot(n,numerical[0]),-dot(n,numerical[1])},exactP{exact.at(face.centre).p-gauge,-2*face.areaVector.y-3*face.areaVector.x};
        for(auto [z,qw]:gauss(errorOrder)){const double s=z-.5;const Point2D p{face.centre.x-s*face.areaVector.y,face.centre.y+s*face.areaVector.x};const auto ex=exact.traction(p,n,gauge);projected[0]=projected[0]+ex*qw;projected[1]=projected[1]+ex*(12*s*qw);}
        w.force=w.force+load.tractionMoments[0];w.exactForce=w.exactForce+projected[0]*length;
        for(auto [z,qw]:gauss(errorOrder)){const double s=z-.5,weight=qw*length;const Point2D p{face.centre.x-s*face.areaVector.y,face.centre.y+s*face.areaVector.x};const auto ex=exact.traction(p,n,gauge),num=numerical[0]+numerical[1]*s,pr=projected[0]+projected[1]*s;
            const double dpc=cellP[0]+s*cellP[1]-exactP[0]-s*exactP[1],dpr=reactionP[0]+s*reactionP[1]-exactP[0]-s*exactP[1];const auto e=num-ex,em=num-pr,proj=pr-ex;
            w.pCell+=weight*dpc*dpc;w.pReaction+=weight*dpr*dpr;w.traction+=weight*dot(e,e);const double shear=-n.y*e.x+n.x*e.y;w.shear+=weight*shear*shear;w.tractionMoment+=weight*dot(em,em);w.projection+=weight*dot(proj,proj);
        }
        for(double s:{-.5,.5}){w.pCellMax=std::max(w.pCellMax,std::abs(cellP[0]+s*cellP[1]-exactP[0]-s*exactP[1]));w.pReactionMax=std::max(w.pReactionMax,std::abs(reactionP[0]+s*reactionP[1]-exactP[0]-s*exactP[1]));}
        wallCsv<<load.face<<','<<face.owner<<','<<embedded<<','<<face.centre.x<<','<<face.centre.y<<','<<length<<','<<n.x<<','<<n.y;
        for(auto a:{cellP,reactionP,exactP})for(double x:a)wallCsv<<','<<x;wallCsv<<','<<numerical[0].x<<','<<numerical[1].x<<','<<numerical[0].y<<','<<numerical[1].y<<','<<projected[0].x<<','<<projected[1].x<<','<<projected[0].y<<','<<projected[1].y<<'\n';
    }checkedClose(wallCsv);
    report<<",\"area\":"<<area<<",\"minimumCellArea\":"<<minArea<<",\"maximumCellDiameter\":"<<maxDiameter<<",\"pressureReferenceAtGaugeCell\":"<<gauge
          <<",\"velocityP2Rms\":"<<std::sqrt(uError/area)<<",\"velocityP1Rms\":"<<std::sqrt(uP1Error/area)<<",\"pressureP1Rms\":"<<std::sqrt(pError/area)
          <<",\"pressureMeanError\":"<<pMeanError/area<<",\"pressureP1MeanRemovedRms\":"<<std::sqrt(std::max(0.,pError/area-pMeanError*pMeanError/area/area))<<",\"pressureVertexMaximum\":"<<pMaximum
          <<",\"velocityGradientRms\":"<<std::sqrt(gradientError/area)<<",\"exactVelocityRms\":"<<std::sqrt(uNorm/area)<<",\"exactGaugedPressureRms\":"<<std::sqrt(pNorm/area)<<",\"embeddedWall\":";walls[0].write(report);report<<",\"outerWall\":";walls[1].write(report);
    report<<",\"momentumImbalance\":["<<loads.momentumImbalance.x<<','<<loads.momentumImbalance.y<<"],\"boundaryVolumeFlux\":"<<loads.boundaryVolumeFlux
          <<",\"totalSeconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";checkedClose(report);
    std::cout<<"cells="<<mesh.cells.size()<<" uP2="<<std::sqrt(uError/area)<<" p="<<std::sqrt(pError/area)<<" inner-wall-p="<<std::sqrt(walls[0].pReaction/walls[0].length)<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}

#endif
