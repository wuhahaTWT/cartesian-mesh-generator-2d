// Research-only lowest-order cell/face hybrid Stokes prototype.
// Uses the native polygon topology and unchanged finite-volume geometry gate.
// This is a Laplacian Stokes discretisation, not a product solver or a claim
// of reproducing the full pressure-robust HHO Navier--Stokes method.
// G_T = sum_F u_F tensor S_TF / |T|, r_T(x)=u_T+G_T(x-c_T).
// a_T = nu |T| G:G + nu sum_F |F|/h_T |u_F-r_T(c_F)|^2.
// B_T u = sum_F S_TF dot u_F; pressure force is exactly -B^T p.
// h_T is the polygon diameter. All boundary velocities are imposed as exact
// face averages. The pressure gauge removes one redundant divergence row.
// For constant cell body forces, optional first-moment load reconstruction
// integrates f dot R(v) = f dot sum_F (c_F-c_T)(v_F dot S_TF).
// This suffices for a global affine pressure force, not general forcing.

#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#define main inherited_native_face_test_main
#include "../../tests/flow_face_test.cpp"
#undef main

using Vec=std::vector<double>;
using namespace cartmesh2d;
using namespace cartmesh2d::fv;

struct Fixture { FvMesh2D mesh; Vec meanY2, diameter; };

Fixture makeFixture(const std::vector<Polygon2D>& polys) {
    Fixture result{makeFvMesh2D(fromPolygons(polys)), {}, {}};
    for (const auto& p:polys) {
        double moment=0.,diameter=0.;
        for(std::size_t j=0;j<p.vertices.size();++j) {
            const auto a=p.vertices[j],b=p.vertices[(j+1)%p.vertices.size()];
            moment+=(a.x*b.y-b.x*a.y)*(a.y*a.y+a.y*b.y+b.y*b.y)/12.;
            for(const auto q:p.vertices) diameter=std::max(diameter,std::hypot(a.x-q.x,a.y-q.y));
        }
        result.meanY2.push_back(moment/p.area());result.diameter.push_back(diameter);
    }
    return result;
}

// Intersections are computed once in canonical endpoint order so adjacent
// polygons retain identical vertex identity, including very short real edges.
Polygon2D clip(const Polygon2D& p,double offset) {
    const auto distance=[offset](Point2D q){return q.y-.31*q.x-offset;};
    Polygon2D r;
    for(std::size_t j=0;j<p.vertices.size();++j) {
        auto a=p.vertices[j],b=p.vertices[(j+1)%p.vertices.size()];
        const double da=distance(a),db=distance(b);
        if(da>=0) r.vertices.push_back(a);
        if((da<0&&db>0)||(da>0&&db<0)) {
            if(std::make_pair(b.x,b.y)<std::make_pair(a.x,a.y)) std::swap(a,b);
            const double t=distance(a)/(distance(a)-distance(b));
            r.vertices.push_back({a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)});
        }
    }
    return r;
}

Fixture grid(int n,double shear,bool cut) {
    std::vector<Polygon2D> polys;
    for(int j=0;j<n;++j) for(int i=0;i<n;++i) {
        const auto pt=[n,shear](int x,int y){return Point2D{(x+shear*y)/n,double(y)/n};};
        Polygon2D p{{pt(i,j),pt(i+1,j),pt(i+1,j+1),pt(i,j+1)}};
        if(cut) p=clip(p,.127);
        if(p.vertices.size()>=3 && p.area()>0) polys.push_back(p);
    }
    return makeFixture(polys);
}

struct DenseLU {
    int n;Vec a;std::vector<int> pivots;
    explicit DenseLU(Vec matrix,int size):n(size),a(std::move(matrix)),pivots(n) {
        for(int k=0;k<n;++k) {
            int p=k;for(int i=k+1;i<n;++i) if(std::abs(a[i*n+k])>std::abs(a[p*n+k]))p=i;
            if(!std::isfinite(a[p*n+k]) || a[p*n+k]==0.)
                throw std::runtime_error("singular research matrix; no solution accepted");
            pivots[k]=p;if(p!=k)for(int j=0;j<n;++j)std::swap(a[k*n+j],a[p*n+j]);
            for(int i=k+1;i<n;++i) {
                const double t=(a[i*n+k]/=a[k*n+k]);
                for(int j=k+1;j<n;++j)a[i*n+j]-=t*a[k*n+j];
            }
        }
    }
    Vec solve(Vec x)const {
        for(int k=0;k<n;++k)if(pivots[k]!=k)std::swap(x[k],x[pivots[k]]);
        for(int i=0;i<n;++i)for(int j=0;j<i;++j)x[i]-=a[i*n+j]*x[j];
        for(int i=n-1;i>=0;--i){for(int j=i+1;j<n;++j)x[i]-=a[i*n+j]*x[j];x[i]/=a[i*n+i];}
        return x;
    }
};

struct Result {
    int cells=0,faces=0,unknowns=0;
    double velocityRms=0,velocityMax=0,pressureRms=0,pressureMax=0;
    double massMax=0,netFlux=0,pressureWork=0,residualMax=0;
    double minimumArea=1e300,area=0,maxAffineMomentDefect=0,seconds=0;
};

Result solveStokes(const Fixture& fixture,const std::string& problem,bool momentLoad, const std::string& outputPrefix = {}) {
    const auto started=std::chrono::steady_clock::now();
    const auto& mesh=fixture.mesh;const int nc=int(mesh.cells.size()),nf=int(mesh.faces.size());
    std::vector<int> map(2*(nc+nf),-1);int nv=2*nc;
    std::iota(map.begin(),map.begin()+nv,0);
    for(int f=0;f<nf;++f)if(mesh.faces[f].neighbour){map[2*(nc+f)]=nv++;map[2*(nc+f)+1]=nv++;}
    const int n=nv+nc-1;
    if(n>1500)throw std::runtime_error("dense research driver limited to 1500 unknowns");
    Vec known(map.size()),matrix(n*n),rhs(n),body(2*nc),exact(2*nc),exactP(nc);
    const double nu=1.;
    const auto velocity=[&](Point2D c,double yy){
        if(problem=="couette")return Vector2D{c.y,0};
        if(problem=="rotation")return Vector2D{-c.y,c.x};
        if(problem=="poiseuille")return Vector2D{4*(c.y-yy),0};
        return Vector2D{};
    };
    for(int t=0;t<nc;++t) {
        const auto c=mesh.cells[t].centre;const auto u=velocity(c,fixture.meanY2[t]);
        exact[2*t]=u.x;exact[2*t+1]=u.y;
        if(problem=="hydrostatic"){exactP[t]=c.x+2*c.y;body[2*t]=1;body[2*t+1]=2;}
        if(problem=="poiseuille")exactP[t]=-8*nu*c.x;
    }
    for(int f=0;f<nf;++f)if(!mesh.faces[f].neighbour) {
        const auto& face=mesh.faces[f];
        const double yy=face.centre.y*face.centre.y+face.areaVector.x*face.areaVector.x/12.;
        const auto u=velocity(face.centre,yy);known[2*(nc+f)]=u.x;known[2*(nc+f)+1]=u.y;
    }
    const auto add=[&](int i,int j,double value){
        const int row=map[i],col=map[j];if(row<0)return;
        if(col>=0)matrix[row*n+col]+=value;else rhs[row]-=value*known[j];
    };
    for(int t=0;t<nc;++t) {
        const auto& cell=mesh.cells[t];const int m=int(cell.faces.size())+1;
        std::vector<int> ids(m);ids[0]=t;
        Vec gx(m),gy(m);
        for(int f=1;f<m;++f) {
            const auto id=cell.faces[f-1];ids[f]=nc+int(id);const auto& face=mesh.faces[id];
            const double sign=face.owner==std::size_t(t)?1:-1;
            gx[f]=sign*face.areaVector.x/cell.area;gy[f]=sign*face.areaVector.y/cell.area;
        }
        for(int c=0;c<2;++c)for(int i=0;i<m;++i)for(int j=0;j<m;++j)
            add(2*ids[i]+c,2*ids[j]+c,nu*cell.area*(gx[i]*gx[j]+gy[i]*gy[j]));
        for(int f=1;f<m;++f) {
            const auto& face=mesh.faces[cell.faces[f-1]];
            const double dx=face.centre.x-cell.centre.x,dy=face.centre.y-cell.centre.y;
            Vec r(m);r[0]=-1;r[f]=1;
            for(int j=1;j<m;++j)r[j]-=gx[j]*dx+gy[j]*dy;
            const double weight=nu*std::hypot(face.areaVector.x,face.areaVector.y)/fixture.diameter[t];
            for(int c=0;c<2;++c)for(int i=0;i<m;++i)for(int j=0;j<m;++j)
                add(2*ids[i]+c,2*ids[j]+c,weight*r[i]*r[j]);
            if(momentLoad)for(int c=0;c<2;++c) {
                const int row=map[2*ids[f]+c];if(row<0)continue;
                rhs[row]+=(body[2*t]*dx+body[2*t+1]*dy)*cell.area*(c==0?gx[f]:gy[f]);
            }
            if(t<nc-1)for(int c=0;c<2;++c) {
                const double b=cell.area*(c==0?gx[f]:gy[f]);const int row=map[2*ids[f]+c],p=nv+t;
                if(row>=0){matrix[row*n+p]-=b;matrix[p*n+row]-=b;}
                else rhs[p]+=b*known[2*ids[f]+c];
            }
        }
        if(!momentLoad)for(int c=0;c<2;++c)rhs[map[2*t+c]]+=cell.area*body[2*t+c];
    }
    const Vec x=DenseLU(matrix,n).solve(rhs);Vec u=known;
    for(const auto value:x)if(!std::isfinite(value))throw std::runtime_error("nonfinite Stokes candidate");
    for(std::size_t i=0;i<map.size();++i)if(map[i]>=0)u[i]=x[map[i]];
    Result result;result.cells=nc;result.faces=nf;result.unknowns=n;
    double gauge=0.;
    for(int t=0;t<nc;++t){result.area+=mesh.cells[t].area;gauge+=mesh.cells[t].area*((t<nc-1?x[nv+t]:0)-exactP[t]);}
    gauge/=result.area;
    for(int t=0;t<nc;++t) {
        const auto& cell=mesh.cells[t];const double err=std::hypot(u[2*t]-exact[2*t],u[2*t+1]-exact[2*t+1]);
        result.velocityRms+=cell.area*err*err;result.velocityMax=std::max(result.velocityMax,err);
        const double p=t<nc-1?x[nv+t]:0,pe=p-exactP[t]-gauge;
        result.pressureRms+=cell.area*pe*pe;result.pressureMax=std::max(result.pressureMax,std::abs(pe));
        double flux=0,xx=0,xy=0,yx=0,yy=0;
        for(const auto f:cell.faces) {
            const auto& face=mesh.faces[f];const double sign=face.owner==std::size_t(t)?1:-1;
            const double sx=sign*face.areaVector.x,sy=sign*face.areaVector.y;
            flux+=sx*u[2*(nc+f)]+sy*u[2*(nc+f)+1];
            const double dx=face.centre.x-cell.centre.x,dy=face.centre.y-cell.centre.y;
            xx+=dx*sx/cell.area;xy+=dx*sy/cell.area;yx+=dy*sx/cell.area;yy+=dy*sy/cell.area;
        }
        result.maxAffineMomentDefect=std::max({result.maxAffineMomentDefect,std::abs(xx-1),std::abs(xy),std::abs(yx),std::abs(yy-1)});
        result.massMax=std::max(result.massMax,std::abs(flux)/cell.area);result.netFlux+=flux;
        result.pressureWork-=p*flux;result.minimumArea=std::min(result.minimumArea,cell.area);
    }
    result.velocityRms=std::sqrt(result.velocityRms/result.area);result.pressureRms=std::sqrt(result.pressureRms/result.area);
    for(int i=0;i<n;++i){double r=-rhs[i];for(int j=0;j<n;++j)r+=matrix[i*n+j]*x[j];result.residualMax=std::max(result.residualMax,std::abs(r));}
    if(!outputPrefix.empty()) {
        std::ofstream field(outputPrefix+".cells.csv"),faces(outputPrefix+".faces.csv");
        field<<std::setprecision(17)<<"cell,x,y,area,u_mean,v_mean,p,p_exact_plus_gauge,u_exact_mean,v_exact_mean\n";
        for(int t=0;t<nc;++t)field<<t<<','<<mesh.cells[t].centre.x<<','<<mesh.cells[t].centre.y<<','<<mesh.cells[t].area<<','
            <<u[2*t]<<','<<u[2*t+1]<<','<<(t<nc-1?x[nv+t]:0)<<','<<exactP[t]+gauge<<','<<exact[2*t]<<','<<exact[2*t+1]<<'\n';
        faces<<std::setprecision(17)<<"face,owner,neighbour,x,y,Sx,Sy,u_mean,v_mean,owner_flux\n";
        for(int f=0;f<nf;++f){const auto& q=mesh.faces[f];faces<<f<<','<<q.owner<<','<<(q.neighbour?int(*q.neighbour):-1)<<','
            <<q.centre.x<<','<<q.centre.y<<','<<q.areaVector.x<<','<<q.areaVector.y<<','<<u[2*(nc+f)]<<','<<u[2*(nc+f)+1]<<','
            <<q.areaVector.x*u[2*(nc+f)]+q.areaVector.y*u[2*(nc+f)+1]<<'\n';}
        field.flush();faces.flush();if(!field||!faces)throw std::runtime_error("research field write failed");
    }
    result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    return result;
}

void print(const std::string& mesh,const std::string& problem,bool moment,const Result& r) {
    std::cout<<std::setprecision(17)<<"{\"mesh\":\""<<mesh<<"\",\"problem\":\""<<problem<<"\",\"moment_load\":"<<(moment?"true":"false")
    <<",\"cells\":"<<r.cells<<",\"faces\":"<<r.faces<<",\"unknowns\":"<<r.unknowns<<",\"velocity_rms\":"<<r.velocityRms<<",\"velocity_max\":"<<r.velocityMax
    <<",\"pressure_rms\":"<<r.pressureRms<<",\"pressure_max\":"<<r.pressureMax<<",\"mass_max\":"<<r.massMax<<",\"net_flux\":"<<r.netFlux
    <<",\"pressure_work\":"<<r.pressureWork<<",\"linear_residual_max\":"<<r.residualMax<<",\"minimum_area\":"<<r.minimumArea<<",\"area\":"<<r.area
    <<",\"affine_moment_defect\":"<<r.maxAffineMomentDefect<<",\"seconds\":"<<r.seconds<<"}\n";
}

#ifndef CARTMESH_HYBRID_STOKES_NO_MAIN
int main(int argc,char** argv)try {
    if(argc>2)throw std::runtime_error("usage: probe [existing-output-directory]");
    const std::string output=argc==2?argv[1]:"";
    if(!output.empty()&&!std::filesystem::is_directory(output))throw std::runtime_error("output directory does not exist");
    const auto prefix=[&](const std::string& name,const std::string& problem,bool moment){
        return output.empty()?std::string{}:output+"/"+name+"-"+problem+(moment?"-moment":"-cell");
    };
    for(const auto& meshName:{std::string("square"),std::string("sheared"),std::string("cut")})for(int n:{4,8,12}) try {
        const auto fixture=grid(n,meshName=="sheared"?.7:0,meshName=="cut");
        for(const std::string problem:{"couette","rotation","hydrostatic","poiseuille"}) {
            print(meshName+std::to_string(n),problem,true,solveStokes(fixture,problem,true,prefix(meshName+std::to_string(n),problem,true)));
            if(problem=="hydrostatic")print(meshName+std::to_string(n),problem,false,solveStokes(fixture,problem,false,prefix(meshName+std::to_string(n),problem,false)));
        }
    } catch(const std::exception& e) {
        if(!std::string(e.what()).starts_with("FVM mesh: existing Solver quality gate failed"))throw;
        std::cout<<"{\"mesh\":\""<<meshName<<n<<"\",\"status\":\"rejected\",\"reason\":\""<<e.what()<<"\"}\n";
    }
    const std::vector<std::pair<std::string,std::vector<Polygon2D>>> special{
        {"split_face",{{{{0,0},{1,0},{1,1},{1,2},{0,2}}},{{{1,0},{2,0},{2,1},{1,1}}},{{{1,1},{2,1},{2,2},{1,2}}}}},
        {"triangular_tip",{{{{0,0},{1,0},{0,1}}},{{{1,0},{1,1},{0,1}}},{{{1,0},{2,0},{2,1},{1,1}}},{{{0,1},{1,1},{1,2},{0,2}}}}}
    };
    for(const auto& [name,polys]:special) {
        const auto fixture=makeFixture(polys);
        for(const std::string problem:{"couette","rotation","hydrostatic","poiseuille"})
            print(name,problem,true,solveStokes(fixture,problem,true,prefix(name,problem,true)));
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}

#endif
