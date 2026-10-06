// Research: nonconstant body-force integration with a local RT0 fan lift.
// Reuses the unchanged k=0 hybrid viscosity and B/-B^T pair. No product change.
// On each centroid fan triangle R(v)=a_j+b(x-c_T), 2b=D_T(v).
// Match every outer normal flux and every internal radial normal trace.
// The one remaining circulation minimises integral |R(v)-v_T|^2.
// This is a 2D k=0 construction, not the complete high-order method in 2203.07180.
#define CARTMESH_HYBRID_STOKES_NO_MAIN
#include "native-laminar-hybrid-stokes.cpp"
#include "cartmesh2d/io/MeshIO2D.hpp"

struct QuadPoint { Point2D p; double w; };
const double nodes[6]={.033765242898423986,.16939530676686776,.38069040695840156,.6193095930415985,.8306046932331322,.966234757101576};
const double weights[6]={.08566224618958517,.1803807865240693,.23395696728634552,.23395696728634552,.1803807865240693,.08566224618958517};
std::vector<QuadPoint> triangleQuad(Point2D c,Point2D a,Point2D b) {
    const double jac=(a.x-c.x)*(b.y-c.y)-(a.y-c.y)*(b.x-c.x);
    if(!(jac>0))throw std::runtime_error("centroid fan is not positive: unsupported research polygon");
    std::vector<QuadPoint> q;
    for(int i=0;i<6;++i)for(int j=0;j<6;++j){double s=nodes[i],t=nodes[j];
        q.push_back({{c.x+s*(a.x-c.x)+(1-s)*t*(b.x-c.x),c.y+s*(a.y-c.y)+(1-s)*t*(b.y-c.y)},weights[i]*weights[j]*(1-s)*jac});}
    return q;
}
struct Fan {
    std::vector<std::size_t> ids;
    std::vector<Vector2D> s,rad;
    std::vector<Point2D> a,b,centroid;
    std::vector<std::vector<QuadPoint>> quad;
    Vec area;
    std::vector<Vec> lift; // one vector [a0x,a0y,...,b] per local velocity basis
    double continuity=0,divergence=0,outer=0;
};
Fan buildFan(const FvMesh2D& mesh,std::size_t ti) {
    const auto& cell=mesh.cells[ti];Fan f;f.ids=cell.faces;
    std::sort(f.ids.begin(),f.ids.end(),[&](auto a,auto b){const auto x=mesh.faces[a].centre-cell.centre,y=mesh.faces[b].centre-cell.centre;return std::atan2(x.y,x.x)<std::atan2(y.y,y.x);});
    const int m=int(f.ids.size()),n=4*m-1;
    Vec k(n*n);
    for(int i=0;i<m;++i){const auto& face=mesh.faces[f.ids[i]];const double sign=face.owner==ti?1:-1;
        Vector2D s{sign*face.areaVector.x,sign*face.areaVector.y};f.s.push_back(s);
        Point2D a{face.centre.x+.5*s.y,face.centre.y-.5*s.x},b{face.centre.x-.5*s.y,face.centre.y+.5*s.x};
        f.a.push_back(a);f.b.push_back(b);f.rad.push_back({a.y-cell.centre.y,-(a.x-cell.centre.x)});
        f.centroid.push_back({(cell.centre.x+a.x+b.x)/3,(cell.centre.y+a.y+b.y)/3});
        auto q=triangleQuad(cell.centre,a,b);double area=0;for(auto z:q)area+=z.w;f.area.push_back(area);f.quad.push_back(q);
        for(int d=0;d<2;++d)k[(2*i+d)*n+2*i+d]=area;
        for(int d=0;d<2;++d){double v=d?s.y:s.x;k[(2*m+i)*n+2*i+d]=v;k[(2*i+d)*n+2*m+i]=v;}
    }
    for(int i=1;i<m;++i)for(int d=0;d<2;++d){const double v=d?f.rad[i].y:f.rad[i].x;const int r=3*m+i-1;
        k[r*n+2*i+d]=v;k[(2*i+d)*n+r]=v;k[r*n+2*(i-1)+d]=-v;k[(2*(i-1)+d)*n+r]=-v;}
    for(int i=0;i<m;++i){const auto d=f.a[i]-f.b[(i+m-1)%m];
        const double scale=1+std::hypot(cell.centre.x,cell.centre.y)+std::sqrt(cell.area);
        if(std::hypot(d.x,d.y)>256*std::numeric_limits<double>::epsilon()*scale)
            throw std::runtime_error("fan endpoints do not join at construction roundoff; no welding performed");
    }
    const auto lu=DenseLU(k,n);
    for(int j=0;j<2*(m+1);++j){Vec rhs(n);double div=0;
        if(j>=2)div=(j%2?f.s[(j-2)/2].y:f.s[(j-2)/2].x)/cell.area;
        const double b=.5*div;
        for(int i=0;i<m;++i){const auto r=f.centroid[i]-cell.centre;
            rhs[2*i]=f.area[i]*((j==0?1:0)-b*r.x);rhs[2*i+1]=f.area[i]*((j==1?1:0)-b*r.y);
            const auto rface=mesh.faces[f.ids[i]].centre-cell.centre;
            rhs[2*m+i]=(j>=2&&(j-2)/2==i?(j%2?f.s[i].y:f.s[i].x):0)-b*dot(f.s[i],rface);
        }
        auto a=lu.solve(rhs);a.resize(2*m);a.push_back(b);f.lift.push_back(a);
        for(int i=0;i<m;++i){int prev=(i+m-1)%m;
            f.continuity=std::max(f.continuity,std::abs(f.rad[i].x*(a[2*i]-a[2*prev])+f.rad[i].y*(a[2*i+1]-a[2*prev+1])));
            auto r=mesh.faces[f.ids[i]].centre-cell.centre;
            const double target=j>=2&&(j-2)/2==i?(j%2?f.s[i].y:f.s[i].x):0;
            f.outer=std::max(f.outer,std::abs(f.s[i].x*(a[2*i]+b*r.x)+f.s[i].y*(a[2*i+1]+b*r.y)-target));
            f.divergence=std::max(f.divergence,std::abs(2*b-div));
        }
    }
    return f;
}
double potential(Point2D p){return p.x*p.x*p.x+p.x*p.y*p.y;}
Vector2D gradPotential(Point2D p){return {3*p.x*p.x+p.y*p.y,2*p.x*p.y};}
Vector2D exactVelocity(Point2D p,const std::string& problem){
    if(problem=="hydrostatic")return {};
    if(problem=="polynomial")return {p.x*p.x,-2*p.x*p.y};
    const double pi=std::acos(-1.);return {std::sin(pi*p.x)*std::cos(pi*p.y),-std::cos(pi*p.x)*std::sin(pi*p.y)};
}
Vector2D force(Point2D p,const std::string& problem,double lambda){
    auto f=gradPotential(p);f.x*=lambda;f.y*=lambda;
    if(problem=="polynomial")f.x-=2.;
    if(problem=="vortex"){auto u=exactVelocity(p,problem);double k=2*std::acos(-1.)*std::acos(-1.);f.x+=k*u.x;f.y+=k*u.y;}
    return f;
}
struct LoadResult {Vec u,p;double urms=0,prms=0,div=0,residual=0,pressureWork=0,seconds=0,trace=0;};
LoadResult solveLoad(const Fixture& fixture,const std::string& problem,double lambda,bool lifted,const std::string& prefix){
    const auto start=std::chrono::steady_clock::now();const auto& mesh=fixture.mesh;int nc=int(mesh.cells.size()),nf=int(mesh.faces.size());
    std::vector<int> map(2*(nc+nf),-1);int nv=2*nc;std::iota(map.begin(),map.begin()+nv,0);
    for(int i=0;i<nf;++i)if(mesh.faces[i].neighbour){map[2*(nc+i)]=nv++;map[2*(nc+i)+1]=nv++;}
    int n=nv+nc-1;if(n>1500)throw std::runtime_error("dense prototype limit 1500");
    Vec mat(n*n),rhs(n),known(map.size()),up(2*nc),pp(nc);LoadResult result;
    for(int i=0;i<nf;++i)if(!mesh.faces[i].neighbour){const auto& f=mesh.faces[i];
        for(int j=0;j<6;++j){Point2D p{f.centre.x-(nodes[j]-.5)*f.areaVector.y,f.centre.y+(nodes[j]-.5)*f.areaVector.x};auto u=exactVelocity(p,problem);known[2*(nc+i)]+=weights[j]*u.x;known[2*(nc+i)+1]+=weights[j]*u.y;}}
    const auto add=[&](int i,int j,double a){if(map[i]<0)return;if(map[j]>=0)mat[map[i]*n+map[j]]+=a;else rhs[map[i]]-=a*known[j];};
    for(int t=0;t<nc;++t){const auto& cell=mesh.cells[t];const auto fan=buildFan(mesh,t);int m=int(fan.ids.size())+1;std::vector<int> ids(m);ids[0]=t;Vec gx(m),gy(m);
        result.trace=std::max({result.trace,fan.continuity,fan.outer});
        for(int i=1;i<m;++i){ids[i]=nc+int(fan.ids[i-1]);gx[i]=fan.s[i-1].x/cell.area;gy[i]=fan.s[i-1].y/cell.area;}
        for(int d=0;d<2;++d)for(int i=0;i<m;++i)for(int j=0;j<m;++j)add(2*ids[i]+d,2*ids[j]+d,cell.area*(gx[i]*gx[j]+gy[i]*gy[j]));
        for(int i=1;i<m;++i){const auto& face=mesh.faces[fan.ids[i-1]];auto r=face.centre-cell.centre;Vec v(m);v[0]=-1;v[i]=1;for(int j=1;j<m;++j)v[j]-=gx[j]*r.x+gy[j]*r.y;
            double w=std::hypot(face.areaVector.x,face.areaVector.y)/fixture.diameter[t];
            for(int d=0;d<2;++d)for(int j=0;j<m;++j)for(int l=0;l<m;++l)add(2*ids[j]+d,2*ids[l]+d,w*v[j]*v[l]);
            if(t<nc-1)for(int d=0;d<2;++d){double b=cell.area*(d?gy[i]:gx[i]);int row=map[2*ids[i]+d],p=nv+t;
                if(row>=0){mat[row*n+p]-=b;mat[p*n+row]-=b;}else rhs[p]+=b*known[2*ids[i]+d];}
        }
        for(int tri=0;tri<m-1;++tri)for(const auto& q:fan.quad[tri]){
            const auto f=force(q.p,problem,lambda),ex=exactVelocity(q.p,problem);auto r=q.p-cell.centre;
            up[2*t]+=q.w*ex.x/cell.area;up[2*t+1]+=q.w*ex.y/cell.area;pp[t]+=q.w*lambda*potential(q.p)/cell.area;
            if(lifted){for(int j=0;j<2*m;++j){int row=map[2*ids[j/2]+j%2];if(row<0)continue;const auto& a=fan.lift[j];double b=a.back();rhs[row]+=q.w*(f.x*(a[2*tri]+b*r.x)+f.y*(a[2*tri+1]+b*r.y));}}
            else {rhs[map[2*t]]+=q.w*f.x;rhs[map[2*t+1]]+=q.w*f.y;}
        }
    }
    const auto lu=DenseLU(mat,n);Vec x=lu.solve(rhs);
    for(int pass=0;pass<2;++pass){Vec r=rhs;for(int i=0;i<n;++i)for(int j=0;j<n;++j)r[i]-=mat[i*n+j]*x[j];auto dx=lu.solve(r);for(int i=0;i<n;++i)x[i]+=dx[i];}
    for(double v:x)if(!std::isfinite(v))throw std::runtime_error("nonfinite Stokes candidate");
    result.u=known;result.p.resize(nc);for(std::size_t i=0;i<map.size();++i)if(map[i]>=0)result.u[i]=x[map[i]];
    double area=0,gauge=0;for(int t=0;t<nc;++t){result.p[t]=t<nc-1?x[nv+t]:0;area+=mesh.cells[t].area;gauge+=mesh.cells[t].area*(result.p[t]-pp[t]);}gauge/=area;
    for(int t=0;t<nc;++t){const auto& cell=mesh.cells[t];double du=std::hypot(result.u[2*t]-up[2*t],result.u[2*t+1]-up[2*t+1]),dp=result.p[t]-pp[t]-gauge;result.urms+=cell.area*du*du;result.prms+=cell.area*dp*dp;
        double flux=0;for(auto fi:cell.faces){const auto& f=mesh.faces[fi];double sign=f.owner==std::size_t(t)?1:-1;flux+=sign*(f.areaVector.x*result.u[2*(nc+fi)]+f.areaVector.y*result.u[2*(nc+fi)+1]);}
        result.div=std::max(result.div,std::abs(flux)/cell.area);result.pressureWork-=result.p[t]*flux;
    }
    result.urms=std::sqrt(result.urms/area);result.prms=std::sqrt(result.prms/area);
    for(int i=0;i<n;++i){double r=-rhs[i];for(int j=0;j<n;++j)r+=mat[i*n+j]*x[j];result.residual=std::max(result.residual,std::abs(r));}
    std::ofstream cells(prefix+".cells.csv"),faces(prefix+".faces.csv");cells<<std::setprecision(17)<<"cell,area,u,v,p,u_exact,v_exact,p_exact_plus_gauge\n";faces<<std::setprecision(17)<<"face,u,v,flux\n";
    for(int i=0;i<nc;++i)cells<<i<<','<<mesh.cells[i].area<<','<<result.u[2*i]<<','<<result.u[2*i+1]<<','<<result.p[i]<<','<<up[2*i]<<','<<up[2*i+1]<<','<<pp[i]+gauge<<'\n';
    for(int i=0;i<nf;++i)faces<<i<<','<<result.u[2*(nc+i)]<<','<<result.u[2*(nc+i)+1]<<','<<mesh.faces[i].areaVector.x*result.u[2*(nc+i)]+mesh.faces[i].areaVector.y*result.u[2*(nc+i)+1]<<'\n';
    cells.flush();faces.flush();if(!cells||!faces)throw std::runtime_error("field write failure");result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return result;
}
int auditMesh(const std::string& path){
    const auto start=std::chrono::steady_clock::now();auto read=readCm2dTopology(path);
    if(!read.valid())throw std::runtime_error(read.error);
    auto mesh=makeFvMesh2D(read.topology);
    double trace=0,identity=0,area=0;std::size_t triangles=0;
    for(std::size_t t=0;t<mesh.cells.size();++t){auto fan=buildFan(mesh,t);const auto& cell=mesh.cells[t];
        trace=std::max({trace,fan.continuity,fan.outer});area+=cell.area;triangles+=fan.ids.size();
        double integral=0;for(const auto& qs:fan.quad)for(auto q:qs)integral+=q.w*potential(q.p);
        for(std::size_t j=0;j<fan.lift.size();++j){const auto& a=fan.lift[j];double lhs=0,rhs=-2*a.back()*integral;
            for(std::size_t i=0;i<fan.ids.size();++i){for(auto q:fan.quad[i]){auto r=q.p-cell.centre,f=gradPotential(q.p);lhs+=q.w*(f.x*(a[2*i]+a.back()*r.x)+f.y*(a[2*i+1]+a.back()*r.y));}
                if(j>=2&&(j-2)/2==i){double avg=0;const auto& face=mesh.faces[fan.ids[i]];
                    for(int k=0;k<6;++k){Point2D p{face.centre.x-(nodes[k]-.5)*fan.s[i].y,face.centre.y+(nodes[k]-.5)*fan.s[i].x};avg+=weights[k]*potential(p);}
                    rhs+=avg*(j%2?fan.s[i].y:fan.s[i].x);}
            }
            identity=std::max(identity,std::abs(lhs-rhs));
        }
    }
    std::cout<<std::setprecision(17)<<"{\"mesh\":\""<<path<<"\",\"cells\":"<<mesh.cells.size()<<",\"faces\":"<<mesh.faces.size()<<",\"fanTriangles\":"<<triangles<<",\"area\":"<<area<<",\"maxNormalTraceIntegralResidual\":"<<trace<<",\"maxGradientLoadIdentityResidual\":"<<identity<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    return 0;
}
int main(int argc,char** argv)try{
    if(argc==3&&std::string(argv[1])=="--audit")return auditMesh(argv[2]);
    if(argc!=2||!std::filesystem::is_directory(argv[1]))throw std::runtime_error("usage: hdiv-load existing-output-directory");
    for(const std::string name:{"square","sheared","cut"})for(int n:{4,8,12})try{
        const auto f=grid(n,name=="sheared"?.7:0,name=="cut");
        for(const std::string problem:{"hydrostatic","polynomial","vortex"})for(bool lifted:{false,true}){LoadResult base;
            for(double lambda:{0.,1.,10000.}){const std::string label=name+std::to_string(n)+"-"+problem+(lifted?"-lift":"-cell")+"-"+std::to_string(int(lambda));
                const auto r=solveLoad(f,problem,lambda,lifted,std::string(argv[1])+"/"+label);if(lambda==0)base=r;double du=0,area=0;
                for(std::size_t i=0;i<f.mesh.cells.size();++i){double a=f.mesh.cells[i].area;area+=a;du+=a*(std::pow(r.u[2*i]-base.u[2*i],2)+std::pow(r.u[2*i+1]-base.u[2*i+1],2));}
                std::cout<<std::setprecision(17)<<"{\"label\":\""<<label<<"\",\"cells\":"<<f.mesh.cells.size()<<",\"velocityRms\":"<<r.urms<<",\"pressureRms\":"<<r.prms<<",\"addedGradientVelocityRms\":"<<std::sqrt(du/area)<<",\"maxDivergence\":"<<r.div<<",\"maxTraceResidual\":"<<r.trace<<",\"pressureWork\":"<<r.pressureWork<<",\"linearResidual\":"<<r.residual<<",\"seconds\":"<<r.seconds<<"}\n";
            }
        }
    }catch(const std::exception& e){if(!std::string(e.what()).starts_with("FVM mesh: existing Solver quality gate failed"))throw;std::cout<<"{\"label\":\""<<name<<n<<"\",\"rejected\":\""<<e.what()<<"\"}\n";}
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
