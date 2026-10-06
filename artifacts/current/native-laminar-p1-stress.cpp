// Research only: P1 hybrid symmetric-stress Stokes, RT1 fan load, static condensation.
// No collocated-product/default/geometry/quality-gate modifications; no convection.
#define CARTMESH_HYBRID_STOKES_P1_NO_MAIN
#include "native-laminar-hybrid-stokes-p1.cpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <cstdint>
#include <limits>
#include <sstream>

struct Tri {
    Point2D c,a,b; double h,det; Vector2D j0,j1;
    Tri(Point2D C,Point2D A,Point2D B,double H):c(C),a(A),b(B),h(H),j0{(A.x-C.x)/H,(A.y-C.y)/H},j1{(B.x-C.x)/H,(B.y-C.y)/H} {
        det=j0.x*j1.y-j0.y*j1.x;
        if(!(det>0))throw std::runtime_error("unsupported nonpositive centroid fan; geometry untouched");
    }
    Vector2D ref(Point2D p)const {auto d=p-c;return {(d.x*j1.y-d.y*j1.x)/(h*det),(j0.x*d.y-j0.y*d.x)/(h*det)};}
    std::array<Vector2D,8> shape(Point2D p)const {
        auto z=ref(p);std::array<Vector2D,8> a{{{1,0},{z.x,0},{z.y,0},{0,1},{0,z.x},{0,z.y},{z.x*z.x,z.x*z.y},{z.x*z.y,z.y*z.y}}};
        for(auto& v:a)v={(j0.x*v.x+j1.x*v.y)/det,(j0.y*v.x+j1.y*v.y)/det};
        return a;
    }
    std::vector<Q> quadrature(int order)const {
        std::vector<Q> out;auto g=gauss(order);
        for(auto [r,w]:g)for(auto [s,z]:g)out.push_back({{c.x+r*(a.x-c.x)+s*(1-r)*(b.x-c.x),c.y+r*(a.y-c.y)+s*(1-r)*(b.y-c.y)},w*z*(1-r)*h*h*det});
        return out;
    }
};
struct Lift {
    std::vector<Tri> tri;std::vector<int> faceLocal;Mat coefficients;
    double traceResidual=0,divResidual=0;
    Lift(const Fixture& f,int t,const Local& a,int order):coefficients(8*int(f.mesh.cells[t].faces.size()),2*a.m) {
        const auto& mesh=f.mesh;const auto& cell=mesh.cells[t];int nf=int(cell.faces.size()),nv=8*nf,ncon=7*nf-1,n=nv+ncon;
        std::vector<int> perm(nf);std::iota(perm.begin(),perm.end(),0);
        std::sort(perm.begin(),perm.end(),[&](int i,int j){auto x=mesh.faces[cell.faces[i]].centre-cell.centre,y=mesh.faces[cell.faces[j]].centre-cell.centre;return std::atan2(x.y,x.x)<std::atan2(y.y,y.x);});
        for(int i:perm){const auto& face=mesh.faces[cell.faces[i]];double s=face.owner==std::size_t(t)?1:-1;auto S=face.areaVector;
            tri.emplace_back(cell.centre,Point2D{face.centre.x+s*S.y/2,face.centre.y-s*S.x/2},Point2D{face.centre.x-s*S.y/2,face.centre.y+s*S.x/2},a.basis.h);faceLocal.push_back(i);}
        for(int i=0;i<nf;++i){auto d=tri[i].a-tri[(i+nf-1)%nf].b;double scale=1+std::hypot(cell.centre.x,cell.centre.y)+std::sqrt(cell.area);
            if(std::hypot(d.x,d.y)>256*std::numeric_limits<double>::epsilon()*scale)throw std::runtime_error("fan endpoints not closed at construction roundoff");}
        Mat k(n,n),r(n,2*a.m),constraints(ncon,nv),targets(ncon,2*a.m);int cr=0;auto g=gauss(order);
        for(int i=0;i<nf;++i){const auto& tr=tri[i];
            for(auto q:tr.quadrature(order)){auto s=tr.shape(q.p);auto phi=a.basis.phi(q.p);double w=q.w/cell.area;
                for(int j=0;j<8;++j){for(int l=0;l<8;++l)k(8*i+j,8*i+l)+=w*dot(s[j],s[l]);for(int c=0;c<2;++c)for(int l=0;l<3;++l)r(8*i+j,c*a.m+l)+=w*(c?s[j].y:s[j].x)*phi[l];}}
            constraints(cr,8*i+1)=1;constraints(cr,8*i+5)=1;constraints(cr+1,8*i+6)=3;constraints(cr+2,8*i+7)=3;
            for(int c=0;c<2;++c)for(int j=0;j<a.m;++j){const auto& d=c?a.gy:a.gx;double hdet=tr.h*tr.det;
                targets(cr,c*a.m+j)=hdet*d(0,j);targets(cr+1,c*a.m+j)=hdet*(d(1,j)*tr.j0.x+d(2,j)*tr.j0.y);targets(cr+2,c*a.m+j)=hdet*(d(1,j)*tr.j1.x+d(2,j)*tr.j1.y);}
            cr+=3;
            const auto& face=mesh.faces[cell.faces[faceLocal[i]]];double sign=face.owner==std::size_t(t)?1:-1;Vector2D S{sign*face.areaVector.x/tr.h,sign*face.areaVector.y/tr.h};
            for(auto [z,w]:g){double s=z-.5;Point2D p{face.centre.x-face.areaVector.y*s,face.centre.y+face.areaVector.x*s};auto shape=tr.shape(p);
                for(int l=0;l<2;++l){double wt=w*(l?s:1);for(int j=0;j<8;++j)constraints(cr+l,8*i+j)+=wt*dot(shape[j],S);
                    for(int c=0;c<2;++c)for(int j=0;j<2;++j)targets(cr+l,c*a.m+3+2*faceLocal[i]+j)+=wt*(j?s:1)*(c?S.y:S.x);}}
            cr+=2;
        }
        // Omit only the redundant mean normal continuity on one radial face;
        // its full trace is checked after reconstruction along with all others.
        for(int i=0;i<nf;++i){int prev=(i+nf-1)%nf;auto d=tri[i].a-cell.centre;Vector2D S{d.y/a.basis.h,-d.x/a.basis.h};
            for(int l=(i==0?1:0);l<2;++l){for(auto [z,w]:g){Point2D p{cell.centre.x+z*d.x,cell.centre.y+z*d.y};auto x=tri[i].shape(p),y=tri[prev].shape(p);double wt=w*(l?z-.5:1);
                for(int j=0;j<8;++j){constraints(cr,8*i+j)+=wt*dot(x[j],S);constraints(cr,8*prev+j)-=wt*dot(y[j],S);}}++cr;}}
        if(cr!=ncon)throw std::runtime_error("RT1 constraint size");
        for(int i=0;i<ncon;++i){double norm=0;for(int j=0;j<nv;++j)norm+=constraints(i,j)*constraints(i,j);norm=std::sqrt(norm);if(!(norm>0))throw std::runtime_error("zero RT1 constraint");
            for(int j=0;j<nv;++j)k(nv+i,j)=k(j,nv+i)=constraints(i,j)/norm;
            for(int j=0;j<2*a.m;++j)r(nv+i,j)=targets(i,j)/norm;}
        DenseLU lu(k.v,n);
        for(int j=0;j<2*a.m;++j){Vec rhs(n);for(int i=0;i<n;++i)rhs[i]=r(i,j);auto x=lu.solve(rhs);
            for(int pass=0;pass<1;++pass){auto rr=rhs;for(int i=0;i<n;++i)for(int l=0;l<n;++l)rr[i]-=k(i,l)*x[l];auto dx=lu.solve(rr);for(int i=0;i<n;++i)x[i]+=dx[i];}
            for(int i=0;i<nv;++i)coefficients(i,j)=x[i];
            for(int i=0;i<ncon;++i){double v=-targets(i,j);for(int l=0;l<nv;++l)v+=constraints(i,l)*x[l];divResidual=std::max(divResidual,std::abs(v));}
            for(int i=0;i<nf;++i){int prev=(i+nf-1)%nf;auto d=tri[i].a-cell.centre;Vector2D S{d.y/a.basis.h,-d.x/a.basis.h};
                for(double z:{0.,.5,1.}){Point2D p{cell.centre.x+z*d.x,cell.centre.y+z*d.y};auto s=tri[i].shape(p),v=tri[prev].shape(p);double jump=0;
                    for(int l=0;l<8;++l)jump+=dot(s[l],S)*x[8*i+l]-dot(v[l],S)*x[8*prev+l];
                    traceResidual=std::max(traceResidual,std::abs(jump));}}
        }
    }
    Vector2D value(int i,int j,Point2D p)const{auto s=tri[i].shape(p);Vector2D v{};for(int l=0;l<8;++l){v.x+=s[l].x*coefficients(8*i+l,j);v.y+=s[l].y*coefficients(8*i+l,j);}return v;}
};

double potentialAt(Point2D p){return p.x*p.x*p.x+p.x*p.y*p.y;}
Exact manufactured(Point2D p,const std::string& problem,double lambda){
    auto e=problem=="cylinder"?Exact{{0,0},{0,0},0}:exactAt(p,problem);
    e.p+=lambda*potentialAt(p);e.f.x+=lambda*(3*p.x*p.x+p.y*p.y);e.f.y+=lambda*2*p.x*p.y;return e;
}
struct Element {
    Local a;Mat matrix;Vec rhs;std::vector<int> inside,outside;Mat eliminated;Vec loadInternal;
    double trace=0,constraint=0,loadIdentity=0;
    Element(const Fixture& f,int t,const std::string& problem,double lambda,bool lifted,bool symmetric,int order):a(f,t,problem=="cylinder"?"couette":problem,order),matrix(2*a.m+3,2*a.m+3),rhs(2*a.m+3),eliminated(8,2*a.m-5) {
        int m=a.m;const auto& cell=f.mesh.cells[t];
        for(int c=0;c<2;++c)for(int i=0;i<m;++i)for(int j=0;j<m;++j){matrix(c*m+i,c*m+j)=a.stiffness(i,j);
            if(symmetric)for(int k=0;k<3;++k)for(int l=0;l<3;++l)matrix(c*m+i,c*m+j)+=a.mass(k,l)*(c?a.gy(k,i)*a.gy(l,j):a.gx(k,i)*a.gx(l,j));}
        if(symmetric)for(int i=0;i<m;++i)for(int j=0;j<m;++j)for(int k=0;k<3;++k)for(int l=0;l<3;++l){double v=a.mass(k,l)*a.gy(k,i)*a.gx(l,j);matrix(i,m+j)+=v;matrix(m+j,i)+=v;}
        for(int c=0;c<2;++c)for(int i=0;i<m;++i)for(int p=0;p<3;++p)matrix(c*m+i,2*m+p)=matrix(2*m+p,c*m+i)=-(c?a.by(p,i):a.bx(p,i));
        if(lifted){Lift lift(f,t,a,order);trace=lift.traceResidual;constraint=lift.divResidual;Vec identity(2*m);
            for(std::size_t tri=0;tri<lift.tri.size();++tri)for(auto q:lift.tri[tri].quadrature(order)){
                auto force=manufactured(q.p,problem,lambda).f;Vector2D grad{3*q.p.x*q.p.x+q.p.y*q.p.y,2*q.p.x*q.p.y};auto phi=a.basis.phi(q.p);
                for(int j=0;j<2*m;++j){auto v=lift.value(int(tri),j,q.p);rhs[j]+=q.w*dot(force,v);int c=j/m,k=j%m;double div=0;for(int l=0;l<3;++l)div+=phi[l]*(c?a.gy(l,k):a.gx(l,k));identity[j]+=q.w*(dot(grad,v)+potentialAt(q.p)*div);}}
            for(std::size_t i=0;i<cell.faces.size();++i){const auto& face=f.mesh.faces[cell.faces[i]];double sign=face.owner==std::size_t(t)?1:-1;
                for(auto [z,w]:gauss(order)){double s=z-.5;Point2D p{face.centre.x-face.areaVector.y*s,face.centre.y+face.areaVector.x*s};
                    for(int c=0;c<2;++c)for(int j=0;j<2;++j)identity[c*m+3+2*i+j]-=w*potentialAt(p)*(j?s:1)*sign*(c?face.areaVector.y:face.areaVector.x);}}
            for(double x:identity)loadIdentity=std::max(loadIdentity,std::abs(x));
        }else for(auto q:a.q){auto e=manufactured(q.p,problem,lambda);auto phi=a.basis.phi(q.p);for(int k=0;k<3;++k){rhs[k]+=q.w*phi[k]*e.f.x;rhs[m+k]+=q.w*phi[k]*e.f.y;}}
        inside={0,1,2,m,m+1,m+2,2*m+1,2*m+2};
        for(int i=0;i<2*m+3;++i)if(std::find(inside.begin(),inside.end(),i)==inside.end())outside.push_back(i);
        Mat ii(8,8);for(int i=0;i<8;++i)for(int j=0;j<8;++j)ii(i,j)=matrix(inside[i],inside[j]);DenseLU lu(ii.v,8);
        for(std::size_t j=0;j<outside.size();++j){Vec v(8);for(int i=0;i<8;++i)v[i]=matrix(inside[i],outside[j]);v=lu.solve(v);for(int i=0;i<8;++i)eliminated(i,int(j))=v[i];}
        Vec r(8);for(int i=0;i<8;++i)r[i]=rhs[inside[i]];loadInternal=lu.solve(r);
    }
    std::pair<Mat,Vec> condensed()const {int n=int(outside.size());Mat k(n,n);Vec b(n);for(int i=0;i<n;++i){b[i]=rhs[outside[i]];for(int l=0;l<8;++l)b[i]-=matrix(outside[i],inside[l])*loadInternal[l];
        for(int j=0;j<n;++j){k(i,j)=matrix(outside[i],outside[j]);for(int l=0;l<8;++l)k(i,j)-=matrix(outside[i],inside[l])*eliminated(l,j);}}return {k,b};}
    Vec recover(const Vec& ext)const{Vec v(rhs.size());for(std::size_t i=0;i<outside.size();++i)v[outside[i]]=ext[i];for(int i=0;i<8;++i){v[inside[i]]=loadInternal[i];for(std::size_t j=0;j<outside.size();++j)v[inside[i]]-=eliminated(i,int(j))*ext[j];}return v;}
};
Fixture readFixture(const std::string& name,int n){
    if(name=="square"||name=="sheared"||name=="cut")return grid(n,name=="sheared"?.7:0,name=="cut");
    if(name=="tip")return makeFixture({{{{0,0},{1,0},{0,1}}},{{{1,0},{1,1},{0,1}}},{{{1,0},{2,0},{2,1},{1,1}}},{{{0,1},{1,1},{1,2},{0,2}}}});
    if(name=="split")return makeFixture({{{{0,0},{1,0},{1,1},{1,2},{0,2}}},{{{1,0},{2,0},{2,1},{1,1}}},{{{1,1},{2,1},{2,2},{1,2}}}});
    auto read=readCm2dTopology(name);if(!read.valid())throw std::runtime_error(read.error);Fixture f{makeFvMesh2D(read.topology),{}, {}};
    for(const auto& c:f.mesh.cells){double h=0;std::vector<Point2D> pts;for(auto fi:c.faces){const auto& face=f.mesh.faces[fi];pts.push_back({face.centre.x+face.areaVector.y/2,face.centre.y-face.areaVector.x/2});pts.push_back({face.centre.x-face.areaVector.y/2,face.centre.y+face.areaVector.x/2});}
        for(auto p:pts)for(auto q:pts)h=std::max(h,std::hypot(p.x-q.x,p.y-q.y));
        f.diameter.push_back(h);f.meanY2.push_back(0);}return f;
}
std::string optionalMetric(double x,bool valid){if(!valid)return "null";std::ostringstream s;s<<std::setprecision(17)<<x;return s.str();}
struct Entry {std::int64_t i,j;double value;};
static_assert(sizeof(Entry)==24);
#ifndef CARTMESH_P1_STRESS_NO_MAIN
int main(int argc,char** argv)try{
    if(argc!=10)throw std::runtime_error("usage: p1-stress assemble|recover mesh n problem lambda lift|cell symmetric|laplace order prefix");
    const std::string mode=argv[1],name=argv[2],problem=argv[4],prefix=argv[9];int n=std::stoi(argv[3]),order=std::stoi(argv[8]);double lambda=std::stod(argv[5]);bool lifted=std::string(argv[6])=="lift",symmetric=std::string(argv[7])=="symmetric";
    if((std::string(argv[6])!="lift"&&std::string(argv[6])!="cell")||(std::string(argv[7])!="symmetric"&&std::string(argv[7])!="laplace")||(mode!="assemble"&&mode!="recover")||order<4||order>12)throw std::runtime_error("invalid mode/order");
    if((problem=="noslip"&&name!="square")||(problem=="noslip-sheared"&&name!="sheared"))throw std::runtime_error("stationary manufactured geometry mismatch");
    auto start=std::chrono::steady_clock::now();auto f=readFixture(name,n);const auto& mesh=f.mesh;
    if(problem=="cylinder")for(const auto& face:mesh.faces)if(!face.neighbour&&face.patch!=BoundaryPatch2D::DomainBoundary&&face.patch!=BoundaryPatch2D::EmbeddedBoundary)throw std::runtime_error("cylinder control requires explicit Domain/Embedded boundary patches");
    int nc=int(mesh.cells.size()),nf=int(mesh.faces.size()),raw=4*nf+nc,count=0;std::vector<int> map(raw,-1);Vec known(raw);
    for(int i=0;i<nf;++i)if(mesh.faces[i].neighbour)for(int j=0;j<4;++j)map[4*i+j]=count++;
    for(int i=0;i<nc-1;++i)map[4*nf+i]=count++;
    for(int i=0;i<nf;++i)if(!mesh.faces[i].neighbour){const auto& face=mesh.faces[i];for(auto [z,w]:gauss(order)){double s=z-.5;Point2D p{face.centre.x-face.areaVector.y*s,face.centre.y+face.areaVector.x*s};
        Vector2D u;if(problem=="noslip"||problem=="noslip-sheared"){u={0,0};}else if(problem=="cylinder"){bool outer=face.patch==BoundaryPatch2D::DomainBoundary;u={outer?1.:0.,0};}else u=manufactured(p,problem,lambda).u;
        known[4*i]+=w*u.x;known[4*i+1]+=12*w*s*u.x;known[4*i+2]+=w*u.y;known[4*i+3]+=12*w*s*u.y;}}
    Vec solution;std::ofstream entries,rhsout,fields,faces;std::uint64_t nnz=0;Vec rhs(count);double trace=0,constraint=0,identity=0,symmetry=0,area=0,gauge=0,urms=0,prms=0,pmax=0,divmax=0,umax=0,pmin=1e300,pmaxfield=-1e300,work=0,energy=0,forceWork=0,internalResidual=0;Vec reaction(4*nf),stressTrace(4*nf);std::vector<Vec> states;
    if(mode=="assemble")entries.open(prefix+".entries",std::ios::binary);
    else {std::ifstream in(prefix+".solution",std::ios::binary);solution.resize(count);in.read(reinterpret_cast<char*>(solution.data()),count*sizeof(double));if(!in||in.peek()!=EOF)throw std::runtime_error("missing/truncated/trailing sparse solution");for(int i=0;i<raw;++i)if(map[i]>=0)known[i]=solution[map[i]];fields.open(prefix+".cells.csv");faces.open(prefix+".faces.csv");fields<<std::setprecision(17)<<"cell,x,y,area,h,mean_X2,mean_XY,mean_Y2,potential0,potentialX,potentialY,u0,uX,uY,v0,vX,vY,p0,pX,pY,ru0,ruX,ruY,ruXX,ruXY,ruYY,rv0,rvX,rvY,rvXX,rvXY,rvYY\n";faces<<std::setprecision(17)<<"face,x,y,Sx,Sy,owner,neighbour,u0,us,v0,vs,traction_x,traction_y,stress_x,stress_y\n";}
    for(int t=0;t<nc;++t){Element e(f,t,problem,lambda,lifted,symmetric,order);const auto& a=e.a;int m=a.m;trace=std::max(trace,e.trace);constraint=std::max(constraint,e.constraint);identity=std::max(identity,e.loadIdentity);area+=mesh.cells[t].area;
        std::vector<int> ids;for(int j:e.outside){if(j==2*m)ids.push_back(4*nf+t);else {int c=j/m,k=j%m;ids.push_back(4*int(mesh.cells[t].faces[(k-3)/2])+2*c+(k-3)%2);}}
        auto [k,b]=e.condensed();for(int i=0;i<k.nr;++i)for(int j=0;j<k.nc;++j)symmetry=std::max(symmetry,std::abs(k(i,j)-k(j,i)));
        if(mode=="assemble"){for(std::size_t i=0;i<ids.size();++i){int row=map[ids[i]];if(row<0)continue;rhs[row]+=b[i];for(std::size_t j=0;j<ids.size();++j){int col=map[ids[j]];double v=k(int(i),int(j));if(col<0)rhs[row]-=v*known[ids[j]];else if(v!=0){Entry z{row,col,v};entries.write(reinterpret_cast<char*>(&z),sizeof z);++nnz;}}}}
        else {Vec ext;for(int id:ids)ext.push_back(known[id]);auto v=e.recover(ext);states.push_back(v);P6 ru{},rv{};
            for(int l=0;l<6;++l)for(int j=0;j<m;++j){ru[l]+=a.potential(l,j)*v[j];rv[l]+=a.potential(l,j)*v[m+j];}
            Vec div(3);for(int p=0;p<3;++p)for(int j=0;j<m;++j)div[p]+=a.bx(p,j)*v[j]+a.by(p,j)*v[m+j];for(int p=0;p<3;++p)work-=v[2*m+p]*div[p];div=DenseLU(a.mass.v,3).solve(div);
            for(auto q:a.q){auto phi=a.basis.phi(q.p);auto theta=a.basis.theta(q.p);auto ex=manufactured(q.p,problem,lambda);double u=0,w=0,p=0,d=0;for(int j=0;j<6;++j){u+=theta[j]*ru[j];w+=theta[j]*rv[j];}for(int j=0;j<3;++j){p+=phi[j]*v[2*m+j];d+=phi[j]*div[j];}
                urms+=q.w*((u-ex.u.x)*(u-ex.u.x)+(w-ex.u.y)*(w-ex.u.y));gauge+=q.w*(p-ex.p);prms+=q.w*(p-ex.p)*(p-ex.p);umax=std::max(umax,std::hypot(u,w));pmin=std::min(pmin,p);pmaxfield=std::max(pmaxfield,p);divmax=std::max(divmax,std::abs(d));}
            for(std::size_t fi=0;fi<mesh.cells[t].faces.size();++fi){int id=int(mesh.cells[t].faces[fi]);const auto& face=mesh.faces[id];if(face.neighbour)continue;
                for(auto [z,w]:gauss(order)){double ss=z-.5;auto S=face.areaVector;auto phi=a.basis.phi({face.centre.x-ss*S.y,face.centre.y+ss*S.x});double ux=0,uy=0,vx=0,vy=0,p=0;
                    for(int l=0;l<3;++l){p+=phi[l]*v[2*m+l];for(int j=0;j<m;++j){ux+=phi[l]*a.gx(l,j)*v[j];uy+=phi[l]*a.gy(l,j)*v[j];vx+=phi[l]*a.gx(l,j)*v[m+j];vy+=phi[l]*a.gy(l,j)*v[m+j];}}
                    double tx=(symmetric?2*ux:ux)*S.x+(symmetric?uy+vx:uy)*S.y-p*S.x,ty=(symmetric?uy+vx:vx)*S.x+(symmetric?2*vy:vy)*S.y-p*S.y;
                    for(int j=0;j<2;++j){stressTrace[4*id+j]+=w*(j?ss:1)*tx;stressTrace[4*id+2+j]+=w*(j?ss:1)*ty;}}}
            for(int i=0;i<2*m;++i){forceWork+=v[i]*e.rhs[i];for(int j=0;j<2*m;++j)energy+=v[i]*e.matrix(i,j)*v[j];}
            for(int i=0;i<2*m+3;++i){double r=-e.rhs[i];for(int j=0;j<2*m+3;++j)r+=e.matrix(i,j)*v[j];if(std::find(e.inside.begin(),e.inside.end(),i)!=e.inside.end())internalResidual=std::max(internalResidual,std::abs(r));
                else if(i<2*m){int c=i/m,j=i%m;reaction[4*mesh.cells[t].faces[(j-3)/2]+2*c+(j-3)%2]+=r;}}
            fields<<t<<','<<mesh.cells[t].centre.x<<','<<mesh.cells[t].centre.y<<','<<mesh.cells[t].area<<','<<a.basis.h;for(double z:a.basis.moments)fields<<','<<z;
            Vec projected(3);for(auto q:a.q){auto phi=a.basis.phi(q.p);for(int j=0;j<3;++j)projected[j]+=q.w*phi[j]*potentialAt(q.p);}projected=DenseLU(a.mass.v,3).solve(projected);for(double z:projected)fields<<','<<z;
            for(int c=0;c<2;++c)for(int j=0;j<3;++j)fields<<','<<v[c*m+j];
            for(int j=0;j<3;++j)fields<<','<<v[2*m+j];
            for(double z:ru)fields<<','<<z;
            for(double z:rv)fields<<','<<z;
            fields<<'\n';
        }
        if((t+1)%2000==0)std::cerr<<mode<<" cells "<<t+1<<'/'<<nc<<'\n';
    }
    double boundaryWork=0,interiorResidual=0,wallFx=0,wallFy=0,stressFx=0,stressFy=0,tractionDifference=0,wallNormal=0;
    if(mode=="assemble"){entries.close();rhsout.open(prefix+".rhs",std::ios::binary);rhsout.write(reinterpret_cast<char*>(rhs.data()),rhs.size()*sizeof(double));rhsout.close();if(!entries||!rhsout)throw std::runtime_error("matrix write failure");}
    else {for(int i=0;i<nf;++i){const auto& face=mesh.faces[i];faces<<i<<','<<face.centre.x<<','<<face.centre.y<<','<<face.areaVector.x<<','<<face.areaVector.y<<','<<face.owner<<','<<(face.neighbour?int(*face.neighbour):-1);for(int j=0;j<4;++j)faces<<','<<known[4*i+j];faces<<','<<reaction[4*i]<<','<<reaction[4*i+2]<<','<<stressTrace[4*i]<<','<<stressTrace[4*i+2]<<'\n';
        for(int j=0;j<4;++j)if(face.neighbour)interiorResidual=std::max(interiorResidual,std::abs(reaction[4*i+j]));else boundaryWork+=known[4*i+j]*reaction[4*i+j];
        if(!face.neighbour){double length=std::hypot(face.areaVector.x,face.areaVector.y);for(int j=0;j<4;++j)tractionDifference=std::max(tractionDifference,std::abs(reaction[4*i+j]-stressTrace[4*i+j])/length);
            if(face.patch==BoundaryPatch2D::EmbeddedBoundary){wallFx-=reaction[4*i];wallFy-=reaction[4*i+2];stressFx-=stressTrace[4*i];stressFy-=stressTrace[4*i+2];
                for(double ss:{-.5,.5})wallNormal=std::max(wallNormal,std::abs((known[4*i]+ss*known[4*i+1])*face.areaVector.x+(known[4*i+2]+ss*known[4*i+3])*face.areaVector.y)/length);}}}
        fields.close();faces.close();if(!fields||!faces)throw std::runtime_error("field write failure");gauge/=area;prms=0;
        for(int t=0;t<nc;++t){const auto& v=states[t];Basis basis{mesh.cells[t].centre,f.diameter[t],{}};int m=(int(v.size())-3)/2;
            for(auto q:cellQuadrature(mesh,t,order)){auto phi=basis.phi(q.p);double p=0;for(int j=0;j<3;++j)p+=phi[j]*v[2*m+j];double d=p-manufactured(q.p,problem,lambda).p-gauge;prms+=q.w*d*d;}}
        prms=std::sqrt(prms/area);pmax=pmaxfield-pmin;}
    std::cout<<std::setprecision(17)<<"{\"mode\":\""<<mode<<"\",\"cells\":"<<nc<<",\"faces\":"<<nf<<",\"unknowns\":"<<count<<",\"entries\":"<<nnz<<",\"area\":"<<area<<",\"traceResidual\":"<<trace<<",\"constraintResidual\":"<<constraint<<",\"loadIdentity\":"<<identity<<",\"condensedAsymmetry\":"<<symmetry<<",\"velocityRms\":"<<optionalMetric(std::sqrt(urms/area),problem!="cylinder")<<",\"pressureRms\":"<<optionalMetric(prms,problem!="cylinder")<<",\"maxSpeed\":"<<umax<<",\"pressureRange\":"<<pmax<<",\"maxDivergence\":"<<divmax<<",\"pressureWork\":"<<work<<",\"viscousEnergy\":"<<energy<<",\"forceWork\":"<<forceWork<<",\"boundaryWork\":"<<boundaryWork<<",\"energyBalance\":"<<energy+work-forceWork-boundaryWork<<",\"internalResidual\":"<<internalResidual<<",\"interiorFaceResidual\":"<<interiorResidual<<",\"innerBodyFx\":"<<wallFx<<",\"innerBodyFy\":"<<wallFy<<",\"stressBodyFx\":"<<stressFx<<",\"stressBodyFy\":"<<stressFy<<",\"tractionMomentDifferencePerLength\":"<<tractionDifference<<",\"wallNormalSpeed\":"<<wallNormal<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}

#endif
