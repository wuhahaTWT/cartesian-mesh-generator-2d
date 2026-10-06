// Degree-one native hybrid Stokes research; native polygon fixtures and dense LU.
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

// P1 cell/face velocity and P1 cell pressure. Reconstructed gradients are P1;
// the potential reconstruction is P2. Its mean is the cell velocity mean.
// Stabilisation is the face projection of r - pi_T^1 r + u_T minus u_F,
// scaled by face length / cell diameter, as in the order-zero comparison.
// The pressure block is the exact transpose of the same weak divergence.
// The load here is the ordinary cell L2 load. This is not the full general
// pressure-robust H(div) reconstruction or a Navier--Stokes solver.

struct Mat {
    int nr=0,nc=0;Vec v;
    Mat(int r,int c):nr(r),nc(c),v(r*c){}
    double& operator()(int r,int c){return v[r*nc+c];}
    double operator()(int r,int c)const{return v[r*nc+c];}
};
struct Q {Point2D p;double w;};
using P3=std::array<double,3>;
using P6=std::array<double,6>;

std::vector<std::pair<double,double>> gauss(int n) {
    std::vector<std::pair<double,double>> q;
    const double pi=std::acos(-1.);
    for(int i=0;i<n;++i) {
        double x=std::cos(pi*(i+.75)/(n+.5)),derivative=0;
        for(int it=0;it<30;++it) {
            double p0=1,p1=x;
            for(int k=2;k<=n;++k){const double p=((2*k-1)*x*p1-(k-1)*p0)/k;p0=p1;p1=p;}
            derivative=n*(x*p1-p0)/(x*x-1);
            const double dx=p1/derivative;x-=dx;
            if(std::abs(dx)<2e-16)break;
        }
        // Re-evaluate derivative at the final root.
        double p0=1,p1=x;
        for(int k=2;k<=n;++k){const double p=((2*k-1)*x*p1-(k-1)*p0)/k;p0=p1;p1=p;}
        derivative=n*(x*p1-p0)/(x*x-1);
        q.push_back({.5*(1+x),1./((1-x*x)*derivative*derivative)});
    }
    return q;
}

std::vector<Q> cellQuadrature(const FvMesh2D& mesh,int t,int order) {
    std::vector<Q> q;const auto c=mesh.cells[t].centre;const auto rule=gauss(order);
    for(const auto f:mesh.cells[t].faces) {
        const auto& face=mesh.faces[f];const double s=face.owner==std::size_t(t)?1:-1;
        const Point2D a{face.centre.x+s*face.areaVector.y/2,face.centre.y-s*face.areaVector.x/2};
        const Point2D b{face.centre.x-s*face.areaVector.y/2,face.centre.y+s*face.areaVector.x/2};
        const double jac=(a.x-c.x)*(b.y-c.y)-(a.y-c.y)*(b.x-c.x);
        for(const auto& [r,wr]:rule)for(const auto& [z,wz]:rule)
            q.push_back({{(1-r)*c.x+r*(1-z)*a.x+r*z*b.x,(1-r)*c.y+r*(1-z)*a.y+r*z*b.y},wr*wz*r*jac});
    }
    return q;
}

struct Basis {
    Point2D c;double h;P3 moments{};
    P3 phi(Point2D p)const{return {1,(p.x-c.x)/h,(p.y-c.y)/h};}
    P6 theta(Point2D p)const{const auto a=phi(p);return {1,a[1],a[2],a[1]*a[1]-moments[0],a[1]*a[2]-moments[1],a[2]*a[2]-moments[2]};}
    std::array<Vector2D,6> grad(Point2D p)const{const auto a=phi(p);return {{{0,0},{1/h,0},{0,1/h},{2*a[1]/h,0},{a[2]/h,a[1]/h},{0,2*a[2]/h}}};}
};

struct Exact {Vector2D u,f;double p;std::array<Vector2D,2> gradient{};};
Exact exactAt(Point2D p,const std::string& problem) {
    if(problem=="couette")return {{p.y,0},{0,0},0,{{{0,1},{0,0}}}};
    if(problem=="rotation")return {{-p.y,p.x},{0,0},0,{{{0,-1},{1,0}}}};
    if(problem=="hydrostatic")return {{0,0},{1,2},p.x+2*p.y};
    if(problem=="poiseuille")return {{4*p.y*(1-p.y),0},{0,0},-8*p.x,{{{0,4-8*p.y},{0,0}}}};
    // Same velocity/gradient, with the static-pressure reference p=0 at x=1.
    if(problem=="outlet-poiseuille")return {{4*p.y*(1-p.y),0},{0,0},8*(1-p.x),{{{0,4-8*p.y},{0,0}}}};
    if(problem=="vortex") {
        const double pi=std::acos(-1.),sx=std::sin(pi*p.x),sy=std::sin(pi*p.y),cx=std::cos(pi*p.x),cy=std::cos(pi*p.y);
        const Vector2D u{sx*cy,-cx*sy};
        return {u,{2*pi*pi*u.x+pi*cx*sy,2*pi*pi*u.y+pi*sx*cy},sx*sy};
    }
    if(problem=="noslip" || problem=="noslip-sheared") {
        // psi=64*g(x-s*y)*g(y), g(t)=t^2(1-t)^2. Both g and g' vanish
        // at each physical wall, giving smooth nonzero interior flow with
        // stationary, impermeable, no-slip boundaries on the parallelogram.
        const double s=problem=="noslip-sheared"?.7:0.,x=p.x-s*p.y,y=p.y;
        const auto g=[](double t){return std::array<double,4>{t*t*(1-t)*(1-t),
            2*t*(1-t)*(1-2*t),2-12*t+12*t*t,-12+24*t};};
        const auto a=g(x),b=g(y);const double pi=std::acos(-1.);
        const double px=pi*std::cos(pi*x)*std::sin(pi*y);
        const double py=pi*std::sin(pi*x)*std::cos(pi*y)-s*px;
        const double dlx=64*((1+s*s)*a[3]*b[0]-2*s*a[2]*b[1]+a[1]*b[2]);
        const double dly=64*(-s*(1+s*s)*a[3]*b[0]+(1+3*s*s)*a[2]*b[1]
                            -3*s*a[1]*b[2]+a[0]*b[3]);
        const double ux=64*(a[1]*b[1]-s*a[2]*b[0]);
        const double uy=64*(a[0]*b[2]-2*s*a[1]*b[1]+s*s*a[2]*b[0]);
        const double vx=-64*a[2]*b[0];
        return {{64*(a[0]*b[1]-s*a[1]*b[0]),-64*a[1]*b[0]},
                {-dly+px,dlx+py},std::sin(pi*x)*std::sin(pi*y),
                {{{ux,uy},{vx,-ux}}}};
    }
    throw std::runtime_error("unknown manufactured Stokes problem");
}

struct Local {
    Basis basis;int m;std::vector<Q> q;
    Mat mass,bx,by,gx,gy,potential,stiffness;std::array<Vec,2> force;
    Local(const Fixture& fixture,int t,const std::string& problem,int order):
        basis{fixture.mesh.cells[t].centre,fixture.diameter[t],{}},
        m(3+2*int(fixture.mesh.cells[t].faces.size())),q(cellQuadrature(fixture.mesh,t,order)),
        mass(3,3),bx(3,m),by(3,m),gx(3,m),gy(3,m),potential(6,m),stiffness(m,m),force{Vec(3),Vec(3)} {
        const auto& mesh=fixture.mesh;const auto& cell=mesh.cells[t];
        for(const auto& v:q) {
            const auto phi=basis.phi(v.p);const auto exact=exactAt(v.p,problem);
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)mass(i,j)+=v.w*phi[i]*phi[j];
            for(int j=0;j<3;++j) {
                bx(1,j)-=v.w*phi[j]/basis.h;by(2,j)-=v.w*phi[j]/basis.h;
                force[0][j]+=v.w*phi[j]*exact.f.x;force[1][j]+=v.w*phi[j]*exact.f.y;
            }
        }
        basis.moments={mass(1,1)/cell.area,mass(1,2)/cell.area,mass(2,2)/cell.area};
        Mat k(5,5),r(5,m),project(3,6);
        for(const auto& v:q) {
            const auto phi=basis.phi(v.p);const auto theta=basis.theta(v.p);const auto grad=basis.grad(v.p);
            for(int i=0;i<5;++i)for(int j=0;j<5;++j)k(i,j)+=v.w*dot(grad[i+1],grad[j+1]);
            for(int i=0;i<3;++i)for(int j=0;j<6;++j)project(i,j)+=v.w*phi[i]*theta[j];
            for(int j=0;j<3;++j){r(2,j)-=2*v.w*phi[j]/(basis.h*basis.h);r(4,j)-=2*v.w*phi[j]/(basis.h*basis.h);}
        }
        const auto rule=gauss(order);
        for(std::size_t f=0;f<cell.faces.size();++f) {
            const auto& face=mesh.faces[cell.faces[f]];const double sign=face.owner==std::size_t(t)?1:-1;
            const auto S=face.areaVector;
            for(const auto& [z,w]:rule) {
                const double s=z-.5;const Point2D p{face.centre.x-S.y*s,face.centre.y+S.x*s};
                const auto phi=basis.phi(p);const auto grad=basis.grad(p);
                for(int j=0;j<2;++j) {
                    const double psi=j?s:1;const int col=3+2*int(f)+j;
                    for(int i=0;i<3;++i){bx(i,col)+=w*psi*sign*S.x*phi[i];by(i,col)+=w*psi*sign*S.y*phi[i];}
                    for(int i=0;i<5;++i)r(i,col)+=w*psi*sign*dot(S,grad[i+1]);
                }
            }
        }
        const DenseLU ml(mass.v,3),kl(k.v,5);
        for(int j=0;j<m;++j) {
            Vec x(3),y(3),z(5);for(int i=0;i<3;++i){x[i]=bx(i,j);y[i]=by(i,j);}for(int i=0;i<5;++i)z[i]=r(i,j);
            x=ml.solve(x);y=ml.solve(y);z=kl.solve(z);
            for(int i=0;i<3;++i){gx(i,j)=x[i];gy(i,j)=y[i];}for(int i=0;i<5;++i)potential(i+1,j)=z[i];
        }
        potential(0,0)=1;
        for(int j=0;j<6;++j){Vec x(3);for(int i=0;i<3;++i)x[i]=project(i,j);x=ml.solve(x);for(int i=0;i<3;++i)project(i,j)=x[i];}
        for(int i=0;i<m;++i)for(int j=0;j<m;++j)for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            stiffness(i,j)+=mass(a,b)*(gx(a,i)*gx(b,j)+gy(a,i)*gy(b,j));
        for(std::size_t f=0;f<cell.faces.size();++f) {
            const auto& face=mesh.faces[cell.faces[f]];const auto S=face.areaVector;const double length=std::hypot(S.x,S.y);
            Mat residual(2,m);
            for(const auto& [z,w]:rule) {
                const double s=z-.5;const Point2D p{face.centre.x-S.y*s,face.centre.y+S.x*s};
                const auto phi=basis.phi(p);auto theta=basis.theta(p);
                for(int a=0;a<6;++a)for(int b=0;b<3;++b)theta[a]-=phi[b]*project(b,a);
                for(int a=0;a<2;++a)for(int j=0;j<m;++j) {
                    double value=j<3?phi[j]:0.;for(int b=0;b<6;++b)value+=theta[b]*potential(b,j);
                    residual(a,j)-=w*(a?12*s:1)*value;
                }
            }
            residual(0,3+2*int(f))+=1;residual(1,4+2*int(f))+=1;
            for(int i=0;i<m;++i)for(int j=0;j<m;++j)
                stiffness(i,j)+=length/basis.h*(residual(0,i)*residual(0,j)+residual(1,i)*residual(1,j)/12.);
        }
    }
};

#ifndef CARTMESH_HYBRID_STOKES_P1_NO_MAIN
int main(int argc,char** argv)try {
    if(argc<4||argc>6)throw std::runtime_error("usage: probe square|sheared|cut|split|tip n problem [quadrature_order] [output_prefix]");
    const std::string name=argv[1],problem=argv[3],output=argc==6?argv[5]:"";const int resolution=std::stoi(argv[2]),order=argc>=5?std::stoi(argv[4]):8;
    if(resolution<1||resolution>12||order<3||order>16)throw std::runtime_error("research size/order outside cost cap");
    const bool stationary=problem=="noslip" || problem=="noslip-sheared";
    if((problem=="noslip" && name!="square") || (problem=="noslip-sheared" && name!="sheared"))
        throw std::runtime_error("stationary manufacture requires its matching square/sheared domain");
    const auto started=std::chrono::steady_clock::now();
    if(name!="square"&&name!="sheared"&&name!="cut"&&name!="split"&&name!="tip")throw std::runtime_error("unknown mesh fixture");
    const Fixture fixture=name=="split"?makeFixture({{{{0,0},{1,0},{1,1},{1,2},{0,2}}},{{{1,0},{2,0},{2,1},{1,1}}},{{{1,1},{2,1},{2,2},{1,2}}}}):
        name=="tip"?makeFixture({{{{0,0},{1,0},{0,1}}},{{{1,0},{1,1},{0,1}}},{{{1,0},{2,0},{2,1},{1,1}}},{{{0,1},{1,1},{1,2},{0,2}}}}):
        grid(resolution,name=="sheared"?.7:0,name=="cut");const auto& mesh=fixture.mesh;
    const int nc=int(mesh.cells.size()),nf=int(mesh.faces.size()),rawNv=6*nc+4*nf;
    std::vector<int> map(rawNv,-1);int nv=6*nc;std::iota(map.begin(),map.begin()+nv,0);
    for(int f=0;f<nf;++f)if(mesh.faces[f].neighbour)for(int j=0;j<4;++j)map[6*nc+4*f+j]=nv++;
    const int n=nv+3*nc-1;if(n>2600)throw std::runtime_error("dense driver cost cap 2600 unknowns");
    const auto pIndex=[&](int t,int k){const int raw=3*t+k,gauge=3*(nc-1);return raw==gauge?-1:nv+raw-(raw>gauge?1:0);};
    const auto id=[&](int t,int comp,int i){if(i<3)return 6*t+3*comp+i;const int f=int(mesh.cells[t].faces[(i-3)/2]);return 6*nc+4*f+2*comp+(i-3)%2;};
    Vec known(rawNv),matrix(n*n),rhs(n);std::vector<Local> local;local.reserve(nc);const auto rule=gauss(order);
    for(int f=0;f<nf;++f)if(!mesh.faces[f].neighbour)for(const auto& [z,w]:rule) {
        const auto& face=mesh.faces[f];const double s=z-.5;const auto u=stationary?Vector2D{}:exactAt({face.centre.x-face.areaVector.y*s,face.centre.y+face.areaVector.x*s},problem).u;
        known[6*nc+4*f]+=w*u.x;known[6*nc+4*f+1]+=12*w*s*u.x;known[6*nc+4*f+2]+=w*u.y;known[6*nc+4*f+3]+=12*w*s*u.y;
    }
    for(int t=0;t<nc;++t) {
        local.emplace_back(fixture,t,problem,order);const auto& a=local.back();
        for(int c=0;c<2;++c)for(int i=0;i<a.m;++i) {
            const int row=map[id(t,c,i)];if(row<0)continue;
            if(i<3)rhs[row]+=a.force[c][i];
            for(int j=0;j<a.m;++j){const int raw=id(t,c,j),col=map[raw];if(col>=0)matrix[row*n+col]+=a.stiffness(i,j);else rhs[row]-=a.stiffness(i,j)*known[raw];}
        }
        for(int p=0;p<3;++p) {
            const int row=pIndex(t,p);if(row<0)continue;
            for(int c=0;c<2;++c)for(int i=0;i<a.m;++i) {
                const double b=c?a.by(p,i):a.bx(p,i);const int raw=id(t,c,i),col=map[raw];
                if(col>=0){matrix[row*n+col]-=b;matrix[col*n+row]-=b;}else rhs[row]+=b*known[raw];
            }
        }
    }
    const Vec x=DenseLU(matrix,n).solve(rhs);Vec u=known;
    for(const double value:x)if(!std::isfinite(value))throw std::runtime_error("nonfinite solution");
    for(int i=0;i<rawNv;++i)if(map[i]>=0)u[i]=x[map[i]];
    double area=0,gauge=0,umean=0,ur=0,pr=0,pmax=0,divmax=0,meanDiv=0,pressureWork=0,residual=0;
    for(int t=0;t<nc;++t)for(const auto& q:local[t].q) {
        const auto phi=local[t].basis.phi(q.p);double p=0;for(int k=0;k<3;++k)if(pIndex(t,k)>=0)p+=phi[k]*x[pIndex(t,k)];
        area+=q.w;gauge+=q.w*(p-exactAt(q.p,problem).p);
    }
    gauge/=area;
    std::ofstream fields,faces;
    if(!output.empty()){fields.open(output+".cells.csv");faces.open(output+".faces.csv");fields<<std::setprecision(17)<<"cell,x,y,area,h,mean_X2,mean_XY,mean_Y2,u0,uX,uY,v0,vX,vY,p0,pX,pY,ru0,ruX,ruY,ruXX,ruXY,ruYY,rv0,rvX,rvY,rvXX,rvXY,rvYY\n";faces<<std::setprecision(17)<<"face,owner,neighbour,x,y,Sx,Sy,u0,us,v0,vs\n";}
    for(int t=0;t<nc;++t) {
        const auto& a=local[t];Vec uc(a.m),vc(a.m),div(3);P6 ru{},rv{};P3 p{};Vector2D ue{};
        for(int j=0;j<a.m;++j){uc[j]=u[id(t,0,j)];vc[j]=u[id(t,1,j)];}
        for(int k=0;k<3;++k) {
            if(pIndex(t,k)>=0)p[k]=x[pIndex(t,k)];
            for(int j=0;j<a.m;++j)div[k]+=a.bx(k,j)*uc[j]+a.by(k,j)*vc[j];
            pressureWork-=p[k]*div[k];
        }
        meanDiv=std::max(meanDiv,std::abs(div[0])/mesh.cells[t].area);div=DenseLU(a.mass.v,3).solve(div);
        for(int k=0;k<6;++k)for(int j=0;j<a.m;++j){ru[k]+=a.potential(k,j)*uc[j];rv[k]+=a.potential(k,j)*vc[j];}
        for(const auto& q:a.q) {
            const auto phi=a.basis.phi(q.p);const auto theta=a.basis.theta(q.p);const auto e=exactAt(q.p,problem);
            double ux=0,vy=0,pv=0;for(int k=0;k<6;++k){ux+=theta[k]*ru[k];vy+=theta[k]*rv[k];}for(int k=0;k<3;++k)pv+=phi[k]*p[k];
            ur+=q.w*((ux-e.u.x)*(ux-e.u.x)+(vy-e.u.y)*(vy-e.u.y));pr+=q.w*(pv-e.p-gauge)*(pv-e.p-gauge);pmax=std::max(pmax,std::abs(pv-e.p-gauge));
            ue.x+=q.w*e.u.x/mesh.cells[t].area;ue.y+=q.w*e.u.y/mesh.cells[t].area;
        }
        umean+=mesh.cells[t].area*((uc[0]-ue.x)*(uc[0]-ue.x)+(vc[0]-ue.y)*(vc[0]-ue.y));
        for(const auto f:mesh.cells[t].faces)for(const double sign:{-1.,1.}) {
            const auto& face=mesh.faces[f];const auto phi=a.basis.phi({face.centre.x+sign*face.areaVector.y/2,face.centre.y-sign*face.areaVector.x/2});
            double d=0;for(int k=0;k<3;++k)d+=phi[k]*div[k];divmax=std::max(divmax,std::abs(d));
        }
        if(fields.is_open()){
            fields<<t<<','<<mesh.cells[t].centre.x<<','<<mesh.cells[t].centre.y<<','<<mesh.cells[t].area<<','<<a.basis.h;for(const double v:a.basis.moments)fields<<','<<v;
            for(int c=0;c<2;++c)for(int k=0;k<3;++k)fields<<','<<u[id(t,c,k)];for(const auto v:p)fields<<','<<v;for(const auto v:ru)fields<<','<<v;for(const auto v:rv)fields<<','<<v;fields<<'\n';
        }
    }
    if(faces.is_open())for(int f=0;f<nf;++f) {
        const auto& face=mesh.faces[f];faces<<f<<','<<face.owner<<','<<(face.neighbour?int(*face.neighbour):-1)<<','<<face.centre.x<<','<<face.centre.y<<','<<face.areaVector.x<<','<<face.areaVector.y;
        for(int i=0;i<4;++i)faces<<','<<u[6*nc+4*f+i];faces<<'\n';
    }
    if(!output.empty()){fields.flush();faces.flush();if(!fields||!faces)throw std::runtime_error("field write failed");}
    for(int i=0;i<n;++i){double r=-rhs[i];for(int j=0;j<n;++j)r+=matrix[i*n+j]*x[j];residual=std::max(residual,std::abs(r));}
    double wallError=0,wallMax=0,wallLength=0,wallSpeed=0;Vector2D momentum{};
    std::ofstream walls;
    if(stationary && !output.empty()) {
        walls.open(output+".walls.csv");
        if(!walls)throw std::runtime_error("wall field open failed");
        walls<<std::setprecision(17)<<"face,length,nx,ny,traction_x0,traction_xs,traction_y0,traction_ys\n";
    }
    if(stationary)for(int t=0;t<nc;++t) {
        const auto& a=local[t];momentum.x+=a.force[0][0];momentum.y+=a.force[1][0];
        P3 p{};for(int k=0;k<3;++k)if(pIndex(t,k)>=0)p[k]=x[pIndex(t,k)];p[0]-=gauge;
        for(std::size_t j=0;j<mesh.cells[t].faces.size();++j) {
            const int f=int(mesh.cells[t].faces[j]);const auto& face=mesh.faces[f];
            if(face.neighbour)continue;
            const double length=std::hypot(face.areaVector.x,face.areaVector.y);
            const Vector2D normal=face.areaVector*(1/length);double reaction[2][2]{};
            // Boundary residual moments define conservative variational
            // traction, including the stabilisation contribution. A raw
            // gradient trace alone is not the method's momentum flux.
            for(int c=0;c<2;++c)for(int k=0;k<2;++k) {
                const int row=3+2*int(j)+k;
                for(int l=0;l<a.m;++l)reaction[c][k]+=a.stiffness(row,l)*u[id(t,c,l)];
                for(int l=0;l<3;++l)reaction[c][k]-=(c?a.by(l,row):a.bx(l,row))*p[l];
            }
            momentum.x+=reaction[0][0];momentum.y+=reaction[1][0];wallLength+=length;
            for(const auto& [z,w]:rule) {
                const double s=z-.5;const Point2D point{face.centre.x-face.areaVector.y*s,face.centre.y+face.areaVector.x*s};
                const auto e=exactAt(point,problem);const auto gu=e.gradient[0],gv=e.gradient[1];
                // Physical symmetric traction agrees with Laplacian traction
                // for this divergence-free, stationary no-slip exact field.
                const Vector2D exactTraction{(2*gu.x-e.p)*normal.x+(gu.y+gv.x)*normal.y,
                    (gu.y+gv.x)*normal.x+(2*gv.y-e.p)*normal.y};
                const double dx=(reaction[0][0]+12*s*reaction[0][1])/length-exactTraction.x;
                const double dy=(reaction[1][0]+12*s*reaction[1][1])/length-exactTraction.y;
                wallError+=w*length*(dx*dx+dy*dy);wallMax=std::max(wallMax,std::hypot(dx,dy));
                const double uface=u[6*nc+4*f]+s*u[6*nc+4*f+1],vface=u[6*nc+4*f+2]+s*u[6*nc+4*f+3];
                wallSpeed=std::max(wallSpeed,std::hypot(uface,vface));
            }
            if(walls.is_open())walls<<f<<','<<length<<','<<normal.x<<','<<normal.y<<','
                <<reaction[0][0]/length<<','<<12*reaction[0][1]/length<<','
                <<reaction[1][0]/length<<','<<12*reaction[1][1]/length<<'\n';
        }
    }
    if(walls.is_open()){walls.flush();if(!walls)throw std::runtime_error("wall field write failed");}

    std::cout<<std::setprecision(17)<<"{\"mesh\":\""<<name<<"\",\"n\":"<<resolution<<",\"problem\":\""<<problem<<"\",\"quadrature_order\":"<<order<<",\"cells\":"<<nc<<",\"unknowns\":"<<n
        <<",\"cell_mean_velocity_rms\":"<<std::sqrt(umean/area)<<",\"reconstructed_velocity_rms\":"<<std::sqrt(ur/area)<<",\"pressure_rms\":"<<std::sqrt(pr/area)<<",\"pressure_max_at_quadrature\":"<<pmax
        <<",\"weak_divergence_max_at_vertices\":"<<divmax<<",\"mean_flux_divergence_max\":"<<meanDiv<<",\"pressure_work\":"<<pressureWork<<",\"linear_residual_max\":"<<residual
        <<",\"area\":"<<area;
    if(stationary)std::cout<<",\"wall_traction_rms\":"<<std::sqrt(wallError/wallLength)
        <<",\"wall_traction_max_at_quadrature\":"<<wallMax<<",\"wall_length\":"<<wallLength
        <<",\"wall_velocity_max\":"<<wallSpeed<<",\"global_momentum_balance_norm\":"<<std::hypot(momentum.x,momentum.y);
    std::cout<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
#endif
