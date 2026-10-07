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
#include "cartmesh2d/fv/detail/CompatibleFlowElement2D.hpp"
#include "../../tests/fixtures/PolygonMesh2D.hpp"

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::fv::detail::compatible;
using cartmesh2d::test::fromPolygons;

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

// Manufactured load adapter; every geometric/discrete coefficient comes
// from the same native kernel used by the higher-order coupled equations.
struct Local : P1Local {
    std::array<Vec,2> force;
    Local(const Fixture& f,int t,const std::string& problem,int order):
        P1Local(f.mesh,t,f.diameter.at(t),order),force{Vec(3),Vec(3)} {
        for(const auto& v:q){const auto phi=basis.phi(v.p);const auto ex=exactAt(v.p,problem);
            for(int j=0;j<3;++j){force[0][j]+=v.w*phi[j]*ex.f.x;force[1][j]+=v.w*phi[j]*ex.f.y;}}
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
        for(int c=0;c<2;++c)for(std::size_t i=0;i<a.m;++i) {
            const int row=map[id(t,c,i)];if(row<0)continue;
            if(i<3)rhs[row]+=a.force[c][i];
            for(std::size_t j=0;j<a.m;++j){const int raw=id(t,c,j),col=map[raw];if(col>=0)matrix[row*n+col]+=a.stiffness(i,j);else rhs[row]-=a.stiffness(i,j)*known[raw];}
        }
        for(int p=0;p<3;++p) {
            const int row=pIndex(t,p);if(row<0)continue;
            for(int c=0;c<2;++c)for(std::size_t i=0;i<a.m;++i) {
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
        for(std::size_t j=0;j<a.m;++j){uc[j]=u[id(t,0,j)];vc[j]=u[id(t,1,j)];}
        for(int k=0;k<3;++k) {
            if(pIndex(t,k)>=0)p[k]=x[pIndex(t,k)];
            for(std::size_t j=0;j<a.m;++j)div[k]+=a.bx(k,j)*uc[j]+a.by(k,j)*vc[j];
            pressureWork-=p[k]*div[k];
        }
        meanDiv=std::max(meanDiv,std::abs(div[0])/mesh.cells[t].area);div=DenseLU(a.mass.v,3).solve(div);
        for(int k=0;k<6;++k)for(std::size_t j=0;j<a.m;++j){ru[k]+=a.potential(k,j)*uc[j];rv[k]+=a.potential(k,j)*vc[j];}
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
                for(std::size_t l=0;l<a.m;++l)reaction[c][k]+=a.stiffness(row,l)*u[id(t,c,l)];
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
