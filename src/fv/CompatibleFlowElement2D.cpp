#include "cartmesh2d/fv/detail/CompatibleFlowElement2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace cartmesh2d::fv::detail::compatible {
DenseLU::DenseLU(Vec matrix,std::size_t size):n(size),a(std::move(matrix)),pivots(n) {
        if(n==0||a.size()!=n*n)throw std::invalid_argument("invalid compatible local matrix dimensions");
        for(std::size_t k=0;k<n;++k) {
            std::size_t p=k;for(std::size_t i=k+1;i<n;++i) if(std::abs(a[i*n+k])>std::abs(a[p*n+k]))p=i;
            if(!std::isfinite(a[p*n+k]) || a[p*n+k]==0.)
                throw std::runtime_error("singular compatible local matrix; no solution accepted");
            pivots[k]=p;if(p!=k)for(std::size_t j=0;j<n;++j)std::swap(a[k*n+j],a[p*n+j]);
            for(std::size_t i=k+1;i<n;++i) {
                const double t=(a[i*n+k]/=a[k*n+k]);
                for(std::size_t j=k+1;j<n;++j)a[i*n+j]-=t*a[k*n+j];
            }
        }
    }
Vec DenseLU::solve(Vec x)const {
        if(x.size()!=n)throw std::invalid_argument("compatible local right-hand side dimension mismatch");
        for(std::size_t k=0;k<n;++k)if(pivots[k]!=k)std::swap(x[k],x[pivots[k]]);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<i;++j)x[i]-=a[i*n+j]*x[j];
        for(std::size_t i=n;i-- > 0;){for(std::size_t j=i+1;j<n;++j)x[i]-=a[i*n+j]*x[j];x[i]/=a[i*n+i];}
        for(double v:x)if(!std::isfinite(v))throw std::runtime_error("nonfinite compatible local solution");
        return x;
    }
std::vector<std::pair<double,double>> gauss(int n) {
    if(n<3||n>16)throw std::invalid_argument("compatible quadrature order must be 3..16");
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
std::vector<Q> cellQuadrature(const FvMesh2D& mesh,std::size_t t,int order) {
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
    P1Local::P1Local(const FvMesh2D& mesh,std::size_t t,double diameter,int order):
        basis{mesh.cells.at(t).centre,diameter,{}},
        m(3+2*mesh.cells[t].faces.size()),q(cellQuadrature(mesh,t,order)),
        mass(3,3),bx(3,m),by(3,m),gx(3,m),gy(3,m),potential(6,m),stiffness(m,m) {
        if(!std::isfinite(diameter)||!(diameter>0))throw std::invalid_argument("invalid compatible cell scale");
        const auto& cell=mesh.cells[t];
        for(const auto& v:q) {
            const auto phi=basis.phi(v.p);
            for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)mass(i,j)+=v.w*phi[i]*phi[j];
            for(std::size_t j=0;j<3;++j) {
                bx(1,j)-=v.w*phi[j]/basis.h;by(2,j)-=v.w*phi[j]/basis.h;
            }
        }
        basis.moments={mass(1,1)/cell.area,mass(1,2)/cell.area,mass(2,2)/cell.area};
        Mat k(5,5),r(5,m),project(3,6);
        for(const auto& v:q) {
            const auto phi=basis.phi(v.p);const auto theta=basis.theta(v.p);const auto grad=basis.grad(v.p);
            for(std::size_t i=0;i<5;++i)for(std::size_t j=0;j<5;++j)k(i,j)+=v.w*dot(grad[i+1],grad[j+1]);
            for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<6;++j)project(i,j)+=v.w*phi[i]*theta[j];
            for(std::size_t j=0;j<3;++j){r(2,j)-=2*v.w*phi[j]/(basis.h*basis.h);r(4,j)-=2*v.w*phi[j]/(basis.h*basis.h);}
        }
        const auto rule=gauss(order);
        for(std::size_t f=0;f<cell.faces.size();++f) {
            const auto& face=mesh.faces[cell.faces[f]];const double sign=face.owner==std::size_t(t)?1:-1;
            const auto S=face.areaVector;
            for(const auto& [z,w]:rule) {
                const double s=z-.5;const Point2D p{face.centre.x-S.y*s,face.centre.y+S.x*s};
                const auto phi=basis.phi(p);const auto grad=basis.grad(p);
                for(std::size_t j=0;j<2;++j) {
                    const double psi=j?s:1;const std::size_t col=3+2*f+j;
                    for(std::size_t i=0;i<3;++i){bx(i,col)+=w*psi*sign*S.x*phi[i];by(i,col)+=w*psi*sign*S.y*phi[i];}
                    for(std::size_t i=0;i<5;++i)r(i,col)+=w*psi*sign*dot(S,grad[i+1]);
                }
            }
        }
        const DenseLU ml(mass.v,3),kl(k.v,5);
        for(std::size_t j=0;j<m;++j) {
            Vec x(3),y(3),z(5);for(std::size_t i=0;i<3;++i){x[i]=bx(i,j);y[i]=by(i,j);}for(std::size_t i=0;i<5;++i)z[i]=r(i,j);
            x=ml.solve(x);y=ml.solve(y);z=kl.solve(z);
            for(std::size_t i=0;i<3;++i){gx(i,j)=x[i];gy(i,j)=y[i];}for(std::size_t i=0;i<5;++i)potential(i+1,j)=z[i];
        }
        potential(0,0)=1;
        for(std::size_t j=0;j<6;++j){Vec x(3);for(std::size_t i=0;i<3;++i)x[i]=project(i,j);x=ml.solve(x);for(std::size_t i=0;i<3;++i)project(i,j)=x[i];}
        for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<m;++j)for(std::size_t a=0;a<3;++a)for(std::size_t b=0;b<3;++b)
            stiffness(i,j)+=mass(a,b)*(gx(a,i)*gx(b,j)+gy(a,i)*gy(b,j));
        for(std::size_t f=0;f<cell.faces.size();++f) {
            const auto& face=mesh.faces[cell.faces[f]];const auto S=face.areaVector;const double length=std::hypot(S.x,S.y);
            Mat residual(2,m);
            for(const auto& [z,w]:rule) {
                const double s=z-.5;const Point2D p{face.centre.x-S.y*s,face.centre.y+S.x*s};
                const auto phi=basis.phi(p);auto theta=basis.theta(p);
                for(std::size_t a=0;a<6;++a)for(std::size_t b=0;b<3;++b)theta[a]-=phi[b]*project(b,a);
                for(std::size_t a=0;a<2;++a)for(std::size_t j=0;j<m;++j) {
                    double value=j<3?phi[j]:0.;for(std::size_t b=0;b<6;++b)value+=theta[b]*potential(b,j);
                    residual(a,j)-=w*(a?12*s:1)*value;
                }
            }
            residual(0,3+2*f)+=1;residual(1,4+2*f)+=1;
            for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<m;++j)
                stiffness(i,j)+=length/basis.h*(residual(0,i)*residual(0,j)+residual(1,i)*residual(1,j)/12.);
        }
    }
    Lift::Lift(const FvMesh2D& mesh,std::size_t t,const P1Local& a,int order):coefficients(8*mesh.cells[t].faces.size(),2*a.m) {
        const auto& cell=mesh.cells[t];const std::size_t nf=cell.faces.size(),nv=8*nf,ncon=7*nf-1,n=nv+ncon;
        std::vector<std::size_t> perm(nf);std::iota(perm.begin(),perm.end(),0);
        std::sort(perm.begin(),perm.end(),[&](std::size_t i,std::size_t j){auto x=mesh.faces[cell.faces[i]].centre-cell.centre,y=mesh.faces[cell.faces[j]].centre-cell.centre;return std::atan2(x.y,x.x)<std::atan2(y.y,y.x);});
        for(std::size_t i:perm){const auto& face=mesh.faces[cell.faces[i]];double s=face.owner==std::size_t(t)?1:-1;auto S=face.areaVector;
            tri.emplace_back(cell.centre,Point2D{face.centre.x+s*S.y/2,face.centre.y-s*S.x/2},Point2D{face.centre.x-s*S.y/2,face.centre.y+s*S.x/2},a.basis.h);faceLocal.push_back(i);}
        for(std::size_t i=0;i<nf;++i){auto d=tri[i].a-tri[(i+nf-1)%nf].b;double scale=1+std::hypot(cell.centre.x,cell.centre.y)+std::sqrt(cell.area);
            if(std::hypot(d.x,d.y)>256*std::numeric_limits<double>::epsilon()*scale)throw std::runtime_error("fan endpoints not closed at construction roundoff");}
        Mat k(n,n),r(n,2*a.m),constraints(ncon,nv),targets(ncon,2*a.m);std::size_t cr=0;auto g=gauss(order);
        for(std::size_t i=0;i<nf;++i){const auto& tr=tri[i];
            for(auto q:tr.quadrature(order)){auto s=tr.shape(q.p);auto phi=a.basis.phi(q.p);double w=q.w/cell.area;
                for(std::size_t j=0;j<8;++j){for(std::size_t l=0;l<8;++l)k(8*i+j,8*i+l)+=w*dot(s[j],s[l]);for(std::size_t c=0;c<2;++c)for(std::size_t l=0;l<3;++l)r(8*i+j,c*a.m+l)+=w*(c?s[j].y:s[j].x)*phi[l];}}
            constraints(cr,8*i+1)=1;constraints(cr,8*i+5)=1;constraints(cr+1,8*i+6)=3;constraints(cr+2,8*i+7)=3;
            for(std::size_t c=0;c<2;++c)for(std::size_t j=0;j<a.m;++j){const auto& d=c?a.gy:a.gx;double hdet=tr.h*tr.det;
                targets(cr,c*a.m+j)=hdet*d(0,j);targets(cr+1,c*a.m+j)=hdet*(d(1,j)*tr.j0.x+d(2,j)*tr.j0.y);targets(cr+2,c*a.m+j)=hdet*(d(1,j)*tr.j1.x+d(2,j)*tr.j1.y);}
            cr+=3;
            const auto& face=mesh.faces[cell.faces[faceLocal[i]]];double sign=face.owner==std::size_t(t)?1:-1;Vector2D S{sign*face.areaVector.x/tr.h,sign*face.areaVector.y/tr.h};
            for(auto [z,w]:g){double s=z-.5;Point2D p{face.centre.x-face.areaVector.y*s,face.centre.y+face.areaVector.x*s};auto shape=tr.shape(p);
                for(std::size_t l=0;l<2;++l){double wt=w*(l?s:1);for(std::size_t j=0;j<8;++j)constraints(cr+l,8*i+j)+=wt*dot(shape[j],S);
                    for(std::size_t c=0;c<2;++c)for(std::size_t j=0;j<2;++j)targets(cr+l,c*a.m+3+2*faceLocal[i]+j)+=wt*(j?s:1)*(c?S.y:S.x);}}
            cr+=2;
        }
        // Omit only the redundant mean normal continuity on one radial face;
        // its full trace is checked after reconstruction along with all others.
        for(std::size_t i=0;i<nf;++i){std::size_t prev=(i+nf-1)%nf;auto d=tri[i].a-cell.centre;Vector2D S{d.y/a.basis.h,-d.x/a.basis.h};
            for(std::size_t l=(i==0?1:0);l<2;++l){for(auto [z,w]:g){Point2D p{cell.centre.x+z*d.x,cell.centre.y+z*d.y};auto x=tri[i].shape(p),y=tri[prev].shape(p);double wt=w*(l?z-.5:1);
                for(std::size_t j=0;j<8;++j){constraints(cr,8*i+j)+=wt*dot(x[j],S);constraints(cr,8*prev+j)-=wt*dot(y[j],S);}}++cr;}}
        if(cr!=ncon)throw std::runtime_error("RT1 constraint size");
        for(std::size_t i=0;i<ncon;++i){double norm=0;for(std::size_t j=0;j<nv;++j)norm+=constraints(i,j)*constraints(i,j);norm=std::sqrt(norm);if(!(norm>0))throw std::runtime_error("zero RT1 constraint");
            for(std::size_t j=0;j<nv;++j)k(nv+i,j)=k(j,nv+i)=constraints(i,j)/norm;
            for(std::size_t j=0;j<2*a.m;++j)r(nv+i,j)=targets(i,j)/norm;}
        DenseLU lu(k.v,n);
        for(std::size_t j=0;j<2*a.m;++j){Vec rhs(n);for(std::size_t i=0;i<n;++i)rhs[i]=r(i,j);auto x=lu.solve(rhs);
            for(std::size_t pass=0;pass<1;++pass){auto rr=rhs;for(std::size_t i=0;i<n;++i)for(std::size_t l=0;l<n;++l)rr[i]-=k(i,l)*x[l];auto dx=lu.solve(rr);for(std::size_t i=0;i<n;++i)x[i]+=dx[i];}
            for(std::size_t i=0;i<nv;++i)coefficients(i,j)=x[i];
            for(std::size_t i=0;i<ncon;++i){double v=-targets(i,j);for(std::size_t l=0;l<nv;++l)v+=constraints(i,l)*x[l];divResidual=std::max(divResidual,std::abs(v));}
            for(std::size_t i=0;i<nf;++i){std::size_t prev=(i+nf-1)%nf;auto d=tri[i].a-cell.centre;Vector2D S{d.y/a.basis.h,-d.x/a.basis.h};
                for(double z:{0.,.5,1.}){Point2D p{cell.centre.x+z*d.x,cell.centre.y+z*d.y};auto s=tri[i].shape(p),v=tri[prev].shape(p);double jump=0;
                    for(std::size_t l=0;l<8;++l)jump+=dot(s[l],S)*x[8*i+l]-dot(v[l],S)*x[8*prev+l];
                    traceResidual=std::max(traceResidual,std::abs(jump));}}
        }
    }
    Vector2D Lift::value(std::size_t i,std::size_t j,Point2D p)const{auto s=tri[i].shape(p);Vector2D v{};for(std::size_t l=0;l<8;++l){v.x+=s[l].x*coefficients(8*i+l,j);v.y+=s[l].y*coefficients(8*i+l,j);}return v;}
P1System::P1System(P1Local local,bool symmetric):a(std::move(local)),matrix(2*a.m+3,2*a.m+3),rhs(2*a.m+3),eliminated(8,2*a.m-5){
    const std::size_t m=a.m;
        for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<m;++j){matrix(c*m+i,c*m+j)=a.stiffness(i,j);
            if(symmetric)for(std::size_t k=0;k<3;++k)for(std::size_t l=0;l<3;++l)matrix(c*m+i,c*m+j)+=a.mass(k,l)*(c?a.gy(k,i)*a.gy(l,j):a.gx(k,i)*a.gx(l,j));}
        if(symmetric)for(std::size_t i=0;i<m;++i)for(std::size_t j=0;j<m;++j)for(std::size_t k=0;k<3;++k)for(std::size_t l=0;l<3;++l){double v=a.mass(k,l)*a.gy(k,i)*a.gx(l,j);matrix(i,m+j)+=v;matrix(m+j,i)+=v;}
        for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<m;++i)for(std::size_t p=0;p<3;++p)matrix(c*m+i,2*m+p)=matrix(2*m+p,c*m+i)=-(c?a.by(p,i):a.bx(p,i));
        inside={0,1,2,m,m+1,m+2,2*m+1,2*m+2};
        for(std::size_t i=0;i<2*m+3;++i)if(std::find(inside.begin(),inside.end(),i)==inside.end())outside.push_back(i);
}
void P1System::condense(){
        // Invalidate the previous factor before any operation that can fail.
        // A failed candidate cannot recover using stale accepted-step data.
        loadInternal.clear();
        Mat ii(8,8);for(std::size_t i=0;i<8;++i)for(std::size_t j=0;j<8;++j)ii(i,j)=matrix(inside[i],inside[j]);DenseLU lu(ii.v,8);
        for(std::size_t j=0;j<outside.size();++j){Vec v(8);for(std::size_t i=0;i<8;++i)v[i]=matrix(inside[i],outside[j]);v=lu.solve(v);for(std::size_t i=0;i<8;++i)eliminated(i,j)=v[i];}
        Vec r(8);for(std::size_t i=0;i<8;++i)r[i]=rhs[inside[i]];loadInternal=lu.solve(r);
}
    std::pair<Mat,Vec> P1System::condensed()const {if(loadInternal.size()!=8)throw std::logic_error("compatible system has not been condensed");const std::size_t n=outside.size();Mat k(n,n);Vec b(n);for(std::size_t i=0;i<n;++i){b[i]=rhs[outside[i]];for(std::size_t l=0;l<8;++l)b[i]-=matrix(outside[i],inside[l])*loadInternal[l];
        for(std::size_t j=0;j<n;++j){k(i,j)=matrix(outside[i],outside[j]);for(std::size_t l=0;l<8;++l)k(i,j)-=matrix(outside[i],inside[l])*eliminated(l,j);}}return {k,b};}
    Vec P1System::recover(const Vec& ext)const{if(loadInternal.size()!=8||ext.size()!=outside.size())throw std::invalid_argument("compatible recovery requires a condensed system and matching retained state");Vec v(rhs.size());for(std::size_t i=0;i<outside.size();++i)v[outside[i]]=ext[i];for(std::size_t i=0;i<8;++i){v[inside[i]]=loadInternal[i];for(std::size_t j=0;j<outside.size();++j)v[inside[i]]-=eliminated(i,j)*ext[j];}return v;}
Vec liftedBodyForce(const P1Local& a,const Lift& lift,int order,
                   const std::function<Vector2D(Point2D)>& force){
    if(!force)throw std::invalid_argument("missing compatible body force");
    Vec rhs(2*a.m+3);
    for(std::size_t k=0;k<lift.tri.size();++k)for(auto q:lift.tri[k].quadrature(order)){
        const auto f=force(q.p);
        if(!std::isfinite(f.x)||!std::isfinite(f.y))throw std::runtime_error("nonfinite compatible body force");
        for(std::size_t j=0;j<2*a.m;++j)rhs[j]+=q.w*dot(f,lift.value(k,j,q.p));
    }
    return rhs;
}
} // namespace cartmesh2d::fv::detail::compatible
