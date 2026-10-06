// Local Stokes coupling experiment, not a product acceptance test.
// A is the positive integrated outward viscous-flux operator. G is the
// integrated product pressure gradient; D H is divergence of the product
// skew-corrected velocity interpolation. C is the Rhie--Chow pressure term,
// with rAU from the unrelaxed scalar viscous diagonal, nu=1 m^2/s.
// The patch Schur operator is S=C-D H A_patch^(-1) G. Values outside the
// selected support are zero; this is not a restriction of the full inverse A.
// We examine symmetric parts and direct quadratic forms, not eigenvalues of
// the nonsymmetric evolution Jacobian. Their signs do not certify a physical
// instability, select a branch, or justify a new mesh/solution threshold.
// All actual polygons remain intact. Homogeneous velocity perturbations and
// unknown pressure traces are used on physical boundaries of this study.
// Optional minimum pressure splitting and midpoint interpolation are bounded
// research controls only. They never change the product source or binary.

#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
#define main inherited_native_face_test_main
#include "../../tests/flow_face_test.cpp"
#undef main

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using Vec=std::vector<double>;

// Matrix-free native viscous face operator.  This research probe uses the
// product gradient/transpose-stress functions and the same scalar face flux
// as momentum(). Homogeneous perturbations are supported on one mesh patch.
// No incompressibility projection: negative work is not a physical instability
// certificate.  No product equations, controls or acceptance limits are changed.
struct Action {
    const FvMesh2D& mesh;
    bool minimumPressure=false;
    std::vector<bool> fixed;
    Vec zero;
    detail::FlowGradientStencil2D stencil;
    explicit Action(const FvMesh2D& m):mesh(m),fixed(m.faces.size(),true),
        zero(m.faces.size()),stencil(detail::buildFlowGradientStencil2D(m,fixed)) {}
    Vec operator()(const Vec& value, bool reconstruct=true, bool symmetric=true) const {
        const auto n=mesh.cells.size();
        Vec u(value.begin(),value.begin()+n),v(value.begin()+n,value.end());
        auto gu=stencil.apply(u,zero),gv=stencil.apply(v,zero);
        if(!reconstruct) {std::fill(gu.begin(),gu.end(),Vector2D{});std::fill(gv.begin(),gv.end(),Vector2D{});}
        const auto correction=detail::symmetricViscousCorrection(mesh,u,v,gu,gv,
            zero,zero,fixed,fixed,fixed,fixed,1.);
        Vec result(2*n);
        for(std::size_t id=0;id<mesh.faces.size();++id) {
            const auto& f=mesh.faces[id];const auto i=f.owner;
            auto a=gu[i],b=gv[i];
            if(f.neighbour) {
                const auto j=*f.neighbour;const double w=f.neighbourWeight;
                a={a.x*(1-w)+gu[j].x*w,a.y*(1-w)+gu[j].y*w};
                b={b.x*(1-w)+gv[j].x*w,b.y*(1-w)+gv[j].y*w};
            } else if(!symmetric) {a={};b={};}
            const double uj=f.neighbour?u[*f.neighbour]:0;
            const double vj=f.neighbour?v[*f.neighbour]:0;
            const double x=-f.transmissibility*(uj-u[i])-dot(a,f.correction)+(symmetric?correction[id].x:0);
            const double y=-f.transmissibility*(vj-v[i])-dot(b,f.correction)+(symmetric?correction[id].y:0);
            result[i]+=x;result[n+i]+=y;
            if(f.neighbour) {result[*f.neighbour]-=x;result[n+*f.neighbour]-=y;}
        }
        return result;
    }
    Vec pressureForce(const Vec& p) const {
        const auto n=mesh.cells.size();
        std::vector<bool> unknown(mesh.faces.size(),false);
        const auto gradient=detail::flowGradient(mesh,p,zero,unknown,true);
        const auto force=detail::conservativePressureGradient(mesh,
            detail::pressureFaceValues(mesh,p,gradient,zero,unknown));
        Vec result(2*n);
        for(std::size_t i=0;i<n;++i){result[i]=mesh.cells[i].area*force[i].x;result[n+i]=mesh.cells[i].area*force[i].y;}
        return result;
    }
    Vec divergence(const Vec& value) const {
        const auto n=mesh.cells.size();
        Vec u(value.begin(),value.begin()+n),v(value.begin()+n,value.end()),result(n);
        const auto gu=stencil.apply(u,zero),gv=stencil.apply(v,zero);
        for(const auto& f:mesh.faces)if(f.neighbour){
            const auto i=f.owner,j=*f.neighbour;const double w=f.neighbourWeight;
            const Point2D point{(1-w)*mesh.cells[i].centre.x+w*mesh.cells[j].centre.x,
                                (1-w)*mesh.cells[i].centre.y+w*mesh.cells[j].centre.y};
            const auto skew=f.centre-point;
            const Vector2D a{(1-w)*gu[i].x+w*gu[j].x,(1-w)*gu[i].y+w*gu[j].y};
            const Vector2D b{(1-w)*gv[i].x+w*gv[j].x,(1-w)*gv[i].y+w*gv[j].y};
            const double flux=((1-w)*u[i]+w*u[j]+dot(a,skew))*f.areaVector.x+
                ((1-w)*v[i]+w*v[j]+dot(b,skew))*f.areaVector.y;
            result[i]+=flux;result[j]-=flux;
        }
        return result;
    }
    Vec pressure(const Vec& p,bool compactOnly=false) const {
        const auto n=mesh.cells.size();
        std::vector<bool> unknown(mesh.faces.size(),false);
        const auto pg=detail::flowGradient(mesh,p,zero,unknown,true);
        const auto force=detail::conservativePressureGradient(mesh,
            detail::pressureFaceValues(mesh,p,pg,zero,unknown));
        Vec diagonal(n),ra(n),result(n);
        for(const auto& f:mesh.faces){diagonal[f.owner]+=f.transmissibility;if(f.neighbour)diagonal[*f.neighbour]+=f.transmissibility;}
        for(std::size_t i=0;i<n;++i)ra[i]=mesh.cells[i].area/diagonal[i];
        for(const auto& f:mesh.faces)if(f.neighbour) {
            const auto i=f.owner,j=*f.neighbour;const double w=f.neighbourWeight;
            const double rf=(1-w)*ra[i]+w*ra[j];
            const Vector2D g{(1-w)*pg[i].x+w*pg[j].x,(1-w)*pg[i].y+w*pg[j].y};
            const Vector2D rg{(1-w)*ra[i]*force[i].x+w*ra[j]*force[j].x,
                              (1-w)*ra[i]*force[i].y+w*ra[j]*force[j].y};
            const auto delta=mesh.cells[j].centre-mesh.cells[i].centre;
            const double t=minimumPressure?dot(f.areaVector,delta)/dot(delta,delta):f.transmissibility;
            const Vector2D correction{f.areaVector.x-t*delta.x,f.areaVector.y-t*delta.y};
            const double flux=(compactOnly?0:dot(rg,f.areaVector))-
                rf*(t*(p[j]-p[i])+dot(g,correction));
            result[i]+=flux;result[j]-=flux;
        }
        return result;
    }
};

Vec solve(Vec a,Vec rhs) {
    const auto n=rhs.size();
    for(std::size_t k=0;k<n;++k){
        std::size_t pivot=k;for(std::size_t i=k+1;i<n;++i)if(std::abs(a[i*n+k])>std::abs(a[pivot*n+k]))pivot=i;
        if(std::abs(a[pivot*n+k])<std::numeric_limits<double>::min())throw std::runtime_error("Singular patch velocity operator");
        if(pivot!=k){for(std::size_t j=k;j<n;++j)std::swap(a[k*n+j],a[pivot*n+j]);std::swap(rhs[k],rhs[pivot]);}
        for(std::size_t i=k+1;i<n;++i){const double f=a[i*n+k]/a[k*n+k];for(std::size_t j=k+1;j<n;++j)a[i*n+j]-=f*a[k*n+j];rhs[i]-=f*rhs[k];}
    }
    for(std::size_t i=n;i-->0;){for(std::size_t j=i+1;j<n;++j)rhs[i]-=a[i*n+j]*rhs[j];rhs[i]/=a[i*n+i];}
    return rhs;
}

struct Eigen {double value,residual;Vec vector;};
Eigen smallest(Vec a,std::size_t n) {
    const Vec original=a;Vec v(n*n);for(std::size_t i=0;i<n;++i)v[i*n+i]=1;
    double scale=0;for(double x:a)scale=std::max(scale,std::abs(x));
    const double stop=64*std::numeric_limits<double>::epsilon()*scale;
    std::size_t step=0;
    for(;step<100*n*n;++step) {
        std::size_t p=0,q=1;double off=0;
        for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j)
            if(std::abs(a[i*n+j])>off){off=std::abs(a[i*n+j]);p=i;q=j;}
        if(off<=stop)break;
        const double theta=.5*std::atan2(2*a[p*n+q],a[q*n+q]-a[p*n+p]);
        const double c=std::cos(theta),s=std::sin(theta);
        const double pp=a[p*n+p],qq=a[q*n+q],pq=a[p*n+q];
        for(std::size_t k=0;k<n;++k)if(k!=p&&k!=q){
            const double kp=a[k*n+p],kq=a[k*n+q];
            a[k*n+p]=a[p*n+k]=c*kp-s*kq;
            a[k*n+q]=a[q*n+k]=s*kp+c*kq;
        }
        a[p*n+p]=c*c*pp-2*c*s*pq+s*s*qq;
        a[q*n+q]=s*s*pp+2*c*s*pq+c*c*qq;
        a[p*n+q]=a[q*n+p]=0;
        for(std::size_t k=0;k<n;++k){const double kp=v[k*n+p],kq=v[k*n+q];v[k*n+p]=c*kp-s*kq;v[k*n+q]=s*kp+c*kq;}
    }
    if(step==100*n*n)throw std::runtime_error("Jacobi eigensolver did not converge");
    std::size_t index=0;for(std::size_t i=1;i<n;++i)if(a[i*n+i]<a[index*n+index])index=i;
    Vec x(n);for(std::size_t i=0;i<n;++i)x[i]=v[i*n+index];
    double residual=0;
    for(std::size_t i=0;i<n;++i){double r=-a[index*n+index]*x[i];for(std::size_t j=0;j<n;++j)r+=original[i*n+j]*x[j];residual=std::max(residual,std::abs(r));}
    return {a[index*n+index],residual,x};
}

int main(int argc,char** argv) {
    if(argc<2||argc>4)throw std::runtime_error("usage: coupling-energy mesh.cm2d|square|shear [rings] [minimum|midpoint]");
    const std::string meshName=argv[1];
    const auto load=[&]{const auto input=readCm2dTopology(meshName);if(!input.valid())throw std::runtime_error(input.error);return makeFvMesh2D(input.topology);};
    auto mesh=meshName=="square"?rectangularMesh(16,16):meshName=="shear"?rectangularMesh(16,16,1.3):load();
    const bool midpoint=argc==4 && std::string(argv[3])=="midpoint";
    if(midpoint)for(auto& f:mesh.faces)if(f.neighbour)f.neighbourWeight=.5;
    Action action(mesh);
    const auto n=mesh.cells.size();const int rings=argc>=3?std::stoi(argv[2]):3;
    if(rings<1||rings>8)throw std::runtime_error("Research patch radius must be between one and eight rings");
    if(argc==4 && std::string(argv[3])!="minimum" && std::string(argv[3])!="midpoint")
        throw std::runtime_error("Unknown research variant");
    action.minimumPressure=argc==4 && std::string(argv[3])=="minimum";
    double worst=0;std::size_t seed=0;
    for(const auto& f:mesh.faces) {
        const double ratio=std::hypot(f.correction.x,f.correction.y)/std::hypot(f.areaVector.x,f.areaVector.y);
        if(ratio>worst){worst=ratio;seed=f.owner;}
    }
    std::vector<std::size_t> patch{seed};
    for(int k=0;k<rings;++k){auto next=patch;for(auto i:patch)for(auto id:mesh.cells[i].faces){const auto& f=mesh.faces[id];if(f.neighbour)next.push_back(f.owner==i?*f.neighbour:f.owner);}std::sort(next.begin(),next.end());next.erase(std::unique(next.begin(),next.end()),next.end());patch=std::move(next);}
    const auto p=patch.size(),d=2*p;Vec matrix(d*d);
    for(std::size_t j=0;j<d;++j){Vec x(2*n);x[(j/p)*n+patch[j%p]]=1;const auto y=action(x);for(std::size_t i=0;i<d;++i)matrix[i*d+j]=y[(i/p)*n+patch[i%p]];}
    Vec symmetric=matrix;double scale=0,asymmetry=0;
    for(std::size_t i=0;i<d;++i)for(std::size_t j=0;j<d;++j){symmetric[i*d+j]=.5*(matrix[i*d+j]+matrix[j*d+i]);scale=std::max(scale,std::abs(matrix[i*d+j]));asymmetry=std::max(asymmetry,std::abs(matrix[i*d+j]-matrix[j*d+i]));}
    const auto eigen=smallest(symmetric,d);Vec x(2*n);for(std::size_t j=0;j<d;++j)x[(j/p)*n+patch[j%p]]=eigen.vector[j];
    const auto work=[&](bool reconstruct,bool sym){const auto y=action(x,reconstruct,sym);return std::inner_product(x.begin(),x.end(),y.begin(),0.);};
    Vec pressureMatrix(p*p);
    for(std::size_t j=0;j<p;++j){Vec q(n);q[patch[j]]=1;const auto y=action.pressure(q);for(std::size_t i=0;i<p;++i)pressureMatrix[i*p+j]=y[patch[i]];}
    Vec pressureSymmetric=pressureMatrix;double pressureScale=0;
    for(std::size_t i=0;i<p;++i)for(std::size_t j=0;j<p;++j){pressureSymmetric[i*p+j]=.5*(pressureMatrix[i*p+j]+pressureMatrix[j*p+i]);pressureScale=std::max(pressureScale,std::abs(pressureMatrix[i*p+j]));}
    const auto pe=smallest(pressureSymmetric,p);Vec q(n);for(std::size_t i=0;i<p;++i)q[patch[i]]=pe.vector[i];
    const auto cy=action.pressure(q),ly=action.pressure(q,true);
    const double pressureWork=std::inner_product(q.begin(),q.end(),cy.begin(),0.);
    const double compactPressureWork=std::inner_product(q.begin(),q.end(),ly.begin(),0.);
    Vec schur=pressureMatrix;
    double maximumVelocitySolveResidual=0;
    for(std::size_t j=0;j<p;++j){
        Vec pressure(n);pressure[patch[j]]=1;const auto force=action.pressureForce(pressure);
        Vec localForce(d);for(std::size_t i=0;i<d;++i)localForce[i]=force[(i/p)*n+patch[i%p]];
        const auto localVelocity=solve(matrix,localForce);Vec velocity(2*n);
        for(std::size_t i=0;i<d;++i){velocity[(i/p)*n+patch[i%p]]=localVelocity[i];double residual=-localForce[i];for(std::size_t k=0;k<d;++k)residual+=matrix[i*d+k]*localVelocity[k];maximumVelocitySolveResidual=std::max(maximumVelocitySolveResidual,std::abs(residual));}
        const auto div=action.divergence(velocity);
        for(std::size_t i=0;i<p;++i)schur[i*p+j]-=div[patch[i]];
    }
    Vec symmetricSchur=schur;double schurScale=0;
    for(std::size_t i=0;i<p;++i)for(std::size_t j=0;j<p;++j){symmetricSchur[i*p+j]=.5*(schur[i*p+j]+schur[j*p+i]);schurScale=std::max(schurScale,std::abs(schur[i*p+j]));}
    const auto se=smallest(symmetricSchur,p);
    std::cout<<std::setprecision(17)
        <<"{\"cells\":"<<n<<",\"seed\":"<<seed<<",\"rings\":"<<rings<<",\"patchCells\":"<<p
        <<",\"maximumNonorthogonalRatio\":"<<worst<<",\"matrixScale\":"<<scale<<",\"asymmetry\":"<<asymmetry
        <<",\"minimumSymmetricEigenvalue\":"<<eigen.value<<",\"eigenResidual\":"<<eigen.residual
        <<",\"nativeViscousWork\":"<<work(true,true)<<",\"withoutReconstructionWork\":"<<work(false,true)
        <<",\"scalarDiffusionWork\":"<<work(true,false)
        <<",\"pressureStabilizationMinimumEigenvalue\":"<<pe.value<<",\"pressureStabilizationWork\":"<<pressureWork
        <<",\"pressureMatrixScale\":"<<pressureScale<<",\"pressureEigenResidual\":"<<pe.residual
        <<",\"compactPressureWork\":"<<compactPressureWork
        <<",\"minimumPressureSplit\":"<<(action.minimumPressure?"true":"false")
        <<",\"midpointInterpolation\":"<<(midpoint?"true":"false")
        <<",\"localSchurMinimumEigenvalue\":"<<se.value<<",\"localSchurMatrixScale\":"<<schurScale
        <<",\"localSchurEigenResidual\":"<<se.residual<<",\"velocitySolveResidual\":"<<maximumVelocitySolveResidual<<"}\n";
}
