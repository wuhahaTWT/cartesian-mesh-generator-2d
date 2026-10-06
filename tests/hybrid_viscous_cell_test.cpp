#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/HybridViscousCell2D.hpp"
#include <iomanip>
#include <iostream>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}require(rejected,"invalid hybrid input accepted");}
Vector2D rotate(Vector2D p,double a){return {std::cos(a)*p.x-std::sin(a)*p.y,std::sin(a)*p.x+std::cos(a)*p.y};}
FvMesh2D meshAt(double length,double scale,double angle) {
    auto mesh=fv_test::rectangle(5,4,length,true);
    for(auto& c:mesh.cells){auto p=rotate({c.centre.x,c.centre.y},angle);c.centre={scale*p.x,scale*p.y};c.area*=scale*scale;}
    for(auto& f:mesh.faces){auto p=rotate({f.centre.x,f.centre.y},angle),s=rotate(f.areaVector,angle),c=rotate(f.correction,angle);f.centre={scale*p.x,scale*p.y};f.areaVector={scale*s.x,scale*s.y};f.correction={scale*c.x,scale*c.y};}
    return mesh;
}
double affineError=0,affineRoundoffError=0,standardAffineError=0,matrixError=0,workError=0,conservationError=0,conservationRoundoffError=0,rotationDissipation=0;
void examine(const FvMesh2D& mesh,double scale,bool standardAspect) {
    const double mu=.37;
    for(const auto a:{std::array<double,4>{0,0,0,0},{0,-2,2,0},{1,0,0,1},{0,3,0,0},{.4,-.7,1.2,-.9}}) {
        auto velocity=[&](Point2D p){return Vector2D{.3+(a[0]*p.x+a[1]*p.y)/scale,-.6+(a[2]*p.x+a[3]*p.y)/scale};};
        std::vector<std::array<double,3>> incidence(mesh.faces.size());
        std::vector<std::array<double,3>> incidenceEnvelope(mesh.faces.size());
        for(std::size_t cell=0;cell<mesh.cells.size();++cell) {
            const HybridViscousCell2D op(mesh,cell,mu);const auto& c=mesh.cells[cell];
            const auto owner=velocity(c.centre);std::vector<Vector2D> trace;
            for(auto id:op.faces())trace.push_back(velocity(mesh.faces[id].centre));
            const auto out=op.evaluate(owner,trace);const double div=(a[0]+a[3])/scale;
            const double xx=mu*(2*a[0]/scale-2./3*div),yy=mu*(2*a[3]/scale-2./3*div),xy=mu*(a[1]+a[2])/scale;
            double heating=0,envelope=0;
            for(std::size_t i=0;i<op.faces().size();++i) {
                const auto id=op.faces()[i];const auto& face=mesh.faces[id];const double sign=face.owner==cell?1:-1;
                const Vector2D s{sign*face.areaVector.x,sign*face.areaVector.y},traction{xx*s.x+xy*s.y,xy*s.x+yy*s.y};
                const std::array<double,3> exact{-traction.x,-traction.y,-dot(trace[i],traction)};
                // Absolute point samples on stretched cells amplify their
                // input rounding. Keep the unscaled error and a standard-
                // aspect gate; additionally bound this amplification with the
                // actual linear traction map, not an arbitrary larger epsilon.
                const auto dimension=2*op.faces().size();const auto& matrix=op.tractionMatrix();
                std::array<double,3> sensitivity{};
                for(std::size_t j=0;j<op.faces().size();++j)for(std::size_t k=0;k<2;++k)
                    sensitivity[k]+=std::abs(matrix[(2*i+k)*dimension+2*j])*(std::abs(trace[j].x)+std::abs(owner.x))+std::abs(matrix[(2*i+k)*dimension+2*j+1])*(std::abs(trace[j].y)+std::abs(owner.y));
                sensitivity[2]=std::abs(trace[i].x)*sensitivity[0]+std::abs(trace[i].y)*sensitivity[1]+std::abs(exact[2]);
                for(std::size_t k=0;k<3;++k){const double error=std::abs(out.outwardFlux[i][k]-exact[k]);affineError=std::max(affineError,error/(1+std::abs(exact[k])));if(standardAspect)standardAffineError=std::max(standardAffineError,error/(1+std::abs(exact[k])));affineRoundoffError=std::max(affineRoundoffError,error/(1+std::abs(exact[k])+sensitivity[k]));if(face.neighbour){incidence[id][k]+=out.outwardFlux[i][k];incidenceEnvelope[id][k]+=sensitivity[k]+std::abs(exact[k]);}}
                // Cell internal-energy production equals total-energy input
                // minus the kinetic-energy input from the momentum flux.
                const auto& flux=out.outwardFlux[i];
                heating+=-flux[2]+owner.x*flux[0]+owner.y*flux[1];
                envelope+=std::abs(flux[2])+std::abs(owner.x*flux[0])+std::abs(owner.y*flux[1]);
            }
            workError=std::max(workError,std::abs(heating-out.dissipation)/(1+envelope+out.dissipation));
            require(out.dissipation>=0,"negative sum-of-squares dissipation");
            if(a==std::array<double,4>{0,-2,2,0})rotationDissipation=std::max(rotationDissipation,out.dissipation/(1+c.area/(scale*scale)));
            if(a==std::array<double,4>{0,0,0,0}){require(out.dissipation==0,"translation dissipates energy");for(const auto& f:out.outwardFlux)for(double x:f)require(x==0,"constant velocity has traction/work");}
        }
        for(std::size_t id=0;id<mesh.faces.size();++id)if(mesh.faces[id].neighbour)for(std::size_t k=0;k<3;++k){conservationError=std::max(conservationError,std::abs(incidence[id][k]));conservationRoundoffError=std::max(conservationRoundoffError,std::abs(incidence[id][k])/(1+incidenceEnvelope[id][k]));}
    }
    // Nonsmooth arbitrary traces exercise the stabilization, including modes
    // invisible to a cell gradient. Check the evaluated physical identity and
    // compare the assembled matrix to traction derivatives of every face DOF.
    for(std::size_t cell=0;cell<mesh.cells.size();++cell) {
        const HybridViscousCell2D op(mesh,cell,mu);const auto count=op.faces().size(),n=2*count;
        const Vector2D owner{.31,-.26};std::vector<Vector2D> trace(count);
        for(std::size_t i=0;i<count;++i)trace[i]={std::sin(double(3*i+cell)),std::cos(double(7*i+cell))};
        const auto out=op.evaluate(owner,trace);double heating=0,envelope=0,quadratic=0,quadraticScale=0;
        std::vector<double> delta(n);for(std::size_t i=0;i<count;++i){delta[2*i]=trace[i].x-owner.x;delta[2*i+1]=trace[i].y-owner.y;const auto& f=out.outwardFlux[i];heating+=-f[2]+owner.x*f[0]+owner.y*f[1];envelope+=std::abs(f[2])+std::abs(owner.x*f[0])+std::abs(owner.y*f[1]);}
        const auto& matrix=op.tractionMatrix();
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){require(matrix[i*n+j]==matrix[j*n+i],"local block is not symmetric");const double value=delta[i]*matrix[i*n+j]*delta[j];quadratic+=value;quadraticScale+=std::abs(value);}
        require(out.dissipation>0,"arbitrary trace mode not dissipated");
        workError=std::max({workError,std::abs(heating-out.dissipation)/(1+envelope+out.dissipation),std::abs(quadratic-out.dissipation)/(1+quadraticScale+out.dissipation)});
        for(std::size_t j=0;j<n;++j) {
            auto plus=trace,minus=trace;(j%2?plus[j/2].y:plus[j/2].x)+=.5;(j%2?minus[j/2].y:minus[j/2].x)-=.5;
            const auto p=op.evaluate(owner,plus),m=op.evaluate(owner,minus);
            for(std::size_t i=0;i<n;++i){const double evaluated=m.outwardFlux[i/2][i%2]-p.outwardFlux[i/2][i%2];double row=0;for(std::size_t k=0;k<n;++k)row+=std::abs(matrix[i*n+k]);matrixError=std::max(matrixError,std::abs(evaluated-matrix[i*n+j])/(1+row));}
        }
        auto bad=trace;bad[0].x=INFINITY;rejects([&]{(void)op.evaluate(owner,bad);});
    }
}
}
int main(){try {
    for(double length:{.03,3.,30.})for(double scale:{1e-3,1.,1e3})for(double angle:{0.,.731})examine(meshAt(length,scale,angle),scale,length==3);
    auto mesh=fv_test::rectangle(2,2,1);
    rejects([&]{HybridViscousCell2D bad(mesh,0,-1);});
    auto wrong=mesh;wrong.cells[0].area*=1.01;rejects([&]{HybridViscousCell2D bad(wrong,0,1);});
    wrong=mesh;wrong.cells[0].centre={-1,-1};rejects([&]{HybridViscousCell2D bad(wrong,0,1);});
    // Existing 4096-epsilon algebra budget; this is not a PDE-error, arbitrary
    // grid stability, trace-solver or complete compressible-flow qualification.
    const double budget=4096*std::numeric_limits<double>::epsilon();
    std::cout<<std::setprecision(17)<<"{\"affineFluxRelativeError\":"<<affineError<<",\"standardAspectAffineError\":"<<standardAffineError<<",\"affineInputRoundoffError\":"<<affineRoundoffError<<",\"tractionMatrixScaledError\":"<<matrixError<<",\"mechanicalWorkIdentityError\":"<<workError<<",\"affineSharedFaceImbalance\":"<<conservationError<<",\"affineSharedFaceRoundoffError\":"<<conservationRoundoffError<<",\"rigidRotationDissipation\":"<<rotationDissipation<<"}\n";
    require(standardAffineError<budget&&affineRoundoffError<budget,"affine hybrid traction/work incorrect");require(matrixError<budget,"matrix differs from evaluated traction");require(workError<budget,"viscous work is not compatible with dissipation");require(conservationRoundoffError<budget,"affine internal traction/work does not cancel within input-roundoff envelope");require(rotationDissipation<budget*budget,"rigid rotation produces dissipation");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
