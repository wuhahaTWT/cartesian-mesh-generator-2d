#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/HeatConduction2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include <iostream>
#include <iomanip>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
struct Evidence {double temperature,flux,gradient,residual,sensitivity;};
Evidence solve(int n,bool warped,WallGradient2D scheme) {
    const auto mesh=fv_test::rectangle(n,n,1,warped);const double pi=std::acos(-1.),k=.37;
    auto exact=[&](Point2D p){return 2+.2*std::sin(pi*p.x)*std::sinh(pi*p.y)/std::sinh(pi);};
    auto gradient=[&](Point2D p){return Vector2D{.2*pi*std::cos(pi*p.x)*std::sinh(pi*p.y)/std::sinh(pi),.2*pi*std::sin(pi*p.x)*std::cosh(pi*p.y)/std::sinh(pi)};};
    std::vector<HeatBoundary2D> bc;for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour)bc.push_back({id,HeatBoundaryKind2D::Temperature,exact(mesh.faces[id].centre),{}});
    const HeatConductionOperator2D op(mesh,bc,k,scheme);const auto count=mesh.cells.size();std::vector<double> t(count,2),capacity(count,1);
    const auto base=op.evaluate(t,capacity).cellResidual;
    std::vector<std::vector<std::pair<std::size_t,double>>> coefficients(count);std::vector<std::pair<std::size_t,std::size_t>> edges;
    // Recover the ACTUAL affine operator by independent unit perturbations.
    // A test-only steady linear solve isolates spatial error from time error.
    for(std::size_t j=0;j<count;++j){t[j]+=1;const auto next=op.evaluate(t,capacity).cellResidual;t[j]-=1;for(std::size_t i=0;i<count;++i){const double a=next[i]-base[i];if(a!=0){coefficients[i].push_back({j,a});if(i!=j)edges.push_back({i,j});}}}
    const detail::SparsePattern2D graph(count,edges);detail::SparseSystem2D system(graph);detail::LinearWorkspace2D workspace(count);std::vector<double> correction(count);
    for(std::size_t i=0;i<count;++i){system.rhs[i]=-base[i];for(const auto& [j,a]:coefficients[i]){if(i==j)system.diag[i]=a;else system.add(i,j,a);}}
    (void)system.solve(correction,workspace,1e-11,1e-11,detail::LinearSolveMethod2D::ILU0,1e-11);
    const auto previous=correction;
    (void)system.solve(correction,workspace,1e-13,1e-13,detail::LinearSolveMethod2D::ILU0,1e-13);
    Evidence evidence{};for(std::size_t i=0;i<count;++i){t[i]+=correction[i];evidence.sensitivity=std::max(evidence.sensitivity,std::abs(previous[i]-correction[i]));evidence.temperature+=mesh.cells[i].area*std::abs(t[i]-exact(mesh.cells[i].centre));}
    // Refine against the evaluated flux operator itself, so finite-difference
    // coefficient rounding cannot become the spatial verification residual.
    for(int iteration=0;iteration<3;++iteration){const auto residual=op.evaluate(t,capacity).cellResidual;for(std::size_t i=0;i<count;++i)system.rhs[i]=-residual[i];std::vector<double> delta(count);(void)system.solve(delta,workspace,1e-14,1e-14,detail::LinearSolveMethod2D::ILU0,1e-13);for(std::size_t i=0;i<count;++i)t[i]+=delta[i];}
    evidence.temperature=0;for(std::size_t i=0;i<count;++i)evidence.temperature+=mesh.cells[i].area*std::abs(t[i]-exact(mesh.cells[i].centre));
    const auto result=op.evaluate(t,capacity);double totalMagnitude=0;
    for(const auto& b:bc){const auto& f=mesh.faces[b.face];const double half=.5*std::hypot(f.areaVector.x,f.areaVector.y),x=f.centre.x,y=f.centre.y;double q;
        if(std::abs(f.areaVector.y)>0){const double integral=(std::cos(pi*(x-half))-std::cos(pi*(x+half)))/pi;q=-k*.2*pi*std::cosh(pi*y)/std::sinh(pi)*integral*std::copysign(1.,f.areaVector.y);}
        else {const double integral=(std::cosh(pi*(y+half))-std::cosh(pi*(y-half)))/pi;q=-k*.2*pi*std::cos(pi*x)/std::sinh(pi)*integral*std::copysign(1.,f.areaVector.x);}
        evidence.flux+=std::abs(result.faceHeatFlux[b.face]-q);totalMagnitude+=std::abs(q);
        evidence.gradient+=std::abs(result.faceHeatFlux[b.face]+k*dot(gradient(f.centre),f.areaVector));
    }
    evidence.flux/=totalMagnitude;evidence.gradient/=totalMagnitude;
    for(double r:result.cellResidual)evidence.residual+=std::abs(r);evidence.residual/=totalMagnitude;
    // Residual and tolerance sensitivity must be negligible compared with the
    // discretization errors being measured, not merely labelled converged.
    require(evidence.residual<1e-8*evidence.flux,"steady algebraic residual contaminates spatial error");
    require(evidence.sensitivity<.01*evidence.temperature,"linear solve tolerance contaminates temperature error");
    return evidence;
}
}
int main(){try {
    std::cout<<std::setprecision(17)<<"{\"harmonic\":[";bool first=true;
    for(bool warped:{false,true}) {
        std::array<Evidence,3> linear{},quadratic{};
        for(std::size_t level=0;level<3;++level)for(const auto scheme:{WallGradient2D::Linear,WallGradient2D::Quadratic}){
            const int n=8*(1<<level);const auto e=solve(n,warped,scheme);(scheme==WallGradient2D::Linear?linear:quadratic)[level]=e;
            if(!first)std::cout<<',';first=false;
            std::cout<<"{\"n\":"<<n<<",\"warped\":"<<(warped?"true":"false")<<",\"scheme\":\""<<(scheme==WallGradient2D::Linear?"linear":"quadratic")<<"\",\"temperatureL1\":"<<e.temperature<<",\"wallFluxRelativeL1\":"<<e.flux<<",\"wallGradientRelativeL1\":"<<e.gradient<<",\"residualRelative\":"<<e.residual<<",\"toleranceSensitivity\":"<<e.sensitivity<<'}';
        }
        // Non-polynomial harmonic heat solution: integrated wall heat L1 error
        // normalized by total absolute exact wall heat. At 32x32 aim below 0.5%,
        // and retain at least the asymptotic second-order trend (>1.7).
        require(quadratic[2].flux<.005,"fine-grid wall heat accuracy exceeds 0.5 percent");
        require(std::log2(quadratic[1].flux/quadratic[2].flux)>1.7,"wall heat refinement lost second order");
        require(quadratic[2].flux<linear[2].flux/2,"quadratic walls do not improve integrated heat accuracy");
    }
    std::cout<<"]}\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
