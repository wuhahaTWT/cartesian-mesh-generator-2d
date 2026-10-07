#include "cartmesh2d/fv/CompatibleIncompressible2D.hpp"
#include "cartmesh2d/fv/detail/CompatibleFlowElement2D.hpp"
#include "cartmesh2d/fv/detail/CompatibleFlowLinear2D.hpp"
#include "fixtures/PolygonMesh2D.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sstream>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::fv::detail::compatible;
namespace {
void check(bool x,const char* message){if(!x)throw std::runtime_error(message);}
FvMesh2D grid(int n,double shear=0,double L=1) {
    std::vector<Polygon2D> cells;auto p=[&](int i,int j){return Point2D{L*(double(i)/n+shear*double(j)/n),L*double(j)/n};};
    for(int j=0;j<n;++j)for(int i=0;i<n;++i)cells.push_back({{p(i,j),p(i+1,j),p(i+1,j+1),p(i,j+1)}});
    return makeFvMesh2D(cartmesh2d::test::fromPolygons(cells));
}
FvMesh2D rotatedGrid(int n,double angle) {
    const double c=std::cos(angle),s=std::sin(angle);std::vector<Polygon2D> cells;
    auto p=[&](int i,int j){const double x=double(i)/n,y=double(j)/n;return Point2D{c*x-s*y,s*x+c*y};};
    for(int j=0;j<n;++j)for(int i=0;i<n;++i)cells.push_back({{p(i,j),p(i+1,j),p(i+1,j+1),p(i,j+1)}});
    return makeFvMesh2D(cartmesh2d::test::fromPolygons(cells));
}
Vector2D velocity(Point2D p,bool rotation){return rotation?Vector2D{-p.y,p.x}:Vector2D{p.x*p.x-p.y*p.y,-2*p.x*p.y};}
double pressure(Point2D p){return 1+2*p.x-3*p.y;}
CompatibleFlowControls2D control(const FvMesh2D& mesh,bool rotation,bool outlet,double U=1,double L=1) {
    CompatibleFlowControls2D c;c.viscosity=.1*U*L;c.referenceLength=L;c.referenceVelocity=U;
    c.equation=rotation?CompatibleEquation2D::NavierStokes:CompatibleEquation2D::Stokes;
    c.globalization=rotation?CompatibleGlobalization2D::PseudoTime:CompatibleGlobalization2D::Backtracking;
    c.acceleration=[=](Point2D p){p.x/=L;p.y/=L;return Vector2D{(2-(rotation?p.x:0))*U*U/L,(-3-(rotation?p.y:0))*U*U/L};};
    for(std::size_t f=0;f<mesh.faces.size();++f)if(!mesh.faces[f].neighbour){const auto& face=mesh.faces[f];CompatibleBoundary2D b;b.face=f;
        if(outlet&&face.areaVector.x>0&&std::abs(face.areaVector.y)<1e-14){b.kind=CompatibleBoundaryKind2D::Traction;auto n=face.areaVector;const double len=std::hypot(n.x,n.y);n.x/=len;n.y/=len;
            b.value=[=](Point2D p){p.x/=L;p.y/=L;const double pr=pressure(p);return rotation?Vector2D{-pr*n.x*U*U,-pr*n.y*U*U}:Vector2D{((.4*p.x-pr)*n.x-.4*p.y*n.y)*U*U,(-.4*p.y*n.x+(-.4*p.x-pr)*n.y)*U*U};};}
        else b.value=[=](Point2D p){p.x/=L;p.y/=L;auto u=velocity(p,rotation);return Vector2D{U*u.x,U*u.y};};
        c.boundaries.push_back(std::move(b));}
    return c;
}
double h(const FvMesh2D& mesh,std::size_t t){std::vector<Point2D> vertices;for(auto id:mesh.cells[t].faces){const auto& f=mesh.faces[id];vertices.push_back({f.centre.x-f.areaVector.y/2,f.centre.y+f.areaVector.x/2});vertices.push_back({f.centre.x+f.areaVector.y/2,f.centre.y-f.areaVector.x/2});}double d=0;for(auto a:vertices)for(auto b:vertices)d=std::max(d,std::hypot(a.x-b.x,a.y-b.y));return d;}
std::pair<double,double> errors(const FvMesh2D& mesh,const CompatibleFlowState2D& state,bool rotation,bool outlet,double U=1,double L=1) {
    double ue=0,pe=0;const auto last=mesh.cells.back().centre;const double gauge=outlet?0:pressure({last.x/L,last.y/L});
    for(std::size_t t=0;t<mesh.cells.size();++t){P1Local a(mesh,t,h(mesh,t),6);Vec local(2*a.m+3);for(std::size_t c=0;c<2;++c)for(std::size_t j=0;j<3;++j)local[c*a.m+j]=state.cells[t][3*c+j];
        for(std::size_t f=0;f<mesh.cells[t].faces.size();++f)for(std::size_t c=0;c<2;++c)for(std::size_t j=0;j<2;++j)local[c*a.m+3+2*f+j]=state.faces[mesh.cells[t].faces[f]][2*c+j];
        for(const auto& q:a.q){const auto theta=a.basis.theta(q.p);const auto phi=a.basis.phi(q.p);Vector2D got{};double pr=0;
            for(std::size_t j=0;j<a.m;++j)for(std::size_t k=0;k<6;++k){got.x+=theta[k]*a.potential(k,j)*local[j];got.y+=theta[k]*a.potential(k,j)*local[a.m+j];}
            for(std::size_t j=0;j<3;++j)pr+=phi[j]*state.cells[t][6+j];
            const auto exact=velocity({q.p.x/L,q.p.y/L},rotation);
            ue=std::max(ue,std::hypot(got.x/U-exact.x,got.y/U-exact.y));pe=std::max(pe,std::abs(pr/U/U-pressure({q.p.x/L,q.p.y/L})+gauge));}
    }
    return {ue,pe};
}
void print(const char* name,const CompatibleFlowResult2D& r,double u=std::numeric_limits<double>::quiet_NaN(),double p=std::numeric_limits<double>::quiet_NaN()){
    const auto number=[](double x){if(!std::isfinite(x))return std::string("null");std::ostringstream s;s<<std::setprecision(17)<<x;return s.str();};
    const auto m=r.iterations.empty()?CompatibleFlowMetrics2D{}:r.iterations.back().metrics.value_or(CompatibleFlowMetrics2D{});
    std::cout<<std::setprecision(17)<<"{\"case\":\""<<name<<"\",\"converged\":"<<(r.converged()?"true":"false")<<",\"iterations\":"<<r.iterations.size()<<",\"velocityMax\":"<<number(u)<<",\"pressureMax\":"<<number(p)<<",\"originalEquationsEvaluated\":"<<(!r.iterations.empty()&&r.iterations.back().metrics?"true":"false")<<",\"cellMomentum\":"<<m.cellMomentum<<",\"faceMomentum\":"<<m.faceMomentum<<",\"divergence\":"<<m.divergence<<",\"stateChange\":"<<m.stateChange<<",\"reason\":\""<<r.reason<<"\"}\n";
}
bool same(const CompatibleFlowState2D& a,const CompatibleFlowState2D& b){return a.cells==b.cells&&a.faces==b.faces;}
void oneLevelFillControl() {
    using detail::SparsePattern2D;using detail::SparseSystem2D;
    using detail::linearNorm;using detail::linearEnsure;
    using cartmesh2d::fv::detail::compatible::linear::oneLevelFillConnections;
    // Nonsymmetric three-row star becomes an exact LU with one fill level.
    // This detects accidental symmetrization as well as missing numeric fill.
    const detail::compatible::linear::Matrix a(3,{{0,0,4},{1,1,6},{2,2,5},{0,1,-1},{1,0,-3},{0,2,-2},{2,0,-1}});
    const SparsePattern2D pattern(3,oneLevelFillConnections(a,3));SparseSystem2D system(pattern);
    for(std::size_t i=0;i<a.n;++i)for(auto j=a.rows[i];j<a.rows[i+1];++j)
        if(a.columns[j]==i)system.diag[i]=a.values[j];else system.off[pattern.slot(i,a.columns[j])]=a.values[j];
    system.factorILU0();const Vec exact{.5,-1,.75},rhs=a.apply(exact);Vec answer(3);
    system.preconditionILU0(rhs,answer);
    linearEnsure(linearNorm(a.residual(rhs,answer))<=256*std::numeric_limits<double>::epsilon()*(1+linearNorm(rhs)),
        "Native one-level ILU failed known nonsymmetric inverse");
    // A generated (1,2) edge must not then create level-2 edge (2,3).
    const detail::compatible::linear::Matrix chain(4,{{0,1,1},{0,2,1},{1,3,1}});
    const auto graph=oneLevelFillConnections(chain,4);
    const std::vector<std::pair<std::size_t,std::size_t>> expected{{0,1},{0,2},{1,2},{1,3}};
    linearEnsure(graph==expected,"ILU1 symbolic fill propagated a higher level");
}

}
int main()try {
    oneLevelFillControl();
    // The ordinary product rounds 1e16+1-1e16 before b-Ax.  The strict
    // residual path must retain the unit term without changing K, b or x.
    const detail::compatible::linear::Matrix cancellation(3,{{0,0,1e16},{0,1,1},{0,2,-1e16}});
    const detail::LinearVector2D cancellationX{1,1,1},cancellationRhs{1,0,0};
    check(cancellation.apply(cancellationX)[0]!=cancellationRhs[0],"Cancellation control unexpectedly exact");
    check(cancellation.residual(cancellationRhs,cancellationX)[0]==0,"Compensated compatible residual lost assembled equation");
    // These analytic polynomial cases exercise actual coupled global assembly,
    // gauges, natural stress, Newton and inhomogeneous trace elimination. The
    // 1e-7 field allowance is a numerical regression target in U,L units,
    // larger than the 1e-9 nonlinear stop; it is not a mesh/physical error gate.
    for(bool outlet:{false,true}){auto mesh=grid(3,outlet?0:.3);auto c=control(mesh,false,outlet);if(outlet)c.pressureInverse=CompatiblePressureInverse2D::DiagonalSchur;
        const auto r=solveCompatibleIncompressible2D(mesh,c);print(outlet?"quadratic-traction":"quadratic-closed",r);check(r.converged()&&r.lastAccepted.has_value(),"Quadratic Stokes global solve failed");const auto [u,p]=errors(mesh,*r.lastAccepted,false,outlet);print("quadratic-errors",r,u,p);check(u<1e-7&&p<1e-7,"Quadratic Stokes analytic field failed");}
    // One full Krylov subspace suffices for this small analytic Stokes
    // system. Early .1 Newton forcing must not consume a linear restart.
    const auto boundedMesh=grid(3);auto bounded=control(boundedMesh,false,true);
    bounded.pressureInverse=CompatiblePressureInverse2D::DiagonalSchur;bounded.maximumLinearRestarts=1;
    const auto strictLinear=solveCompatibleIncompressible2D(boundedMesh,bounded);
    print("one-restart-linear-target",strictLinear);
    check(strictLinear.converged()&&strictLinear.lastAccepted&&strictLinear.iterations.size()==1&&strictLinear.iterations[0].linearRestarts==1,"Linear Newton forcing exhausted a sufficient Krylov subspace");
    const auto [bu,bp]=errors(boundedMesh,*strictLinear.lastAccepted,false,true);
    check(bu<1e-7&&bp<1e-7,"One-restart analytic Stokes field failed");
    // More than 32 pressure rows exercises coarsening, not just the dense
    // terminal solve. Keep the same polynomial field regression allowance.
    const auto aggregationMesh=grid(8);auto aggregationControls=control(aggregationMesh,false,true);
    aggregationControls.pressureInverse=CompatiblePressureInverse2D::DiagonalSchurAggregation;
    const auto aggregated=solveCompatibleIncompressible2D(aggregationMesh,aggregationControls);
    check(aggregated.converged()&&aggregated.lastAccepted,"Aggregation Schur global solve failed");
    const auto [au,ap]=errors(aggregationMesh,*aggregated.lastAccepted,false,true);
    print("aggregation-schur-polynomial",aggregated,au,ap);
    check(au<1e-7&&ap<1e-7,"Aggregation Schur changed the analytic Stokes field");
    auto filledControls=aggregationControls;filledControls.velocityInverse=CompatibleVelocityInverse2D::ILU1;
    const auto filled=solveCompatibleIncompressible2D(aggregationMesh,filledControls);
    check(filled.converged()&&filled.lastAccepted,"ILU1 velocity global solve failed");
    const auto [fu,fp]=errors(aggregationMesh,*filled.lastAccepted,false,true);
    print("ilu1-velocity-polynomial",filled,fu,fp);
    check(fu<1e-7&&fp<1e-7,"ILU1 velocity changed the analytic Stokes field");
    auto invalidVelocity=filledControls;invalidVelocity.velocityInverse=static_cast<CompatibleVelocityInverse2D>(99);
    bool velocityRejected=false;try{(void)solveCompatibleIncompressible2D(aggregationMesh,invalidVelocity);}catch(const std::invalid_argument&){velocityRejected=true;}
    check(velocityRejected,"Unknown compatible velocity inverse accepted");
    auto closedAggregation=control(boundedMesh,false,false);closedAggregation.pressureInverse=CompatiblePressureInverse2D::DiagonalSchurAggregation;
    bool closedRejected=false;try{(void)solveCompatibleIncompressible2D(boundedMesh,closedAggregation);}catch(const std::invalid_argument&){closedRejected=true;}
    check(closedRejected,"Aggregation silently expanded the API pressure gauge scope");
    const auto polygonPatch=[&](const std::vector<Polygon2D>& polygons){const auto mesh=makeFvMesh2D(cartmesh2d::test::fromPolygons(polygons));const auto r=solveCompatibleIncompressible2D(mesh,control(mesh,false,false));check(r.converged(),"Actual polygon global solve failed");const auto [u,p]=errors(mesh,*r.lastAccepted,false,false);print("actual-polygon-patch",r,u,p);check(u<1e-7&&p<1e-7,"Actual polygon polynomial consistency lost");};
    polygonPatch({{{{0,0},{1,0},{.75,1},{0,1}}}});
    polygonPatch({{{{0,0},{1,0},{1,1},{1,2},{0,2}}},{{{1,0},{2,0},{2,1},{1,1}}},{{{1,1},{2,1},{2,2},{1,2}}}});
    auto mesh=grid(3,.3);auto c=control(mesh,true,false);const auto nonlinear=solveCompatibleIncompressible2D(mesh,c);print("nonlinear-rotation",nonlinear);check(nonlinear.converged(),"Nonlinear rotation failed");auto [u,p]=errors(mesh,*nonlinear.lastAccepted,true,false);print("nonlinear-errors",nonlinear,u,p);check(u<1e-7&&p<1e-7,"Nonlinear analytic field failed");
    auto warmControls=c;warmControls.linearInitialGuess=CompatibleLinearInitialGuess2D::CurrentState;const auto warm=solveCompatibleIncompressible2D(mesh,warmControls);check(warm.converged()&&warm.lastAccepted,"Current-state linear initial guess failed");
    check(warm.iterations.size()>1&&warm.iterations[0].linearInitialRelativeResidual==1&&warm.iterations[1].linearInitialRelativeResidual<1,"Accepted-state initial residual was not applied after the first step");
    auto [wu,wp]=errors(mesh,*warm.lastAccepted,true,false);check(wu<1e-7&&wp<1e-7,"Accepted-state linear initial guess changed analytic accuracy");
    auto scaledMesh=grid(3,.3,2.5);auto scaledControls=control(scaledMesh,true,false,3,2.5);const auto scaled=solveCompatibleIncompressible2D(scaledMesh,scaledControls);check(scaled.converged(),"Physical scaling solve failed");auto [su,sp]=errors(scaledMesh,*scaled.lastAccepted,true,false,3,2.5);print("scaled-nonlinear",scaled,su,sp);check(su<1e-7&&sp<1e-7,"Physical unit conversion failed");
    auto one=c;one.maximumIterations=1;const auto budget=solveCompatibleIncompressible2D(mesh,one);print("nonlinear-budget",budget);check(budget.stop==CompatibleFlowStop2D::NonlinearBudget&&budget.lastAccepted&&budget.seed,"Iteration limit mislabeled");check(!same(*budget.seed,*budget.lastAccepted),"Seed substituted for accepted iteration");
    const auto resumed=solveCompatibleIncompressible2D(mesh,c,budget.lastAccepted);check(resumed.converged(),"Algebraic seed continuation failed");auto [ru,rp]=errors(mesh,*resumed.lastAccepted,true,false);check(ru<1e-7&&rp<1e-7,"Resumed algebraic solve differs");
    bool stop=false;auto cancelled=c;cancelled.stopRequested=[&]{return stop;};cancelled.iterationAccepted=[&](const auto&){stop=true;};const auto stopped=solveCompatibleIncompressible2D(mesh,cancelled);check(stopped.stop==CompatibleFlowStop2D::Cancelled&&stopped.lastAccepted&&same(*stopped.lastAccepted,*budget.lastAccepted),"Cancellation replaced last accepted field");print("cancelled-after-accepted",stopped);
    auto early=c;early.stopRequested=[]{return true;};const auto empty=solveCompatibleIncompressible2D(mesh,early);check(empty.stop==CompatibleFlowStop2D::Cancelled&&!empty.seed&&!empty.lastAccepted,"Cancelled preparation created accepted state");
    std::size_t polls=0;auto during=c;during.stopRequested=[&]{return ++polls==60;};const auto interrupted=solveCompatibleIncompressible2D(mesh,during);
    check(interrupted.stop==CompatibleFlowStop2D::Cancelled&&!interrupted.lastAccepted&&interrupted.iterations.size()==1&&interrupted.iterations[0].matrixProducts>0&&!interrupted.iterations[0].metrics,"Krylov cancellation published an unchecked field");
    print("cancelled-inside-krylov",interrupted);
    const auto oneCell=grid(1);CompatibleFlowControls2D incompatible;incompatible.equation=CompatibleEquation2D::Stokes;incompatible.globalization=CompatibleGlobalization2D::Backtracking;incompatible.maximumIterations=1;
    for(std::size_t f=0;f<oneCell.faces.size();++f){CompatibleBoundary2D b;b.face=f;if(oneCell.faces[f].areaVector.x>0)b.value=[](Point2D){return Vector2D{1,0};};incompatible.boundaries.push_back(b);}
    const auto massFailure=solveCompatibleIncompressible2D(oneCell,incompatible);print("incompatible-prescribed-flux",massFailure);
    check(!massFailure.converged()&&!massFailure.iterations.empty()&&massFailure.iterations.back().metrics&&massFailure.iterations.back().metrics->divergence>.5,"Pressure gauge omission hid a global continuity failure");
    auto linear=c;linear.maximumLinearRestarts=1;linear.krylovDirections=1;const auto failed=solveCompatibleIncompressible2D(mesh,linear);print("linear-budget",failed);check(failed.stop==CompatibleFlowStop2D::LinearBudget&&!failed.lastAccepted,"Failed linear trial replaced initial state");
    auto bad=c;bad.acceleration=[](Point2D){return Vector2D{std::numeric_limits<double>::quiet_NaN(),0};};const auto nan=solveCompatibleIncompressible2D(mesh,bad);check(nan.stop==CompatibleFlowStop2D::NumericalFailure&&!nan.lastAccepted,"NaN accepted");
    struct Sentinel{};bad=c;bad.acceleration=[](Point2D)->Vector2D{throw Sentinel{};};bool propagated=false;try{(void)solveCompatibleIncompressible2D(mesh,bad);}catch(const Sentinel&){propagated=true;}check(propagated,"User callback exception swallowed");
    bad=c;bad.boundaries.pop_back();bool rejected=false;try{(void)solveCompatibleIncompressible2D(mesh,bad);}catch(const std::invalid_argument&){rejected=true;}check(rejected,"Missing boundary silently inferred");
    // Pseudo-traction Couette: true symmetric correction nu*G^T*n at the
    // outlet is nonzero. Omitting it changes the actual velocity/pressure.
    const auto square=grid(3);CompatibleFlowControls2D couette;couette.viscosity=.1;
    for(std::size_t f=0;f<square.faces.size();++f)if(!square.faces[f].neighbour){CompatibleBoundary2D b;b.face=f;
        if(square.faces[f].areaVector.x>0&&square.faces[f].areaVector.y==0)b.kind=CompatibleBoundaryKind2D::PseudoTraction;
        else b.value=[](Point2D q){return Vector2D{q.y,0};};
        couette.boundaries.push_back(b);}
    auto shear=solveCompatibleIncompressible2D(square,couette);print("pseudo-traction-couette",shear);check(shear.converged(),"Pseudo-traction nonlinear solve failed");double error=0;
    for(std::size_t t=0;t<square.cells.size();++t){error=std::max({error,std::abs(shear.lastAccepted->cells[t][0]-square.cells[t].centre.y),std::abs(shear.lastAccepted->cells[t][3]),std::abs(shear.lastAccepted->cells[t][6])});}check(error<1e-7,"Pseudo-traction transpose-gradient correction lost");
    // A rotated straight slip wall constrains only its normal trace moments.
    // Uniform tangential flow must remain exact in physical x/y coefficients.
    const double angle=.37,ca=std::cos(angle),sa=std::sin(angle);auto slipMesh=rotatedGrid(3,angle);CompatibleFlowControls2D slip;slip.viscosity=.1;slip.equation=CompatibleEquation2D::Stokes;slip.globalization=CompatibleGlobalization2D::Backtracking;slip.pressureInverse=CompatiblePressureInverse2D::DiagonalSchur;
    for(std::size_t f=0;f<slipMesh.faces.size();++f)if(!slipMesh.faces[f].neighbour){CompatibleBoundary2D b;b.face=f;const auto& face=slipMesh.faces[f];const double ex=face.areaVector.x*ca+face.areaVector.y*sa,ey=-face.areaVector.x*sa+face.areaVector.y*ca;
        if(ex>0&&std::abs(ey)<1e-13)b.kind=CompatibleBoundaryKind2D::PseudoTraction;
        else if(std::abs(ex)<1e-13)b.kind=CompatibleBoundaryKind2D::NormalVelocity;
        else b.value=[=](Point2D){return Vector2D{ca,sa};};
        slip.boundaries.push_back(std::move(b));}
    const auto slipResult=solveCompatibleIncompressible2D(slipMesh,slip);print("rotated-symmetry-uniform",slipResult);check(slipResult.converged()&&slipResult.lastAccepted,"Rotated symmetry solve failed");double slipError=0;
    for(const auto& cell:slipResult.lastAccepted->cells)slipError=std::max({slipError,std::abs(cell[0]-ca),std::abs(cell[3]-sa),std::abs(cell[1]),std::abs(cell[2]),std::abs(cell[4]),std::abs(cell[5]),std::abs(cell[6]),std::abs(cell[7]),std::abs(cell[8])});
    for(const auto& face:slipResult.lastAccepted->faces)slipError=std::max({slipError,std::abs(face[0]-ca),std::abs(face[2]-sa),std::abs(face[1]),std::abs(face[3])});
    check(slipError<1e-10,"Rotated symmetry changed uniform tangential flow");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
