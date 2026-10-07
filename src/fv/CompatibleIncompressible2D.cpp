#include "cartmesh2d/fv/CompatibleIncompressible2D.hpp"
#include "cartmesh2d/fv/detail/CompatibleFlowTransport2D.hpp"
#include "cartmesh2d/fv/detail/CompatibleFlowLinear2D.hpp"
#include <exception>
#include <limits>
#include <numeric>
namespace cartmesh2d::fv {
namespace {
using namespace detail::compatible;
using detail::linearFinite;
using detail::linearNorm;
namespace sparse=detail::compatible::linear;
constexpr auto absent=std::numeric_limits<std::size_t>::max();
struct Cancelled {};
struct UserException {std::exception_ptr exception;};
template<class F> auto userCall(F&& f) {
    try{return f();}catch(...){throw UserException{std::current_exception()};}
}
void poll(const CompatibleFlowControls2D& c){if(c.stopRequested&&userCall([&]{return c.stopRequested();}))throw Cancelled{};}
Vector2D value(const std::function<Vector2D(Point2D)>& f,Point2D p){
    const auto v=f?userCall([&]{return f(p);}):Vector2D{};
    return {linearFinite(v.x),linearFinite(v.y)};
}
void require(bool b,const char* s){if(!b)throw std::invalid_argument(s);}
void validate(const FvMesh2D& mesh,const CompatibleFlowControls2D& c,const std::optional<CompatibleFlowState2D>& initial) {
    validateFvMesh2D(mesh);
    for(double x:{c.viscosity,c.referenceLength,c.referenceVelocity,c.equationTolerance,c.stateTolerance,c.linearTolerance})require(std::isfinite(x)&&x>0,"Compatible scales, viscosity and tolerances must be positive finite");
    require(std::isfinite(c.referenceVelocity*c.referenceVelocity)&&c.referenceVelocity*c.referenceVelocity>0&&std::isfinite(c.referenceLength*c.referenceLength)&&c.referenceLength*c.referenceLength>0,"Compatible reference scale range exceeded");
    require(c.equation==CompatibleEquation2D::Stokes||c.equation==CompatibleEquation2D::NavierStokes,"Invalid compatible equation");
    require(c.globalization==CompatibleGlobalization2D::Backtracking||c.globalization==CompatibleGlobalization2D::PseudoTime,"Invalid compatible globalization");
    require(c.pressureInverse==CompatiblePressureInverse2D::ViscousMass||c.pressureInverse==CompatiblePressureInverse2D::DiagonalSchur,"Invalid compatible pressure inverse");
    require(c.linearTolerance<=.01&&c.maximumIterations>0&&c.maximumLinearRestarts>0&&c.krylovDirections>0&&c.maximumBacktracks>0&&c.maximumBacktracks<=static_cast<std::size_t>(std::numeric_limits<double>::max_exponent-std::numeric_limits<double>::min_exponent),"Invalid compatible iteration budget");
    require(c.quadratureOrder>=4&&c.quadratureOrder<=12,"Compatible quadrature must be 4..12");
    require(std::isfinite(c.armijo)&&c.armijo>0&&c.armijo<.5,"Invalid compatible Armijo control");
    require(std::isfinite(c.initialPseudoStep)&&c.initialPseudoStep>0&&std::isfinite(c.maximumPseudoStep)&&c.maximumPseudoStep>=c.initialPseudoStep,"Invalid compatible pseudo-time interval");
    std::vector<bool> seen(mesh.faces.size());std::size_t velocityFaces=0,openFaces=0;
    for(const auto& b:c.boundaries){require(b.face<mesh.faces.size()&&!mesh.faces[b.face].neighbour&&!seen[b.face],"Compatible boundary coverage/ownership conflict");seen[b.face]=true;
        require(b.kind==CompatibleBoundaryKind2D::Velocity||b.kind==CompatibleBoundaryKind2D::Traction||b.kind==CompatibleBoundaryKind2D::PseudoTraction,"Invalid compatible boundary kind");
        (b.kind==CompatibleBoundaryKind2D::Velocity?velocityFaces:openFaces)++;}
    for(std::size_t f=0;f<mesh.faces.size();++f)require(mesh.faces[f].neighbour.has_value()||seen[f],"Compatible boundary missing");
    require(velocityFaces>0,"Compatible solver requires a velocity boundary to fix rigid motion");
    require(c.pressureInverse!=CompatiblePressureInverse2D::DiagonalSchur||openFaces>0,"Diagonal Schur currently requires a traction pressure reference");
    // A single pressure gauge cannot anchor disconnected components. Report
    // this unsupported input; do not solve an unanchored saddle system.
    std::vector<bool> visited(mesh.cells.size());std::vector<std::size_t> pending{0};visited[0]=true;
    for(std::size_t i=0;i<pending.size();++i)for(auto f:mesh.cells[pending[i]].faces){const auto& face=mesh.faces[f];if(!face.neighbour)continue;
        const auto other=face.owner==pending[i]?*face.neighbour:face.owner;if(!visited[other]){visited[other]=true;pending.push_back(other);}}
    require(pending.size()==mesh.cells.size(),"Compatible solver requires a connected fluid component");
    if(initial){require(initial->cells.size()==mesh.cells.size()&&initial->faces.size()==mesh.faces.size(),"Compatible seed dimensions differ from mesh");
        for(const auto& a:initial->cells)for(double x:a)require(std::isfinite(x),"Nonfinite compatible seed");
        for(const auto& a:initial->faces)for(double x:a)require(std::isfinite(x),"Nonfinite compatible seed");}
}
double diameter(const FvMesh2D& mesh,std::size_t t){std::vector<Point2D> v;
    for(auto id:mesh.cells[t].faces){const auto& f=mesh.faces[id];v.push_back({f.centre.x+f.areaVector.y/2,f.centre.y-f.areaVector.x/2});v.push_back({f.centre.x-f.areaVector.y/2,f.centre.y+f.areaVector.x/2});}
    double h=0;for(auto a:v)for(auto b:v)h=std::max(h,std::hypot(a.x-b.x,a.y-b.y));return h;
}
struct CellData {
    P1System base;Lift lift;
    CellData(const FvMesh2D& mesh,std::size_t t,int order):base(P1Local(mesh,t,diameter(mesh,t),order)),lift(mesh,t,base.a,order){}
};
struct Problem {
    FvMesh2D mesh;const CompatibleFlowControls2D& control;
    std::vector<CellData> cells;std::vector<bool> open;std::vector<std::size_t> boundary,map;
    Vec known,areas;std::size_t count=0;bool outlet=false;double nu;
    Problem(const FvMesh2D& original,const CompatibleFlowControls2D& c):mesh(original),control(c),open(mesh.faces.size()),boundary(mesh.faces.size(),absent),map(4*mesh.faces.size()+mesh.cells.size(),absent),known(map.size()),nu(c.viscosity/c.referenceVelocity/c.referenceLength) {
        require(std::isfinite(nu)&&nu>0,"Compatible dimensionless viscosity out of range");
        const double L=c.referenceLength,U=c.referenceVelocity;
        for(auto& x:mesh.cells){x.centre.x/=L;x.centre.y/=L;x.area/=L*L;areas.push_back(x.area);}
        for(auto& f:mesh.faces){f.centre.x/=L;f.centre.y/=L;f.areaVector.x/=L;f.areaVector.y/=L;f.correction.x/=L;f.correction.y/=L;}
        validateFvMesh2D(mesh);
        for(std::size_t i=0;i<c.boundaries.size();++i){const auto& b=c.boundaries[i];boundary[b.face]=i;open[b.face]=b.kind!=CompatibleBoundaryKind2D::Velocity;outlet=outlet||open[b.face];}
        for(std::size_t id=0;id<mesh.faces.size();++id){poll(c);const auto& f=mesh.faces[id];
            if(f.neighbour||open[id]){for(std::size_t j=0;j<4;++j)map[4*id+j]=count++;}
            else for(auto [z,w]:gauss(c.quadratureOrder)){const double s=z-.5;const auto u=value(c.boundaries[boundary[id]].value,{(f.centre.x-s*f.areaVector.y)*L,(f.centre.y+s*f.areaVector.x)*L});
                known[4*id]+=w*u.x/U;known[4*id+1]+=12*w*s*u.x/U;known[4*id+2]+=w*u.y/U;known[4*id+3]+=12*w*s*u.y/U;}}
        for(std::size_t t=0;t<mesh.cells.size()-(outlet?0:1);++t)map[4*mesh.faces.size()+t]=count++;
        cells.reserve(mesh.cells.size());
        for(std::size_t t=0;t<mesh.cells.size();++t){poll(c);cells.emplace_back(mesh,t,c.quadratureOrder);auto& e=cells.back().base;const auto& a=e.a;const auto m=a.m;
            for(std::size_t i=0;i<2*m;++i)for(std::size_t j=0;j<2*m;++j)e.matrix(i,j)*=nu;
            e.rhs=liftedBodyForce(a,cells.back().lift,c.quadratureOrder,[&](Point2D p){const auto force=value(c.acceleration,{p.x*L,p.y*L});return Vector2D{force.x/U/U*L,force.y/U/U*L};});
            for(std::size_t l=0;l<mesh.cells[t].faces.size();++l){const auto id=mesh.cells[t].faces[l];if(!open[id])continue;const auto& f=mesh.faces[id];const auto& b=c.boundaries[boundary[id]];const auto S=f.areaVector;const double length=std::hypot(S.x,S.y);
                for(auto [z,w]:gauss(c.quadratureOrder)){const double s=z-.5;const Point2D p{f.centre.x-s*S.y,f.centre.y+s*S.x};const auto traction=value(b.value,{p.x*L,p.y*L});const auto phi=a.basis.phi(p);
                    for(std::size_t k=0;k<2;++k){const double test=w*(k?s:1);const auto rx=3+2*l+k,ry=m+rx;
                        e.rhs[rx]+=test*length*traction.x/U/U;e.rhs[ry]+=test*length*traction.y/U/U;
                        if(b.kind==CompatibleBoundaryKind2D::PseudoTraction)for(std::size_t j=0;j<m;++j){double gx=0,gy=0;for(std::size_t q=0;q<3;++q){gx+=phi[q]*a.gx(q,j);gy+=phi[q]*a.gy(q,j);}
                            e.matrix(rx,j)-=test*nu*gx*S.x;e.matrix(rx,m+j)-=test*nu*gx*S.y;e.matrix(ry,j)-=test*nu*gy*S.x;e.matrix(ry,m+j)-=test*nu*gy*S.y;}
                    }
                }
            }
        }
    }
    Vec seed(const std::optional<CompatibleFlowState2D>& initial)const {
        const auto nc=mesh.cells.size();Vec state(9*nc+4*mesh.faces.size());const double U=control.referenceVelocity;
        if(initial){for(std::size_t t=0;t<nc;++t)for(std::size_t j=0;j<9;++j)state[9*t+j]=initial->cells[t][j]/(j<6?U:U*U);
            for(std::size_t f=0;f<mesh.faces.size();++f)for(std::size_t j=0;j<4;++j)state[9*nc+4*f+j]=initial->faces[f][j]/U;}
        for(std::size_t i=0;i<4*mesh.faces.size();++i)if(map[i]==absent)state[9*nc+i]=known[i];
        if(!outlet){const double p=state[9*(nc-1)+6];for(std::size_t t=0;t<nc;++t)state[9*t+6]-=p;}
        for(auto x:state)linearFinite(x);return state;
    }
    CompatibleFlowState2D physical(const Vec& state)const {
        CompatibleFlowState2D out;out.cells.resize(mesh.cells.size());out.faces.resize(mesh.faces.size());const double U=control.referenceVelocity;
        for(std::size_t t=0;t<out.cells.size();++t)for(std::size_t j=0;j<9;++j)out.cells[t][j]=linearFinite(state[9*t+j]*(j<6?U:U*U));
        for(std::size_t f=0;f<out.faces.size();++f)for(std::size_t j=0;j<4;++j)out.faces[f][j]=linearFinite(state[9*out.cells.size()+4*f+j]*U);return out;
    }
    std::size_t index(std::size_t t,std::size_t m,std::size_t j)const {
        if(j>=2*m)return 9*t+6+j-2*m;const auto c=j/m,k=j%m;
        return k<3?9*t+3*c+k:9*mesh.cells.size()+4*mesh.cells[t].faces[(k-3)/2]+2*c+(k-3)%2;
    }
    Vec local(std::size_t t,const Vec& state)const {const auto m=cells[t].base.a.m;Vec out(2*m+3);for(std::size_t j=0;j<out.size();++j)out[j]=state[index(t,m,j)];return out;}
    std::vector<std::size_t> retained(std::size_t t,const P1System& e)const {std::vector<std::size_t> ids;const auto m=e.a.m;
        for(auto j:e.outside)if(j==2*m)ids.push_back(4*mesh.faces.size()+t);else{const auto c=j/m,k=j%m;ids.push_back(4*mesh.cells[t].faces[(k-3)/2]+2*c+(k-3)%2);}return ids;}
    std::vector<bool> localOpen(std::size_t t)const {std::vector<bool> flags;for(auto f:mesh.cells[t].faces)flags.push_back(open[f]);return flags;}
    P1System equations(std::size_t t,const Vec& state)const {auto e=cells[t].base;
        if(control.equation==CompatibleEquation2D::NavierStokes){const auto C=transportMatrix(mesh,t,e.a,cells[t].lift,local(t,state),localOpen(t),control.quadratureOrder);
            for(std::size_t i=0;i<C.nr;++i)for(std::size_t j=0;j<C.nc;++j)e.matrix(i,j)+=C(i,j);}return e;}
    CompatibleFlowMetrics2D metrics(const Vec& state)const {
        CompatibleFlowMetrics2D out;Vec residual(state.size());
        for(std::size_t t=0;t<cells.size();++t){poll(control);const auto e=equations(t,state);const auto v=local(t,state);const auto m=e.a.m;
            for(std::size_t i=0;i<v.size();++i){double r=-e.rhs[i];for(std::size_t j=0;j<v.size();++j)r+=e.matrix(i,j)*v[j];residual[index(t,m,i)]+=linearFinite(r);}
            for(const auto& tri:cells[t].lift.tri)for(auto point:{tri.c,tri.a,tri.b}){const auto phi=e.a.basis.phi(point);double d=0;
                for(std::size_t j=0;j<3;++j)for(std::size_t l=0;l<m;++l)d+=phi[j]*(e.a.gx(j,l)*v[l]+e.a.gy(j,l)*v[m+l]);out.divergence=std::max(out.divergence,std::abs(linearFinite(d)));}
        }
        for(std::size_t t=0;t<cells.size();++t)for(std::size_t j=0;j<6;++j)out.cellMomentum=std::max(out.cellMomentum,std::abs(residual[9*t+j])/mesh.cells[t].area);
        for(std::size_t f=0;f<mesh.faces.size();++f)for(std::size_t j=0;j<4;++j){auto& r=residual[9*cells.size()+4*f+j];if(map[4*f+j]==absent)r=0;else out.faceMomentum=std::max(out.faceMomentum,std::abs(r)/std::hypot(mesh.faces[f].areaVector.x,mesh.faces[f].areaVector.y));}
        if(!outlet)residual[9*(cells.size()-1)+6]=0;out.residualNorm=linearNorm(residual);
        linearFinite(out.cellMomentum);linearFinite(out.faceMomentum);return out;
    }
};
struct Assembly {std::vector<P1System> local;std::vector<sparse::Entry> entries,picard;Vec rhs;};
Assembly assemble(const Problem& p,const Vec& state,double inverseStep) {
    Assembly out;out.rhs.resize(p.count);
    for(std::size_t t=0;t<p.cells.size();++t){poll(p.control);auto e=p.equations(t,state);const auto m=e.a.m;const auto old=p.local(t,state);
        if(inverseStep>0)for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j){const double mass=inverseStep*e.a.mass(i,j);e.matrix(c*m+i,c*m+j)+=mass;e.rhs[c*m+i]+=mass*old[c*m+j];}
        e.condense();const auto pk=e.condensed().first;
        if(p.control.equation==CompatibleEquation2D::NavierStokes){const auto d=transportAdvectorDerivative(p.mesh,t,e.a,p.cells[t].lift,old,p.localOpen(t),p.control.quadratureOrder);
            for(std::size_t i=0;i<2*m;++i)for(std::size_t j=0;j<2*m;++j){e.matrix(i,j)+=d(i,j);e.rhs[i]+=d(i,j)*old[j];}e.condense();}
        const auto [k,b]=e.condensed();const auto ids=p.retained(t,e);
        for(std::size_t i=0;i<ids.size();++i){const auto row=p.map[ids[i]];if(row==absent)continue;out.rhs[row]+=b[i];
            for(std::size_t j=0;j<ids.size();++j){const auto col=p.map[ids[j]];
                if(col==absent)out.rhs[row]-=k(i,j)*p.known[ids[j]];
                else {if(k(i,j)!=0)out.entries.push_back({static_cast<std::int64_t>(row),static_cast<std::int64_t>(col),k(i,j)});
                    if(pk(i,j)!=0)out.picard.push_back({static_cast<std::int64_t>(row),static_cast<std::int64_t>(col),pk(i,j)});}}
        }
        out.local.push_back(std::move(e));
    }
    for(auto x:out.rhs)linearFinite(x);return out;
}
std::optional<Vec> linearSolve(const Problem& p,const Assembly& a,CompatibleFlowIteration2D& record) {
    const double initial=linearNorm(a.rhs);Vec x(p.count),r=a.rhs;record.linearRelativeResidual=initial==0?0:1;if(initial==0)return x;
    const sparse::Matrix k(p.count,a.entries),pk(p.count,a.picard);
    sparse::Block block(k,p.areas,"ilu0",p.outlet?(p.control.pressureInverse==CompatiblePressureInverse2D::DiagonalSchur?"outlet-schur-diag":"outlet"):"gauge",p.nu,0,pk);
    for(std::size_t it=0;it<p.control.maximumLinearRestarts;++it){poll(p.control);
        const auto d=detail::newtonKrylovDirection2D([&](const Vec& v){poll(p.control);++record.matrixProducts;return k.apply(block.apply(v));},r,p.control.krylovDirections);
        if(!d)return std::nullopt;const auto step=block.apply(*d);for(std::size_t i=0;i<x.size();++i)x[i]+=step[i];
        const auto ax=k.apply(x);for(std::size_t i=0;i<x.size();++i)r[i]=a.rhs[i]-ax[i];++record.linearRestarts;
        record.linearRelativeResidual=linearNorm(r)/initial;if(record.linearRelativeResidual<=p.control.linearTolerance)return x;
    }
    return std::nullopt;
}
Vec recover(const Problem& p,const Assembly& a,const Vec& solution){Vec retained=p.known,state(9*p.cells.size()+4*p.mesh.faces.size());
    for(std::size_t i=0;i<retained.size();++i)if(p.map[i]!=absent)retained[i]=solution[p.map[i]];
    for(std::size_t t=0;t<a.local.size();++t){poll(p.control);const auto& e=a.local[t];Vec ext;for(auto id:p.retained(t,e))ext.push_back(retained[id]);const auto v=e.recover(ext);
        for(std::size_t j=0;j<v.size();++j)state[p.index(t,e.a.m,j)]=v[j];}
    for(auto x:state)linearFinite(x);return state;
}
}
CompatibleFlowResult2D solveCompatibleIncompressible2D(const FvMesh2D& mesh,const CompatibleFlowControls2D& c,const std::optional<CompatibleFlowState2D>& initial) {
    validate(mesh,c,initial);CompatibleFlowResult2D result;
    try {
        poll(c);Problem problem(mesh,c);auto current=problem.seed(initial);result.seed=problem.physical(current);auto old=problem.metrics(current);
        const bool pseudo=c.globalization==CompatibleGlobalization2D::PseudoTime;double step=c.initialPseudoStep;
        for(std::size_t it=0;it<c.maximumIterations;++it){poll(c);result.iterations.emplace_back();auto& record=result.iterations.back();record.iteration=it;record.pseudoStep=pseudo?step:0;
            const auto assembly=assemble(problem,current,pseudo?1/step:0);const auto solved=linearSolve(problem,assembly,record);
            if(!solved){result.stop=CompatibleFlowStop2D::LinearBudget;result.reason="Sparse linear solve did not meet its true residual target; no candidate recovery";return result;}
            const auto candidate=recover(problem,assembly,*solved);bool accepted=false;
            for(std::size_t bt=0;bt<(pseudo?1:c.maximumBacktracks);++bt){poll(c);const double alpha=std::ldexp(1.,-static_cast<int>(bt));auto trial=current;double change=0;
                for(std::size_t j=0;j<trial.size();++j){const double delta=alpha*(candidate[j]-current[j]);trial[j]=linearFinite(current[j]+delta);change=std::max(change,std::abs(delta));}
                auto metrics=problem.metrics(trial);metrics.stateChange=linearFinite(change);record.metrics=metrics;record.alpha=alpha;++record.trials;
                const bool equations=std::max({metrics.cellMomentum,metrics.faceMomentum,metrics.divergence})<=c.equationTolerance;
                const bool converged=equations&&((c.equation==CompatibleEquation2D::Stokes&&!pseudo)||change<=c.stateTolerance);
                const bool sufficient=metrics.residualNorm<=std::sqrt(1-2*c.armijo*alpha)*old.residualNorm;
                if(pseudo||sufficient||converged){double next=step;
                    if(pseudo){next=metrics.residualNorm==0||old.residualNorm==0||converged?c.maximumPseudoStep:std::min(c.maximumPseudoStep,step*(old.residualNorm/metrics.residualNorm));
                        if(!std::isfinite(next)||next<=0)throw std::runtime_error("Compatible pseudo-time step invalid; trial not accepted");}
                    auto physical=problem.physical(trial);poll(c);result.lastAccepted=std::move(physical);current=std::move(trial);old=metrics;step=next;accepted=true;record.accepted=true;
                    if(c.iterationAccepted)userCall([&]{c.iterationAccepted(record);});
                    if(converged){result.stop=CompatibleFlowStop2D::Converged;result.reason="Original momentum, divergence, linear and applicable state-change targets satisfied";return result;}
                    break;
                }
                result.lastRejected=problem.physical(trial);
            }
            if(!accepted){result.stop=CompatibleFlowStop2D::BacktrackingBudget;result.reason="All backtracking trials rejected; last accepted iterate retained";return result;}
        }
        result.stop=CompatibleFlowStop2D::NonlinearBudget;result.reason="Nonlinear iteration budget exhausted; no converged flow";
    }catch(const UserException& e){std::rethrow_exception(e.exception);}
    catch(const Cancelled&){result.stop=CompatibleFlowStop2D::Cancelled;result.reason="Cancelled; unfinished trial not accepted";}
    catch(const std::exception& e){result.stop=CompatibleFlowStop2D::NumericalFailure;result.reason=e.what();}
    return result;
}
}
