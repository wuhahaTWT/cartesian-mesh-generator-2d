#include "cartmesh2d/fv/CompatibleIncompressible2D.hpp"
#include "cartmesh2d/fv/detail/CompatibleFlowTransport2D.hpp"
#include "cartmesh2d/fv/detail/CompatibleFlowLinear2D.hpp"
#include "cartmesh2d/fv/detail/CompatibleFlowCheckpoint2D.hpp"
#include <bit>
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
struct IncompatibleCheckpoint:std::invalid_argument {using std::invalid_argument::invalid_argument;};
struct UserException {std::exception_ptr exception;};
template<class F> auto userCall(F&& f) {
    try{return f();}catch(...){throw UserException{std::current_exception()};}
}
void poll(const CompatibleFlowControls2D& c){if(c.stopRequested&&userCall([&]{return c.stopRequested();}))throw Cancelled{};}
Vector2D value(const std::function<Vector2D(Point2D)>& f,Point2D p){
    const auto v=f?userCall([&]{return f(p);}):Vector2D{};
    return {linearFinite(v.x),linearFinite(v.y)};
}
double scalarValue(const std::function<double(Point2D)>& f,Point2D p){return f?linearFinite(userCall([&]{return f(p);})):0;}
void require(bool b,const char* s){if(!b)throw std::invalid_argument(s);}
void validate(const FvMesh2D& mesh,const CompatibleFlowControls2D& c,const std::optional<CompatibleFlowState2D>& initial) {
    validateFvMesh2D(mesh);
    for(double x:{c.viscosity,c.referenceLength,c.referenceVelocity,c.equationTolerance,c.stateTolerance,c.linearTolerance})require(std::isfinite(x)&&x>0,"Compatible scales, viscosity and tolerances must be positive finite");
    require(std::isfinite(c.referenceVelocity*c.referenceVelocity)&&c.referenceVelocity*c.referenceVelocity>0&&std::isfinite(c.referenceLength*c.referenceLength)&&c.referenceLength*c.referenceLength>0,"Compatible reference scale range exceeded");
    require(c.equation==CompatibleEquation2D::Stokes||c.equation==CompatibleEquation2D::NavierStokes,"Invalid compatible equation");
    require(c.globalization==CompatibleGlobalization2D::Backtracking||c.globalization==CompatibleGlobalization2D::PseudoTime,"Invalid compatible globalization");
    require(c.pressureInverse==CompatiblePressureInverse2D::ViscousMass||c.pressureInverse==CompatiblePressureInverse2D::DiagonalSchur||
        c.pressureInverse==CompatiblePressureInverse2D::DiagonalSchurAggregation,"Invalid compatible pressure inverse");
    require(c.velocityInverse==CompatibleVelocityInverse2D::ILU0||c.velocityInverse==CompatibleVelocityInverse2D::ILU1,"Invalid compatible velocity inverse");
    require(c.linearInitialGuess==CompatibleLinearInitialGuess2D::Zero||c.linearInitialGuess==CompatibleLinearInitialGuess2D::CurrentState,"Invalid compatible linear initial guess");
    require(c.linearTolerance<=.01&&c.maximumIterations>0&&c.maximumLinearRestarts>0&&c.krylovDirections>0&&c.maximumBacktracks>0&&c.maximumBacktracks<=static_cast<std::size_t>(std::numeric_limits<double>::max_exponent-std::numeric_limits<double>::min_exponent),"Invalid compatible iteration budget");
    require(c.quadratureOrder>=4&&c.quadratureOrder<=12,"Compatible quadrature must be 4..12");
    require(std::isfinite(c.armijo)&&c.armijo>0&&c.armijo<.5,"Invalid compatible Armijo control");
    require(std::isfinite(c.initialPseudoStep)&&c.initialPseudoStep>0&&std::isfinite(c.maximumPseudoStep)&&c.maximumPseudoStep>=c.initialPseudoStep,"Invalid compatible pseudo-time interval");
    std::vector<bool> seen(mesh.faces.size());std::size_t velocityFaces=0,openFaces=0;
    for(const auto& b:c.boundaries){require(b.face<mesh.faces.size()&&!mesh.faces[b.face].neighbour&&!seen[b.face],"Compatible boundary coverage/ownership conflict");seen[b.face]=true;
        require(b.kind==CompatibleBoundaryKind2D::Velocity||b.kind==CompatibleBoundaryKind2D::Traction||b.kind==CompatibleBoundaryKind2D::PseudoTraction||b.kind==CompatibleBoundaryKind2D::NormalVelocity,"Invalid compatible boundary kind");
        const bool mixed=b.kind==CompatibleBoundaryKind2D::NormalVelocity;
        require(mixed?!b.value:(!b.normalVelocity&&!b.tangentialTraction),"Conflicting compatible boundary data");
        require(!b.rejectBackflow||(b.kind==CompatibleBoundaryKind2D::Traction||b.kind==CompatibleBoundaryKind2D::PseudoTraction),"Backflow rejection requires a traction boundary");
        if(b.kind==CompatibleBoundaryKind2D::Velocity)++velocityFaces;else if(!mixed)++openFaces;}
    for(std::size_t f=0;f<mesh.faces.size();++f)require(mesh.faces[f].neighbour.has_value()||seen[f],"Compatible boundary missing");
    require(velocityFaces>0,"Compatible solver requires a velocity boundary to fix rigid motion");
    require(c.pressureInverse==CompatiblePressureInverse2D::ViscousMass||openFaces>0,"Diagonal Schur currently requires a traction pressure reference");
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
struct BoundaryLoadData {std::array<double,4> prescribed{};Mat transpose{0,0};};
struct Problem {
    FvMesh2D mesh;const CompatibleFlowControls2D& control;
    std::vector<CellData> cells;std::vector<bool> open,mixed;std::vector<Vector2D> normals;std::vector<std::size_t> boundary,map;
    Vec known,areas;std::size_t count=0;bool outlet=false;double nu;
    std::vector<BoundaryLoadData> loadData;Vector2D volumeForce{};
    Problem(const FvMesh2D& original,const CompatibleFlowControls2D& c,bool prepareLoads=false):mesh(original),control(c),open(mesh.faces.size()),mixed(mesh.faces.size()),normals(mesh.faces.size()),boundary(mesh.faces.size(),absent),map(4*mesh.faces.size()+mesh.cells.size(),absent),known(map.size()),nu(c.viscosity/c.referenceVelocity/c.referenceLength) {
        require(std::isfinite(nu)&&nu>0,"Compatible dimensionless viscosity out of range");
        if(prepareLoads)loadData.resize(c.boundaries.size());
        const double L=c.referenceLength,U=c.referenceVelocity;
        for(auto& x:mesh.cells){x.centre.x/=L;x.centre.y/=L;x.area/=L*L;areas.push_back(x.area);}
        for(auto& f:mesh.faces){f.centre.x/=L;f.centre.y/=L;f.areaVector.x/=L;f.areaVector.y/=L;f.correction.x/=L;f.correction.y/=L;}
        validateFvMesh2D(mesh);
        for(std::size_t i=0;i<c.boundaries.size();++i){const auto& b=c.boundaries[i];boundary[b.face]=i;mixed[b.face]=b.kind==CompatibleBoundaryKind2D::NormalVelocity;open[b.face]=b.kind!=CompatibleBoundaryKind2D::Velocity;outlet=outlet||(open[b.face]&&!mixed[b.face]);
            const auto sf=mesh.faces[b.face].areaVector;const double length=std::hypot(sf.x,sf.y);normals[b.face]={sf.x/length,sf.y/length};}
        for(std::size_t id=0;id<mesh.faces.size();++id){poll(c);const auto& f=mesh.faces[id];
            if(mixed[id]){map[4*id+2]=count++;map[4*id+3]=count++;
                for(auto [z,w]:gauss(c.quadratureOrder)){const double s=z-.5;const double un=scalarValue(c.boundaries[boundary[id]].normalVelocity,{(f.centre.x-s*f.areaVector.y)*L,(f.centre.y+s*f.areaVector.x)*L})/U;
                    known[4*id]+=w*un;known[4*id+1]+=12*w*s*un;}}
            else if(f.neighbour||open[id]){for(std::size_t j=0;j<4;++j)map[4*id+j]=count++;}
            else for(auto [z,w]:gauss(c.quadratureOrder)){const double s=z-.5;const auto u=value(c.boundaries[boundary[id]].value,{(f.centre.x-s*f.areaVector.y)*L,(f.centre.y+s*f.areaVector.x)*L});
                known[4*id]+=w*u.x/U;known[4*id+1]+=12*w*s*u.x/U;known[4*id+2]+=w*u.y/U;known[4*id+3]+=12*w*s*u.y/U;}}
        for(std::size_t t=0;t<mesh.cells.size()-(outlet?0:1);++t)map[4*mesh.faces.size()+t]=count++;
        cells.reserve(mesh.cells.size());
        for(std::size_t t=0;t<mesh.cells.size();++t){poll(c);cells.emplace_back(mesh,t,c.quadratureOrder);auto& e=cells.back().base;const auto& a=e.a;const auto m=a.m;
            for(std::size_t i=0;i<2*m;++i)for(std::size_t j=0;j<2*m;++j)e.matrix(i,j)*=nu;
            e.rhs=liftedBodyForce(a,cells.back().lift,c.quadratureOrder,[&](Point2D p){const auto force=value(c.acceleration,{p.x*L,p.y*L});return Vector2D{force.x/U/U*L,force.y/U/U*L};});
            if(prepareLoads){
                Vector2D total{e.rhs[0],e.rhs[m]};
                for(std::size_t l=0;l<mesh.cells[t].faces.size();++l){total.x+=e.rhs[3+2*l];total.y+=e.rhs[m+3+2*l];}
                volumeForce.x+=total.x;volumeForce.y+=total.y;
            }
            for(std::size_t l=0;l<mesh.cells[t].faces.size();++l){const auto id=mesh.cells[t].faces[l];if(!open[id])continue;const auto& f=mesh.faces[id];const auto& b=c.boundaries[boundary[id]];const auto S=f.areaVector;const double length=std::hypot(S.x,S.y);
                if(prepareLoads&&b.kind==CompatibleBoundaryKind2D::PseudoTraction)loadData[boundary[id]].transpose=Mat(4,2*m);
                for(auto [z,w]:gauss(c.quadratureOrder)){const double s=z-.5;const Point2D p{f.centre.x-s*S.y,f.centre.y+s*S.x};Vector2D traction{};
                    if(mixed[id]){const auto n=normals[id];const double shear=scalarValue(b.tangentialTraction,{p.x*L,p.y*L});traction={-n.y*shear,n.x*shear};}
                    else traction=value(b.value,{p.x*L,p.y*L});
                    const auto phi=a.basis.phi(p);
                    for(std::size_t k=0;k<2;++k){const double test=w*(k?s:1);const auto rx=3+2*l+k,ry=m+rx;
                        e.rhs[rx]+=test*length*traction.x/U/U;e.rhs[ry]+=test*length*traction.y/U/U;
                        if(prepareLoads){auto& load=loadData[boundary[id]].prescribed;load[k]+=test*length*traction.x/U/U;load[2+k]+=test*length*traction.y/U/U;}
                        if(b.kind==CompatibleBoundaryKind2D::PseudoTraction)for(std::size_t j=0;j<m;++j){double gx=0,gy=0;for(std::size_t q=0;q<3;++q){gx+=phi[q]*a.gx(q,j);gy+=phi[q]*a.gy(q,j);}
                            e.matrix(rx,j)-=test*nu*gx*S.x;e.matrix(rx,m+j)-=test*nu*gx*S.y;e.matrix(ry,j)-=test*nu*gy*S.x;e.matrix(ry,m+j)-=test*nu*gy*S.y;
                            if(prepareLoads){auto& correction=loadData[boundary[id]].transpose;correction(k,j)+=test*nu*gx*S.x;correction(k,m+j)+=test*nu*gx*S.y;correction(2+k,j)+=test*nu*gy*S.x;correction(2+k,m+j)+=test*nu*gy*S.y;}}
                    }
                }
            }
        }
    }
    std::shared_ptr<const std::vector<std::uint64_t>> checkpointContext(const FvMesh2D& original)const {
        auto context=std::make_shared<std::vector<std::uint64_t>>();
        context->reserve(detail::compatibleCheckpointContextWords2D(original));
        const auto word=[&](std::uint64_t x){context->push_back(x);};
        const auto real=[&](double x){word(std::bit_cast<std::uint64_t>(x));};
        // Scheme version must change when the underlying discrete problem changes.
        word(1);word(original.cells.size());word(original.faces.size());
        real(control.viscosity);real(control.referenceVelocity);real(control.referenceLength);
        word(static_cast<std::uint64_t>(control.equation));word(static_cast<std::uint64_t>(control.quadratureOrder));
        word(static_cast<std::uint64_t>(control.globalization));
        for(const auto& cell:original.cells){real(cell.centre.x);real(cell.centre.y);real(cell.area);word(cell.faces.size());for(auto f:cell.faces)word(f);}
        for(std::size_t f=0;f<original.faces.size();++f){const auto& face=original.faces[f];
            word(face.owner);word(face.neighbour?static_cast<std::uint64_t>(*face.neighbour):std::numeric_limits<std::uint64_t>::max());
            word(static_cast<std::uint64_t>(face.patch));real(face.centre.x);real(face.centre.y);
            real(face.areaVector.x);real(face.areaVector.y);real(face.correction.x);real(face.correction.y);
            real(face.transmissibility);real(face.neighbourWeight);
            if(face.neighbour){word(std::numeric_limits<std::uint64_t>::max());word(0);}
            else {const auto& b=control.boundaries[boundary[f]];word(static_cast<std::uint64_t>(b.kind));word(b.rejectBackflow?1:0);}
        }
        for(double x:known)real(x);
        // These are the loads actually integrated by the native preparation,
        // including arbitrary callbacks. No second quadrature/PDE implementation.
        for(const auto& cell:cells)for(double x:cell.base.rhs)real(x);
        require(context->size()==detail::compatibleCheckpointContextWords2D(original),"Internal compatible checkpoint context extent");
        return context;
    }
    Vec seed(const std::optional<CompatibleFlowState2D>& initial)const {
        const auto nc=mesh.cells.size();Vec state(9*nc+4*mesh.faces.size());const double U=control.referenceVelocity;
        if(initial){for(std::size_t t=0;t<nc;++t)for(std::size_t j=0;j<9;++j)state[9*t+j]=initial->cells[t][j]/(j<6?U:U*U);
            for(std::size_t f=0;f<mesh.faces.size();++f)for(std::size_t j=0;j<4;++j)state[9*nc+4*f+j]=initial->faces[f][j]/U;}
        for(std::size_t id=0;id<mesh.faces.size();++id){
            if(mixed[id]){const auto n=normals[id];for(std::size_t j=0;j<2;++j){const auto x=9*nc+4*id+j,y=x+2;const double tangent=-n.y*state[x]+n.x*state[y];
                state[x]=n.x*known[4*id+j]-n.y*tangent;state[y]=n.y*known[4*id+j]+n.x*tangent;}}
            else for(std::size_t j=0;j<4;++j)if(map[4*id+j]==absent)state[9*nc+4*id+j]=known[4*id+j];}
        if(!outlet){const double p=state[9*(nc-1)+6];for(std::size_t t=0;t<nc;++t)state[9*t+6]-=p;}
        for(auto x:state)linearFinite(x);
        return state;
    }
    CompatibleFlowState2D physical(const Vec& state)const {
        CompatibleFlowState2D out;out.cells.resize(mesh.cells.size());out.faces.resize(mesh.faces.size());const double U=control.referenceVelocity;
        for(std::size_t t=0;t<out.cells.size();++t)for(std::size_t j=0;j<9;++j)out.cells[t][j]=linearFinite(state[9*t+j]*(j<6?U:U*U));
        for(std::size_t f=0;f<out.faces.size();++f)for(std::size_t j=0;j<4;++j)out.faces[f][j]=linearFinite(state[9*out.cells.size()+4*f+j]*U);
        return out;
    }
    std::size_t index(std::size_t t,std::size_t m,std::size_t j)const {
        if(j>=2*m)return 9*t+6+j-2*m;
        const auto c=j/m,k=j%m;
        return k<3?9*t+3*c+k:9*mesh.cells.size()+4*mesh.cells[t].faces[(k-3)/2]+2*c+(k-3)%2;
    }
    Vec local(std::size_t t,const Vec& state)const {const auto m=cells[t].base.a.m;Vec out(2*m+3);for(std::size_t j=0;j<out.size();++j)out[j]=state[index(t,m,j)];return out;}
    std::vector<std::size_t> retained(std::size_t t,const P1System& e)const {std::vector<std::size_t> ids;const auto m=e.a.m;
        for(auto j:e.outside)if(j==2*m)ids.push_back(4*mesh.faces.size()+t);else{const auto c=j/m,k=j%m;ids.push_back(4*mesh.cells[t].faces[(k-3)/2]+2*c+(k-3)%2);}
        return ids;}
    Vec reduced(const Vec& state)const {
        Vec out(count);const auto offset=9*cells.size();
        for(std::size_t f=0;f<mesh.faces.size();++f){
            if(mixed[f]){const auto n=normals[f];for(std::size_t j=0;j<2;++j){
                const auto id=map[4*f+2+j];require(id!=absent,"Internal compatible mixed map");
                out[id]=-n.y*state[offset+4*f+j]+n.x*state[offset+4*f+2+j];
            }}else for(std::size_t j=0;j<4;++j)if(map[4*f+j]!=absent)out[map[4*f+j]]=state[offset+4*f+j];
        }
        for(std::size_t t=0;t<cells.size();++t)if(map[4*mesh.faces.size()+t]!=absent)out[map[4*mesh.faces.size()+t]]=state[9*t+6];
        for(auto x:out)linearFinite(x);
        return out;
    }
    // Orthogonal face coordinates affect both test and trial functions:
    // K_frame=W^T K_cart W, f_frame=W^T f_cart. Stored fields stay Cartesian.
    std::array<std::pair<std::size_t,double>,2> project(std::size_t raw)const {
        if(raw>=4*mesh.faces.size()||!mixed[raw/4])return {{{raw,1},{raw,0}}};
        const auto f=raw/4,c=(raw%4)/2,j=raw%2;const auto n=normals[f];
        return {{{4*f+j,c?n.y:n.x},{4*f+2+j,c?n.x:-n.y}}};
    }
    std::vector<bool> localOpen(std::size_t t)const {std::vector<bool> flags;for(auto f:mesh.cells[t].faces)flags.push_back(open[f]);return flags;}
    P1System equations(std::size_t t,const Vec& state)const {auto e=cells[t].base;
        if(control.equation==CompatibleEquation2D::NavierStokes){const auto C=transportMatrix(mesh,t,e.a,cells[t].lift,local(t,state),localOpen(t),control.quadratureOrder);
            for(std::size_t i=0;i<C.nr;++i)for(std::size_t j=0;j<C.nc;++j)e.matrix(i,j)+=C(i,j);}return e;}
    std::optional<std::string> boundaryFailure(const Vec& state)const {
        for(const auto& b:control.boundaries)if(b.rejectBackflow){const auto& face=mesh.faces[b.face];const auto offset=9*cells.size()+4*b.face;const double length=std::hypot(face.areaVector.x,face.areaVector.y);
            for(double s:{-.5,.5}){const double normal=((state[offset]+s*state[offset+1])*face.areaVector.x+(state[offset+2]+s*state[offset+3])*face.areaVector.y)/length;
                if(normal < -1e-12)return "Backflow at ordinary compatible pressure outlet face "+std::to_string(b.face)+"; trial not accepted";}}
        return std::nullopt;
    }
    CompatibleFlowLoads2D loads(const Vec& state,Point2D origin)const {
        CompatibleFlowLoads2D out;out.momentOrigin=origin;out.absolutePressureReference=outlet;
        const double U=control.referenceVelocity,L=control.referenceLength,forceScale=U*U*L;
        require(std::isfinite(forceScale)&&forceScale>0,"Compatible load scale out of range");
        const auto scale=[&](Vector2D v){return Vector2D{linearFinite(v.x*forceScale),linearFinite(v.y*forceScale)};};
        out.bodyForce=scale(volumeForce);
        for(std::size_t t=0;t<cells.size();++t){poll(control);
            // Context and volume forcing cover every cell; only physical
            // boundary owners need an additional transport/traction assembly.
            if(std::none_of(mesh.cells[t].faces.begin(),mesh.cells[t].faces.end(),[&](auto f){return !mesh.faces[f].neighbour;}))continue;
            const auto e=equations(t,state);const auto v=local(t,state);const auto m=e.a.m;
            for(std::size_t l=0;l<mesh.cells[t].faces.size();++l){const auto id=mesh.cells[t].faces[l];const auto& face=mesh.faces[id];if(face.neighbour)continue;
                CompatibleBoundaryLoad2D load;load.face=id;const auto& data=loadData[boundary[id]];
                for(std::size_t k=0;k<2;++k)for(std::size_t c=0;c<2;++c){const auto row=c*m+3+2*l+k;double value=-e.rhs[row];
                    for(std::size_t j=0;j<v.size();++j)value+=e.matrix(row,j)*v[j];
                    value+=data.prescribed[2*c+k];
                    for(std::size_t j=0;j<data.transpose.nc;++j)value+=data.transpose(2*c+k,j)*v[j];
                    (c?load.tractionMoments[k].y:load.tractionMoments[k].x)=linearFinite(value);
                }
                // The conservative transport form omits physical advective
                // boundary flux only on prescribed full-velocity faces. Put
                // that flux back before reporting a physical traction.
                const auto S=face.areaVector;
                for(auto [z,w]:gauss(control.quadratureOrder)){const double s=z-.5;const Vector2D u{v[3+2*l]+s*v[4+2*l],v[m+3+2*l]+s*v[m+4+2*l]};
                    const double volume=dot(u,S),flux=control.equation==CompatibleEquation2D::NavierStokes?volume:0;
                    load.volumeFlux+=w*volume;load.momentumFlux.x+=w*flux*u.x;load.momentumFlux.y+=w*flux*u.y;
                    const auto phi=e.a.basis.phi({face.centre.x-s*S.y,face.centre.y+s*S.x});double p=0;
                    for(std::size_t j=0;j<3;++j)p+=phi[j]*v[2*m+j];
                    for(std::size_t k=0;k<2;++k){const double test=w*(k?s:1);
                        load.pressureMoments[k].x-=test*p*S.x;load.pressureMoments[k].y-=test*p*S.y;
                        if(!open[id]){load.tractionMoments[k].x+=test*flux*u.x;load.tractionMoments[k].y+=test*flux*u.y;}
                    }
                }
                for(auto& moment:load.tractionMoments)moment=scale(moment);
                for(auto& moment:load.pressureMoments)moment=scale(moment);
                load.momentumFlux=scale(load.momentumFlux);load.volumeFlux=linearFinite(load.volumeFlux*U*L);
                const Vector2D arm{face.centre.x*L-origin.x,face.centre.y*L-origin.y},tangent{-S.y*L,S.x*L};
                const auto mean=load.tractionMoments[0],first=load.tractionMoments[1];
                load.torqueOnFluid=linearFinite(arm.x*mean.y-arm.y*mean.x+tangent.x*first.y-tangent.y*first.x);
                out.boundaryTraction.x+=mean.x;out.boundaryTraction.y+=mean.y;out.boundaryTorqueOnFluid+=load.torqueOnFluid;
                out.boundaryMomentumFlux.x+=load.momentumFlux.x;out.boundaryMomentumFlux.y+=load.momentumFlux.y;out.boundaryVolumeFlux+=load.volumeFlux;
                out.boundaries.push_back(load);
            }
        }
        std::sort(out.boundaries.begin(),out.boundaries.end(),[](const auto& a,const auto& b){return a.face<b.face;});
        for(double x:{out.boundaryTraction.x,out.boundaryTraction.y,out.boundaryMomentumFlux.x,out.boundaryMomentumFlux.y,out.boundaryVolumeFlux,out.boundaryTorqueOnFluid})linearFinite(x);
        out.momentumImbalance={linearFinite(out.boundaryTraction.x+out.bodyForce.x-out.boundaryMomentumFlux.x),linearFinite(out.boundaryTraction.y+out.bodyForce.y-out.boundaryMomentumFlux.y)};
        return out;
    }
    CompatibleFlowMetrics2D metrics(const Vec& state)const {
        CompatibleFlowMetrics2D out;Vec residual(state.size());
        for(std::size_t t=0;t<cells.size();++t){poll(control);const auto e=equations(t,state);const auto v=local(t,state);const auto m=e.a.m;
            for(std::size_t i=0;i<v.size();++i){double r=-e.rhs[i];for(std::size_t j=0;j<v.size();++j)r+=e.matrix(i,j)*v[j];residual[index(t,m,i)]+=linearFinite(r);}
            for(const auto& tri:cells[t].lift.tri)for(auto point:{tri.c,tri.a,tri.b}){const auto phi=e.a.basis.phi(point);double d=0;
                for(std::size_t j=0;j<3;++j)for(std::size_t l=0;l<m;++l)d+=phi[j]*(e.a.gx(j,l)*v[l]+e.a.gy(j,l)*v[m+l]);
                out.divergence=std::max(out.divergence,std::abs(linearFinite(d)));}
        }
        for(std::size_t t=0;t<cells.size();++t)for(std::size_t j=0;j<6;++j)out.cellMomentum=std::max(out.cellMomentum,std::abs(residual[9*t+j])/mesh.cells[t].area);
        for(std::size_t f=0;f<mesh.faces.size();++f){
            if(mixed[f]){const auto n=normals[f];for(std::size_t j=0;j<2;++j){auto& rx=residual[9*cells.size()+4*f+j];auto& ry=residual[9*cells.size()+4*f+2+j];const double tangent=-n.y*rx+n.x*ry;rx=n.x*rx+n.y*ry;ry=tangent;}}
            for(std::size_t j=0;j<4;++j){auto& q=residual[9*cells.size()+4*f+j];if(map[4*f+j]==absent)q=0;else out.faceMomentum=std::max(out.faceMomentum,std::abs(q)/std::hypot(mesh.faces[f].areaVector.x,mesh.faces[f].areaVector.y));}}
        if(!outlet)residual[9*(cells.size()-1)+6]=0;
        out.residualNorm=linearNorm(residual);
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
            for(std::size_t i=0;i<2*m;++i)for(std::size_t j=0;j<2*m;++j){e.matrix(i,j)+=d(i,j);e.rhs[i]+=d(i,j)*old[j];}
            e.condense();}
        const auto [k,b]=e.condensed();const auto ids=p.retained(t,e);
        for(std::size_t i=0;i<ids.size();++i)for(const auto& [rawRow,rowWeight]:p.project(ids[i])){
            if(rowWeight==0)continue;
            const auto row=p.map[rawRow];if(row==absent)continue;out.rhs[row]+=rowWeight*b[i];
            for(std::size_t j=0;j<ids.size();++j)for(const auto& [rawCol,colWeight]:p.project(ids[j])){
                if(colWeight==0)continue;
                const auto col=p.map[rawCol];const double value=rowWeight*k(i,j)*colWeight,preconditioner=rowWeight*pk(i,j)*colWeight;
                if(col==absent)out.rhs[row]-=value*p.known[rawCol];
                else {if(value!=0)out.entries.push_back({static_cast<std::int64_t>(row),static_cast<std::int64_t>(col),value});
                    if(preconditioner!=0)out.picard.push_back({static_cast<std::int64_t>(row),static_cast<std::int64_t>(col),preconditioner});}}
        }
        out.local.push_back(std::move(e));
    }
    for(auto x:out.rhs)linearFinite(x);
    return out;
}
std::optional<Vec> linearSolve(const Problem& p,const Assembly& a,CompatibleFlowIteration2D& record,const Vec* guess) {
    const double initial=linearNorm(a.rhs);Vec x=guess?*guess:Vec(p.count),r=a.rhs;if(initial==0)return Vec(p.count);
    const sparse::Matrix k(p.count,a.entries),pk(p.count,a.picard);
    if(guess)r=k.residual(a.rhs,x);
    record.linearInitialRelativeResidual=record.linearRelativeResidual=linearNorm(r)/initial;
    if(record.linearRelativeResidual<=p.control.linearTolerance)return x;
    const char* pressure="gauge";
    if(p.outlet) {
        pressure=p.control.pressureInverse==CompatiblePressureInverse2D::DiagonalSchurAggregation?"outlet-schur-aggregation":
            p.control.pressureInverse==CompatiblePressureInverse2D::DiagonalSchur?"outlet-schur-diag":"outlet";
    }
    sparse::Block block(k,p.areas,p.control.velocityInverse==CompatibleVelocityInverse2D::ILU1?"ilu1":"ilu0",pressure,p.nu,0,pk);
    for(std::size_t it=0;it<p.control.maximumLinearRestarts;++it){poll(p.control);
        // This is a strict linear solve, not an inexact Newton direction.
        // Keep the available subspace until the original relative target or
        // direction cap, then always recompute b-K*x below before accepting.
        const double remaining=p.control.linearTolerance/record.linearRelativeResidual;
        const auto d=detail::newtonKrylovDirection2D([&](const Vec& v){poll(p.control);++record.matrixProducts;return k.apply(block.apply(v));},r,p.control.krylovDirections,remaining);
        if(!d)return std::nullopt;
        const auto step=block.apply(*d);for(std::size_t i=0;i<x.size();++i)x[i]+=step[i];
        r=k.residual(a.rhs,x);++record.linearRestarts;
        record.linearRelativeResidual=linearNorm(r)/initial;if(record.linearRelativeResidual<=p.control.linearTolerance)return x;
    }
    return std::nullopt;
}
Vec recover(const Problem& p,const Assembly& a,const Vec& solution){Vec retained=p.known,state(9*p.cells.size()+4*p.mesh.faces.size());
    for(std::size_t i=0;i<retained.size();++i)if(p.map[i]!=absent)retained[i]=solution[p.map[i]];
    for(std::size_t t=0;t<a.local.size();++t){poll(p.control);const auto& e=a.local[t];Vec ext;for(auto id:p.retained(t,e)){double value=0;for(const auto& [raw,weight]:p.project(id))if(weight!=0)value+=weight*retained[raw];ext.push_back(value);}const auto v=e.recover(ext);
        for(std::size_t j=0;j<v.size();++j)state[p.index(t,e.a.m,j)]=v[j];}
    for(auto x:state)linearFinite(x);
    return state;
}
CompatibleFlowResult2D run(const FvMesh2D& mesh,const CompatibleFlowControls2D& c,const std::optional<CompatibleFlowState2D>& initial,const CompatibleFlowCheckpoint2D* restart) {
    validate(mesh,c,initial);CompatibleFlowResult2D result;
    using CheckpointAccess=detail::CompatibleCheckpointAccess2D;
    try {
        poll(c);Problem problem(mesh,c);const auto context=problem.checkpointContext(mesh);
        Vec current;double step=c.initialPseudoStep;std::size_t acceptedCount=0;
        if(restart){const auto& saved=CheckpointAccess::get(*restart);
            if(!saved.context||*saved.context!=*context||saved.normalizedState.size()!=9*mesh.cells.size()+4*mesh.faces.size())
                throw IncompatibleCheckpoint("Compatible checkpoint mesh, scales, equation, boundary or native load differs");
            if(saved.nextPseudoStep>c.maximumPseudoStep)throw IncompatibleCheckpoint("Checkpoint next pseudo-step exceeds the requested maximum");
            if(saved.acceptedIterations>std::numeric_limits<std::size_t>::max()-c.maximumIterations)throw IncompatibleCheckpoint("Compatible checkpoint iteration range exceeded");
            current=saved.normalizedState;step=saved.nextPseudoStep;acceptedCount=saved.acceptedIterations;
            result.resumed=true;result.acceptedIterationsBefore=acceptedCount;result.lastAccepted=problem.physical(current);result.checkpoint=*restart;
        }else {current=problem.seed(initial);result.seed=problem.physical(current);}
        auto old=problem.metrics(current);
        const bool pseudo=c.globalization==CompatibleGlobalization2D::PseudoTime;
        for(std::size_t it=0;it<c.maximumIterations;++it){poll(c);result.iterations.emplace_back();auto& record=result.iterations.back();record.iteration=result.acceptedIterationsBefore+it;record.pseudoStep=pseudo?step:0;
            const auto assembly=assemble(problem,current,pseudo?1/step:0);const auto guess=c.linearInitialGuess==CompatibleLinearInitialGuess2D::CurrentState?std::optional<Vec>(problem.reduced(current)):std::nullopt;
            const auto solved=linearSolve(problem,assembly,record,guess?&*guess:nullptr);
            if(!solved){result.stop=CompatibleFlowStop2D::LinearBudget;result.reason="Sparse linear solve did not meet its true residual target; no candidate recovery";return result;}
            const auto candidate=recover(problem,assembly,*solved);bool accepted=false;
            for(std::size_t bt=0;bt<(pseudo?1:c.maximumBacktracks);++bt){poll(c);const double alpha=std::ldexp(1.,-static_cast<int>(bt));auto trial=current;double change=0;
                for(std::size_t j=0;j<trial.size();++j){const double delta=alpha*(candidate[j]-current[j]);trial[j]=linearFinite(current[j]+delta);change=std::max(change,std::abs(delta));}
                auto metrics=problem.metrics(trial);metrics.stateChange=linearFinite(change);record.metrics=metrics;record.alpha=alpha;++record.trials;
                if(const auto failure=problem.boundaryFailure(trial)){result.lastRejected=problem.physical(trial);
                    if(pseudo||bt+1==c.maximumBacktracks){result.stop=CompatibleFlowStop2D::BoundaryFailure;result.reason=*failure;return result;}continue;}
                const bool equations=std::max({metrics.cellMomentum,metrics.faceMomentum,metrics.divergence})<=c.equationTolerance;
                const bool converged=equations&&((c.equation==CompatibleEquation2D::Stokes&&!pseudo)||change<=c.stateTolerance);
                const bool sufficient=metrics.residualNorm<=std::sqrt(1-2*c.armijo*alpha)*old.residualNorm;
                if(pseudo||sufficient||converged){double next=step;
                    if(pseudo){next=metrics.residualNorm==0||old.residualNorm==0||converged?c.maximumPseudoStep:std::min(c.maximumPseudoStep,step*(old.residualNorm/metrics.residualNorm));
                        if(!std::isfinite(next)||next<=0)throw std::runtime_error("Compatible pseudo-time step invalid; trial not accepted");}
                    auto physical=problem.physical(trial);poll(c);result.lastAccepted=std::move(physical);current=std::move(trial);old=metrics;step=next;accepted=true;record.accepted=true;
                    ++acceptedCount;
                    result.checkpoint=CheckpointAccess::make({context,current,mesh.cells.size(),mesh.faces.size(),acceptedCount,step,metrics});
                    if(c.checkpointAccepted)userCall([&]{c.checkpointAccepted(*result.checkpoint);});
                    if(c.iterationAccepted)userCall([&]{c.iterationAccepted(record);});
                    if(converged){result.stop=CompatibleFlowStop2D::Converged;result.reason="Original momentum, divergence, linear and applicable state-change targets satisfied";return result;}
                    break;
                }
                result.lastRejected=problem.physical(trial);
            }
            if(!accepted){result.stop=CompatibleFlowStop2D::BacktrackingBudget;result.reason="All backtracking trials rejected; last accepted iterate retained";return result;}
        }
        result.stop=CompatibleFlowStop2D::NonlinearBudget;result.reason="Nonlinear iteration budget exhausted; no converged flow";
    }catch(const IncompatibleCheckpoint&){throw;}
    catch(const UserException& e){std::rethrow_exception(e.exception);}
    catch(const Cancelled&){result.stop=CompatibleFlowStop2D::Cancelled;result.reason="Cancelled; unfinished trial not accepted";}
    catch(const std::exception& e){result.stop=CompatibleFlowStop2D::NumericalFailure;result.reason=e.what();}
    return result;
}
}
CompatibleFlowLoads2D evaluateCompatibleFlowLoads2D(const FvMesh2D& mesh,const CompatibleFlowControls2D& c,const CompatibleFlowCheckpoint2D& checkpoint,Point2D origin) {
    validate(mesh,c,std::nullopt);require(std::isfinite(origin.x)&&std::isfinite(origin.y),"Nonfinite compatible load moment origin");
    try {
        poll(c);Problem problem(mesh,c,true);const auto& saved=detail::CompatibleCheckpointAccess2D::get(checkpoint);
        if(!saved.context||*saved.context!=*problem.checkpointContext(mesh)||saved.normalizedState.size()!=9*mesh.cells.size()+4*mesh.faces.size())
            throw IncompatibleCheckpoint("Compatible load checkpoint mesh, scales, equation, boundary or native load differs");
        auto result=problem.loads(saved.normalizedState,origin);poll(c);result.acceptedIterations=saved.acceptedIterations;return result;
    }catch(const UserException& e){std::rethrow_exception(e.exception);}
    catch(const Cancelled&){throw std::runtime_error("Compatible load evaluation cancelled; no report published");}
}
CompatibleFlowResult2D solveCompatibleIncompressible2D(const FvMesh2D& mesh,const CompatibleFlowControls2D& c,const std::optional<CompatibleFlowState2D>& initial) {
    return run(mesh,c,initial,nullptr);
}
CompatibleFlowResult2D resumeCompatibleIncompressible2D(const FvMesh2D& mesh,const CompatibleFlowControls2D& c,const CompatibleFlowCheckpoint2D& checkpoint) {
    return run(mesh,c,std::nullopt,&checkpoint);
}
}
