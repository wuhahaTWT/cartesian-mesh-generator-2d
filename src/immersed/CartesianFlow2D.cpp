#include "cartmesh2d/immersed/CartesianFlow2D.hpp"
#include "cartmesh2d/spatial/BoundarySegmentIndex2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>

namespace cartmesh2d::immersed {
namespace {
using Clock = std::chrono::steady_clock;
using namespace fv::detail;
double seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }
void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
bool positive(double x) { return std::isfinite(x) && x > 0; }
double referenceSpeed(const Grid& g) { const auto& c=g.controls; return c.drive*c.height*c.height/(12*c.viscosity); }
double mask(double signedDistance, double width) {
    if (signedDistance <= -width) return 1;
    if (signedDistance >= width) return 0;
    const double r=signedDistance/width;
    return .5*(1-r-std::sin(std::numbers::pi*r)/std::numbers::pi);
}
std::size_t previous(std::size_t i, std::size_t n) { return i ? i-1 : n-1; }
std::size_t next(std::size_t i, std::size_t n) { return i+1==n ? 0 : i+1; }
void validState(const Grid& g,const State& s) {
    require(s.u.size()==g.maskU.size() && s.v.size()==g.maskV.size() && s.p.size()==s.u.size(),"immersed: state dimensions do not match grid");
    for(const auto* field:{&s.u,&s.v,&s.p}) for(double a:*field)
        require(std::isfinite(a),"immersed: nonfinite field");
    for(std::size_t i=0;i<g.controls.nx;++i)
        require(s.v[g.index(i,0)]==0 && s.v[g.index(i,g.controls.ny)]==0,"immersed: wall normal velocity is not zero");
}
double upwind(double flux,double left,double right) { return flux*(flux>=0 ? left : right); }
// Conservative donor-cell advection on the staggered velocity control volumes.
// No solid cells are deleted. Diffusion uses an odd tangential ghost value at
// the two domain walls. The original input geometry only defines the drag mask.
void unforced(const Grid& g,const State& s,std::vector<double>& ru,std::vector<double>& rv) {
    const auto& c=g.controls; const double hx=g.dx(),hy=g.dy();
    ru.assign(s.u.size(),0);rv.assign(s.v.size(),0);
    for(std::size_t j=0;j<c.ny;++j) for(std::size_t i=0;i<c.nx;++i) {
        const auto k=g.index(i,j),im=previous(i,c.nx),ip=next(i,c.nx);
        const double u=s.u[k],w=s.u[g.index(im,j)],e=s.u[g.index(ip,j)];
        const double n=j+1<c.ny ? s.u[g.index(i,j+1)] : -u;
        const double b=j ? s.u[g.index(i,j-1)] : -u;
        const double fn=.5*(s.v[g.index(im,j+1)]+s.v[g.index(i,j+1)]);
        const double fs=.5*(s.v[g.index(im,j)]+s.v[g.index(i,j)]);
        const double adv=(upwind(.5*(u+e),u,e)-upwind(.5*(w+u),w,u))/hx
                        +(upwind(fn,u,n)-upwind(fs,b,u))/hy;
        ru[k]=c.viscosity*((w-2*u+e)/(hx*hx)+(b-2*u+n)/(hy*hy))-adv+c.drive;
    }
    for(std::size_t j=1;j<c.ny;++j) for(std::size_t i=0;i<c.nx;++i) {
        const auto k=g.index(i,j),im=previous(i,c.nx),ip=next(i,c.nx);
        const double v=s.v[k],w=s.v[g.index(im,j)],e=s.v[g.index(ip,j)];
        const double n=s.v[g.index(i,j+1)],b=s.v[g.index(i,j-1)];
        const double fe=.5*(s.u[g.index(ip,j-1)]+s.u[g.index(ip,j)]);
        const double fw=.5*(s.u[g.index(i,j-1)]+s.u[g.index(i,j)]);
        const double adv=(upwind(fe,v,e)-upwind(fw,w,v))/hx
                        +(upwind(.5*(v+n),v,n)-upwind(.5*(b+v),b,v))/hy;
        rv[k]=c.viscosity*((w-2*v+e)/(hx*hx)+(b-2*v+n)/(hy*hy))-adv;
    }
}
double divergence(const Grid& g,const State& s,std::size_t i,std::size_t j) {
    return (s.u[g.index(next(i,g.controls.nx),j)]-s.u[g.index(i,j)])/g.dx()
          +(s.v[g.index(i,j+1)]-s.v[g.index(i,j)])/g.dy();
}
// Bilinear interpolation at the unmodified polygon segments, including the
// periodic seam and tangential no-slip ghost samples at domain walls.
double sample(const Grid& g,const std::vector<double>& f,Point2D point,bool uComponent) {
    const int nx=static_cast<int>(g.controls.nx),ny=static_cast<int>(g.controls.ny);
    double x=point.x/g.dx()-(uComponent?0:.5), y=point.y/g.dy()-(uComponent?.5:0);
    const int i=static_cast<int>(std::floor(x)),j=static_cast<int>(std::floor(y));
    x-=i;y-=j;
    const auto at=[&](int a,int b) {
        a=(a%nx+nx)%nx;
        if(uComponent && b<0)return -f[g.index(static_cast<std::size_t>(a),0)];
        if(uComponent && b>=ny)return -f[g.index(static_cast<std::size_t>(a),g.controls.ny-1)];
        require(b>=0 && b<=ny,"immersed: interpolation outside domain");
        return f[g.index(static_cast<std::size_t>(a),static_cast<std::size_t>(b))];
    };
    return (1-y)*((1-x)*at(i,j)+x*at(i+1,j))+y*((1-x)*at(i,j+1)+x*at(i+1,j+1));
}
Metrics equations(const Grid& g,const State& s) {
    Metrics m;std::vector<double> ru,rv;unforced(g,s,ru,rv);
    const auto& c=g.controls;const double ref=referenceSpeed(g);
    const double scale=c.drive; // relative to the imposed driving acceleration
    for(std::size_t j=0;j<c.ny;++j) for(std::size_t i=0;i<c.nx;++i) {
        const auto k=g.index(i,j),im=previous(i,c.nx);
        const double residual=ru[k]-(s.p[k]-s.p[g.index(im,j)])/g.dx()-g.maskU[k]*s.u[k]/c.penaltyTime;
        require(std::isfinite(residual),"immersed: nonfinite momentum residual");
        m.momentum=std::max(m.momentum,std::abs(residual)/scale);
        m.continuity=std::max(m.continuity,std::abs(divergence(g,s,i,j))*c.height/ref);
    }
    for(std::size_t j=1;j<c.ny;++j) for(std::size_t i=0;i<c.nx;++i) {
        const auto k=g.index(i,j);
        const double residual=rv[k]-(s.p[k]-s.p[g.index(i,j-1)])/g.dy()-g.maskV[k]*s.v[k]/c.penaltyTime;
        require(std::isfinite(residual),"immersed: nonfinite momentum residual");
        m.momentum=std::max(m.momentum,std::abs(residual)/scale);
    }
    require(std::isfinite(m.momentum) && std::isfinite(m.continuity),"immersed: nonfinite equation residual");
    return m;
}
}

Grid makeGrid(const Controls& c,const std::vector<BoundaryLoop>& solids) {
    const auto start=Clock::now();
    require(c.nx>=4 && c.ny>=4 && c.nx<=2048 && c.ny<=2048 && c.nx*c.ny<=1048576,"immersed: grid size must be at least 4x4 and at most 1048576 cells");
    require(c.maxSteps>0 && c.maxSteps<=10000000,"immersed: invalid step budget");
    require(positive(c.length)&&positive(c.height)&&positive(c.viscosity)&&positive(c.drive)&&positive(c.penaltyTime)&&positive(c.maxTimeStep),"immersed: invalid dimensional control");
    require(positive(c.maskHalfWidthCells)&&c.maskHalfWidthCells<=2,"immersed: mask half width must be in (0,2] cells");
    require(positive(c.steadyTolerance)&&c.steadyTolerance<1&&positive(c.continuityTolerance)&&c.continuityTolerance<1&&positive(c.linearTolerance)&&c.linearTolerance<=1e-2,"immersed: invalid residual control");
    require(!c.systemCholesky || fv::detail::systemCholeskyAvailable2D(),"immersed: system Cholesky requires macOS");
    Grid g;g.controls=c;g.solids=solids;
    require(positive(referenceSpeed(g)) && positive(g.dx()) && positive(g.dy()),"immersed: numerical scale overflow");
    g.maskU.assign(c.nx*c.ny,0);g.maskV.assign(c.nx*(c.ny+1),0);
    g.maskCell.assign(c.nx*c.ny,0);g.classification.assign(c.nx*c.ny,0);
    g.gridSeconds=seconds(start);const auto boundaryStart=Clock::now();
    if(!solids.empty()) {
        const BoundaryRegion2D region(solids);
        require(region.diagnose().valid(),"immersed: invalid solid geometry (zero area, duplicate edges, self-intersection or touching loops)");
        g.solidArea=region.area();
        const auto bounds=region.bounds();const double width=c.maskHalfWidthCells*std::min(g.dx(),g.dy());
        const double guard=2*std::max(g.dx(),g.dy())+width;
        require(bounds.min.x>guard && bounds.max.x<c.length-guard && bounds.min.y>guard && bounds.max.y<c.height-guard,"immersed: solid requires two-cell clearance from domain walls and periodic seam");
        const BoundarySegmentIndex2D index(region);require(index.valid(),"immersed: invalid boundary index");
        // Point-to-segment distance is evaluated on original segments. No
        // polygon clipping or smoothing of the input boundary is performed.
        const auto chi=[&](Point2D p) {
            double distance=std::numeric_limits<double>::infinity();
            for(const auto& loop:solids) for(std::size_t k=0;k<loop.vertices().size();++k) {
                const auto a=loop.vertices()[k],b=loop.vertices()[(k+1)%loop.vertices().size()];
                const auto ab=b-a,ap=p-a;const double t=std::clamp(dot(ap,ab)/squaredNorm(ab),0.,1.);
                distance=std::min(distance,std::hypot(ap.x-t*ab.x,ap.y-t*ab.y));
            }
            if(index.classifyPoint(p)!=PointInPolygon::Outside)distance=-distance;
            return mask(distance,width);
        };
        for(std::size_t j=0;j<c.ny;++j) for(std::size_t i=0;i<c.nx;++i) {
            const auto k=g.index(i,j);const double x=static_cast<double>(i)*g.dx(),y=static_cast<double>(j)*g.dy();
            g.maskU[k]=chi({x,y+.5*g.dy()});g.maskV[k]=chi({x+.5*g.dx(),y});
            const Point2D centre{x+.5*g.dx(),y+.5*g.dy()};g.maskCell[k]=chi(centre);
            g.classification[k]=index.intersects({{x,y},{x+g.dx(),y+g.dy()}}) ? 2U :
                (index.classifyPoint(centre)==PointInPolygon::Outside ? 0U : 1U);
        }
        require(std::count(g.classification.begin(),g.classification.end(),1U)>0,"immersed: solid is under-resolved (no fully interior cell)");
    }
    g.boundarySeconds=seconds(boundaryStart);return g;
}
State zeroState(const Grid& g) { State s;s.u.assign(g.maskU.size(),0);s.v.assign(g.maskV.size(),0);s.p.assign(g.maskU.size(),0);return s; }
std::vector<WallSample> sampleWalls(const Grid& g,const State& s) {
    validState(g,s);std::vector<WallSample> result;
    const auto depths=g.solids.empty() ? std::vector<std::size_t>{} : BoundaryRegion2D(g.solids).nestingDepths();
    for(std::size_t l=0;l<g.solids.size();++l) {
      const auto& loop=g.solids[l];
      for(std::size_t k=0;k<loop.vertices().size();++k) {
        const auto a=loop.vertices()[k],b=loop.vertices()[(k+1)%loop.vertices().size()];
        const double length=std::hypot(b.x-a.x,b.y-a.y);
        const auto n=static_cast<std::size_t>(std::max(1.,std::ceil(length/(.5*std::min(g.dx(),g.dy())))));
        const double sign=(loop.polygon().signedArea()>0 ? 1 : -1)*(depths[l]%2 ? -1 : 1);
        for(std::size_t q=0;q<n;++q) {
            const double t=(static_cast<double>(q)+.5)/static_cast<double>(n);
            const Point2D p{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)};
            result.push_back({p,{sign*(b.y-a.y)/length,sign*(a.x-b.x)/length},
                {sample(g,s.u,p,true),sample(g,s.v,p,false)},length/static_cast<double>(n)});
        }
    }
    }
    return result;
}
Metrics evaluate(const Grid& g,const State& s) {
    validState(g,s);auto m=equations(g,s);const auto& c=g.controls;
    double minFlux=std::numeric_limits<double>::infinity(),maxFlux=-minFlux,error2=0,exact2=0;
    for(std::size_t i=0;i<c.nx;++i) {
        double flux=0;
        for(std::size_t j=0;j<c.ny;++j) {
            const auto k=g.index(i,j);flux+=s.u[k]*g.dy();
            const double u=.5*(s.u[k]+s.u[g.index(next(i,c.nx),j)]),v=.5*(s.v[k]+s.v[g.index(i,j+1)]);
            const double speed=std::hypot(u,v);m.maxSpeed=std::max(m.maxSpeed,speed);
            if(g.maskCell[k]==1)m.deepSolidSpeed=std::max(m.deepSolidSpeed,speed);
            m.meanVelocity+=s.u[k]/static_cast<double>(c.nx*c.ny);
            m.penaltyDrag+=g.maskU[k]*s.u[k]*g.dx()*g.dy()/c.penaltyTime;
            m.penaltyLift+=g.maskV[k]*s.v[k]*g.dx()*g.dy()/c.penaltyTime;
            const double y=(static_cast<double>(j)+.5)*g.dy(),exact=c.drive*y*(c.height-y)/(2*c.viscosity);
            error2+=(u-exact)*(u-exact)+v*v;exact2+=exact*exact;
        }
        minFlux=std::min(minFlux,flux);maxFlux=std::max(maxFlux,flux);
    }
    m.fluxSpread=(maxFlux-minFlux)/(referenceSpeed(g)*c.height);
    if(g.solids.empty())m.channelRelativeL2=std::sqrt(error2/exact2);
    for(const auto& w:sampleWalls(g,s)) {
        m.wallSpeed=std::max(m.wallSpeed,std::hypot(w.velocity.x,w.velocity.y));
        m.wallNormalSpeed=std::max(m.wallNormalSpeed,std::abs(dot(w.normal,w.velocity)));
        m.wallTangentialSpeed=std::max(m.wallTangentialSpeed,std::abs(w.normal.x*w.velocity.y-w.normal.y*w.velocity.x));
    }
    return m;
}
Result solve(const Grid& g,const std::function<bool(const State&,const Metrics&)>& progress) {
    const auto start=Clock::now();const auto& c=g.controls;Result result;result.state=zeroState(g);
    const double ref=referenceSpeed(g),hx=g.dx(),hy=g.dy();
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    for(std::size_t j=0;j<c.ny;++j)for(std::size_t i=0;i<c.nx;++i) {
        edges.emplace_back(g.index(i,j),g.index(next(i,c.nx),j));
        if(j+1<c.ny)edges.emplace_back(g.index(i,j),g.index(i,j+1));
    }
    SparsePattern2D pattern(c.nx*c.ny,edges);SparseSystem2D pressure(pattern);
    LinearWorkspace2D workspace(c.nx*c.ny);
    std::vector<double> bu(g.maskU.size()),bv(g.maskV.size()),ru,rv,correction(c.nx*c.ny),applied;
    double matrixDt=-1;
    // Monotone explicit transport/diffusion bound; implicit drag is NOT included
    // in this CFL restriction. Decreasing dt never alters the requested residuals.
    double dt=std::min(c.maxTimeStep,.45/(2*c.viscosity*(1/(hx*hx)+1/(hy*hy))+2*ref/hx+ref/hy));
    try {
        for(std::size_t step=0;step<c.maxSteps;++step) {
            const auto& old=result.state;
            double umax=0,vmax=0;
            for(double u:old.u)umax=std::max(umax,std::abs(u));
            for(double v:old.v)vmax=std::max(vmax,std::abs(v));
            dt=std::min(dt,.45/(2*c.viscosity*(1/(hx*hx)+1/(hy*hy))+umax/hx+vmax/hy));
            require(positive(dt),"immersed: invalid time step");
            if(dt!=matrixDt) {
                pressure.reset();
                for(std::size_t k=0;k<bu.size();++k)bu[k]=1/(1+dt*g.maskU[k]/c.penaltyTime);
                for(std::size_t k=0;k<bv.size();++k)bv[k]=1/(1+dt*g.maskV[k]/c.penaltyTime);
                const auto addEdge=[&](std::size_t a,std::size_t b,double weight) {
                    pressure.diag[a]+=weight;pressure.diag[b]+=weight;
                    pressure.add(a,b,-weight);pressure.add(b,a,-weight);
                };
                for(std::size_t j=0;j<c.ny;++j)for(std::size_t i=0;i<c.nx;++i) {
                    const auto ip=next(i,c.nx);addEdge(g.index(i,j),g.index(ip,j),bu[g.index(ip,j)]/(hx*hx));
                    if(j+1<c.ny)addEdge(g.index(i,j),g.index(i,j+1),bv[g.index(i,j+1)]/(hy*hy));
                }
                pressure.pin(0);matrixDt=dt;
            }
            unforced(g,old,ru,rv);State candidate=old;
            for(std::size_t j=0;j<c.ny;++j)for(std::size_t i=0;i<c.nx;++i) {
                const auto k=g.index(i,j);
                candidate.u[k]=bu[k]*(old.u[k]+dt*(ru[k]-(old.p[k]-old.p[g.index(previous(i,c.nx),j)])/hx));
                if(j)candidate.v[k]=bv[k]*(old.v[k]+dt*(rv[k]-(old.p[k]-old.p[g.index(i,j-1)])/hy));
            }
            for(std::size_t j=0;j<c.ny;++j)for(std::size_t i=0;i<c.nx;++i)
                pressure.rhs[g.index(i,j)]=-divergence(g,candidate,i,j)/dt;
            pressure.rhs[0]=0;std::fill(correction.begin(),correction.end(),0);
            const auto ps=Clock::now();
            const auto iterations=pressure.solvePressure(correction,workspace,c.systemCholesky?LinearPressureMethod2D::SystemCholesky:LinearPressureMethod2D::IC0,c.linearTolerance);
            result.pressureSeconds+=seconds(ps);result.pressureIterations+=iterations;
            pressure.apply(correction,applied);double residual2=0,rhs2=0;
            for(std::size_t k=0;k<applied.size();++k) {
                residual2+=(applied[k]-pressure.rhs[k])*(applied[k]-pressure.rhs[k]);rhs2+=pressure.rhs[k]*pressure.rhs[k];
                candidate.p[k]+=correction[k];
            }
            for(std::size_t j=0;j<c.ny;++j)for(std::size_t i=0;i<c.nx;++i) {
                const auto k=g.index(i,j);
                candidate.u[k]-=dt*bu[k]*(correction[k]-correction[g.index(previous(i,c.nx),j)])/hx;
                if(j)candidate.v[k]-=dt*bv[k]*(correction[k]-correction[g.index(i,j-1)])/hy;
            }
            validState(g,candidate);Metrics m=equations(g,candidate);
            for(std::size_t k=0;k<old.u.size();++k)m.fieldChange=std::max(m.fieldChange,std::abs(candidate.u[k]-old.u[k])/ref);
            for(std::size_t k=0;k<old.v.size();++k)m.fieldChange=std::max(m.fieldChange,std::abs(candidate.v[k]-old.v[k])/ref);
            m.dt=dt;m.pressureIterations=iterations;m.pressureLinearResidual=std::sqrt(residual2);m.pressureLinearRhsNorm=std::sqrt(rhs2);
            // Conservation is checked after the velocity correction, including
            // the pressure reference cell omitted by the linear solve.
            require(m.continuity<=c.continuityTolerance,"immersed: corrected candidate fails continuity tolerance");
            candidate.steps=old.steps+1;candidate.pseudoTime=old.pseudoTime+dt;
            result.state=std::move(candidate);result.metrics=m;
            if(progress && !progress(result.state,m)) { result.stopReason="cancelled";break; }
            if(result.state.steps>=20 && m.momentum<=c.steadyTolerance && m.fieldChange<=c.steadyTolerance) {
                result.converged=true;result.stopReason="steady-converged";break;
            }
        }
    } catch(const std::exception& e) { result.stopReason="candidate-failed";result.error=e.what(); }
    const auto iterationMetrics=result.metrics;result.metrics=evaluate(g,result.state);
    result.metrics.fieldChange=iterationMetrics.fieldChange;result.metrics.dt=iterationMetrics.dt;
    result.metrics.pressureIterations=iterationMetrics.pressureIterations;
    result.metrics.pressureLinearResidual=iterationMetrics.pressureLinearResidual;
    result.metrics.pressureLinearRhsNorm=iterationMetrics.pressureLinearRhsNorm;
    result.solveSeconds=seconds(start);return result;
}
}
