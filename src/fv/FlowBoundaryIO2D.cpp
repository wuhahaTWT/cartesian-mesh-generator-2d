#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(std::string("Flow boundaries: ")+message);
}
void token(std::istream& input,const char* expected) {
    std::string value;
    require(bool(input>>value) && value==expected,"unexpected token or truncated file");
}
std::size_t boundaryCount(const FvMesh2D& mesh) {
    return static_cast<std::size_t>(std::count_if(mesh.faces.begin(),mesh.faces.end(),
        [](const Face& face){return !face.neighbour;}));
}
}

const char* flowBoundaryKindName2D(FlowBoundaryKind2D kind) {
    switch(kind) {
    case FlowBoundaryKind2D::VelocityInlet:return "velocity-inlet";
    case FlowBoundaryKind2D::PressureOutlet:return "pressure-outlet";
    case FlowBoundaryKind2D::PressureOpening:return "pressure-opening";
    case FlowBoundaryKind2D::Wall:return "wall";
    case FlowBoundaryKind2D::MovingWall:return "moving-wall";
    case FlowBoundaryKind2D::SmoothMovingWall:return "smooth-moving-wall";
    }
    throw std::runtime_error("Flow boundaries: unknown condition type");
}
FlowBoundaryKind2D flowBoundaryKindFromName2D(const std::string& name) {
    for(auto kind:{FlowBoundaryKind2D::VelocityInlet,FlowBoundaryKind2D::PressureOutlet,FlowBoundaryKind2D::PressureOpening,
                  FlowBoundaryKind2D::Wall,FlowBoundaryKind2D::MovingWall,FlowBoundaryKind2D::SmoothMovingWall})
        if(name==flowBoundaryKindName2D(kind))return kind;
    throw std::runtime_error("Flow boundaries: unknown condition type");
}

std::vector<FlowBoundaryCondition2D> rotatingAnnulusBoundaryPreset2D(
    const FvMesh2D& mesh, double speed) {
    validateFvMesh2D(mesh);
    require(std::isfinite(speed) && speed > 0, "inner surface speed must be positive");
    double xmin=std::numeric_limits<double>::infinity(),ymin=xmin,xmax=-xmin,ymax=-xmin;
    for (const auto& f : mesh.faces) if (!f.neighbour) {
        xmin=std::min(xmin,f.centre.x-.5*std::abs(f.areaVector.y));
        xmax=std::max(xmax,f.centre.x+.5*std::abs(f.areaVector.y));
        ymin=std::min(ymin,f.centre.y-.5*std::abs(f.areaVector.x));
        ymax=std::max(ymax,f.centre.y+.5*std::abs(f.areaVector.x));
    }
    const Point2D centre{.5*(xmin+xmax),.5*(ymin+ymax)};
    const double eps=16*TolerancePolicy{}.scale(std::max(xmax-xmin,ymax-ymin));
    const double angleEps=128*TolerancePolicy{}.scale(1.);
    const double pi=std::acos(-1.);
    struct Ring { double distance=0,length=0; std::vector<double> angles; std::vector<std::size_t> faces; };
    Ring rings[2]; // inner, outer
    std::vector<FlowBoundaryCondition2D> result;
    for (std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id]; if (f.neighbour) continue;
        const double length=std::hypot(f.areaVector.x,f.areaVector.y);
        const auto n=f.areaVector*(1/length);
        const auto d=f.centre-centre;
        const double offset=d.x*n.x+d.y*n.y;
        require(std::abs(offset)>eps,"annulus boundary crosses its inferred centre");
        const bool inner=offset<0;
        auto& ring=rings[inner?0:1];
        if (ring.faces.empty()) ring.distance=std::abs(offset);
        require(std::abs(std::abs(offset)-ring.distance)<=eps,
                "annulus template requires two concentric regular circular polygons");
        double angle=std::atan2(n.y,n.x); if (angle<0) angle+=2*pi;
        ring.angles.push_back(angle);ring.faces.push_back(id);ring.length+=length;
        result.push_back({id,inner?FlowBoundaryKind2D::SmoothMovingWall:FlowBoundaryKind2D::Wall,
            inner?Vector2D{speed*n.y,-speed*n.x}:Vector2D{},0,inner?"rotor":"housing"});
    }
    double radii[2]{};
    for (int side=0;side<2;++side) {
        auto& ring=rings[side];auto& angles=ring.angles;
        require(!angles.empty(),"annulus template needs both inner and outer walls");
        std::sort(angles.begin(),angles.end());
        angles.erase(std::unique(angles.begin(),angles.end(),[&](double a,double b){return b-a<=angleEps;}),angles.end());
        if (angles.size()>1 && angles.front()+2*pi-angles.back()<=angleEps) angles.pop_back();
        const auto count=angles.size();
        require(count>=16 && count%2==0,"annulus template requires even regular polygons with at least 16 sides");
        for (std::size_t i=0;i<count;++i) {
            const double next=i+1<count?angles[i+1]:angles[0]+2*pi;
            require(std::abs(next-angles[i]-2*pi/count)<=angleEps,"annulus facet directions are not regularly spaced");
        }
        radii[side]=ring.distance/std::cos(pi/count);
        const double perimeter=2*count*ring.distance*std::tan(pi/count);
        require(std::abs(perimeter-ring.length)<=count*eps,"annulus perimeter is incomplete or duplicated");
        for (auto id:ring.faces) {
            const auto& f=mesh.faces[id];
            for (double sign:{-1.,1.}) {
                const double x=f.centre.x-centre.x+sign*.5*f.areaVector.y;
                const double y=f.centre.y-centre.y-sign*.5*f.areaVector.x;
                require(std::hypot(x,y)<=radii[side]+eps,"annulus face extends outside its regular polygon");
            }
        }
    }
    require(radii[0]+eps<rings[1].distance,"annulus inner and outer polygons overlap");
    FlowControls2D controls;controls.scenario="custom";controls.speed=speed;controls.boundaryConditions=result;
    validateFlowBoundaryConditions2D(mesh,controls);
    return result;
}

std::vector<FlowBoundaryCondition2D> readFlowBoundaryConditions2D(
    std::istream& input,const FvMesh2D& mesh,const FlowControls2D& controls) {
    validateFvMesh2D(mesh);
    require(controls.scenario=="custom","explicit file requires custom scenario");
    token(input,"CARTMESH2D_FLOW_BOUNDARIES");token(input,"1");token(input,"COUNTS");
    std::size_t nc=0,nf=0,nb=0;
    require(bool(input>>nc>>nf>>nb),"truncated mesh counts");
    require(nc==mesh.cells.size() && nf==mesh.faces.size() && nb==boundaryCount(mesh),
            "configuration mesh counts differ from final mesh");
    double xmin=std::numeric_limits<double>::infinity(),ymin=xmin,xmax=-xmin,ymax=-xmin;
    for(const auto& f:mesh.faces) {
        xmin=std::min(xmin,f.centre.x);xmax=std::max(xmax,f.centre.x);
        ymin=std::min(ymin,f.centre.y);ymax=std::max(ymax,f.centre.y);
    }
    const auto tolerance=TolerancePolicy{};
    const double positionTolerance=tolerance.scale(std::max(xmax-xmin,ymax-ymin));
    auto same=[](double a,double b,double eps){return std::isfinite(a)&&std::abs(a-b)<=eps;};
    std::vector<FlowBoundaryCondition2D> result;
    result.reserve(nb);
    for(std::size_t i=0;i<nb;++i) {
        token(input,"BOUNDARY");
        FlowBoundaryCondition2D condition;
        std::size_t owner=0;double x=0,y=0,sx=0,sy=0;std::string kind;
        require(bool(input>>condition.face>>owner>>x>>y>>sx>>sy>>kind>>std::quoted(condition.name)
                    >>condition.velocity.x>>condition.velocity.y>>condition.pressure),"truncated boundary record");
        require(condition.face<nf,"boundary face ID out of range");
        const auto& face=mesh.faces[condition.face];
        require(!face.neighbour && owner==face.owner,"face owner or boundary identity differs from mesh");
        const double vectorTolerance=tolerance.scale(std::hypot(face.areaVector.x,face.areaVector.y));
        require(same(x,face.centre.x,positionTolerance)&&same(y,face.centre.y,positionTolerance)&&
                same(sx,face.areaVector.x,vectorTolerance)&&same(sy,face.areaVector.y,vectorTolerance),
                "boundary face geometry differs from final mesh");
        condition.kind=flowBoundaryKindFromName2D(kind);
        result.push_back(std::move(condition));
    }
    token(input,"END");input>>std::ws;
    require(input.peek()==std::char_traits<char>::eof(),"trailing configuration data");
    auto check=controls;check.boundaryConditions=result;
    validateFlowBoundaryConditions2D(mesh,check);
    std::sort(result.begin(),result.end(),[](const auto& a,const auto& b){return a.face<b.face;});
    return result;
}

void writeFlowBoundaryConditions2D(std::ostream& output,const FvMesh2D& mesh,
                                  const FlowControls2D& controls) {
    validateFvMesh2D(mesh);
    require(controls.scenario=="custom","explicit writer requires custom scenario");
    validateFlowBoundaryConditions2D(mesh,controls);
    const auto oldFlags=output.flags();const auto oldPrecision=output.precision();
    struct Restore {
        std::ostream& out;std::ios_base::fmtflags flags;std::streamsize precision;
        ~Restore(){out.flags(flags);out.precision(precision);}
    } restore{output,oldFlags,oldPrecision};
    output<<std::defaultfloat<<std::dec<<std::noshowpos<<std::noshowbase<<std::setprecision(17);
    output<<"CARTMESH2D_FLOW_BOUNDARIES 1\nCOUNTS "<<mesh.cells.size()<<' '<<mesh.faces.size()
          <<' '<<controls.boundaryConditions.size()<<'\n';
    auto entries=controls.boundaryConditions;
    std::sort(entries.begin(),entries.end(),[](const auto& a,const auto& b){return a.face<b.face;});
    for(const auto& b:entries) {
        const auto& f=mesh.faces[b.face];
        output<<"BOUNDARY "<<b.face<<' '<<f.owner<<' '<<f.centre.x<<' '<<f.centre.y<<' '
              <<f.areaVector.x<<' '<<f.areaVector.y<<' '<<flowBoundaryKindName2D(b.kind)<<' '
              <<std::quoted(b.name)<<' '<<b.velocity.x<<' '<<b.velocity.y<<' '<<b.pressure<<'\n';
    }
    output<<"END\n";
    require(bool(output),"could not write configuration");
}
}
