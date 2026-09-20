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
    case FlowBoundaryKind2D::Wall:return "wall";
    case FlowBoundaryKind2D::MovingWall:return "moving-wall";
    }
    throw std::runtime_error("Flow boundaries: unknown condition type");
}
FlowBoundaryKind2D flowBoundaryKindFromName2D(const std::string& name) {
    for(auto kind:{FlowBoundaryKind2D::VelocityInlet,FlowBoundaryKind2D::PressureOutlet,
                  FlowBoundaryKind2D::Wall,FlowBoundaryKind2D::MovingWall})
        if(name==flowBoundaryKindName2D(kind))return kind;
    throw std::runtime_error("Flow boundaries: unknown condition type");
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
