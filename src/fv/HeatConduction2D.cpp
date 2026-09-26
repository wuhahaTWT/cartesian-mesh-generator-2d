#include "cartmesh2d/fv/HeatConduction2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(std::string("heat conduction: ")+message);}
double checked(double x){require(std::isfinite(x),"non-finite arithmetic");return x;}
}
HeatConductionOperator2D::HeatConductionOperator2D(const FvMesh2D& mesh,
    const std::vector<HeatBoundary2D>& boundaries,double conductivity):conductivity_(conductivity) {
    validateFvMesh2D(mesh);require(std::isfinite(conductivity)&&conductivity>0,"conductivity must be positive");
    std::vector<const HeatBoundary2D*> lookup(mesh.faces.size());
    for(const auto& b:boundaries) {
        require(b.face<lookup.size()&&!mesh.faces[b.face].neighbour&&!lookup[b.face],"invalid/duplicate thermal boundary");
        require(std::isfinite(b.value),"non-finite boundary value");
        require(b.kind==HeatBoundaryKind2D::Insulated||b.kind==HeatBoundaryKind2D::Temperature||b.kind==HeatBoundaryKind2D::OutwardFlux||b.kind==HeatBoundaryKind2D::Periodic,"unknown boundary kind");
        require(b.kind!=HeatBoundaryKind2D::Temperature||b.value>0,"wall temperature must be positive Kelvin");
        require((b.kind==HeatBoundaryKind2D::Periodic)==b.partner.has_value(),"thermal periodic pairing mismatch");
        require((b.kind!=HeatBoundaryKind2D::Insulated&&b.kind!=HeatBoundaryKind2D::Periodic)||b.value==0,"inactive boundary value must be zero");
        lookup[b.face]=&b;
    }
    faces_.resize(mesh.faces.size());gradients_.resize(mesh.cells.size());
    for(const auto& c:mesh.cells)areas_.push_back(c.area);
    for(std::size_t id=0;id<faces_.size();++id) {
        const auto& f=mesh.faces[id];auto& out=faces_[id];const auto* b=lookup[id];
        require(f.neighbour||b,"missing thermal boundary");
        out.owner=f.owner;out.neighbour=f.neighbour;out.length=std::hypot(f.areaVector.x,f.areaVector.y);
        out.transmissibility=f.transmissibility;out.correction=f.correction;out.weight=f.neighbourWeight;
        if(!b)continue;
        out.kind=b->kind;out.value=b->value;out.partner=b->partner;
        if(!b->partner)continue;
        const auto partner=*b->partner;
        require(partner<lookup.size()&&partner!=id&&lookup[partner]&&lookup[partner]->partner==id,"nonreciprocal thermal periodic pairing");
        const auto& other=mesh.faces[partner];
        require(std::hypot(f.areaVector.x+other.areaVector.x,f.areaVector.y+other.areaVector.y)<=1e-10*out.length,"periodic thermal normals differ");
        out.neighbour=other.owner;
        const Vector2D d{mesh.cells[other.owner].centre.x+f.centre.x-other.centre.x-mesh.cells[f.owner].centre.x,
                         mesh.cells[other.owner].centre.y+f.centre.y-other.centre.y-mesh.cells[f.owner].centre.y};
        const double sd=dot(f.areaVector,d);require(sd>0,"nonpositive periodic normal distance");
        out.transmissibility=checked(out.length*out.length/sd);
        out.correction={f.areaVector.x-out.transmissibility*d.x,f.areaVector.y-out.transmissibility*d.y};
        out.weight=dot(f.areaVector,{f.centre.x-mesh.cells[f.owner].centre.x,f.centre.y-mesh.cells[f.owner].centre.y})/sd;
        require(out.weight>0&&out.weight<1,"periodic face outside centre bracket");
    }
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        struct Constraint {Vector2D direction;std::optional<std::size_t> cell;double value;bool derivative;};
        std::vector<Constraint> constraints;double xx=0,xy=0,yy=0;
        const auto centre=mesh.cells[i].centre;
        for(auto id:mesh.cells[i].faces) {
            const auto& f=mesh.faces[id];const auto& h=faces_[id];
            std::optional<std::size_t> j;Point2D point=f.centre;
            if(f.neighbour){j=f.owner==i?*f.neighbour:f.owner;point=mesh.cells[*j].centre;}
            else if(h.partner){j=h.neighbour;const auto& other=mesh.faces[*h.partner];point={mesh.cells[*j].centre.x+f.centre.x-other.centre.x,mesh.cells[*j].centre.y+f.centre.y-other.centre.y};}
            const bool derivative=!j&&h.kind!=HeatBoundaryKind2D::Temperature;
            Vector2D d{point.x-centre.x,point.y-centre.y};double value=h.value;
            if(derivative){d={f.areaVector.x/h.length,f.areaVector.y/h.length};value=-value/conductivity_;}
            const double norm2=dot(d,d);require(norm2>0,"zero thermal stencil distance");
            xx+=d.x*d.x/norm2;xy+=d.x*d.y/norm2;yy+=d.y*d.y/norm2;
            constraints.push_back({{d.x/norm2,d.y/norm2},j,value,derivative});
        }
        const double det=checked(xx*yy-xy*xy),trace=xx+yy;
        require(det>64*std::numeric_limits<double>::epsilon()*trace*trace,"rank-deficient temperature gradient stencil");
        auto& gradient=gradients_[i];
        for(const auto& c:constraints) {
            Vector2D w{checked((yy*c.direction.x-xy*c.direction.y)/det),checked((xx*c.direction.y-xy*c.direction.x)/det)};
            if(c.derivative){gradient.constant.x+=w.x*c.value;gradient.constant.y+=w.y*c.value;}
            else gradient.samples.push_back({c.cell,c.value,w});
        }
    }
    // Assemble only the constant Jacobian for a dimensional time-step estimate.
    // Runtime flux evaluation below uses temperature DIFFERENCES, preserving
    // a constant field exactly without cancellation of large absolute Kelvin.
    std::vector<std::map<std::size_t,double>> rows(mesh.cells.size());
    for(std::size_t id=0;id<faces_.size();++id) {
        const auto& f=faces_[id];if(f.partner&&id>*f.partner)continue;
        if(!f.neighbour&&f.kind!=HeatBoundaryKind2D::Temperature)continue;
        std::map<std::size_t,double> coefficients;
        const double a=conductivity_*f.transmissibility;
        coefficients[f.owner]+=a;if(f.neighbour)coefficients[*f.neighbour]-=a;
        const auto addGradient=[&](std::size_t i,double weight){
            for(const auto& sample:gradients_[i].samples) {
                const double c=-conductivity_*weight*dot(sample.weight,f.correction);
                coefficients[i]-=c;if(sample.cell)coefficients[*sample.cell]+=c;
            }
        };
        addGradient(f.owner,f.neighbour?1-f.weight:1);
        if(f.neighbour)addGradient(*f.neighbour,f.weight);
        for(const auto& [j,c]:coefficients){rows[f.owner][j]+=c;if(f.neighbour)rows[*f.neighbour][j]-=c;}
    }
    rowNorm_.resize(rows.size());
    for(std::size_t i=0;i<rows.size();++i) {
        double sum=0;for(const auto& [j,c]:rows[i])sum+=std::abs(c);
        rowNorm_[i]=checked(.5*sum);
        bool nonMonotone=false;
        for(const auto& [j,c]:rows[i])if(j!=i&&c>64*std::numeric_limits<double>::epsilon()*sum)nonMonotone=true;
        if(nonMonotone)++nonMonotoneRows_;
    }
}
HeatConductionResult2D HeatConductionOperator2D::evaluate(const std::vector<double>& t,const std::vector<double>& capacity) const {
    require(t.size()==areas_.size()&&capacity.size()==t.size(),"temperature/capacity size mismatch");
    HeatConductionResult2D out;out.faceHeatFlux.resize(faces_.size());out.cellResidual.resize(t.size());out.rate.resize(t.size());
    std::vector<Vector2D> gradients(t.size());
    for(std::size_t i=0;i<t.size();++i) {
        require(std::isfinite(t[i])&&t[i]>0&&std::isfinite(capacity[i])&&capacity[i]>0,"nonpositive/nonfinite temperature or heat capacity");
        auto g=gradients_[i].constant;
        for(const auto& s:gradients_[i].samples){const double d=(s.cell?t[*s.cell]:s.value)-t[i];g.x+=s.weight.x*d;g.y+=s.weight.y*d;}
        gradients[i]={checked(g.x),checked(g.y)};out.rate[i]=checked(rowNorm_[i]/(areas_[i]*capacity[i]));
    }
    for(std::size_t id=0;id<faces_.size();++id) {
        const auto& f=faces_[id];if(f.partner&&id>*f.partner)continue;
        double q=f.kind==HeatBoundaryKind2D::OutwardFlux?f.value*f.length:0;
        if(f.neighbour||f.kind==HeatBoundaryKind2D::Temperature) {
            auto g=gradients[f.owner];
            if(f.neighbour){const auto h=gradients[*f.neighbour];g={g.x*(1-f.weight)+h.x*f.weight,g.y*(1-f.weight)+h.y*f.weight};}
            q=-conductivity_*(f.transmissibility*((f.neighbour?t[*f.neighbour]:f.value)-t[f.owner])+dot(g,f.correction));
        }
        q=checked(q);out.faceHeatFlux[id]=q;out.cellResidual[f.owner]+=q;
        if(f.neighbour)out.cellResidual[*f.neighbour]-=q;
        if(f.partner)out.faceHeatFlux[*f.partner]=-q;
    }
    for(auto q:out.cellResidual)checked(q);
    return out;
}
} // namespace cartmesh2d::fv
