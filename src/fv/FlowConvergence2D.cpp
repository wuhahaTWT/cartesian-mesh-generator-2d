#include "FlowSolverDetail2D.hpp"
#include <map>
namespace cartmesh2d::fv::solver_detail {
void prepareMonitors(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,FlowResult2D& r,
    std::vector<std::size_t>& monitorGroup,Vec& monitorLengths){
    const auto nf=m.faces.size();
    if(c.convergence==FlowConvergence2D::Engineering) {
        r.monitorNames.push_back("kinetic-energy/(area*Uref^2)");
        std::map<std::string,std::size_t> groups;
        std::map<std::size_t,std::string> names;
        for(const auto& condition:c.boundaryConditions)names[condition.face]=condition.name;
        for(std::size_t id=0;id<nf;++id)if(!m.faces[id].neighbour) {
            const std::string key=c.scenario=="custom" ? names.at(id)
                : "role-"+std::to_string(static_cast<int>(b.role[id]))+"-patch-"+std::to_string(static_cast<int>(m.faces[id].patch));
            const auto [it,added]=groups.emplace(key,groups.size());
            if(added) {
                monitorLengths.push_back(0);
                for(const auto* label:{"flux/(Uref*length)","forceX/(pressureScale*length)","forceY/(pressureScale*length)","mean-pressure/pressureScale","torque/(pressureScale*length^2)"})
                    r.monitorNames.push_back(key+":"+label);
            }
            monitorGroup[id]=it->second;
            monitorLengths[it->second]+=std::hypot(m.faces[id].areaVector.x,m.faces[id].areaVector.y);
        }
    }
}
Vec physicalMonitors(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,const FlowResult2D& r,
    const std::vector<std::size_t>& monitorGroup,const Vec& monitorLengths,double pressureScale,
    const Vec& pressureFaces,const std::vector<Vector2D>& gu,const std::vector<Vector2D>& gv,
    const std::vector<Vector2D>& stressCorrection){
    const auto n=m.cells.size(),nf=m.faces.size();
        Vec values(r.monitorNames.size());
        if(values.empty())return values;
        double totalArea=0;
        for(std::size_t i=0;i<n;++i) {
            totalArea+=m.cells[i].area;
            values[0]+=.5*m.cells[i].area*(r.u[i]*r.u[i]+r.v[i]*r.v[i]);
        }
        values[0]/=totalArea*c.speed*c.speed;
        for(std::size_t id=0;id<nf;++id)if(!m.faces[id].neighbour) {
            const auto& face=m.faces[id];const auto i=face.owner,g=monitorGroup[id],offset=1+5*g;
            values[offset]+=r.flux[id]/(c.speed*monitorLengths[g]);
            values[offset+3]+=pressureFaces[id]*std::hypot(face.areaVector.x,face.areaVector.y)/(pressureScale*monitorLengths[g]);
            if(b.role[id]!=Role::Wall && b.role[id]!=Role::Lid)continue;
            const double dx=-faceNu(c,id)*(face.transmissibility*(b.u[id]-r.u[i])+dot(gu[i],face.correction))
                +(stressCorrection.empty()?0:stressCorrection[id].x);
            const double dy=-faceNu(c,id)*(face.transmissibility*(b.v[id]-r.v[i])+dot(gv[i],face.correction))
                +(stressCorrection.empty()?0:stressCorrection[id].y);
            values[offset+1]+=(pressureFaces[id]*face.areaVector.x+dx)/(pressureScale*monitorLengths[g]);
            values[offset+2]+=(pressureFaces[id]*face.areaVector.y+dy)/(pressureScale*monitorLengths[g]);
            const double x=face.centre.x-.5*(b.xmin+b.xmax),y=face.centre.y-.5*(b.ymin+b.ymax);
            values[offset+4]+=(x*(pressureFaces[id]*face.areaVector.y+dy)-y*(pressureFaces[id]*face.areaVector.x+dx))
                /(pressureScale*monitorLengths[g]*monitorLengths[g]);
        }
        for(double value:values)finite(value);
        return values;
}
}
