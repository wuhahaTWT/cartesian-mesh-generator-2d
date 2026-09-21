#pragma once
#include "cartmesh2d/fv/Euler2D.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace cartmesh2d::fv {
namespace euler_checkpoint_detail {
inline std::string binding(const FvMesh2D& mesh,const std::vector<EulerBoundary2D>& boundaries,
                           const IdealGas2D& gas,const std::string& problem) {
    validateFvMesh2D(mesh);validateEulerBoundaries2D(mesh,boundaries,gas);
    if(problem.find_first_of("\r\n")!=std::string::npos)throw std::runtime_error("Euler checkpoint: invalid problem label");
    std::ostringstream out;out<<std::setprecision(17)<<"CM2D_EULER_CHECKPOINT 1\nGAS "<<gas.gamma<<' '<<gas.gasConstant
        <<"\nPROBLEM "<<std::quoted(problem)<<"\nCELLS "<<mesh.cells.size()<<'\n';
    for(const auto& cell:mesh.cells) {
        out<<cell.centre.x<<' '<<cell.centre.y<<' '<<cell.area<<' '<<cell.faces.size();
        for(auto face:cell.faces)out<<' '<<face;
        out<<'\n';
    }
    out<<"FACES "<<mesh.faces.size()<<'\n';
    for(const auto& f:mesh.faces) {
        out<<f.owner<<' ';if(f.neighbour)out<<*f.neighbour;else out<<'-';
        out<<' '<<static_cast<int>(f.patch)<<' '<<f.centre.x<<' '<<f.centre.y<<' '<<f.areaVector.x<<' '
            <<f.areaVector.y<<' '<<f.transmissibility<<' '<<f.neighbourWeight<<' '<<f.correction.x<<' '<<f.correction.y<<'\n';
    }
    std::vector<EulerBoundary2D> ordered=boundaries;
    std::sort(ordered.begin(),ordered.end(),[](const auto& a,const auto& b){return a.face<b.face;});
    out<<"BOUNDARIES "<<ordered.size()<<'\n';
    for(const auto& b:ordered) {
        if(b.name.find_first_of("\r\n")!=std::string::npos)throw std::runtime_error("Euler checkpoint: invalid boundary label");
        out<<b.face<<' '<<static_cast<int>(b.kind)<<' '<<std::quoted(b.name)<<' ';
        if(b.partner)out<<*b.partner;else out<<'-';
        out<<' '<<b.reference.density<<' '<<b.reference.u<<' '<<b.reference.v<<' '<<b.reference.pressure<<'\n';
    }
    return out.str();
}
inline void validate(const EulerState2D& state,std::size_t count,const IdealGas2D& gas) {
    if(!std::isfinite(state.time)||state.time<0||state.cells.size()!=count||state.steps==std::numeric_limits<std::size_t>::max())
        throw std::runtime_error("Euler checkpoint: invalid state size, time or step counter");
    for(const auto& q:state.cells)(void)eulerPrimitive2D(q,gas);
}
}
// Full geometry and physical boundary binding, not a trusted summary/hash.
// Numerical dt/CFL and requested final time may change on restart.
inline void writeEulerCheckpoint2D(std::ostream& out,const FvMesh2D& mesh,
    const std::vector<EulerBoundary2D>& boundaries,const IdealGas2D& gas,const EulerState2D& state,
    const std::string& problem="") {
    euler_checkpoint_detail::validate(state,mesh.cells.size(),gas);
    std::ostringstream data;data<<std::setprecision(17)<<euler_checkpoint_detail::binding(mesh,boundaries,gas,problem)
        <<"STATE "<<state.time<<' '<<state.steps<<'\n';
    for(const auto& q:state.cells)data<<q[0]<<' '<<q[1]<<' '<<q[2]<<' '<<q[3]<<'\n';
    data<<"END\n";out<<data.str();
    if(!out)throw std::runtime_error("Euler checkpoint: cannot write state");
}
inline EulerState2D readEulerCheckpoint2D(std::istream& in,const FvMesh2D& mesh,
    const std::vector<EulerBoundary2D>& boundaries,const IdealGas2D& gas,const std::string& problem="") {
    std::istringstream expected(euler_checkpoint_detail::binding(mesh,boundaries,gas,problem));
    std::string line,actual;
    while(std::getline(expected,line))if(!std::getline(in,actual)||actual!=line)
        throw std::runtime_error("Euler checkpoint: incompatible geometry, gas, boundary or problem");
    EulerState2D state;std::string token,steps;
    if(!(in>>token>>state.time>>steps)||token!="STATE"||steps.empty()||steps.find_first_not_of("0123456789")!=std::string::npos)
        throw std::runtime_error("Euler checkpoint: invalid state header");
    const auto count=std::stoull(steps);
    if(count>=std::numeric_limits<std::size_t>::max())throw std::runtime_error("Euler checkpoint: step counter overflow");
    state.steps=static_cast<std::size_t>(count);state.cells.resize(mesh.cells.size());
    for(auto& q:state.cells)for(double& value:q)if(!(in>>value))throw std::runtime_error("Euler checkpoint: truncated state");
    if(!(in>>token)||token!="END"||(in>>token))throw std::runtime_error("Euler checkpoint: invalid terminator or trailing data");
    euler_checkpoint_detail::validate(state,mesh.cells.size(),gas);return state;
}
}
