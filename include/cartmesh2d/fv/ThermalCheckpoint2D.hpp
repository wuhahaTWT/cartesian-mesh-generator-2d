#pragma once
#include "cartmesh2d/fv/ThermalFlow2D.hpp"
#include "cartmesh2d/fv/FlowCheckpoint2D.hpp"
#include <istream>
#include <ostream>

namespace cartmesh2d::fv {
// One stream bundles scalar, carrier, clock, full geometry and physical setup.
// The caller writes to a temporary file then atomically replaces the accepted
// checkpoint. A separate carrier file is diagnostic, never the restart record.
inline void writeThermalCheckpoint2D(std::ostream& out,const FvMesh2D& mesh,
    const FlowControls2D& flow,const ThermalSetup2D& setup,
    const ScalarTransportControls2D& controls,const ThermalFlowState2D& state) {
    validateThermalSetup2D(mesh,setup);
    if (state.scalar.size()!=mesh.cells.size())
        flow_checkpoint_detail::fail("thermal scalar dimensions invalid");
    for (double s:state.scalar) flow_checkpoint_detail::finite(s,"thermal scalar");
    struct RestoreFormat {
        std::ostream& s; std::ios_base::fmtflags flags; std::streamsize precision;
        ~RestoreFormat(){s.flags(flags);s.precision(precision);}
    } restore{out,out.flags(),out.precision()};
    out<<std::defaultfloat<<std::dec<<std::noshowpos<<std::noshowbase<<std::setprecision(17);
    out<<"CARTMESH2D_THERMAL_CHECKPOINT 1\nCOUPLING new-time-flux-Euler-v1\nTHERMAL_CONFIG "
       <<setup.diffusivity<<' '<<flow_checkpoint_detail::convectionName(controls.convection)<<'\n';
    out<<"SOURCES "<<setup.sourceDensity.size();
    for (double s:setup.sourceDensity) out<<' '<<s;
    std::size_t count=0;for(const auto& face:mesh.faces) if(!face.neighbour) ++count;
    out<<"\nBOUNDARIES "<<count<<'\n';
    for (std::size_t id=0;id<mesh.faces.size();++id) {
        if (mesh.faces[id].neighbour) continue;
        const auto& b=setup.boundary[id];
        out<<"BC "<<id<<' '<<(b.kind==ScalarBoundaryKind2D::Value?"value":"flux")<<' '
           <<b.value<<' '<<(b.inflowValue?1:0);
        if (b.inflowValue) out<<' '<<*b.inflowValue;
        out<<'\n';
    }
    out<<"SCALAR "<<state.scalar.size();for(double s:state.scalar)out<<' '<<s;
    out<<"\nFLOW\n";
    writeFlowCheckpoint2D(out,mesh,flow,state.flow);
    if(!out)flow_checkpoint_detail::fail("thermal write failed");
}
inline ThermalFlowState2D readThermalCheckpoint2D(std::istream& in,const FvMesh2D& mesh,
    const FlowControls2D& flow,const ThermalSetup2D& setup,const ScalarTransportControls2D& controls) {
    validateThermalSetup2D(mesh,setup);
    using namespace flow_checkpoint_detail;
    token(in,"CARTMESH2D_THERMAL_CHECKPOINT");token(in,"1");
    token(in,"COUPLING");token(in,"new-time-flux-Euler-v1");token(in,"THERMAL_CONFIG");
    double diffusivity=0;std::string scheme;
    if(!(in>>diffusivity>>scheme))fail("truncated thermal configuration");
    exact(diffusivity,setup.diffusivity,"thermal diffusivity");
    if(scheme!=convectionName(controls.convection))fail("thermal convection mismatch");
    token(in,"SOURCES");count(in,setup.sourceDensity.size(),"thermal source");
    for(double expected:setup.sourceDensity) {
        double s=0;if(!(in>>s))fail("truncated thermal sources");exact(s,expected,"thermal source");
    }
    std::size_t boundaryCount=0;for(const auto& f:mesh.faces)if(!f.neighbour)++boundaryCount;
    token(in,"BOUNDARIES");count(in,boundaryCount,"thermal boundary");
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        if(mesh.faces[id].neighbour)continue;
        token(in,"BC");std::size_t index=0;int hasInflow=0;double value=0;std::string kind;
        if(!(in>>index>>kind>>value>>hasInflow))fail("truncated thermal boundary");
        exact(index,id,"thermal boundary face");
        const auto& expected=setup.boundary[id];
        if(kind!=(expected.kind==ScalarBoundaryKind2D::Value?"value":"flux") ||
           hasInflow!=(expected.inflowValue?1:0))fail("thermal boundary type/inflow mismatch");
        exact(value,expected.value,"thermal boundary value");
        if(hasInflow) {double v=0;if(!(in>>v))fail("truncated thermal inflow");exact(v,*expected.inflowValue,"thermal inflow");}
    }
    ThermalFlowState2D state;
    token(in,"SCALAR");count(in,mesh.cells.size(),"thermal scalar");state.scalar.resize(mesh.cells.size());
    for(double& s:state.scalar) {
        if(!(in>>s))fail("truncated thermal scalar");finite(s,"thermal scalar");
    }
    token(in,"FLOW");state.flow=readFlowCheckpoint2D(in,mesh,flow);
    return state;
}
}
