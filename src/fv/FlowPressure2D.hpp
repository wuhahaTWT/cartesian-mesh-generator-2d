#pragma once
#include "FlowSolverDetail2D.hpp"
namespace cartmesh2d::fv::solver_detail {
template<class LinearSolve>
Vec solvePressureCorrection(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,
    FlowResult2D& r,System& ap,Vec& pc,const Vec& ra,const Vec& df,const Vec& predicted,
    const detail::FlowGradientStencil2D& pressureGradientStencil,const Vec& zeros,LinearSolve&& linearSolve){
    const auto nf=m.faces.size();
        ap.reset();
        for (std::size_t id=0;id<nf;++id) {
            const auto& f=m.faces[id]; const auto i=f.owner;
            if (f.neighbour) {
                const auto j=*f.neighbour;
                ap.diag[i]+=df[id]; ap.diag[j]+=df[id];
                ap.add(i,j,-df[id]); ap.add(j,i,-df[id]);
            } else if (b.fixedP[id]) ap.diag[i]+=df[id];
        }
        if (b.closed) ap.pin(0);
        std::fill(pc.begin(),pc.end(),0);Vec correction(nf);
        for(std::size_t pass=0;pass<c.pressureCorrectionPasses;++pass){
            std::fill(ap.rhs.begin(),ap.rhs.end(),0.);
            const auto gc=pressureGradientStencil.apply(pc,zeros);
            for(std::size_t id=0;id<nf;++id){const auto&f=m.faces[id];const auto i=f.owner;
                correction[id]=(f.neighbour||b.fixedP[id])?-interpolate(f,ra)*dot(interpolateGradient(f,gc),f.correction):0;
                ap.rhs[i]-=predicted[id]+correction[id];
                if(f.neighbour)ap.rhs[*f.neighbour]+=predicted[id]+correction[id];
            }
            if(b.closed){ap.rhs[0]=0;pc[0]=0;}
            if (linearSolve(ap, pc, true)==0) {
                // A zero-iteration solve does not modify pc. The next pass
                // would reconstruct the exact same gradient, correction, RHS
                // and matrix, and therefore return zero again. Keep this
                // pass's correction for the accepted flux and omit only those
                // identical repeats; no geometric approximation or new stop.
                if (c.profile) r.performance.pressureCorrectionPassesSkipped+=c.pressureCorrectionPasses-1-pass;
                break;
            }
        }
    return correction;
}
}
