#pragma once

#include "cartmesh2d/fv/Incompressible2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh2d::fv::detail {

// Project-specific engineering criterion, inspired by residual reduction and
// physical-monitor stability guidance, not Fluent-equivalent residuals.
// The max-cell residual cap remains mandatory even after a large reduction.
struct EngineeringWindow2D {
    static constexpr std::size_t length = 50;
    static constexpr double stability = 1e-3;
    double residualMax = 0, velocityPath = 0, pressurePath = 0, monitorRange = 0;
    bool accepted = false;
};
inline EngineeringWindow2D engineeringWindow2D(const std::vector<FlowIteration2D>& history,
                                               double reference, double tolerance) {
    EngineeringWindow2D result;
    if(history.size()<EngineeringWindow2D::length || !(reference>=0) || !std::isfinite(reference)
       || !(tolerance>0) || !std::isfinite(tolerance))return result;
    const auto start=history.size()-EngineeringWindow2D::length;
    const auto count=history.back().monitors.size();
    if(count==0)return result;
    auto low=history[start].monitors, high=low;
    if(low.size()!=count)return result;
    bool valid=true;
    for(std::size_t i=start;i<history.size();++i) {
        const auto& h=history[i];
        if(h.monitors.size()!=count)return result;
        if(i>start && h.iteration!=history[i-1].iteration+1)return result;
        for(double x:{h.momentumResidual,h.velocityChange,h.pressureChange,h.continuity,h.globalRelativeImbalance})
            valid=valid && std::isfinite(x) && x>=0;
        result.residualMax=std::max(result.residualMax,h.momentumResidual);
        result.velocityPath+=h.velocityChange;
        result.pressurePath+=h.pressureChange;
        valid=valid && h.continuity<1e-8 && h.globalRelativeImbalance<1e-8;
        for(std::size_t j=0;j<count;++j) {
            valid=valid && std::isfinite(h.monitors[j]);
            low[j]=std::min(low[j],h.monitors[j]);high[j]=std::max(high[j],h.monitors[j]);
        }
    }
    for(std::size_t j=0;j<count;++j)result.monitorRange=std::max(result.monitorRange,high[j]-low[j]);
    const bool reduced=result.residualMax<=reference*1e-3;
    const bool smallAbsolute=result.residualMax<std::min(tolerance,1e-5);
    result.accepted=valid && result.residualMax<tolerance && (reduced || smallAbsolute)
        && result.velocityPath<EngineeringWindow2D::stability
        && result.pressurePath<EngineeringWindow2D::stability
        && result.monitorRange<EngineeringWindow2D::stability;
    return result;
}

// Bounded forcing sequence. It tightens as the best nonlinear defect falls;
// prolonged stagnation tightens it further. An eventual stopping candidate
// must be repeated using the original strict linear tolerances.
struct AdaptiveLinear2D {
    double best=1, multiplier=1;
    std::size_t stagnant=0;
    void observe(double defect) {
        if(!std::isfinite(defect) || defect<0)throw std::runtime_error("Invalid adaptive linear defect");
        if(defect<best*.95)stagnant=0;
        else if(++stagnant>=20){multiplier=std::max(1e-8,multiplier*.1);stagnant=0;}
        best=std::min(best,defect);
    }
    double relative(bool strict) const {return strict?1e-11:std::clamp(.01*best*multiplier,1e-11,1e-3);}
    double row(double original,bool strict,double scale) const {
        return strict?original:std::max(original,std::min(1e-3,.01*best*multiplier)*scale);
    }
};
} // namespace cartmesh2d::fv::detail
