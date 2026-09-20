#include "cartmesh2d/fv/FlowInitialization2D.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
int main(){
    const FlowInitialVortex2D vortex{{2,-1},.8,.3};
    const auto peak=sampleInitialVortex2D({2+.8/std::sqrt(5.),-1},vortex);
    require(std::abs(peak.velocity.y-.3)<1e-15 && peak.velocity.x==0,"peak speed definition");
    const auto exterior=sampleInitialVortex2D({4,0},vortex);
    require(exterior.velocity.x==0 && exterior.velocity.y==0 && exterior.streamfunction==0,"compact support");
    const Point2D p{2.2,-.85};const double h=1e-6;
    const auto xp=sampleInitialVortex2D({p.x+h,p.y},vortex),xm=sampleInitialVortex2D({p.x-h,p.y},vortex);
    const auto yp=sampleInitialVortex2D({p.x,p.y+h},vortex),ym=sampleInitialVortex2D({p.x,p.y-h},vortex);
    const auto here=sampleInitialVortex2D(p,vortex);
    require(std::abs((yp.streamfunction-ym.streamfunction)/(2*h)-here.velocity.x)<1e-10,"u equals dpsi/dy");
    require(std::abs((xp.streamfunction-xm.streamfunction)/(2*h)+here.velocity.y)<1e-10,"v equals -dpsi/dx");
    require(std::abs((xp.velocity.x-xm.velocity.x+yp.velocity.y-ym.velocity.y)/(2*h))<1e-9,"analytic divergence");
    const double angle=.63,c=std::cos(angle),s=std::sin(angle);
    const Point2D rotated{2+c*(p.x-2)-s*(p.y+1),-1+s*(p.x-2)+c*(p.y+1)};
    const auto result=sampleInitialVortex2D(rotated,vortex);
    require(std::hypot(result.velocity.x-(c*here.velocity.x-s*here.velocity.y),
                       result.velocity.y-(s*here.velocity.x+c*here.velocity.y))<1e-15,"rotational covariance");
    for(int kind=0;kind<3;++kind){
        bool threw=false;auto bad=vortex;auto point=p;
        if(kind==0)bad.radius=0;
        if(kind==1)bad.peakSpeed=std::numeric_limits<double>::infinity();
        if(kind==2)point.x=std::numeric_limits<double>::quiet_NaN();
        try{(void)sampleInitialVortex2D(point,bad);}catch(const std::runtime_error&){threw=true;}
        require(threw,"invalid initial condition must fail");
    }
    std::cout<<"initial vortex analytic definition verified\n";
}
