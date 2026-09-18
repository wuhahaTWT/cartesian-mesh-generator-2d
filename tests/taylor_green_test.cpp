#include "cartmesh2d/fv/TaylorGreen2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

void check(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
int main() {
    try {
        for (double nu : {.01,.1}) for (double speed : {.3,1.}) {
            constexpr double t=.23,h=1e-4;
            for (Point2D x : {Point2D{.17,.29},Point2D{.63,.78},Point2D{.41,.51}}) {
                auto q=taylorGreen2D(x,t,speed,nu);
                auto xp=taylorGreen2D({x.x+h,x.y},t,speed,nu),xm=taylorGreen2D({x.x-h,x.y},t,speed,nu);
                auto yp=taylorGreen2D({x.x,x.y+h},t,speed,nu),ym=taylorGreen2D({x.x,x.y-h},t,speed,nu);
                auto tp=taylorGreen2D(x,t+h,speed,nu),tm=taylorGreen2D(x,t-h,speed,nu);
                double ux=(xp.velocity.x-xm.velocity.x)/(2*h),uy=(yp.velocity.x-ym.velocity.x)/(2*h);
                double vx=(xp.velocity.y-xm.velocity.y)/(2*h),vy=(yp.velocity.y-ym.velocity.y)/(2*h);
                double rx=(tp.velocity.x-tm.velocity.x)/(2*h)+q.velocity.x*ux+q.velocity.y*uy
                    +(xp.pressure-xm.pressure)/(2*h)-nu*(xp.velocity.x+xm.velocity.x+yp.velocity.x+ym.velocity.x-4*q.velocity.x)/(h*h);
                double ry=(tp.velocity.y-tm.velocity.y)/(2*h)+q.velocity.x*vx+q.velocity.y*vy
                    +(yp.pressure-ym.pressure)/(2*h)-nu*(xp.velocity.y+xm.velocity.y+yp.velocity.y+ym.velocity.y-4*q.velocity.y)/(h*h);
                check(std::hypot(rx,ry)<2e-6,"unforced Navier-Stokes derivative balance");
                check(std::abs(ux+vy)<1e-10,"analytic incompressibility");
                check(std::abs((yp.streamfunction-ym.streamfunction)/(2*h)-q.velocity.x)<1e-7,"streamfunction u sign");
                check(std::abs(-(xp.streamfunction-xm.streamfunction)/(2*h)-q.velocity.y)<1e-7,"streamfunction v sign");
            }
            double energy=0;
            for (int j=0;j<32;++j) for (int i=0;i<32;++i) {
                auto q=taylorGreen2D({(i+.5)/32,(j+.5)/32},t,speed,nu);
                energy+=.5*(q.velocity.x*q.velocity.x+q.velocity.y*q.velocity.y)/(32*32);
            }
            const double expected=.25*speed*speed*std::exp(-4*nu*std::numbers::pi*std::numbers::pi*t);
            check(std::abs(energy-expected)<1e-14,"kinetic energy decay integral");
            for (double a : {.0,.3,.7,1.}) for (double wall : {0.,1.}) {
                check(std::abs(taylorGreen2D({wall,a},t,speed,nu).velocity.x)<1e-14,"x wall impermeability");
                check(std::abs(taylorGreen2D({a,wall},t,speed,nu).velocity.y)<1e-14,"y wall impermeability");
                const double shearX=(taylorGreen2D({wall+h,a},t,speed,nu).velocity.y-taylorGreen2D({wall-h,a},t,speed,nu).velocity.y)/(2*h);
                const double shearY=(taylorGreen2D({a,wall+h},t,speed,nu).velocity.x-taylorGreen2D({a,wall-h},t,speed,nu).velocity.x)/(2*h);
                check(std::abs(shearX)+std::abs(shearY)<1e-10,"free slip tangential normal derivative");
            }
        }
        std::cout << "Taylor-Green unforced equations, slip walls, energy and streamfunction verified\n";
        return 0;
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
