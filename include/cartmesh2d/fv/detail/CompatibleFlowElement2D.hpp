#pragma once
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

// Shared local discretization on the actual validated 2D polygon mesh.
// There are no case names, manufactured fields, boundary policies, file I/O
// or global solver decisions in this kernel. The full solver remains separate.
namespace cartmesh2d::fv::detail::compatible {
using Vec=std::vector<double>;
struct DenseLU {
    std::size_t n;Vec a;std::vector<std::size_t> pivots;
    explicit DenseLU(Vec matrix,std::size_t size);
    Vec solve(Vec x)const;
};
struct Mat {
    std::size_t nr=0,nc=0;Vec v;
    Mat(std::size_t r,std::size_t c):nr(r),nc(c),v(r*c){}
    double& operator()(std::size_t r,std::size_t c){return v[r*nc+c];}
    double operator()(std::size_t r,std::size_t c)const{return v[r*nc+c];}
};
struct Q {Point2D p;double w;};
using P3=std::array<double,3>;
using P6=std::array<double,6>;
std::vector<std::pair<double,double>> gauss(int order);
std::vector<Q> cellQuadrature(const FvMesh2D& mesh,std::size_t cell,int order);
struct Basis {
    Point2D c;double h;P3 moments{};
    P3 phi(Point2D p)const{return {1,(p.x-c.x)/h,(p.y-c.y)/h};}
    P6 theta(Point2D p)const{const auto a=phi(p);return {1,a[1],a[2],a[1]*a[1]-moments[0],a[1]*a[2]-moments[1],a[2]*a[2]-moments[2]};}
    std::array<Vector2D,6> grad(Point2D p)const{const auto a=phi(p);return {{{0,0},{1/h,0},{0,1/h},{2*a[1]/h,0},{a[2]/h,a[1]/h},{0,2*a[2]/h}}};}
};

// P1 cell/face velocity, P1 pressure and P2 potential. diameter is the
// geometric coordinate scale from the original polygon, never a cut-cell
// area replacement. mesh must come from makeFvMesh2D's unchanged quality gate.
struct P1Local {
    Basis basis;std::size_t m;std::vector<Q> q;
    Mat mass,bx,by,gx,gy,potential,stiffness;
    P1Local(const FvMesh2D& mesh,std::size_t cell,double diameter,int order);
};
struct Tri {
    Point2D c,a,b; double h,det; Vector2D j0,j1;
    Tri(Point2D C,Point2D A,Point2D B,double H):c(C),a(A),b(B),h(H),j0{(A.x-C.x)/H,(A.y-C.y)/H},j1{(B.x-C.x)/H,(B.y-C.y)/H} {
        det=j0.x*j1.y-j0.y*j1.x;
        if(!(det>0))throw std::runtime_error("unsupported nonpositive centroid fan; geometry untouched");
    }
    Vector2D ref(Point2D p)const {auto d=p-c;return {(d.x*j1.y-d.y*j1.x)/(h*det),(j0.x*d.y-j0.y*d.x)/(h*det)};}
    std::array<Vector2D,8> shape(Point2D p)const {
        auto z=ref(p);std::array<Vector2D,8> values{{{1,0},{z.x,0},{z.y,0},{0,1},{0,z.x},{0,z.y},{z.x*z.x,z.x*z.y},{z.x*z.y,z.y*z.y}}};
        for(auto& v:values)v={(j0.x*v.x+j1.x*v.y)/det,(j0.y*v.x+j1.y*v.y)/det};
        return values;
    }
    std::vector<Q> quadrature(int order)const {
        std::vector<Q> out;auto g=gauss(order);
        for(auto [r,w]:g)for(auto [s,z]:g)out.push_back({{c.x+r*(a.x-c.x)+s*(1-r)*(b.x-c.x),c.y+r*(a.y-c.y)+s*(1-r)*(b.y-c.y)},w*z*(1-r)*h*h*det});
        return out;
    }
};

// RT1 normal-trace/divergence constraints use every original atomic face.
struct Lift {
    std::vector<Tri> tri;std::vector<std::size_t> faceLocal;Mat coefficients;
    double traceResidual=0,divResidual=0;
    Lift(const FvMesh2D& mesh,std::size_t cell,const P1Local& local,int order);
    Vector2D value(std::size_t triangle,std::size_t dof,Point2D point)const;
};

// Mutable local weak equations for subsequent convection/boundary terms.
// Call condense() after each matrix/RHS change and before recover/condensed.
// Retain all face velocities and cell mean pressure; eliminate six cell
// velocities and two pressure slopes. No pressure penalty or pivot fallback.
struct P1System {
    P1Local a;Mat matrix;Vec rhs;std::vector<std::size_t> inside,outside;
    Mat eliminated;Vec loadInternal;
    explicit P1System(P1Local local,bool symmetric=true);
    void condense();
    std::pair<Mat,Vec> condensed()const;
    Vec recover(const Vec& retained)const;
};
// General physical acceleration, integrated against the same H(div) lift.
// The caller supplies force; the kernel contains no analytic reference field.
Vec liftedBodyForce(const P1Local& local,const Lift& lift,int order,
                   const std::function<Vector2D(Point2D)>& force);
} // namespace cartmesh2d::fv::detail::compatible
