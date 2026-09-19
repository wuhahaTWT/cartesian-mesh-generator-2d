#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using cartmesh2d::fv::detail::AggregationHierarchy2D;
using cartmesh2d::fv::detail::LinearPressureMethod2D;
using cartmesh2d::fv::detail::LinearWorkspace2D;
using cartmesh2d::fv::detail::SparsePattern2D;
using cartmesh2d::fv::detail::SparseSystem2D;
int failures = 0;
void check(bool ok, const std::string& msg) { if (!ok) { ++failures; std::cerr << "FAIL: " << msg << '\n'; } }
template<class F> void rejects(F&& f, const std::string& msg) { try { f(); check(false, msg); } catch (const std::exception&) { check(true, msg); } }
double dot(const std::vector<double>& a, const std::vector<double>& b) { double r=0; for (std::size_t i=0;i<a.size();++i) r+=a[i]*b[i]; return r; }
double norm(const std::vector<double>& a) { return std::sqrt(dot(a,a)); }

struct Grid {
    static constexpr std::size_t w=17, n=w*w;
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    std::vector<double> d, h, v;
    SparsePattern2D pattern;
    static double hw(std::size_t r,std::size_t c) { return 1.+.01*static_cast<double>(r+2*c+1); }
    static double vw(std::size_t r,std::size_t c) { return 1.2+.008*static_cast<double>(2*r+c+1); }
    static std::vector<std::pair<std::size_t,std::size_t>> makeEdges() {
        std::vector<std::pair<std::size_t,std::size_t>> e;
        for (std::size_t r=0;r<w;++r) for (std::size_t c=0;c<w;++c) { const auto i=r*w+c;
            if (c+1<w) { e.push_back({i,i+1}); e.push_back({i+1,i}); }
            if (r+1<w) { e.push_back({i,i+w}); e.push_back({i+w,i}); }
        } return e;
    }
    Grid():edges(makeEdges()),d(n,1.),h((w-1)*w),v(w*(w-1)),pattern(n,edges) {
        for (std::size_t r=0;r<w;++r) for (std::size_t c=0;c<w;++c) { const auto i=r*w+c;
            if (c+1<w) { const auto z=hw(r,c); h[r*(w-1)+c]=z; d[i]+=z; d[i+1]+=z; }
            if (r+1<w) { const auto z=vw(r,c); v[r*w+c]=z; d[i]+=z; d[i+w]+=z; }
        }
    }
    std::vector<double> reference(const std::vector<double>& x) const {
        std::vector<double> y(n,0.);
        for (std::size_t r=0;r<w;++r) for (std::size_t c=0;c<w;++c) { const auto i=r*w+c; y[i]+=d[i]*x[i];
            if(c) y[i]-=h[r*(w-1)+c-1]*x[i-1]; if(c+1<w) y[i]-=h[r*(w-1)+c]*x[i+1];
            if(r) y[i]-=v[(r-1)*w+c]*x[i-w]; if(r+1<w) y[i]-=v[r*w+c]*x[i+w];
        } return y;
    }
    void assemble(SparseSystem2D& s) const { for(std::size_t i=0;i<n;++i)s.diag[i]=d[i];
        for(std::size_t r=0;r<w;++r) for(std::size_t c=0;c<w;++c) { const auto i=r*w+c;
            if(c+1<w){ const auto z=-h[r*(w-1)+c]; s.add(i,i+1,z); s.add(i+1,i,z); }
            if(r+1<w){ const auto z=-v[r*w+c]; s.add(i,i+w,z); s.add(i+w,i,z); }
        }
    }
};

void known(const Grid& g,const std::vector<double>& actual,const std::vector<double>& exact,const std::string& label) {
    const auto rhs=g.reference(exact), got=g.reference(actual); std::vector<double> e(actual.size()),r(actual.size());
    for(std::size_t i=0;i<actual.size();++i){e[i]=actual[i]-exact[i];r[i]=got[i]-rhs[i];}
    check(norm(r)<=1.05*(1e-13+1e-11*norm(rhs)),label+" reaches independent residual target");
    check(norm(e)<=5e-10*(1.+norm(exact)),label+" matches known solution");
}

void regression() {
    Grid g; SparseSystem2D s(g.pattern); g.assemble(s);
    AggregationHierarchy2D h(g.pattern.rows,g.pattern.columns,s.diag,s.off);
    check(h.levels()>3 && h.coarseCells()<=32,"17x17 grid exercises multiple bounded levels");
    std::vector<double> x(g.n),y(g.n); for(std::size_t i=0;i<g.n;++i){x[i]=std::sin(.17*(i+1));y[i]=std::cos(.11*(i+2));}
    std::vector<double> bx,by; h.apply(x,bx); h.apply(y,by); const auto scale=1.+std::abs(dot(x,by))+std::abs(dot(y,bx));
    check(std::abs(dot(x,by)-dot(y,bx))<=2e-11*scale,"V-cycle bilinear form is symmetric");
    check(dot(x,bx)>0.&&dot(y,by)>0.,"V-cycle quadratic form is positive");
    std::vector<double> sum=x, linearResult, expected(bx.size());
    for(std::size_t i=0;i<g.n;++i) { sum[i]+=y[i]; expected[i]=bx[i]+by[i]; }
    h.apply(sum,linearResult);
    double linearError=0.; for(std::size_t i=0;i<g.n;++i) linearError=std::max(linearError,std::abs(linearResult[i]-expected[i]));
    check(linearError<=2e-13,"V-cycle is linear for fixed coefficients");
    LinearWorkspace2D w(g.n); std::vector<double> sol(g.n); s.rhs=g.reference(x);
    (void)s.solvePressure(sol,w,LinearPressureMethod2D::Aggregation); known(g,sol,x,"first aggregation solve");
    check(s.hierarchyBuilds()==1&&s.hierarchyReuses()==0,"first solve builds one hierarchy");
    s.rhs=g.reference(y); sol.assign(g.n,0.); (void)s.solvePressure(sol,w,LinearPressureMethod2D::Aggregation);
    check(s.hierarchyBuilds()==1&&s.hierarchyReuses()==1,"changed RHS reuses hierarchy"); known(g,sol,y,"changed-RHS aggregation solve");
    s.diag[10]+=.25; g.d[10]+=.25; s.rhs=g.reference(y); sol.assign(g.n,0.);
    (void)s.solvePressure(sol,w,LinearPressureMethod2D::Aggregation); check(s.hierarchyBuilds()==2,"direct diagonal mutation rebuilds hierarchy"); known(g,sol,y,"diagonal mutation solve");
    const auto a=g.pattern.slot(10,9),at=g.pattern.slot(9,10); s.off[a]-=.1;s.off[at]-=.1;g.h[9]+=.1;
    s.rhs=g.reference(y);sol.assign(g.n,0.);(void)s.solvePressure(sol,w,LinearPressureMethod2D::Aggregation);check(s.hierarchyBuilds()==3,"direct off mutation rebuilds hierarchy");known(g,sol,y,"off mutation solve");
    s.off[a]+=.1; const auto before=s.hierarchyBuilds(); rejects([&]{std::vector<double> z(g.n);(void)s.solvePressure(z,w,LinearPressureMethod2D::Aggregation);},"invalid asymmetric hierarchy rejected");
    check(s.hierarchyBuilds()==before,"failed hierarchy build is not counted"); s.off[a]-=.1; s.rhs=g.reference(y); sol.assign(g.n,0.); (void)s.solvePressure(sol,w,LinearPressureMethod2D::Aggregation);
    check(s.hierarchyBuilds()==before+1,"repaired hierarchy rebuilds after failure"); known(g,sol,y,"repaired hierarchy solve");
    const std::vector<std::size_t> rows{0,1,2},cols{1,0};
    rejects([&]{AggregationHierarchy2D q(rows,cols,{0.,2.},{-1.,-1.});},"zero diagonal rejected");
    rejects([&]{AggregationHierarchy2D q(rows,cols,{std::numeric_limits<double>::infinity(),2.},{-1.,-1.});},"nonfinite diagonal rejected");
    rejects([&]{AggregationHierarchy2D q(rows,cols,{2.,2.},{1.,1.});},"positive off-diagonal rejected");
    rejects([&]{AggregationHierarchy2D q({0,2,1},{1},{2.,2.},{-1.});},"malformed CSR rows rejected");

    std::vector<std::size_t> diagonalRows(41,0);
    AggregationHierarchy2D diagonalHierarchy(diagonalRows,{},std::vector<double>(40,1.),{});
    check(diagonalHierarchy.levels()==1 && diagonalHierarchy.coarseCells()==40,
          "isolated diagonal matrix uses the exact diagonal terminal");

    std::vector<std::size_t> starRows(41,0), starColumns;
    std::vector<double> starOff;
    for(std::size_t row=0;row<40;++row) {
        if(row==0) { for(std::size_t leaf=1;leaf<40;++leaf) { starColumns.push_back(leaf); starOff.push_back(-1.); } }
        else { starColumns.push_back(0); starOff.push_back(-1.); }
        starRows[row+1]=starColumns.size();
    }
    // Minimal failure of pure pair matching: 40 -> 39 coarse vertices used
    // to hit the unchanged memory guard. Attaching connected singletons keeps
    // all matrix rows/couplings while reducing this star to one coarse group.
    std::vector<double> starDiagonal(40,2.); starDiagonal[0]=40.;
    AggregationHierarchy2D starHierarchy(starRows,starColumns,starDiagonal,starOff);
    check(starHierarchy.levels()==2 && starHierarchy.coarseCells()==1,
          "connected unmatched star vertices join the neighbouring aggregate");
    std::vector<std::pair<std::size_t,std::size_t>> starEdges;
    for (std::size_t leaf=1;leaf<40;++leaf) starEdges.emplace_back(0,leaf);
    SparsePattern2D starPattern(40,starEdges); SparseSystem2D star(starPattern);
    star.diag=starDiagonal; star.off=starOff;
    auto starReference=[](const std::vector<double>& z) {
        std::vector<double> value(40); value[0]=40*z[0];
        for(std::size_t i=1;i<40;++i) {value[0]-=z[i]; value[i]=2*z[i]-z[0];}
        return value;
    };
    std::vector<double> exactStar(40),solutionStar(40);
    for(std::size_t i=0;i<40;++i) exactStar[i]=std::sin(.2*(i+1));
    star.rhs=starReference(exactStar); LinearWorkspace2D starWorkspace(40);
    (void)star.solvePressure(solutionStar,starWorkspace,LinearPressureMethod2D::Aggregation);
    auto starResidual=starReference(solutionStar);
    for(std::size_t i=0;i<40;++i) starResidual[i]-=star.rhs[i];
    check(norm(starResidual)<=1.05*(1e-13+1e-11*norm(star.rhs)),
          "repaired star hierarchy reaches the independent true-residual target");
    std::vector<double> diagonalRhs(40,3.),diagonalResult;
    diagonalHierarchy.apply(diagonalRhs,diagonalResult);
    check(diagonalResult==diagonalRhs,"isolated diagonal terminal solves exactly");
    std::vector<std::size_t> weakRows(41,2);weakRows[0]=0;weakRows[1]=1;
    rejects([&]{AggregationHierarchy2D q(weakRows,{1,0},std::vector<double>(40,2.),{-1.,-1.});},
            "mostly isolated graph retains the explicit hierarchy memory guard");

}
}
int main(){try{regression();}catch(const std::exception& e){std::cerr<<"UNEXPECTED EXCEPTION: "<<e.what()<<'\n';return 1;}std::cout<<"flow aggregation regression failures: "<<failures<<'\n';return failures?1:0;}
