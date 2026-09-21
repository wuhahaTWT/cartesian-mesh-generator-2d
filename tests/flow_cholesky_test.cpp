#include <cstdlib>
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include <iostream>
using namespace cartmesh2d::fv::detail;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F&& f){bool rejected=false;try{f();}catch(const std::runtime_error&){rejected=true;}require(rejected,"Invalid Cholesky input accepted");}
int main()try {
#ifdef __APPLE__
    setenv("VECLIB_MAXIMUM_THREADS", "1", 1);
#endif
    constexpr std::size_t nx=13,ny=9,n=nx*ny;
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    for(std::size_t i=0;i<n;++i){if(i%nx+1<nx)edges.emplace_back(i,i+1);if(i+nx<n)edges.emplace_back(i,i+nx);}
    SparsePattern2D pattern(n,edges);SparseSystem2D system(pattern);LinearWorkspace2D work(n);
    for(auto [i,j]:edges){const double value=j==i+1?.1+static_cast<double>(i%7):.01;system.off[pattern.slot(i,j)]=system.off[pattern.slot(j,i)]=-value;system.diag[i]+=value;system.diag[j]+=value;}
    for(std::size_t i=0;i<n;++i)if(i%nx==nx-1)system.diag[i]+=.25;
    LinearVector2D expected(n),x(n),ax(n);
    for(std::size_t i=0;i<n;++i)expected[i]=std::sin(.13*static_cast<double>(i));
    system.apply(expected,system.rhs);
    if(!systemCholeskyAvailable2D()) {
        rejects([&]{system.solvePressure(x,work,LinearPressureMethod2D::SystemCholesky);});
        system.rhs.assign(n,0);rejects([&]{system.solvePressure(x,work,LinearPressureMethod2D::SystemCholesky);});
        std::cout<<"Unsupported system Cholesky is explicitly refused, including zero RHS\n";return 0;
    }
    auto solve=[&](SparseSystem2D& a) {
        a.apply(expected,a.rhs);std::fill(x.begin(),x.end(),0);
        const auto steps=a.solvePressure(x,work,LinearPressureMethod2D::SystemCholesky);
        a.apply(x,ax);double error=0;for(std::size_t i=0;i<n;++i){ax[i]=a.rhs[i]-ax[i];error=std::max(error,std::abs(x[i]-expected[i]));}
        require(steps<=2 && error<1e-10,"Known solution was not recovered");
        require(linearNorm(ax)<=1e-13+1e-11*linearNorm(a.rhs),"Original true-residual gate failed");
    };
    solve(system);require(system.choleskyBuilds()==1,"Initial factor not recorded");
    for(auto& value:expected)value*=.7;solve(system);require(system.choleskyReuses()==1,"Exact coefficient reuse missing");
    system.diag[17]+=.03;solve(system);require(system.choleskyRefactors()==1,"Changed coefficients did not refactor");
    auto copy=system;copy.diag[5]+=.01;solve(copy);
    require(copy.choleskyBuilds()==2,"Copied matrix mutated a shared factor");
    solve(system);require(system.choleskyBuilds()==1,"Original shared factor was invalidated by copy");
    const auto off=system.off;system.off[0]*=1.1;std::fill(x.begin(),x.end(),0);
    rejects([&]{system.solvePressure(x,work,LinearPressureMethod2D::SystemCholesky);});
    system.off=off;solve(system);require(system.choleskyBuilds()==2,"Failed refactor cache was reused");
    rejects([&]{SystemCholesky2D invalid({0,1,1},{1},{1,1},{-.1});});
    rejects([&]{SystemCholesky2D invalid({0,1,2},{1,0},{1,1},{2,2});});
    SystemCholesky2D factor(pattern.rows,pattern.columns,system.diag,system.off);
    auto bad=expected;bad[0]=std::numeric_limits<double>::quiet_NaN();rejects([&]{factor.apply(bad,x);});
    bad.pop_back();rejects([&]{factor.apply(bad,x);});
    std::cout<<"Cholesky exact solutions, true residual, reuse, refactor, copy isolation and invalid-input rejection pass\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
