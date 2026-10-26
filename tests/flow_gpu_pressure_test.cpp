#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include <iostream>

using namespace cartmesh2d::fv::detail;
namespace {
void require(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
double residual(SparseSystem2D& matrix,const LinearVector2D& x) {
    LinearVector2D ax;matrix.apply(x,ax);
    for(std::size_t i=0;i<x.size();++i)ax[i]=matrix.rhs[i]-ax[i];
    return linearNorm(ax);
}
}
int main() {
    try {
        constexpr std::size_t n=64;
        std::vector<std::pair<std::size_t,std::size_t>> edges;
        for(std::size_t i=1;i<n;++i)edges.emplace_back(i-1,i);
        SparsePattern2D graph(n,edges);SparseSystem2D matrix(graph);
        matrix.diag.assign(n,4.);matrix.off.assign(graph.columns.size(),-1.);
        LinearVector2D exact(n);for(std::size_t i=0;i<n;++i)exact[i]=std::sin(.27*static_cast<double>(i))+.1;
        matrix.apply(exact,matrix.rhs);
        LinearWorkspace2D workspace(n);
        for(auto method:{LinearPressureMethod2D::Jacobi,LinearPressureMethod2D::IC0}) {
            LinearVector2D cpu(n),gpu(n);
            matrix.solvePressure(cpu,workspace,method);
            matrix.solvePressure(gpu,workspace,method,1e-11,true);
            const auto cpuResidual=residual(matrix,cpu),gpuResidual=residual(matrix,gpu);
            const double stop=1e-13+1e-11*linearNorm(matrix.rhs);
            require(cpuResidual<=stop && gpuResidual<=stop,"GPU differs from original residual acceptance");
            LinearVector2D delta(n);for(std::size_t i=0;i<n;++i)delta[i]=cpu[i]-gpu[i];
            // diag=4, at most two off-diagonals=-1 => lambda_min >= 2.
            // Compare to the residual-derived solution bound plus a 64-ulp
            // arithmetic allowance on this finite, explicitly scaled test.
            require(linearNorm(delta)<=(cpuResidual+gpuResidual)/2+
                64*std::numeric_limits<double>::epsilon()*(1+linearNorm(cpu)),"GPU solution disagrees beyond residual/roundoff bound");
            const auto uploads=matrix.gpuStats().matrixUploads;
            std::fill(gpu.begin(),gpu.end(),0.);
            matrix.solvePressure(gpu,workspace,method,1e-11,true);
            require(matrix.gpuStats().matrixUploads==uploads,"RHS-only replay rebuilt GPU matrix");
            const auto saved=matrix.rhs;matrix.rhs[0]+=.1;
            std::fill(gpu.begin(),gpu.end(),0.);matrix.solvePressure(gpu,workspace,method,1e-11,true);
            require(matrix.gpuStats().matrixUploads==uploads && residual(matrix,gpu)<=1e-13+1e-11*linearNorm(matrix.rhs),
                    "GPU RHS-only update was stale");matrix.rhs=saved;
        }
        const auto uploads=matrix.gpuStats().matrixUploads;matrix.diag[5]+=.5;
        LinearVector2D changed(n);matrix.solvePressure(changed,workspace,LinearPressureMethod2D::IC0,1e-11,true);
        require(matrix.gpuStats().matrixUploads==uploads+1,"GPU coefficient update reused stale matrix/factors");
        bool rejected=false;
        try{matrix.solvePressure(changed,workspace,LinearPressureMethod2D::Aggregation,1e-11,true);}
        catch(const std::exception&){rejected=true;}
        require(rejected,"Unsupported GPU preconditioner silently fell back");
        GpuPressure2D backend;GpuPressureStats2D stats;LinearVector2D trial(n),empty;
        const auto original=trial;
        rejected=false;
        try {backend.solve(graph.rows,graph.columns,matrix.diag,matrix.off,matrix.rhs,trial,empty,empty,false,
                           1e-13+1e-11*linearNorm(matrix.rhs),1,stats);}
        catch(const std::exception&){rejected=true;}
        require(rejected && trial==original,"GPU budget failure changed retained input");
        backend.solve(graph.rows,graph.columns,matrix.diag,matrix.off,matrix.rhs,trial,empty,empty,false,
                      1e-13+1e-11*linearNorm(matrix.rhs),3000,stats);
        require(residual(matrix,trial)<=1e-13+1e-11*linearNorm(matrix.rhs),"Failed GPU context was not reusable");
        std::cout<<"GPU FP64 Jacobi/IC0, native residual gates, cache updates, rejection and budget recovery passed\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
