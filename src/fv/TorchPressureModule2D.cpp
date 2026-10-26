// Optional Windows /MD module. All arithmetic uses existing precompiled ATen
// CUDA operators; no CUDA compiler, copied Torch installation or Python worker.
#include "FlowGpuPressureAbi2D.hpp"
#include <ATen/ATen.h>
#include <ATen/ops/_sparse_mm.h>
#include <ATen/ops/triangular_solve.h>
#include <c10/core/InferenceMode.h>
#include <windows.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using Clock=std::chrono::steady_clock;
using Tensor=at::Tensor;
double seconds(Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();}
void require(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
void errorText(char* output,std::uint64_t capacity,const char* text) {
    if(!output || !capacity)return;
    const auto size=std::min<std::size_t>(std::strlen(text),static_cast<std::size_t>(capacity-1));
    std::memcpy(output,text,size);output[size]='\0';
}
template<class T> Tensor device(const T* data,std::size_t size,at::ScalarType dtype) {
    return at::from_blob(const_cast<T*>(data),{static_cast<std::int64_t>(size)},
        at::TensorOptions().dtype(dtype)).to(at::kCUDA);
}
Tensor product(const Tensor& matrix,const Tensor& x){return at::_sparse_mm(matrix,x.unsqueeze(1)).squeeze(1);}
double gpuNorm(const Tensor& x){return at::linalg_vector_norm(x).item<double>();}

struct Context {
    std::vector<std::uint64_t> oldRows,oldColumns;
    std::vector<double> oldDiagonal,oldOff;
    Tensor matrix,diagonal,lower,factorDiagonal;
    bool ic0=false;

    void prepare(std::uint64_t n,const std::uint64_t* rows,const std::uint64_t* columns,
        const double* d,const double* off,const double* fd,const double* fl,bool useIc0,
        CartmeshGpuPressureStatsV1& stats) {
        const auto nnz=static_cast<std::size_t>(rows[n]);
        require(n>0 && n<static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()),"GPU matrix size unsupported");
        if(oldDiagonal.size()==n && oldOff.size()==nnz && ic0==useIc0 &&
           std::equal(oldRows.begin(),oldRows.end(),rows) &&
           std::equal(oldColumns.begin(),oldColumns.end(),columns) &&
           std::equal(oldDiagonal.begin(),oldDiagonal.end(),d) && std::equal(oldOff.begin(),oldOff.end(),off))return;
        // A failed upload/factor construction must not authorize a stale hit.
        oldDiagonal.clear();
        const auto start=Clock::now();
        std::vector<std::int64_t> fullRows{0},fullColumns,lowerRows{0},lowerColumns;
        std::vector<double> values,lowerValues;
        fullColumns.reserve(nnz+static_cast<std::size_t>(n));values.reserve(fullColumns.capacity());
        for(std::uint64_t i=0;i<n;++i) {
            require(rows[i]<=rows[i+1] && rows[i+1]<=nnz && std::isfinite(d[i]) && d[i]>0,"GPU pressure matrix invalid");
            bool added=false;
            for(auto k=rows[i];k<rows[i+1];++k) {
                const auto j=columns[k];
                require(j<n && j!=i && std::isfinite(off[k]),"GPU pressure sparse entry invalid");
                if(!added && j>i){fullColumns.push_back(static_cast<std::int64_t>(i));values.push_back(d[i]);added=true;}
                fullColumns.push_back(static_cast<std::int64_t>(j));values.push_back(off[k]);
                if(useIc0 && j<i){lowerColumns.push_back(static_cast<std::int64_t>(j));lowerValues.push_back(fl[k]);}
            }
            if(!added){fullColumns.push_back(static_cast<std::int64_t>(i));values.push_back(d[i]);}
            fullRows.push_back(static_cast<std::int64_t>(values.size()));
            if(useIc0){lowerColumns.push_back(static_cast<std::int64_t>(i));lowerValues.push_back(1.);lowerRows.push_back(static_cast<std::int64_t>(lowerValues.size()));}
        }
        const auto options=at::TensorOptions().dtype(at::kDouble).device(at::kCUDA);
        matrix=at::sparse_csr_tensor(device(fullRows.data(),fullRows.size(),at::kLong),
            device(fullColumns.data(),fullColumns.size(),at::kLong),device(values.data(),values.size(),at::kDouble),
            {static_cast<std::int64_t>(n),static_cast<std::int64_t>(n)},options);
        diagonal=device(d,static_cast<std::size_t>(n),at::kDouble);
        if(useIc0) {
            lower=at::sparse_csr_tensor(device(lowerRows.data(),lowerRows.size(),at::kLong),
                device(lowerColumns.data(),lowerColumns.size(),at::kLong),device(lowerValues.data(),lowerValues.size(),at::kDouble),
                {static_cast<std::int64_t>(n),static_cast<std::int64_t>(n)},options);
            factorDiagonal=device(fd,static_cast<std::size_t>(n),at::kDouble);
        } else {lower=Tensor{};factorDiagonal=Tensor{};}
        oldRows.assign(rows,rows+n+1);oldColumns.assign(columns,columns+nnz);
        oldDiagonal.assign(d,d+n);oldOff.assign(off,off+nnz);ic0=useIc0;
        // Blocking copy synchronizes asynchronous matrix construction before timing.
        (void)diagonal[0].item<double>();
        ++stats.matrixUploads;stats.uploadSeconds+=seconds(start);
    }
    Tensor precondition(const Tensor& r) const {
        if(!ic0)return r/diagonal;
        auto y=std::get<0>(at::triangular_solve(r.unsqueeze(1),lower,false,false,true));
        y=y/factorDiagonal.unsqueeze(1);
        return std::get<0>(at::triangular_solve(y,lower,false,true,true)).squeeze(1);
    }
};
}

extern "C" __declspec(dllexport) int cartmesh2d_pressure_create_v1(void** result,char* error,std::uint64_t capacity) {
    try {
        require(result!=nullptr,"GPU context output is null");*result=nullptr;
        // Keep CUDA's registration library loaded for the process lifetime.
        // The native loader has already added the selected existing Torch directory.
        static HMODULE cuda=LoadLibraryExW(L"torch_cuda.dll",nullptr,LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
        require(cuda!=nullptr,"Cannot load existing torch_cuda.dll or a dependency");
        c10::InferenceMode inference;
        (void)at::ones({1},at::TensorOptions().dtype(at::kDouble).device(at::kCUDA)).item<double>();
        *result=new Context;return 0;
    }catch(const std::exception& e){errorText(error,capacity,e.what());return 1;}
}
extern "C" __declspec(dllexport) void cartmesh2d_pressure_destroy_v1(void* context) {delete static_cast<Context*>(context);}

extern "C" __declspec(dllexport) int cartmesh2d_pressure_solve_v1(void* pointer,std::uint64_t n,
    const std::uint64_t* rows,const std::uint64_t* columns,const double* diagonal,const double* off,
    const double* rhs,double* output,const double* fd,const double* fl,int ic0,double stop,
    std::uint64_t maximumIterations,CartmeshGpuPressureStatsV1* stats,char* error,std::uint64_t capacity) {
    try {
        require(pointer && stats && rows && diagonal && rhs && output,"GPU pressure buffer missing");
        require(std::isfinite(stop) && stop>0 && maximumIterations>0,"GPU stopping control invalid");
        require(ic0==0 || (ic0==1 && fd && (rows[n]==0 || fl)),"GPU IC0 factors missing");
        c10::InferenceMode inference;
        auto& context=*static_cast<Context*>(pointer);
        context.prepare(n,rows,columns,diagonal,off,fd,fl,ic0!=0,*stats);
        auto phase=Clock::now();
        auto b=device(rhs,static_cast<std::size_t>(n),at::kDouble);
        auto x=device(output,static_cast<std::size_t>(n),at::kDouble);
        (void)x[0].item<double>();stats->uploadSeconds+=seconds(phase);phase=Clock::now();
        auto r=b-product(context.matrix,x);
        auto z=context.precondition(r),p=z.clone();
        double rz=at::dot(r,z).item<double>();
        bool accepted=gpuNorm(r)<=stop;
        for(std::uint64_t iteration=0;!accepted && iteration<maximumIterations;++iteration) {
            const auto ap=product(context.matrix,p);
            const double denominator=at::dot(p,ap).item<double>();
            require(std::isfinite(denominator) && denominator>0 && std::isfinite(rz) && rz>0,
                    "PCG lost positive definiteness");
            const double alpha=rz/denominator;
            x.add_(p,alpha);r.add_(ap,-alpha);++stats->iterations;
            const double residual=gpuNorm(r);require(std::isfinite(residual),"PCG residual is nonfinite");
            if(residual<=stop) {
                r=b-product(context.matrix,x);
                accepted=gpuNorm(r)<=stop;
                if(accepted)break;
                z=context.precondition(r);p=z.clone();rz=at::dot(r,z).item<double>();
                continue;
            }
            z=context.precondition(r);
            const double next=at::dot(r,z).item<double>();
            require(std::isfinite(next) && next>0,"PCG preconditioned residual invalid");
            p.mul_(next/rz).add_(z);rz=next;
        }
        stats->iterationSeconds+=seconds(phase);
        require(accepted,"PCG iteration budget exhausted; final state not returned");
        phase=Clock::now();
        const auto host=x.to(at::kCPU).contiguous();
        std::memcpy(output,host.data_ptr<double>(),static_cast<std::size_t>(n)*sizeof(double));
        stats->downloadSeconds+=seconds(phase);return 0;
    }catch(const std::exception& e){
        if(pointer)static_cast<Context*>(pointer)->oldDiagonal.clear();
        errorText(error,capacity,e.what());return 1;
    }
}
