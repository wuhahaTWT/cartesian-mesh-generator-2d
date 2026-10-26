#include "cartmesh2d/fv/detail/FlowGpuPressure2D.hpp"
#include "FlowGpuPressureAbi2D.hpp"
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#ifdef _WIN64
#include <windows.h>
#endif

namespace cartmesh2d::fv::detail {
namespace {
using Clock=std::chrono::steady_clock;
thread_local std::string capturePath;

template<class T> void write(std::ofstream& out,const std::vector<T>& data) {
    out.write(reinterpret_cast<const char*>(data.data()),static_cast<std::streamsize>(data.size()*sizeof(T)));
}
}

void setPressureCapturePath2D(std::string path) {capturePath=std::move(path);}
void capturePressureSystem2D(const std::vector<std::size_t>& rows,const std::vector<std::size_t>& columns,
    const std::vector<double>& diagonal,const std::vector<double>& off,
    const std::vector<double>& rhs,const std::vector<double>& x,double tolerance) {
    if(capturePath.empty())return;
    static_assert(sizeof(double)==8);
    const std::filesystem::path path=std::filesystem::u8path(capturePath);
    if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path,std::ios::binary);
    if(!out)throw std::runtime_error("Cannot create explicit pressure snapshot");
    const std::uint64_t n=diagonal.size(),nnz=off.size();
    out.write("CM2DPRS1",8);
    out.write(reinterpret_cast<const char*>(&n),8);out.write(reinterpret_cast<const char*>(&nnz),8);
    out.write(reinterpret_cast<const char*>(&tolerance),8);
    const std::vector<std::uint64_t> r(rows.begin(),rows.end()),c(columns.begin(),columns.end());
    write(out,r);write(out,c);write(out,diagonal);write(out,off);write(out,rhs);write(out,x);
    if(!out)throw std::runtime_error("Cannot write explicit pressure snapshot");
    capturePath.clear();
}

struct GpuPressure2D::Impl {
#ifdef _WIN64
    HMODULE module=nullptr;
    DLL_DIRECTORY_COOKIE directory=nullptr;
    void* context=nullptr;
    CartmeshGpuPressureDestroyV1 destroy=nullptr;
    CartmeshGpuPressureSolveV1 solve=nullptr;
    double startupSeconds=0;
    Impl() {
        const auto start=Clock::now();
        try {
            const wchar_t* root=_wgetenv(L"CARTMESH2D_TORCH_ROOT");
            if(!root || !*root)throw std::runtime_error("torch-cuda requires CARTMESH2D_TORCH_ROOT pointing to the existing torch package");
            const auto library=std::filesystem::absolute(std::filesystem::path(root)/L"lib");
            directory=AddDllDirectory(library.c_str());
            if(!directory)throw std::runtime_error("Cannot add existing Torch DLL directory");
            std::filesystem::path plugin;
            if(const auto* specified=_wgetenv(L"CARTMESH2D_PRESSURE_PLUGIN");specified && *specified)plugin=specified;
            else {
                std::wstring executable(32768,L'\0');
                const auto count=GetModuleFileNameW(nullptr,executable.data(),static_cast<DWORD>(executable.size()));
                if(!count || count==executable.size())throw std::runtime_error("Cannot locate native executable for GPU module");
                executable.resize(count);plugin=std::filesystem::path(executable).parent_path()/L"cartmesh2d_torch_pressure.dll";
            }
            plugin=std::filesystem::absolute(plugin);
            module=LoadLibraryExW(plugin.c_str(),nullptr,LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR|LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
            if(!module)throw std::runtime_error("Cannot load optional Torch pressure module; Windows error "+std::to_string(GetLastError()));
            const auto create=reinterpret_cast<CartmeshGpuPressureCreateV1>(GetProcAddress(module,"cartmesh2d_pressure_create_v1"));
            destroy=reinterpret_cast<CartmeshGpuPressureDestroyV1>(GetProcAddress(module,"cartmesh2d_pressure_destroy_v1"));
            solve=reinterpret_cast<CartmeshGpuPressureSolveV1>(GetProcAddress(module,"cartmesh2d_pressure_solve_v1"));
            if(!create || !destroy || !solve)throw std::runtime_error("GPU pressure module ABI v1 is incomplete");
            char error[2048]{};
            if(create(&context,error,sizeof(error)))throw std::runtime_error(error);
            startupSeconds=std::chrono::duration<double>(Clock::now()-start).count();
        } catch(...) {release();throw;}
    }
    void release() {
        if(context && destroy)destroy(context);
        context=nullptr;
        if(module)FreeLibrary(module);
        module=nullptr;
        if(directory)RemoveDllDirectory(directory);
        directory=nullptr;
    }
    ~Impl(){release();}
#else
    Impl(){throw std::runtime_error("Experimental torch-cuda pressure module currently supports Windows x64 only");}
#endif
};

GpuPressure2D::GpuPressure2D():impl_(std::make_shared<Impl>()){}
std::size_t GpuPressure2D::solve(const std::vector<std::size_t>& rows,const std::vector<std::size_t>& columns,
    const std::vector<double>& diagonal,const std::vector<double>& off,
    const std::vector<double>& rhs,std::vector<double>& x,
    const std::vector<double>& factorDiagonal,const std::vector<double>& factorLower,
    bool ic0,double stop,std::size_t maximumIterations,GpuPressureStats2D& stats) {
#ifdef _WIN64
    static_assert(sizeof(std::size_t)==sizeof(std::uint64_t));
    CartmeshGpuPressureStatsV1 result;char error[2048]{};
    stats.startupSeconds+=impl_->startupSeconds;impl_->startupSeconds=0;
    const auto status=impl_->solve(impl_->context,diagonal.size(),
        reinterpret_cast<const std::uint64_t*>(rows.data()),reinterpret_cast<const std::uint64_t*>(columns.data()),
        diagonal.data(),off.data(),rhs.data(),x.data(),factorDiagonal.data(),factorLower.data(),ic0?1:0,
        stop,maximumIterations,&result,error,sizeof(error));
    ++stats.solves;stats.iterations+=result.iterations;stats.matrixUploads+=result.matrixUploads;
    stats.startupSeconds+=result.startupSeconds;stats.uploadSeconds+=result.uploadSeconds;
    stats.iterationSeconds+=result.iterationSeconds;stats.downloadSeconds+=result.downloadSeconds;
    if(status)throw std::runtime_error(std::string("torch-cuda pressure: ")+error);
    return static_cast<std::size_t>(result.iterations);
#else
    (void)rows;(void)columns;(void)diagonal;(void)off;(void)rhs;(void)x;(void)factorDiagonal;(void)factorLower;
    (void)ic0;(void)stop;(void)maximumIterations;(void)stats;
    throw std::runtime_error("Experimental torch-cuda pressure module currently supports Windows x64 only");
#endif
}
} // namespace cartmesh2d::fv::detail
