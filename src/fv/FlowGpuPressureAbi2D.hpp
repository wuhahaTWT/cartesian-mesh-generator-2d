#pragma once
#include <cstdint>

// ABI v1 uses only POD buffers and explicit error text; no exception, STL value,
// or allocated memory crosses the DLL boundary. All inputs are host arrays.
struct CartmeshGpuPressureStatsV1 {
    std::uint64_t iterations=0, matrixUploads=0;
    double startupSeconds=0, uploadSeconds=0, iterationSeconds=0, downloadSeconds=0;
};
using CartmeshGpuPressureCreateV1 = int (*)(void**,char*,std::uint64_t);
using CartmeshGpuPressureDestroyV1 = void (*)(void*);
using CartmeshGpuPressureSolveV1 = int (*)(void*,std::uint64_t,const std::uint64_t*,
    const std::uint64_t*,const double*,const double*,const double*,double*,
    const double*,const double*,int,double,std::uint64_t,CartmeshGpuPressureStatsV1*,char*,std::uint64_t);
