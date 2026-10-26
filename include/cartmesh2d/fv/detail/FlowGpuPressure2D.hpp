#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace cartmesh2d::fv::detail {

// The optional module has a C boundary: the standalone solver keeps its static
// CRT and never exchanges C++ objects/allocations with the external Torch CRT.
struct GpuPressureStats2D {
    std::size_t solves=0, iterations=0, matrixUploads=0;
    double startupSeconds=0, uploadSeconds=0, iterationSeconds=0, downloadSeconds=0;
    GpuPressureStats2D& operator+=(const GpuPressureStats2D& other) {
        solves+=other.solves;iterations+=other.iterations;matrixUploads+=other.matrixUploads;
        startupSeconds+=other.startupSeconds;uploadSeconds+=other.uploadSeconds;
        iterationSeconds+=other.iterationSeconds;downloadSeconds+=other.downloadSeconds;
        return *this;
    }
};
inline GpuPressureStats2D operator-(const GpuPressureStats2D& a,const GpuPressureStats2D& b) {
    return {a.solves-b.solves,a.iterations-b.iterations,a.matrixUploads-b.matrixUploads,
        a.startupSeconds-b.startupSeconds,a.uploadSeconds-b.uploadSeconds,
        a.iterationSeconds-b.iterationSeconds,a.downloadSeconds-b.downloadSeconds};
}

class GpuPressure2D {
    struct Impl;
    std::shared_ptr<Impl> impl_;
public:
    GpuPressure2D();
    std::size_t solve(const std::vector<std::size_t>& rows,const std::vector<std::size_t>& columns,
        const std::vector<double>& diagonal,const std::vector<double>& off,
        const std::vector<double>& rhs,std::vector<double>& x,
        const std::vector<double>& factorDiagonal,const std::vector<double>& factorLower,
        bool ic0,double absoluteStop,std::size_t maximumIterations,GpuPressureStats2D& stats);
};

// Explicit research capture, once per calling thread, before a nonzero pressure
// solve. Replay calls the production native/GPU linear solvers, not a second PDE.
void setPressureCapturePath2D(std::string path);
void capturePressureSystem2D(const std::vector<std::size_t>& rows,const std::vector<std::size_t>& columns,
    const std::vector<double>& diagonal,const std::vector<double>& off,
    const std::vector<double>& rhs,const std::vector<double>& x,double relativeTolerance);

} // namespace cartmesh2d::fv::detail
