#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string_view>
#include <limits>
#include <stdexcept>
#include <vector>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif

namespace cartmesh2d::fv::detail {

inline constexpr bool systemCholeskyAvailable2D() {
#ifdef __APPLE__
    return true;
#else
    return false;
#endif
}

// Optional macOS system sparse Cholesky factor, used only as a PCG
// preconditioner. The original CSR matrix and true-residual acceptance remain
// in SparseSystem2D. No coefficient shift, dropped coupling or solver fallback.
// API contract: Apple Accelerate Sparse/Solve.h (CSC symmetric lower triangle,
// SparseFactor, SparseRefactor, SparseSolve, SparseCleanup).
class SystemCholesky2D {
    using Vector=std::vector<double>;
    Vector diagonal_, off_;
    std::vector<std::size_t> mirror_;
#ifdef __APPLE__
    std::vector<std::size_t> slots_;
    std::vector<long> starts_;
    std::vector<int> indices_;
    Vector values_;
    SparseOpaqueFactorization_Double factor_{};
    bool ownsFactor_=false;
    SparseMatrix_Double matrix() {
        SparseAttributes_t attributes{};
        attributes.kind=SparseSymmetric;attributes.triangle=SparseLowerTriangle;
        const auto n=static_cast<int>(diagonal_.size());
        return {{n,n,starts_.data(),indices_.data(),attributes,1},values_.data()};
    }
    void fill(const Vector& diagonal,const Vector& off) {
        for(std::size_t i=0;i<diagonal.size();++i) {
            values_[static_cast<std::size_t>(starts_[i])]=diagonal[i];
            for(long k=starts_[i]+1;k<starts_[i+1];++k)
                values_[static_cast<std::size_t>(k)]=off[slots_[static_cast<std::size_t>(k)]];
        }
    }
#endif
    static void require(bool valid,const char* message) {
        if(!valid)throw std::runtime_error(message);
    }
    void validateValues(const Vector& diagonal,const Vector& off) const {
        require(diagonal.size()==diagonal_.size() && off.size()==off_.size(),"System Cholesky coefficient dimensions differ");
        for(double value:diagonal)require(std::isfinite(value)&&value>0,"System Cholesky diagonal invalid");
        for(std::size_t k=0;k<off.size();++k)
            require(std::isfinite(off[k]) && off[k]==off[mirror_[k]],"System Cholesky requires an exactly symmetric finite matrix");
    }
public:
    SystemCholesky2D(const std::vector<std::size_t>& rows,const std::vector<std::size_t>& columns,
                    const Vector& diagonal,const Vector& off) :diagonal_(diagonal),off_(off) {
        require(systemCholeskyAvailable2D(),"System sparse Cholesky is available only on macOS");
        const char* threads=std::getenv("VECLIB_MAXIMUM_THREADS");
        require(threads && std::string_view(threads)=="1",
            "System Cholesky requires VECLIB_MAXIMUM_THREADS=1 before any Accelerate call for deterministic results");
        const auto n=diagonal.size();
        require(n>0 && n<=static_cast<std::size_t>(std::numeric_limits<int>::max()) &&
            rows.size()==n+1 && rows.front()==0 && rows.back()==off.size() && columns.size()==off.size(),
            "System Cholesky matrix dimensions invalid");
        for(std::size_t i=0;i<n;++i) {
            require(rows[i]<=rows[i+1] && rows[i+1]<=off.size(),"System Cholesky rows invalid");
            for(std::size_t k=rows[i];k<rows[i+1];++k)
                require(columns[k]<n && columns[k]!=i && (k==rows[i] || columns[k-1]<columns[k]),
                    "System Cholesky columns invalid");
        }
        mirror_.resize(off.size());
        for(std::size_t i=0;i<n;++i)for(std::size_t k=rows[i];k<rows[i+1];++k) {
            const auto j=columns[k];
            const auto first=columns.begin()+static_cast<std::ptrdiff_t>(rows[j]);
            const auto last=columns.begin()+static_cast<std::ptrdiff_t>(rows[j+1]);
            const auto found=std::lower_bound(first,last,i);
            require(found!=last && *found==i,"System Cholesky transpose entry missing");
            mirror_[k]=static_cast<std::size_t>(found-columns.begin());
        }
        validateValues(diagonal,off);
#ifdef __APPLE__
        // Exact symmetry lets CSR upper entries represent CSC lower entries.
        // Include the diagonal explicitly; keep the same zero pattern on reuse.
        for(std::size_t i=0;i<n;++i) {
            starts_.push_back(static_cast<long>(values_.size()));
            indices_.push_back(static_cast<int>(i));values_.push_back(diagonal[i]);slots_.push_back(off.size());
            for(std::size_t k=rows[i];k<rows[i+1];++k)if(columns[k]>i) {
                indices_.push_back(static_cast<int>(columns[k]));values_.push_back(off[k]);slots_.push_back(k);
            }
        }
        starts_.push_back(static_cast<long>(values_.size()));
        factor_=SparseFactor(SparseFactorizationCholesky,matrix());ownsFactor_=true;
        if(factor_.status!=SparseStatusOK) {
            SparseCleanup(factor_);ownsFactor_=false;
            throw std::runtime_error("System Cholesky factorization failed; matrix must be positive definite");
        }
#endif
    }
    SystemCholesky2D(const SystemCholesky2D&)=delete;
    SystemCholesky2D& operator=(const SystemCholesky2D&)=delete;
    ~SystemCholesky2D() {
#ifdef __APPLE__
        if(ownsFactor_)SparseCleanup(factor_);
#endif
    }
    bool matches(const Vector& diagonal,const Vector& off) const {
#ifdef __APPLE__
        return ownsFactor_ && factor_.status==SparseStatusOK && diagonal_==diagonal && off_==off;
#else
        (void)diagonal;(void)off;return false;
#endif
    }
    void refresh(const Vector& diagonal,const Vector& off) {
        validateValues(diagonal,off);
#ifdef __APPLE__
        fill(diagonal,off);SparseRefactor(matrix(),&factor_);
        require(factor_.status==SparseStatusOK,"System Cholesky refactorization failed; matrix must be positive definite");
        diagonal_=diagonal;off_=off;
#else
        throw std::runtime_error("System sparse Cholesky is available only on macOS");
#endif
    }
    void apply(const Vector& input,Vector& output) const {
        require(input.size()==diagonal_.size() && output.size()==input.size(),"System Cholesky vector dimensions invalid");
#ifdef __APPLE__
        require(ownsFactor_ && factor_.status==SparseStatusOK,"System Cholesky factor unavailable");
        for(double value:input)require(std::isfinite(value),"System Cholesky right-hand side nonfinite");
        output=input;SparseSolve(factor_,DenseVector_Double{static_cast<int>(output.size()),output.data()});
        for(double value:output)require(std::isfinite(value),"System Cholesky solution nonfinite");
#else
        throw std::runtime_error("System sparse Cholesky is available only on macOS");
#endif
    }
};
} // namespace cartmesh2d::fv::detail
