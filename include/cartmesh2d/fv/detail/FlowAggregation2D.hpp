#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cartmesh2d::fv::detail {

// Experimental unsmoothed aggregation: piecewise-constant P, R=P^T,
// Galerkin coarse matrices, paired forward/backward Gauss-Seidel sweeps and
// an exact small coarse LDL^T solve. One fixed symmetric V-cycle per PCG step.
// Standard multigrid principles; independently implemented, no hypre code/library.
class AggregationHierarchy2D {
    struct Level {
        std::vector<std::size_t> rows, columns, aggregate;
        std::vector<double> diagonal, off;
        mutable std::vector<double> rhs, x;
    };
    std::vector<Level> levels_;
    std::vector<double> coarseLower_, coarseDiagonal_;
    bool diagonalCoarse_ = false;

    static void require(bool condition, const char* message) {
        if (!condition) throw std::runtime_error(message);
    }
    static double finite(double value) {
        require(std::isfinite(value), "Flow aggregation numerical range exceeded");
        return value;
    }
    static void scratch(Level& level) {
        level.rhs.resize(level.diagonal.size());
        level.x.resize(level.diagonal.size());
    }
    void factorCoarse() {
        const auto& level=levels_.back();
        const auto n=level.diagonal.size();
        if (diagonalCoarse_) return;
        require(n<=32,"Flow aggregation coarse system exceeds dense budget");
        coarseLower_.assign(n*n,0.); coarseDiagonal_.resize(n);
        for (std::size_t i=0;i<n;++i)
            for (std::size_t k=level.rows[i];k<level.rows[i+1];++k)
                coarseLower_[i*n+level.columns[k]]=level.off[k];
        for (std::size_t i=0;i<n;++i) {
            double pivot=level.diagonal[i];
            for (std::size_t j=0;j<i;++j) {
                double value=coarseLower_[i*n+j];
                for (std::size_t k=0;k<j;++k)
                    value-=coarseLower_[i*n+k]*coarseDiagonal_[k]*coarseLower_[j*n+k];
                const double lij=finite(value/coarseDiagonal_[j]);
                coarseLower_[i*n+j]=lij; pivot-=lij*lij*coarseDiagonal_[j];
            }
            require(std::isfinite(pivot)&&pivot>0,"Flow aggregation nonpositive coarse pivot");
            coarseDiagonal_[i]=pivot;
        }
    }
    void cycle(std::size_t index) const {
        const auto& level=levels_[index]; const auto n=level.diagonal.size();
        auto& x=level.x;
        if (index+1==levels_.size()) {
            if (diagonalCoarse_) {
                for (std::size_t i=0;i<n;++i) x[i]=finite(level.rhs[i]/level.diagonal[i]);
                return;
            }
            for (std::size_t i=0;i<n;++i) {
                double value=level.rhs[i];
                for (std::size_t j=0;j<i;++j) value-=coarseLower_[i*n+j]*x[j];
                x[i]=finite(value);
            }
            for (std::size_t i=0;i<n;++i) x[i]=finite(x[i]/coarseDiagonal_[i]);
            for (std::size_t i=n;i-- >0;)
                for (std::size_t j=0;j<i;++j) x[j]-=coarseLower_[i*n+j]*x[i];
            return;
        }
        // Forward smoothing from zero: upper unknowns are mathematically zero,
        // so only the already overwritten lower entries are needed. This also
        // avoids clearing/re-reading the full scratch vector each V-cycle.
        for (std::size_t i=0;i<n;++i) {
            double value=level.rhs[i];
            for (std::size_t k=level.rows[i];k<level.rows[i+1] && level.columns[k]<i;++k)
                value-=level.off[k]*x[level.columns[k]];
            x[i]=finite(value/level.diagonal[i]);
        }
        const auto& coarse=levels_[index+1];
        std::fill(coarse.rhs.begin(),coarse.rhs.end(),0.);
        for (std::size_t i=0;i<n;++i) {
            double residual=level.rhs[i]-level.diagonal[i]*x[i];
            for (std::size_t k=level.rows[i];k<level.rows[i+1];++k)
                residual-=level.off[k]*x[level.columns[k]];
            coarse.rhs[level.aggregate[i]]+=finite(residual);
        }
        cycle(index+1);
        for (std::size_t i=0;i<n;++i) x[i]+=coarse.x[level.aggregate[i]];
        // Reverse sweep is the transpose partner of the forward sweep.
        for (std::size_t i=n;i-- >0;) {
            double value=level.rhs[i];
            for (std::size_t k=level.rows[i];k<level.rows[i+1];++k)
                value-=level.off[k]*x[level.columns[k]];
            x[i]=finite(value/level.diagonal[i]);
        }
    }
public:
    AggregationHierarchy2D(const std::vector<std::size_t>& rows,
        const std::vector<std::size_t>& columns,
        const std::vector<double>& diagonal,const std::vector<double>& off) {
        require(!diagonal.empty() && rows.size()==diagonal.size()+1 &&
            rows.front()==0 && rows.back()==off.size() && columns.size()==off.size(),
            "Flow aggregation matrix dimensions invalid");
        Level fine; fine.rows=rows; fine.columns=columns; fine.diagonal=diagonal; fine.off=off;
        scratch(fine); levels_.push_back(std::move(fine));
        for (;;) {
            auto& level=levels_.back(); const auto n=level.diagonal.size();
            for (std::size_t i=0;i<n;++i)
                require(level.rows[i]<=level.rows[i+1] && level.rows[i+1]<=level.off.size(),
                    "Flow aggregation rows invalid");
            for (std::size_t i=0;i<n;++i) {
                require(std::isfinite(level.diagonal[i])&&level.diagonal[i]>0,
                    "Flow aggregation diagonal invalid");
                std::size_t previous=n;
                for (std::size_t k=level.rows[i];k<level.rows[i+1];++k) {
                    const auto j=level.columns[k];
                    require(j<n && j!=i && (previous==n || j>previous),
                        "Flow aggregation column ordering invalid"); previous=j;
                    require(std::isfinite(level.off[k])&&level.off[k]<=0,
                        "Flow aggregation requires nonpositive off-diagonal coefficients");
                }
            }
            // Validate the entire CSR layout before following transpose rows.
            for (std::size_t i=0;i<n;++i) {
                for (std::size_t k=level.rows[i];k<level.rows[i+1];++k) {
                    const auto j=level.columns[k];
                    const auto begin=level.columns.begin()+level.rows[j];
                    const auto end=level.columns.begin()+level.rows[j+1];
                    const auto transpose=std::lower_bound(begin,end,i);
                    require(transpose!=end && *transpose==i &&
                        level.off[static_cast<std::size_t>(transpose-level.columns.begin())]==level.off[k],
                        "Flow aggregation requires an exactly symmetric matrix");
                }
            }
            if (n<=32) break;
            level.aggregate.assign(n,n); std::size_t count=0;
            for (std::size_t i=0;i<n;++i) {
                if (level.aggregate[i]!=n) continue;
                std::size_t partner=n; double strongest=0.;
                for (std::size_t k=level.rows[i];k<level.rows[i+1];++k) {
                    const auto j=level.columns[k];
                    if (level.aggregate[j]==n && -level.off[k]>strongest) {
                        strongest=-level.off[k]; partner=j;
                    }
                }
                level.aggregate[i]=count;
                if (partner!=n) level.aggregate[partner]=count;
                ++count;
            }
            // Maximal matching can leave many unpaired vertices around a hub.
            // Attach each connected singleton to a neighbouring matched group
            // instead of carrying all those vertices to the next level. These
            // vertices remain in P and in every Galerkin contribution; no row
            // or coupling is dropped. Isolated diagonal unknowns stay separate.
            std::vector<std::size_t> groupSize(count,0);
            for (auto group:level.aggregate) ++groupSize[group];
            const auto matchedSize=groupSize;
            for (std::size_t i=0;i<n;++i) {
                const auto original=level.aggregate[i];
                if (matchedSize[original]!=1) continue;
                std::size_t target=count; double strongest=0.;
                for (std::size_t k=level.rows[i];k<level.rows[i+1];++k) {
                    const auto group=level.aggregate[level.columns[k]];
                    if (matchedSize[group]>=2 && -level.off[k]>strongest) {
                        strongest=-level.off[k]; target=group;
                    }
                }
                if (target!=count) {
                    level.aggregate[i]=target; --groupSize[original]; ++groupSize[target];
                }
            }
            std::vector<std::size_t> compact(count,0); count=0;
            for (std::size_t group=0;group<groupSize.size();++group)
                if (groupSize[group]) compact[group]=count++;
            for (auto& group:level.aggregate) group=compact[group];
            if (count==n) { diagonalCoarse_=true; break; } // All remaining off-diagonals are zero.
            require(count<=n-n/5,"Flow aggregation coarsening too weak for memory budget; select IC0");
            require(levels_.size()<32,"Flow aggregation depth budget exceeded");
            Level coarse;
            std::vector<long double> diagonalSum(count,0.);
            // Bucket upper-triangle fine contributions by their coarse row.
            // Unlike a global tuple sort, bounded-degree rows only sort their
            // local columns. Stable order within each coarse pair preserves the
            // original fine-face summation order and long-double reduction.
            std::vector<std::size_t> upperRows(count+1,0);
            for (std::size_t i=0;i<n;++i) {
                const auto ci=level.aggregate[i]; diagonalSum[ci]+=level.diagonal[i];
                for (std::size_t k=level.rows[i];k<level.rows[i+1];++k) {
                    const auto j=level.columns[k]; if (j<=i || level.off[k]==0) continue;
                    const auto cj=level.aggregate[j];
                    if (ci==cj) diagonalSum[ci]+=2*static_cast<long double>(level.off[k]);
                    else ++upperRows[std::min(ci,cj)+1];
                }
            }
            for (std::size_t i=0;i<count;++i) upperRows[i+1]+=upperRows[i];
            using Entry=std::pair<std::size_t,double>;
            std::vector<Entry> upper(upperRows.back());
            auto cursor=upperRows;
            for (std::size_t i=0;i<n;++i) {
                const auto ci=level.aggregate[i];
                for (std::size_t k=level.rows[i];k<level.rows[i+1];++k) {
                    const auto j=level.columns[k]; if (j<=i || level.off[k]==0) continue;
                    const auto cj=level.aggregate[j];
                    if (ci!=cj) upper[cursor[std::min(ci,cj)]++]={std::max(ci,cj),level.off[k]};
                }
            }
            std::vector<std::size_t> uniqueRows(count+1,0);
            std::size_t used=0;
            coarse.rows.assign(count+1,0);
            for (std::size_t i=0;i<count;++i) {
                const auto end=upperRows[i+1];
                std::stable_sort(upper.begin()+static_cast<std::ptrdiff_t>(upperRows[i]),
                    upper.begin()+static_cast<std::ptrdiff_t>(end),
                    [](const Entry& a,const Entry& b) { return a.first<b.first; });
                for (auto k=upperRows[i];k<end;) {
                    const auto j=upper[k].first; long double value=0;
                    do { value+=upper[k++].second; } while (k<end && upper[k].first==j);
                    upper[used++]={j,finite(static_cast<double>(value))};
                    ++coarse.rows[i+1]; ++coarse.rows[j+1];
                }
                uniqueRows[i+1]=used;
            }
            for (std::size_t i=0;i<count;++i) coarse.rows[i+1]+=coarse.rows[i];
            coarse.columns.resize(coarse.rows.back()); coarse.off.resize(coarse.rows.back());
            cursor=coarse.rows;
            for (std::size_t i=0;i<count;++i) {
                for (auto k=uniqueRows[i];k<uniqueRows[i+1];++k) {
                    const auto [j,value]=upper[k];
                    const auto ij=cursor[i]++, ji=cursor[j]++;
                    coarse.columns[ij]=j; coarse.off[ij]=value;
                    coarse.columns[ji]=i; coarse.off[ji]=value;
                }
            }
            // Transposed lower entries arrive in increasing row order, before
            // each row's sorted upper entries: CSR is sorted and exactly
            // symmetric without a second global sort or a doubled tuple list.
            coarse.diagonal.resize(count);
            for (std::size_t i=0;i<count;++i) coarse.diagonal[i]=finite(static_cast<double>(diagonalSum[i]));
            scratch(coarse);levels_.push_back(std::move(coarse));
        }
        factorCoarse();
    }
    bool matches(const std::vector<double>& diagonal,const std::vector<double>& off) const {
        return levels_.front().diagonal==diagonal && levels_.front().off==off;
    }
    void apply(const std::vector<double>& rhs,std::vector<double>& result) const {
        require(rhs.size()==levels_.front().diagonal.size(),"Flow aggregation RHS size invalid");
        levels_.front().rhs=rhs; cycle(0); result=levels_.front().x;
        for (double value:result) finite(value);
    }
    // Read-only matrix diagnostics for checking the Galerkin hierarchy.
    const auto& level(std::size_t index) const { return levels_.at(index); }
    std::size_t levels() const { return levels_.size(); }
    std::size_t coarseCells() const { return levels_.back().diagonal.size(); }
};
}
