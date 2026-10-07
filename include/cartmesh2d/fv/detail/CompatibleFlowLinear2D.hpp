#pragma once
// Shared sparse saddle solve building blocks. K is never modified by its
// approximate inverse; no geometry, boundary or nonlinear acceptance here.
#include "cartmesh2d/fv/detail/NewtonKrylov2D.hpp"
#include <cstdint>
#include <map>
#include <numeric>
#include <tuple>
namespace cartmesh2d::fv::detail::compatible::linear {
using Vec=LinearVector2D;
struct Entry {std::int64_t row,column;double value;};
static_assert(sizeof(Entry)==24);
struct Matrix {
    std::size_t n;std::vector<std::size_t> rows,columns;Vec values;
    Matrix(std::size_t size,std::vector<Entry> entries):n(size),rows(n+1) {
        for(const auto& e:entries)linearEnsure(e.row>=0 && e.column>=0 && static_cast<std::size_t>(e.row)<n &&
            static_cast<std::size_t>(e.column)<n && std::isfinite(e.value),"invalid native entry");
        std::stable_sort(entries.begin(),entries.end(),[](const auto& a,const auto& b){return std::tie(a.row,a.column)<std::tie(b.row,b.column);});
        for(std::size_t k=0;k<entries.size();) {
            const auto i=entries[k].row,j=entries[k].column;double v=0;
            do {v=linearFinite(v+entries[k++].value);}while(k<entries.size() && entries[k].row==i && entries[k].column==j);
            columns.push_back(static_cast<std::size_t>(j));values.push_back(v);++rows[static_cast<std::size_t>(i)+1];
        }
        std::partial_sum(rows.begin(),rows.end(),rows.begin());
    }
    Vec apply(const Vec& x)const {
        linearEnsure(x.size()==n,"native matrix size mismatch");Vec y(n);
        for(std::size_t i=0;i<n;++i)for(auto k=rows[i];k<rows[i+1];++k)y[i]+=values[k]*x[columns[k]];
        for(auto v:y)linearFinite(v);
        return y;
    }
    // Recompute the true residual of the original assembled matrix without
    // first rounding A*x.  Near a strict target, separately rounded products
    // and sums can otherwise make iterative refinement stagnate above the
    // requested gate even though the retained double solution is accurate.
    // This changes neither K nor b nor the target; it mirrors the portable
    // FMA/TwoSum residual already used by the product sparse flow solver.
    Vec residual(const Vec& rhs,const Vec& x)const {
        linearEnsure(rhs.size()==n&&x.size()==n,"native residual size mismatch");Vec r(n);
        for(std::size_t i=0;i<n;++i) {
            LinearTwofold2D value{rhs[i],0};
            for(auto k=rows[i];k<rows[i+1];++k) {
                const double product=linearFinite(-values[k]*x[columns[k]]);
                value.add(product);
                value.add(std::fma(-values[k],x[columns[k]],-product));
            }
            r[i]=value.high;
        }
        return r;
    }
};
inline std::vector<std::pair<std::size_t,std::size_t>> connections(const Matrix& k,std::size_t nv) {
    std::vector<std::pair<std::size_t,std::size_t>> result;
    for(std::size_t i=0;i<nv;++i)for(auto p=k.rows[i];p<k.rows[i+1];++p)
        if(k.columns[p]<nv && i!=k.columns[p])result.emplace_back(std::min(i,k.columns[p]),std::max(i,k.columns[p]));
    // Condensation can round one side of a mathematically symmetric entry to
    // zero. Keep the union graph in the preconditioner; retain the original K.
    std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());
    return result;
}
struct PressureSchurDiagonal {
    struct Data {
        Vec diag;
        std::map<std::pair<std::size_t,std::size_t>,double> off;
        std::vector<std::pair<std::size_t,std::size_t>> connections;
        double symmetryError=0,coefficientMaximum=0;
    } data;
    SparsePattern2D pattern;
    SparseSystem2D system;
    LinearWorkspace2D workspace;
    std::unique_ptr<AggregationHierarchy2D> aggregation;

    static Data build(const Matrix& matrix,const Matrix& velocityMatrix,std::size_t nv) {
        linearEnsure(matrix.n==velocityMatrix.n&&nv<matrix.n,"Schur input dimensions invalid");
        const auto np=matrix.n-nv;Vec diagonalA(nv);
        std::vector<std::vector<std::pair<std::size_t,double>>> divergence(nv),gradient(nv);
        for(std::size_t i=0;i<velocityMatrix.n;++i)for(auto k=velocityMatrix.rows[i];k<velocityMatrix.rows[i+1];++k)
            if(i<nv&&velocityMatrix.columns[k]==i)diagonalA[i]=velocityMatrix.values[k];
        for(double value:diagonalA)linearEnsure(std::isfinite(value)&&value>0,"Schur velocity diagonal invalid");
        for(std::size_t i=0;i<matrix.n;++i)for(auto k=matrix.rows[i];k<matrix.rows[i+1];++k) {
            const auto j=matrix.columns[k];const auto value=matrix.values[k];
            if(i>=nv&&j<nv)divergence[j].emplace_back(i-nv,value);
            if(i<nv&&j>=nv)gradient[i].emplace_back(j-nv,value);
        }
        std::map<std::pair<std::size_t,std::size_t>,long double> coefficients;
        for(std::size_t k=0;k<nv;++k)for(const auto& [i,d]:divergence[k])for(const auto& [j,g]:gradient[k])
            coefficients[{i,j}]+=static_cast<long double>(d)*g/diagonalA[k];
        Data result;result.diag.assign(np,0);
        for(const auto& [key,sum]:coefficients) {
            const double value=linearFinite(static_cast<double>(sum));
            result.coefficientMaximum=std::max(result.coefficientMaximum,std::abs(value));
            if(key.first==key.second)result.diag[key.first]=value;
            else if(value!=0)result.off[key]=value;
        }
        for(std::size_t i=0;i<np;++i)linearEnsure(std::isfinite(result.diag[i])&&result.diag[i]>0,"Schur diagonal is not positive");
        for(const auto& [key,value]:result.off) {
            const auto reverse=result.off.find({key.second,key.first});
            const double other=reverse==result.off.end()?0:reverse->second;
            result.symmetryError=std::max(result.symmetryError,std::abs(value-other));
            if(key.first<key.second)result.connections.push_back(key);
        }
        linearEnsure(result.symmetryError<=1e-12*std::max(1.,result.coefficientMaximum),"Schur diagonal approximation is not symmetric");
        return result;
    }
    PressureSchurDiagonal(const Matrix& matrix,const Matrix& velocityMatrix,std::size_t nv,
        LinearPressureMethod2D inverse=LinearPressureMethod2D::IC0)
        :PressureSchurDiagonal(build(matrix,velocityMatrix,nv),inverse) {}
    explicit PressureSchurDiagonal(Data assembled,LinearPressureMethod2D inverse=LinearPressureMethod2D::IC0)
        :data(std::move(assembled)),pattern(data.diag.size(),data.connections),system(pattern),workspace(data.diag.size()) {
        linearEnsure(inverse==LinearPressureMethod2D::IC0||inverse==LinearPressureMethod2D::Aggregation,"Unsupported Schur inverse");
        system.diag=data.diag;
        for(const auto& [key,value]:data.off)system.off[pattern.slot(key.first,key.second)]=value;
        if(inverse==LinearPressureMethod2D::Aggregation) {
            // A fixed linear V-cycle is compatible with the ordinary outer
            // GMRES. Preserve every coefficient; the hierarchy rejects a
            // nonsymmetry or positive off diagonals rather than repairing them.
            aggregation=std::make_unique<AggregationHierarchy2D>(pattern.rows,pattern.columns,system.diag,system.off);
        } else system.factorIC0();
    }
    void apply(const Vec& input,std::size_t offset,Vec& output) {
        linearEnsure(input.size()==output.size()&&offset+data.diag.size()==input.size(),"Schur pressure vector size mismatch");
        for(std::size_t i=0;i<data.diag.size();++i)workspace.r[i]=input[offset+i];
        if(aggregation)aggregation->apply(workspace.r,workspace.z);
        else system.precondition(workspace,LinearPressureMethod2D::IC0);
        for(std::size_t i=0;i<data.diag.size();++i)output[offset+i]=-workspace.z[i];
    }
};
struct Block {
    const Matrix& matrix;const Vec& volume;std::size_t nv,np;bool corrected,ilu0,oseenScale,schurDiagonal;double viscosity,transportSpeed;
    SparsePattern2D pattern;SparseSystem2D velocity;LinearWorkspace2D workspace;
    std::unique_ptr<PressureSchurDiagonal> pressureSchur;
    LinearPressureMethod2D method;double asymmetry=0,matrixAsymmetry=0,symmetryCorrection=0;std::size_t applications=0;
    Block(const Matrix& k,const Vec& area,const std::string& mode,const std::string& pressure,double nu,double speed,const Matrix& pk):matrix(k),volume(area),
        nv(k.n-(area.size()-((pressure=="outlet"||pressure=="outlet-oseen"||pressure=="outlet-schur-diag"||pressure=="outlet-schur-aggregation")?0:1))),
        np(area.size()-((pressure=="outlet"||pressure=="outlet-oseen"||pressure=="outlet-schur-diag"||pressure=="outlet-schur-aggregation")?0:1)),
        corrected(pressure=="gauge"),ilu0(mode=="ilu0"),oseenScale(pressure=="outlet-oseen"),schurDiagonal(pressure=="outlet-schur-diag"||pressure=="outlet-schur-aggregation"),viscosity(nu),transportSpeed(speed),
        pattern(nv,connections(pk,nv)),velocity(pattern),workspace(nv),method(mode=="ic0"?LinearPressureMethod2D::IC0:LinearPressureMethod2D::Jacobi) {
        linearEnsure(nv>0&&pk.n==k.n,"expected matching velocity/pressure preconditioner dimensions");
        for(std::size_t i=0;i<k.n;++i)for(auto p=k.rows[i];p<k.rows[i+1];++p)if(i>=nv&&k.columns[p]>=nv)linearEnsure(k.values[p]==0,"only the zero retained-pressure block format is supported");
        for(std::size_t i=0;i<pk.n;++i)for(auto p=pk.rows[i];p<pk.rows[i+1];++p) {
            const auto j=pk.columns[p];const double v=pk.values[p];
            if(i<nv && j<nv){if(i==j)velocity.diag[i]=v;else velocity.off[pattern.slot(i,j)]=v;}
            if(i>=nv && j>=nv)linearEnsure(v==0,"only the zero retained-pressure block format is supported");
        }
        for(std::size_t i=0;i<nv;++i)for(auto p=k.rows[i];p<k.rows[i+1];++p){const auto j=k.columns[p];if(j>=nv||j==i)continue;
            const auto first=k.columns.begin()+static_cast<std::ptrdiff_t>(k.rows[j]),last=k.columns.begin()+static_cast<std::ptrdiff_t>(k.rows[j+1]),at=std::lower_bound(first,last,i);
            const double reverse=at==last||*at!=i?0:k.values[std::size_t(at-k.columns.begin())];matrixAsymmetry=std::max(matrixAsymmetry,std::abs(k.values[p]-reverse));}
        // IC0 uses the symmetric part only in its approximate inverse. ILU0
        // preserves all nonsymmetric velocity entries. Original K is unchanged.
        for(std::size_t i=0;i<nv;++i) {
            linearEnsure(std::isfinite(velocity.diag[i])&&velocity.diag[i]>0,"invalid velocity diagonal");
            for(auto p=pattern.rows[i];p<pattern.lowerEnd[i];++p) {
                const auto q=pattern.transpose[p];asymmetry=std::max(asymmetry,std::abs(velocity.off[p]-velocity.off[q]));
                if(mode=="ic0"){const double value=.5*velocity.off[p]+.5*velocity.off[q];
                    symmetryCorrection=std::max(symmetryCorrection,std::max(std::abs(value-velocity.off[p]),std::abs(value-velocity.off[q])));
                    velocity.off[p]=velocity.off[q]=value;}
            }
        }
        if(mode=="ic0")velocity.factorIC0();
        if(ilu0)velocity.factorILU0();
        if(schurDiagonal)pressureSchur=std::make_unique<PressureSchurDiagonal>(matrix,pk,nv,
            pressure=="outlet-schur-aggregation"?LinearPressureMethod2D::Aggregation:LinearPressureMethod2D::IC0);
    }
    Vec apply(const Vec& r) {
        ++applications;Vec z(r.size());long double sum=0;
        for(std::size_t i=nv;i<r.size();++i)sum+=r[i];
        // In pressure coordinates p_last=0, the mean-free mass is
        // M=diag(V_i)-V_i*V_j/Vtotal. Its inverse is exactly
        // diag(1/V_i)+1*1^T/V_last (Sherman-Morrison), not diag(1/V_i).
        // An outlet fixes the pressure level: all cell pressures are retained
        // and the mass inverse is diagonal (no mean removal). The nu scale is
        // a viscous Schur approximation, not an exact Oseen Schur inverse.
        // Explicit outlet-oseen uses nu+|U|*sqrt(V_i), with the same L^2/T
        // dimensions as viscosity in 2D. Only the approximate inverse changes.
        const double shift=corrected?static_cast<double>(sum/volume.back()):0;
        if(schurDiagonal) {
            pressureSchur->apply(r,nv,z);
        } else for(std::size_t i=nv;i<r.size();++i) {
            const double scale=viscosity+(oseenScale?transportSpeed*std::sqrt(volume[i-nv]):0);
            z[i]=scale*(-r[i]/volume[i-nv]-shift);
        }
        for(std::size_t i=0;i<nv;++i) {
            workspace.r[i]=r[i];
            for(auto p=matrix.rows[i];p<matrix.rows[i+1];++p)if(matrix.columns[p]>=nv)
                workspace.r[i]-=matrix.values[p]*z[matrix.columns[p]];
        }
        if(ilu0)velocity.preconditionILU0(workspace.r,workspace.z);
        else velocity.precondition(workspace,method);
        std::copy(workspace.z.begin(),workspace.z.end(),z.begin());
        return z;
    }
};

}
