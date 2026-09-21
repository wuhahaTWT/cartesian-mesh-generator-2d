#pragma once

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <vector>

namespace cartmesh2d::fv::detail {

// Type-II Anderson mixing of a scaled fixed-point map. This only proposes
// an iterate; the caller must check its actual equations before accepting it.
class Anderson2D {
    using Vector = std::vector<double>;
    Vector previousF_, previousG_;
    std::vector<Vector> differencesF_, differencesG_;
    static double product(const Vector& a,const Vector& b) {
        double value=0;
        for(std::size_t i=0;i<a.size();++i)value+=a[i]*b[i];
        return value;
    }
public:
    void clear() { previousF_.clear();previousG_.clear();differencesF_.clear();differencesG_.clear(); }
    std::optional<Vector> propose(const Vector& x,const Vector& g) {
        if(x.empty() || x.size()!=g.size() || (!previousF_.empty() && previousF_.size()!=g.size()))
            throw std::invalid_argument("Anderson state size mismatch");
        Vector f(g.size());
        for(std::size_t i=0;i<f.size();++i)f[i]=g[i]-x[i];
        if(!previousF_.empty()) {
            Vector df(f.size()),dg(f.size());
            for(std::size_t i=0;i<f.size();++i) {df[i]=f[i]-previousF_[i];dg[i]=g[i]-previousG_[i];}
            differencesF_.push_back(std::move(df));differencesG_.push_back(std::move(dg));
            if(differencesF_.size()>4) {differencesF_.erase(differencesF_.begin());differencesG_.erase(differencesG_.begin());}
        }
        previousF_=f;previousG_=g;
        if(differencesF_.empty())return {};
        // Twice-orthogonalized QR; reject rank-deficient windows instead of
        // dividing by a vanishing normal-equation pivot.
        const auto count=differencesF_.size();
        std::vector<Vector> q;std::vector<Vector> upper(count,Vector(count));
        Vector rhs(count),gamma(count);
        double largest=0;
        for(const auto& v:differencesF_)largest=std::max(largest,std::sqrt(product(v,v)));
        if(!(largest>0) || !std::isfinite(largest))return {};
        for(std::size_t j=0;j<count;++j) {
            auto v=differencesF_[j];
            for(int pass=0;pass<2;++pass)for(std::size_t i=0;i<j;++i) {
                const double projection=product(q[i],v);upper[i][j]+=projection;
                for(std::size_t k=0;k<v.size();++k)v[k]-=projection*q[i][k];
            }
            const double norm=std::sqrt(product(v,v));
            if(!std::isfinite(norm) || norm<1e-10*largest)return {};
            upper[j][j]=norm;for(auto& value:v)value/=norm;
            rhs[j]=product(v,f);q.push_back(std::move(v));
        }
        double coefficientNorm=0;
        for(std::size_t i=count;i-- >0;) {
            double value=rhs[i];
            for(std::size_t j=i+1;j<count;++j)value-=upper[i][j]*gamma[j];
            gamma[i]=value/upper[i][i];coefficientNorm+=std::abs(gamma[i]);
        }
        if(!std::isfinite(coefficientNorm) || coefficientNorm>1e4)return {};
        auto candidate=g;
        for(std::size_t j=0;j<count;++j)for(std::size_t i=0;i<g.size();++i)
            candidate[i]-=gamma[j]*differencesG_[j][i];
        for(double value:candidate)if(!std::isfinite(value))return {};
        return candidate;
    }
};
} // namespace cartmesh2d::fv::detail
