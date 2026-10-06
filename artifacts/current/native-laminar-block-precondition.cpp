// Research linear solver for native P1 Stokes/Oseen condensed matrices.
// The matrix and RHS are kept unchanged. No geometry, PDE or acceptance gate
// is implemented here. Reuse product IC0/ILU0 and bounded GMRES building blocks.
#include "cartmesh2d/fv/detail/NewtonKrylov2D.hpp"
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>

using namespace cartmesh2d::fv::detail;
using Vec=LinearVector2D;
struct Entry {std::int64_t row,column;double value;};
static_assert(sizeof(Entry)==24);

template<class T> std::vector<T> binary(const std::string& name) {
    std::ifstream in(name,std::ios::binary|std::ios::ate);
    if(!in)throw std::runtime_error("cannot open native matrix input");
    const auto bytes=in.tellg();
    if(bytes<0 || static_cast<std::uint64_t>(bytes)%sizeof(T))throw std::runtime_error("truncated native matrix input");
    std::vector<T> values(static_cast<std::size_t>(bytes)/sizeof(T));in.seekg(0);
    in.read(reinterpret_cast<char*>(values.data()),bytes);
    if(!in)throw std::runtime_error("native matrix read failed");
    return values;
}
std::vector<std::string> split(const std::string& line) {
    std::istringstream in(line);std::vector<std::string> values;std::string value;
    while(std::getline(in,value,','))values.push_back(value);
    return values;
}
Vec areas(const std::string& name) {
    std::ifstream in(name);std::string line;
    if(!std::getline(in,line))throw std::runtime_error("missing native cell fields");
    const auto header=split(line);
    const auto ci=std::find(header.begin(),header.end(),"cell"),ai=std::find(header.begin(),header.end(),"area");
    if(ci==header.end() || ai==header.end())throw std::runtime_error("cell/area columns missing");
    const auto c=static_cast<std::size_t>(ci-header.begin()),a=static_cast<std::size_t>(ai-header.begin());Vec result;
    while(std::getline(in,line)) {
        const auto values=split(line);if(values.size()!=header.size())throw std::runtime_error("invalid native cell row");
        if(std::stoull(values[c])!=result.size())throw std::runtime_error("native cell order mismatch");
        const double area=std::stod(values[a]);linearEnsure(std::isfinite(area)&&area>0,"invalid native area");result.push_back(area);
    }
    linearEnsure(!result.empty(),"empty cell field");return result;
}
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
};
std::vector<std::pair<std::size_t,std::size_t>> connections(const Matrix& k,std::size_t nv) {
    std::vector<std::pair<std::size_t,std::size_t>> result;
    for(std::size_t i=0;i<nv;++i)for(auto p=k.rows[i];p<k.rows[i+1];++p)
        if(k.columns[p]<nv && i!=k.columns[p])result.emplace_back(std::min(i,k.columns[p]),std::max(i,k.columns[p]));
    // Condensation can round one side of a mathematically symmetric entry to
    // zero. Keep the union graph in the preconditioner; retain the original K.
    std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());
    return result;
}
struct Block {
    const Matrix& matrix;const Vec& volume;std::size_t nv,np;bool corrected,ilu0;double viscosity;
    SparsePattern2D pattern;SparseSystem2D velocity;LinearWorkspace2D workspace;
    LinearPressureMethod2D method;double asymmetry=0,symmetryCorrection=0;std::size_t applications=0;
    Block(const Matrix& k,const Vec& area,const std::string& mode,const std::string& pressure,double nu):matrix(k),volume(area),
        nv(k.n-(area.size()-(pressure=="outlet"?0:1))),np(area.size()-(pressure=="outlet"?0:1)),
        corrected(pressure=="gauge"),ilu0(mode=="ilu0"),viscosity(nu),
        pattern(nv,connections(k,nv)),velocity(pattern),workspace(nv),method(mode=="ic0"?LinearPressureMethod2D::IC0:LinearPressureMethod2D::Jacobi) {
        linearEnsure(nv>0,"expected retained velocity block followed by cell pressures");
        for(std::size_t i=0;i<k.n;++i)for(auto p=k.rows[i];p<k.rows[i+1];++p) {
            const auto j=k.columns[p];const double v=k.values[p];
            if(i<nv && j<nv){if(i==j)velocity.diag[i]=v;else velocity.off[pattern.slot(i,j)]=v;}
            if(i>=nv && j>=nv)linearEnsure(v==0,"only the zero retained-pressure block format is supported");
        }
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
        const double shift=corrected?static_cast<double>(sum/volume.back()):0;
        for(std::size_t i=nv;i<r.size();++i)z[i]=viscosity*(-r[i]/volume[i-nv]-shift);
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

int main(int argc,char** argv)try {
    if(argc<6 || argc>9)throw std::runtime_error("usage: probe input_prefix cell_csv ic0|ilu0|jacobi gauge|plain|outlet output_prefix [relative_tolerance=1e-11] [restarts=50] [viscosity=1]");
    const std::string prefix=argv[1],mode=argv[3],gauge=argv[4],output=argv[5];
    linearEnsure((mode=="ic0"||mode=="ilu0"||mode=="jacobi")&&(gauge=="gauge"||gauge=="plain"||gauge=="outlet"),"invalid preconditioner choice");
    const double tolerance=argc>=7?std::stod(argv[6]):1e-11;const auto maximum=argc>=8?std::stoull(argv[7]):50;
    const double nu=argc>=9?std::stod(argv[8]):1.;
    linearEnsure(std::isfinite(nu)&&nu>0,"invalid viscous Schur scale");
    linearEnsure(std::isfinite(tolerance)&&tolerance>0&&tolerance<=.01&&maximum>0,"invalid research solve control");
    linearEnsure(!std::filesystem::exists(output+".json") && !std::filesystem::exists(output+".solution") &&
        !std::filesystem::exists(output+".candidate"),"output exists; retain prior evidence and choose another prefix");
    const auto start=std::chrono::steady_clock::now();const auto rhs=binary<double>(prefix+".rhs"),area=areas(argv[2]);
    linearEnsure(rhs.size()>area.size()-(gauge=="outlet"?0:1),"invalid retained pressure dimension");
    for(double v:rhs)linearFinite(v);
    const Matrix matrix(rhs.size(),binary<Entry>(prefix+".entries"));Block block(matrix,area,mode,gauge,nu);
    const double setup=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    const double initial=linearNorm(rhs);Vec x(rhs.size()),r=rhs;std::size_t restarts=0,products=0;std::vector<double> history{initial};
    // Use the existing relative linear tolerance; do not add the product's
    // absolute residual floor to mixed physical units. Recheck the true
    // unpreconditioned native residual after every bounded GMRES restart.
    bool converged=initial==0;
    while(!converged && restarts<maximum) {
        const auto apply=[&](const Vec& v){++products;return matrix.apply(block.apply(v));};
        const auto direction=newtonKrylovDirection2D(apply,r,60);
        if(!direction)break;
        const auto step=block.apply(*direction);for(std::size_t i=0;i<x.size();++i)x[i]+=step[i];
        const auto ax=matrix.apply(x);for(std::size_t i=0;i<x.size();++i)r[i]=rhs[i]-ax[i];
        ++restarts;history.push_back(linearNorm(r));converged=history.back()<=tolerance*initial;
    }
    double maximumResidual=0;for(double v:r)maximumResidual=std::max(maximumResidual,std::abs(v));
    const std::string field=output+(converged?".solution":".candidate");std::ofstream result(field,std::ios::binary);
    result.write(reinterpret_cast<const char*>(x.data()),static_cast<std::streamsize>(x.size()*sizeof(double)));result.close();
    linearEnsure(bool(result),"candidate write failed");
    std::ostringstream record;record<<std::setprecision(17)<<"{\"converged\":"<<(converged?"true":"false")<<",\"unknowns\":"<<rhs.size()
        <<",\"velocity_unknowns\":"<<block.nv<<",\"pressure_unknowns\":"<<block.np<<",\"matrix_nonzeros\":"<<matrix.values.size()
        <<",\"velocity_off_nonzeros\":"<<block.velocity.off.size()<<",\"method\":\""<<mode<<"\",\"mass\":\""<<gauge
        <<"\",\"relative_tolerance\":"<<tolerance<<",\"rhs_norm\":"<<initial<<",\"residual_norm\":"<<history.back()
        <<",\"relative_residual\":"<<(initial?history.back()/initial:0)<<",\"maximum_residual\":"<<maximumResidual
        <<",\"restarts\":"<<restarts<<",\"krylov_products\":"<<products<<",\"preconditioner_applications\":"<<block.applications
        <<",\"velocity_preconditioner_symmetry_correction_max\":"<<block.asymmetry
        <<",\"velocity_matrix_asymmetry_max\":"<<block.asymmetry<<",\"velocity_preconditioner_entry_change_max\":"<<block.symmetryCorrection
        <<",\"viscous_pressure_mass_scale\":"<<nu<<",\"ilu0_factor_builds\":"<<block.velocity.ilu0Builds()<<",\"setup_seconds\":"<<setup
        <<",\"total_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"residual_history\":[";
    for(std::size_t i=0;i<history.size();++i){if(i)record<<',';record<<history[i];}record<<"]}\n";
    std::ofstream log(output+".json");log<<record.str();log.close();linearEnsure(bool(log),"result write failed");
    std::cout<<record.str();return converged?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
