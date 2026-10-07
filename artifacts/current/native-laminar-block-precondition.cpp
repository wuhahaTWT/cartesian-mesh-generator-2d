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
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>

using namespace cartmesh2d::fv::detail;
using Vec=LinearVector2D;
#include "cartmesh2d/fv/detail/CompatibleFlowLinear2D.hpp"
using cartmesh2d::fv::detail::compatible::linear::Entry;
using cartmesh2d::fv::detail::compatible::linear::Matrix;
using cartmesh2d::fv::detail::compatible::linear::Block;

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

int main(int argc,char** argv)try {
    if(argc<6 || argc>12)throw std::runtime_error("usage: probe input_prefix cell_csv ic0|ilu0|jacobi gauge|plain|outlet|outlet-oseen|outlet-schur-diag output_prefix [relative_tolerance=1e-11] [restarts=50] [viscosity=1] [transport_speed=0] [--velocity-preconditioner prefix]");
    const std::string prefix=argv[1],mode=argv[3],gauge=argv[4],output=argv[5];
    linearEnsure((mode=="ic0"||mode=="ilu0"||mode=="jacobi")&&(gauge=="gauge"||gauge=="plain"||gauge=="outlet"||gauge=="outlet-oseen"||gauge=="outlet-schur-diag"),"invalid preconditioner choice");
    const double tolerance=argc>=7?std::stod(argv[6]):1e-11;const auto maximum=argc>=8?std::stoull(argv[7]):50;
    const double nu=argc>=9?std::stod(argv[8]):1.;
    // Preserve the published positional transport-speed argument. The
    // independent velocity matrix uses an explicit flag to avoid ambiguity.
    int tail=9;double speed=0;std::string preconditionerPrefix;
    if(tail<argc&&std::string(argv[tail])!="--velocity-preconditioner") {
        const std::string value=argv[tail++];std::size_t used=0;speed=std::stod(value,&used);
        linearEnsure(used==value.size(),"invalid transport speed");
    }
    if(tail<argc) {
        linearEnsure(tail+2==argc&&std::string(argv[tail])=="--velocity-preconditioner","invalid velocity preconditioner option");
        preconditionerPrefix=argv[tail+1];linearEnsure(!preconditionerPrefix.empty(),"empty velocity preconditioner prefix");
    }
    linearEnsure(std::isfinite(nu)&&nu>0,"invalid viscous Schur scale");
    linearEnsure(std::isfinite(speed)&&speed>=0&&((gauge=="outlet-oseen")==(speed>0)),"outlet-oseen requires a positive transport speed; other modes require zero");
    linearEnsure(std::isfinite(tolerance)&&tolerance>0&&tolerance<=.01&&maximum>0,"invalid research solve control");
    linearEnsure(!std::filesystem::exists(output+".json") && !std::filesystem::exists(output+".solution") &&
        !std::filesystem::exists(output+".candidate"),"output exists; retain prior evidence and choose another prefix");
    const auto start=std::chrono::steady_clock::now();const auto rhs=binary<double>(prefix+".rhs"),area=areas(argv[2]);
    const bool outlet=gauge=="outlet"||gauge=="outlet-oseen"||gauge=="outlet-schur-diag";
    linearEnsure(rhs.size()>area.size()-(outlet?0:1),"invalid retained pressure dimension");
    for(double v:rhs)linearFinite(v);
    const Matrix matrix(rhs.size(),binary<Entry>(prefix+".entries"));
    std::optional<Matrix> preconditionerMatrix;if(!preconditionerPrefix.empty())preconditionerMatrix.emplace(rhs.size(),binary<Entry>(preconditionerPrefix+".entries"));
    Block block(matrix,area,mode,gauge,nu,speed,preconditionerMatrix?*preconditionerMatrix:matrix);
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
        <<",\"separate_velocity_preconditioner\":"<<(preconditionerMatrix?"true":"false")<<",\"velocity_off_nonzeros\":"<<block.velocity.off.size()<<",\"method\":\""<<mode<<"\",\"mass\":\""<<gauge
        <<"\",\"relative_tolerance\":"<<tolerance<<",\"rhs_norm\":"<<initial<<",\"residual_norm\":"<<history.back()
        <<",\"relative_residual\":"<<(initial?history.back()/initial:0)<<",\"maximum_residual\":"<<maximumResidual
        <<",\"restarts\":"<<restarts<<",\"krylov_products\":"<<products<<",\"preconditioner_applications\":"<<block.applications
        <<",\"velocity_preconditioner_symmetry_correction_max\":"<<block.asymmetry
        <<",\"velocity_matrix_asymmetry_max\":"<<block.matrixAsymmetry<<",\"velocity_preconditioner_input_asymmetry_max\":"<<block.asymmetry<<",\"velocity_preconditioner_entry_change_max\":"<<block.symmetryCorrection
        <<",\"viscous_pressure_mass_scale\":"<<nu<<",\"transport_speed\":"<<speed
        <<",\"pressure_mass_scale_min\":"<<(nu+(gauge=="outlet-oseen"?speed*std::sqrt(*std::min_element(area.begin(),area.end())):0))
        <<",\"pressure_mass_scale_max\":"<<(nu+(gauge=="outlet-oseen"?speed*std::sqrt(*std::max_element(area.begin(),area.end())):0))
        <<",\"pressure_schur_nonzeros\":"<<(block.pressureSchur?block.pressureSchur->system.diag.size()+block.pressureSchur->system.off.size():0)
        <<",\"pressure_schur_symmetry_error\":"<<(block.pressureSchur?block.pressureSchur->data.symmetryError:0)
        <<",\"pressure_schur_ic0_builds\":"<<(block.pressureSchur?block.pressureSchur->system.ic0Builds():0)
        <<",\"ilu0_factor_builds\":"<<block.velocity.ilu0Builds()<<",\"setup_seconds\":"<<setup
        <<",\"total_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"residual_history\":[";
    for(std::size_t i=0;i<history.size();++i){if(i)record<<',';record<<history[i];}record<<"]}\n";
    std::ofstream log(output+".json");log<<record.str();log.close();linearEnsure(bool(log),"result write failed");
    std::cout<<record.str();return converged?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
