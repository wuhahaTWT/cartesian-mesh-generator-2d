// Replay a frozen production pressure matrix using the existing native/GPU
// solvers. No geometric/PDE reconstruction and no alternative acceptance gate.
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/grid/RectilinearMesh2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace cartmesh2d::fv::detail;
namespace {
using Clock=std::chrono::steady_clock;
struct Input {std::vector<std::size_t> rows,columns;LinearVector2D diagonal,off,rhs,initial;double tolerance=1e-11;};
template<class T> void read(std::ifstream& in,std::vector<T>& values,std::size_t count) {
    values.resize(count);in.read(reinterpret_cast<char*>(values.data()),static_cast<std::streamsize>(count*sizeof(T)));
}
Input snapshot(const std::filesystem::path& path) {
    std::ifstream in(path,std::ios::binary);in.exceptions(std::ios::badbit|std::ios::failbit);
    char magic[8];in.read(magic,8);linearEnsure(std::string(magic,8)=="CM2DPRS1","Pressure snapshot format mismatch");
    std::uint64_t n=0,nnz=0;Input input;
    in.read(reinterpret_cast<char*>(&n),8);in.read(reinterpret_cast<char*>(&nnz),8);
    in.read(reinterpret_cast<char*>(&input.tolerance),8);
    const auto bytes=std::filesystem::file_size(path);
    linearEnsure(n>0 && n<=bytes/32 && nnz<=bytes/16,"Pressure snapshot dimensions exceed file size");
    linearEnsure(bytes==32+8*(4*n+1+2*nnz),"Pressure snapshot length mismatch");
    std::vector<std::uint64_t> rows,columns;read(in,rows,static_cast<std::size_t>(n+1));read(in,columns,static_cast<std::size_t>(nnz));
    input.rows.assign(rows.begin(),rows.end());input.columns.assign(columns.begin(),columns.end());
    read(in,input.diagonal,static_cast<std::size_t>(n));read(in,input.off,static_cast<std::size_t>(nnz));
    read(in,input.rhs,static_cast<std::size_t>(n));read(in,input.initial,static_cast<std::size_t>(n));
    linearEnsure(input.rows.front()==0 && input.rows.back()==nnz,"Pressure snapshot CSR invalid");
    return input;
}
Input grid(std::size_t nx,std::size_t ny) {
    linearEnsure(nx>1 && ny>1 && nx<=std::numeric_limits<std::size_t>::max()/ny,"Grid dimensions invalid");
    Input input;const auto n=nx*ny;
    input.diagonal.assign(n,4.);input.initial.assign(n,0.);input.rhs.resize(n);input.rows.push_back(0);
    for(std::size_t i=0;i<n;++i) {
        std::vector<std::size_t> columns;
        if(i>=nx)columns.push_back(i-nx);if(i%nx)columns.push_back(i-1);
        if(i%nx+1<nx)columns.push_back(i+1);if(i+nx<n)columns.push_back(i+nx);
        for(auto j:columns){input.columns.push_back(j);input.off.push_back(-1.);}
        input.rows.push_back(input.columns.size());
        // Deterministic multi-frequency forcing, not a single CG eigenvector.
        input.rhs[i]=1.+.2*std::sin(.371*static_cast<double>(i))+.1*std::cos(.193*static_cast<double>(i/nx));
    }
    return input;
}
void captureChannel(std::size_t nx,std::size_t ny,const std::string& path) {
    using namespace cartmesh2d;
    linearEnsure(nx>1 && ny>1 && nx<std::numeric_limits<std::size_t>::max() &&
        ny<std::numeric_limits<std::size_t>::max() && nx<=std::numeric_limits<std::size_t>::max()/ny,
        "Channel grid dimensions invalid");
    std::vector<double> x(nx+1),y(ny+1);
    for(std::size_t i=0;i<=nx;++i)x[i]=2.*static_cast<double>(i)/static_cast<double>(nx);
    for(std::size_t j=0;j<=ny;++j)y[j]=static_cast<double>(j)/static_cast<double>(ny);
    const auto topology=makeRectilinearMesh2D(x,y);
    const auto quality=evaluateSolverQuality2D(topology);
    linearEnsure(quality.valid(),"Production channel grid failed original Solver quality gate");
    const auto mesh=fv::makeFvMesh2D(topology);
    fv::FlowControls2D controls;controls.scenario="channel";controls.nu=.1;controls.speed=.1;
    controls.maxIterations=1;controls.tolerance=1e-6;
    setPressureCapturePath2D(path);
    const auto result=fv::solveIncompressible2D(mesh,controls);
    linearEnsure(std::filesystem::exists(std::filesystem::u8path(path)),"Production channel did not capture a pressure solve");
    std::cout<<"captured production channel pressure: cells="<<mesh.cells.size()
        <<" CFD outer steps="<<result.history.back().iteration<<" CFD converged="<<result.converged
        <<" (matrix source only, not a converged large CFD case)\n"<<std::flush;
}
void string(std::ostream& out,const std::string& text) {
    out<<'"';for(unsigned char c:text){if(c=='"'||c=='\\')out<<'\\'<<static_cast<char>(c);
        else if(c=='\n')out<<"\\n";else if(c=='\r')out<<"\\r";else if(c=='\t')out<<"\\t";
        else if(c>=32)out<<static_cast<char>(c);}out<<'"';
}
}
int main(int argc,char** argv) {
    try {
        std::string file,output;std::size_t nx=0,ny=0,repeats=3;bool gpu=false,channel=false;
        const auto count=[](const std::string& text) {
            linearEnsure(!text.empty() && std::all_of(text.begin(),text.end(),[](char c){return c>='0' && c<='9';}),
                         "Replay counts must be unsigned decimal integers");
            const auto value=std::stoull(text);
            linearEnsure(value<=std::numeric_limits<std::size_t>::max(),"Replay count exceeds host index range");
            return static_cast<std::size_t>(value);
        };
        for(int i=1;i<argc;++i) {
            const std::string a=argv[i];if(a=="--gpu"){gpu=true;continue;}
            if(a=="--channel"){channel=true;continue;}
            if(i+1>=argc)throw std::runtime_error("Missing replay argument");const std::string v=argv[++i];
            if(a=="--snapshot")file=v;else if(a=="--output")output=v;else if(a=="--nx")nx=count(v);
            else if(a=="--ny")ny=count(v);else if(a=="--repeats")repeats=count(v);
            else throw std::runtime_error("Unknown replay option");
        }
        linearEnsure(!output.empty() && repeats>0,"Replay output/repeat control missing");
        if(channel){linearEnsure(file.empty(),"Choose snapshot or channel, not both");file=output+".pressure.bin";captureChannel(nx,ny,file);}
        const auto input=file.empty()?grid(nx,ny):snapshot(std::filesystem::u8path(file));
        const auto n=input.diagonal.size();std::vector<std::pair<std::size_t,std::size_t>> connections;
        for(std::size_t i=0;i<n;++i) {
            linearEnsure(input.rows[i]<=input.rows[i+1] && input.rows[i+1]<=input.off.size(),"Replay CSR row invalid");
            for(auto k=input.rows[i];k<input.rows[i+1];++k)if(input.columns[k]>i)connections.emplace_back(i,input.columns[k]);
        }
        SparsePattern2D pattern(n,connections);
        linearEnsure(pattern.rows==input.rows && pattern.columns==input.columns,"Replay CSR not canonical/symmetric");
        const double target=1e-13+input.tolerance*linearNorm(input.rhs);
        std::ofstream out(std::filesystem::u8path(output));out.exceptions(std::ios::badbit|std::ios::failbit);out<<std::setprecision(17);
        out<<"{\"format\":\"cartmesh2d-pressure-replay-v1\",\"source\":";string(out,file.empty()?"synthetic-dirichlet-poisson":file);
        out<<",\"cells\":"<<n<<",\"offDiagonalNnz\":"<<input.off.size()<<",\"relativeTolerance\":"<<input.tolerance
            <<",\"nativeAbsoluteTarget\":"<<target<<",\"scope\":\"frozen linear system; not full CFD speed or physical validation\",\"runs\":[\n";
        bool comma=false;LinearVector2D reference;
        for(bool device:{false,true}) {
            if(device && !gpu)continue;
            for(auto method:{LinearPressureMethod2D::IC0,LinearPressureMethod2D::Jacobi,LinearPressureMethod2D::Aggregation}) {
                if(device && method==LinearPressureMethod2D::Aggregation)continue;
                if(comma)out<<",\n";comma=true;
                const char* name=method==LinearPressureMethod2D::IC0?"ic0":method==LinearPressureMethod2D::Jacobi?"jacobi":"aggregation";
                out<<"{\"backend\":\""<<(device?"torch-cuda":"cpu")<<"\",\"preconditioner\":\""<<name<<'"';
                try {
                    auto start=Clock::now();SparseSystem2D system(pattern);system.diag=input.diagonal;system.off=input.off;system.rhs=input.rhs;
                    LinearWorkspace2D workspace(n);auto field=input.initial;
                    const auto coldIterations=system.solvePressure(field,workspace,method,input.tolerance,device);
                    const double coldSeconds=std::chrono::duration<double>(Clock::now()-start).count();
                    std::vector<double> times;std::size_t iterations=0;
                    for(std::size_t repeat=0;repeat<repeats;++repeat) {
                        field=input.initial;start=Clock::now();
                        iterations=system.solvePressure(field,workspace,method,input.tolerance,device);
                        times.push_back(std::chrono::duration<double>(Clock::now()-start).count());
                    }
                    system.apply(field,workspace.ax);for(std::size_t i=0;i<n;++i)workspace.r[i]=input.rhs[i]-workspace.ax[i];
                    const double residual=linearNorm(workspace.r);linearEnsure(residual<=target,"Replay failed original native residual gate");
                    if(reference.empty())reference=field;double difference=0;
                    for(std::size_t i=0;i<n;++i)difference=std::max(difference,std::abs(field[i]-reference[i]));
                    auto sorted=times;std::sort(sorted.begin(),sorted.end());const auto& stats=system.gpuStats();
                    out<<",\"accepted\":true,\"coldSeconds\":"<<coldSeconds<<",\"warmMedianSeconds\":"<<sorted[sorted.size()/2]
                        <<",\"coldIterations\":"<<coldIterations<<",\"iterations\":"<<iterations<<",\"trueResidualNorm\":"<<residual
                        <<",\"maxAbsDifferenceFromNativeIC0\":"<<difference<<",\"gpuMatrixUploads\":"<<stats.matrixUploads
                        <<",\"gpuStartupSeconds\":"<<stats.startupSeconds<<",\"gpuUploadSeconds\":"<<stats.uploadSeconds
                        <<",\"gpuIterationSeconds\":"<<stats.iterationSeconds<<",\"gpuDownloadSeconds\":"<<stats.downloadSeconds<<",\"warmSeconds\":[";
                    for(std::size_t i=0;i<times.size();++i){if(i)out<<',';out<<times[i];}out<<"]}";
                    std::cout<<(device?"torch-cuda":"cpu")<<' '<<name<<" iterations="<<iterations<<" median="<<sorted[sorted.size()/2]<<'\n'<<std::flush;
                } catch(const std::exception& e){out<<",\"accepted\":false,\"error\":";string(out,e.what());out<<'}';std::cout<<"FAILED "<<name<<": "<<e.what()<<'\n'<<std::flush;}
                out.flush();
            }
        }
        out<<"\n]}\n";
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
