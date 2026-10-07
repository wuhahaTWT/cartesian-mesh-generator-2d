// Same native K, b and approximate inverse; compare restart targets only.
// Reuse binary/native-area readers and the shared GMRES/block implementation.
#define main archivedBlockProbeMain
#include "native-laminar-block-precondition.cpp"
#undef main

int main(int argc,char** argv)try {
    if(argc!=8&&argc!=9)throw std::runtime_error("usage: krylov-target matrix-prefix cell-csv gauge|outlet|outlet-schur-diag viscosity relative-target maximum-restarts output-prefix [velocity-prefix]");
    const std::string input=argv[1],pressure=argv[3],output=argv[7];
    const double nu=std::stod(argv[4]),tolerance=std::stod(argv[5]);const auto maximum=std::stoull(argv[6]);
    linearEnsure((pressure=="gauge"||pressure=="outlet"||pressure=="outlet-schur-diag")&&std::isfinite(nu)&&nu>0&&std::isfinite(tolerance)&&tolerance>0&&tolerance<1&&maximum>0,"invalid comparison controls");
    const auto readStart=std::chrono::steady_clock::now();
    const auto rhs=binary<double>(input+".rhs"),area=areas(argv[2]);
    const Matrix k(rhs.size(),binary<Entry>(input+".entries"));
    std::optional<Matrix> pk;if(argc==9)pk.emplace(rhs.size(),binary<Entry>(std::string(argv[8])+".entries"));
    const double readSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-readStart).count();
    for(double v:rhs)linearFinite(v);
    const double initial=linearNorm(rhs);
    for(bool strict:{false,true}) {
        const std::string prefix=output+(strict?".linear-target":".newton-forcing");
        linearEnsure(!std::filesystem::exists(prefix+".json")&&!std::filesystem::exists(prefix+".solution")&&!std::filesystem::exists(prefix+".candidate"),"comparison output exists");
        const auto start=std::chrono::steady_clock::now();Block block(k,area,"ilu0",pressure,nu,0,pk?*pk:k);
        const double setup=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        Vec x(rhs.size()),r=rhs;std::vector<double> history{initial};
        std::size_t products=0,restarts=0;bool converged=initial==0;
        while(!converged&&restarts<maximum) {
            const double remaining=strict?tolerance/(linearNorm(r)/initial):.1;
            const auto direction=newtonKrylovDirection2D([&](const Vec& v){++products;return k.apply(block.apply(v));},r,60,remaining);
            if(!direction)break;
            const auto step=block.apply(*direction);for(std::size_t i=0;i<x.size();++i)x[i]+=step[i];
            const auto applied=k.apply(x);for(std::size_t i=0;i<r.size();++i)r[i]=rhs[i]-applied[i];
            ++restarts;history.push_back(linearNorm(r));converged=history.back()/initial<=tolerance;
        }
        const double total=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::ofstream field(prefix+(converged?".solution":".candidate"),std::ios::binary);
        field.write(reinterpret_cast<const char*>(x.data()),static_cast<std::streamsize>(x.size()*sizeof(double)));field.close();linearEnsure(bool(field),"comparison field write failed");
        std::ostringstream record;record<<std::setprecision(17)<<"{\"mode\":\""<<(strict?"linear-target":"newton-forcing")<<"\",\"converged\":"<<(converged?"true":"false")
            <<",\"unknowns\":"<<rhs.size()<<",\"relativeTarget\":"<<tolerance<<",\"directionsPerRestart\":60,\"maximumRestarts\":"<<maximum
            <<",\"products\":"<<products<<",\"restarts\":"<<restarts<<",\"preconditionerApplications\":"<<block.applications
            <<",\"relativeTrueResidual\":"<<(initial?history.back()/initial:0)<<",\"matrixReadSeconds\":"<<readSeconds
            <<",\"setupSeconds\":"<<setup<<",\"solveAndSetupSeconds\":"<<total<<",\"trueResidualHistory\":[";
        for(std::size_t i=0;i<history.size();++i){if(i)record<<',';record<<history[i];}record<<"]}\n";
        std::ofstream log(prefix+".json");log<<record.str();log.close();linearEnsure(bool(log),"comparison log write failed");std::cout<<record.str();
    }
    // A comparison completes even when an individual linear attempt failed.
    // Its .candidate is never a recovered/accepted flow state.
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
