// Research comparison on unchanged native saddle matrices. Standard scaled
// least-squares commutator, independently expressed in existing B/G/A blocks:
// L=B D^-1 G,  -S^-1 ~= -L^-1 B D^-1 A D^-1 G L^-1.
// https://petsc.org/release/manualpages/PC/PCLSC/
// Full versus upper block factorization is a separate comparison:
// https://petsc.org/release/manualpages/PC/PCFieldSplitSetSchurFactType/
// No new PDE, boundary condition, nonlinear acceptance or branch selection.
#define main archivedBlockProbeMain
#include "native-laminar-block-precondition.cpp"
#undef main
#include "cartmesh2d/fv/detail/CompatibleFlowElement2D.hpp"
using cartmesh2d::fv::detail::compatible::DenseLU;
using cartmesh2d::fv::detail::compatible::linear::PressureSchurDiagonal;

struct PressureApproximation {
    const Matrix& k;Block velocity;
    PressureSchurDiagonal laplacian;
    std::optional<DenseLU> exact;
    Vec inverseDiagonal;
    bool lsc,full;
    std::size_t applications=0,pressureSolves=0,velocityProducts=0,velocitySolves=0;
    PressureApproximation(const Matrix& matrix,const Matrix& pk,const Vec& area,bool outlet,double nu,bool commutator,bool dense,bool fullFactorization=false)
        :k(matrix),velocity(k,area,"ilu0",outlet?"outlet":"gauge",nu,0,pk),laplacian(k,pk,velocity.nv),inverseDiagonal(velocity.nv),lsc(commutator),full(fullFactorization) {
        for(std::size_t i=0;i<velocity.nv;++i)inverseDiagonal[i]=1/velocity.velocity.diag[i];
        if(dense){
            // Explicit small-reference control only: 256^2 doubles = 0.5 MiB
            // before factor copies; cubic factorization is not a scale backend.
            linearEnsure(velocity.np<=256,"Dense pressure reference limited to 256 unknowns");
            Vec a(velocity.np*velocity.np);
            for(std::size_t i=0;i<velocity.np;++i)a[i*velocity.np+i]=laplacian.data.diag[i];
            for(const auto& [ij,v]:laplacian.data.off)a[ij.first*velocity.np+ij.second]=v;
            exact.emplace(std::move(a),velocity.np);
        }
    }
    Vec solvePressure(const Vec& r){
        ++pressureSolves;if(exact)return exact->solve(r);
        laplacian.workspace.r=r;laplacian.system.precondition(laplacian.workspace,LinearPressureMethod2D::IC0);return laplacian.workspace.z;
    }
    Vec pressure(const Vec& r){
        const auto nv=velocity.nv;Vec y=solvePressure(r);
        if(lsc){
            Vec g(nv),a(nv),b(velocity.np);
            for(std::size_t i=0;i<nv;++i)for(auto p=k.rows[i];p<k.rows[i+1];++p)if(k.columns[p]>=nv)g[i]+=k.values[p]*y[k.columns[p]-nv];
            for(std::size_t i=0;i<nv;++i)g[i]*=inverseDiagonal[i];
            // Include the full nonsymmetric Newton velocity block. D and the
            // velocity ILU still come from the original Picard preconditioner.
            ++velocityProducts;
            for(std::size_t i=0;i<nv;++i)for(auto p=k.rows[i];p<k.rows[i+1];++p)if(k.columns[p]<nv)a[i]+=k.values[p]*g[k.columns[p]];
            for(std::size_t i=0;i<nv;++i)a[i]*=inverseDiagonal[i];
            for(std::size_t i=nv;i<k.n;++i)for(auto p=k.rows[i];p<k.rows[i+1];++p)if(k.columns[p]<nv)b[i-nv]+=k.values[p]*a[k.columns[p]];
            y=solvePressure(b);
        }
        for(auto& v:y)v=linearFinite(-v);
        return y;
    }
    Vec apply(const Vec& r){
        ++applications;const auto nv=velocity.nv;Vec rp(r.begin()+static_cast<std::ptrdiff_t>(nv),r.end());
        if(full){
            // Approximate LDU rather than upper block triangular: pressure
            // sees r_p - B Ahat^-1 r_u before the same pressure inverse.
            std::copy(r.begin(),r.begin()+static_cast<std::ptrdiff_t>(nv),velocity.workspace.r.begin());
            ++velocitySolves;velocity.velocity.preconditionILU0(velocity.workspace.r,velocity.workspace.z);
            for(std::size_t i=nv;i<k.n;++i)for(auto j=k.rows[i];j<k.rows[i+1];++j)if(k.columns[j]<nv)rp[i-nv]-=k.values[j]*velocity.workspace.z[k.columns[j]];
        }
        const auto p=pressure(rp);Vec z(r.size());
        std::copy(p.begin(),p.end(),z.begin()+static_cast<std::ptrdiff_t>(nv));
        for(std::size_t i=0;i<nv;++i){velocity.workspace.r[i]=r[i];for(auto j=k.rows[i];j<k.rows[i+1];++j)if(k.columns[j]>=nv)velocity.workspace.r[i]-=k.values[j]*z[k.columns[j]];}
        ++velocitySolves;velocity.velocity.preconditionILU0(velocity.workspace.r,velocity.workspace.z);std::copy(velocity.workspace.z.begin(),velocity.workspace.z.end(),z.begin());return z;
    }
};

void algebraControl(){
    std::vector<Entry> a{{0,0,2},{1,1,3},{2,2,5},{0,3,1},{1,3,1},{1,4,1},{2,4,1},{3,0,1},{3,1,1},{4,1,1},{4,2,1}};
    Matrix k(5,a);Vec area{1,1};PressureApproximation diagonal(k,k,area,true,1,false,true),lsc(k,k,area,true,1,true,true);
    const Vec r{.3,-.7};const auto x=diagonal.pressure(r),y=lsc.pressure(r);double difference=0;
    for(std::size_t i=0;i<x.size();++i)difference=std::max(difference,std::abs(x[i]-y[i]));
    linearEnsure(difference<=256*std::numeric_limits<double>::epsilon()*(1+linearNorm(x)),"LSC lost exact diagonal-velocity identity or pressure sign");
    PressureApproximation full(k,k,area,true,1,false,true,true);const Vec rhs{.1,.3,-.2,.2,-.1};
    const auto reconstructed=k.apply(full.apply(rhs));double residual=0;
    for(std::size_t i=0;i<rhs.size();++i)residual=std::max(residual,std::abs(reconstructed[i]-rhs[i]));
    linearEnsure(residual<=256*std::numeric_limits<double>::epsilon()*(1+linearNorm(rhs)),"Full block factorization lost the exact coupled inverse");
}

int main(int argc,char** argv)try{
    if(argc!=10)throw std::runtime_error("usage: pressure-lsc MATRIX CELL_AREAS gauge|outlet nu schur|schur-full|lsc ic0|dense relative-target maximum-restarts OUTPUT [Picard entries read from MATRIX.picard.entries]");
    algebraControl();const std::string prefix=argv[1],gauge=argv[3],mode=argv[5],inverse=argv[6],output=argv[9];
    const double nu=std::stod(argv[4]),target=std::stod(argv[7]);const auto maximum=std::stoull(argv[8]);
    linearEnsure((gauge=="gauge"||gauge=="outlet")&&(mode=="schur"||mode=="schur-full"||mode=="lsc")&&(inverse=="ic0"||inverse=="dense")&&nu>0&&std::isfinite(nu)&&target>0&&target<1&&std::isfinite(target)&&maximum>0,"invalid LSC controls");
    linearEnsure(!std::filesystem::exists(output+".json")&&!std::filesystem::exists(output+".solution")&&!std::filesystem::exists(output+".candidate"),"LSC output exists");
    const auto started=std::chrono::steady_clock::now();const auto rhs=binary<double>(prefix+".rhs"),area=areas(argv[2]);
    const Matrix k(rhs.size(),binary<Entry>(prefix+".entries")),pk(rhs.size(),binary<Entry>(prefix+".picard.entries"));
    PressureApproximation block(k,pk,area,gauge=="outlet",nu,mode=="lsc",inverse=="dense",mode=="schur-full");
    const double setup=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    const double initial=linearNorm(rhs);Vec x(rhs.size()),r=rhs;std::vector<double> history{initial};std::size_t products=0,restarts=0;bool converged=initial==0;
    while(!converged&&restarts<maximum){
        const double remaining=target/(linearNorm(r)/initial);
        const auto d=newtonKrylovDirection2D([&](const Vec& v){++products;return k.apply(block.apply(v));},r,60,remaining);if(!d)break;
        const auto step=block.apply(*d);for(std::size_t i=0;i<x.size();++i)x[i]+=step[i];
        const auto ax=k.apply(x);for(std::size_t i=0;i<r.size();++i)r[i]=rhs[i]-ax[i];history.push_back(linearNorm(r));++restarts;converged=history.back()/initial<=target;
    }
    const double total=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::ofstream field(output+(converged?".solution":".candidate"),std::ios::binary);field.write(reinterpret_cast<const char*>(x.data()),static_cast<std::streamsize>(x.size()*sizeof(double)));field.close();linearEnsure(bool(field),"LSC field write failed");
    std::ostringstream out;out<<std::setprecision(17)<<"{\"pressureApproximation\":\""<<mode<<"\",\"pressureInverse\":\""<<inverse<<"\",\"converged\":"<<(converged?"true":"false")
        <<",\"unknowns\":"<<k.n<<",\"pressureUnknowns\":"<<block.velocity.np<<",\"target\":"<<target<<",\"maximumRestarts\":"<<maximum<<",\"krylovDirections\":60,\"restarts\":"<<restarts<<",\"outerMatrixProducts\":"<<products<<",\"preconditionerApplications\":"<<block.applications
        <<",\"velocitySolves\":"<<block.velocitySolves<<",\"internalVelocityProducts\":"<<block.velocityProducts<<",\"pressureSolves\":"<<block.pressureSolves<<",\"relativeTrueResidual\":"<<(initial?history.back()/initial:0)<<",\"setupSecondsIncludingRead\":"<<setup<<",\"totalSeconds\":"<<total<<",\"trueResidualHistory\":[";
    for(std::size_t i=0;i<history.size();++i){if(i)out<<',';out<<history[i];}out<<"]}\n";std::ofstream log(output+".json");log<<out.str();log.close();linearEnsure(bool(log),"LSC record write failed");std::cout<<out.str();return converged?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
