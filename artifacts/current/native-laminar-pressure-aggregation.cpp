// Bounded comparison of one native fixed aggregation V-cycle with pressure IC0.
// Original saddle matrix K, RHS, velocity ILU0, GMRES and stopping gate are shared.
#define main archivedBlockProbeMain
#include "native-laminar-block-precondition.cpp"
#undef main
using cartmesh2d::fv::detail::compatible::linear::PressureSchurDiagonal;

void aggregationWiringControl() {
    const Matrix k(5,{{0,0,2},{1,1,3},{2,2,5},
        {0,3,1},{1,3,-1},{1,4,1},{2,4,-1},{3,0,1},{3,1,-1},{4,1,1},{4,2,-1}});
    const Vec r{.3,-.7};
    PressureSchurDiagonal pressure(k,k,3,LinearPressureMethod2D::Aggregation);
    Vec p(2);pressure.apply(r,0,p);
    const Vec residual{(1./2+1./3)*p[0]-p[1]/3+r[0],-p[0]/3+(1./3+1./5)*p[1]+r[1]};
    linearEnsure(linearNorm(residual)<=256*std::numeric_limits<double>::epsilon()*(1+linearNorm(r)),
        "Aggregation pressure inverse lost diagonal-velocity identity or saddle sign");
}

int main(int argc,char** argv)try {
    if(argc!=9)throw std::runtime_error("usage: pressure-aggregation MATRIX CELL_AREAS gauge|outlet nu ic0|aggregation relative-target maximum-restarts OUTPUT [MATRIX.picard.entries required]");
    aggregationWiringControl();
    const std::string prefix=argv[1],gauge=argv[3],method=argv[5],output=argv[8];
    const double nu=std::stod(argv[4]),target=std::stod(argv[6]);
    const auto maximum=std::stoull(argv[7]);
    linearEnsure((gauge=="gauge"||gauge=="outlet")&&(method=="ic0"||method=="aggregation")&&
        std::isfinite(nu)&&nu>0&&std::isfinite(target)&&target>0&&target<1&&maximum>0,"Invalid aggregation controls");
    linearEnsure(!std::filesystem::exists(output+".json")&&!std::filesystem::exists(output+".solution")&&
        !std::filesystem::exists(output+".candidate"),"Aggregation output already exists");
    const auto started=std::chrono::steady_clock::now();
    const auto rhs=binary<double>(prefix+".rhs"),area=areas(argv[2]);
    const Matrix k(rhs.size(),binary<Entry>(prefix+".entries")),pk(rhs.size(),binary<Entry>(prefix+".picard.entries"));
    // For gauge matrices the final pressure coordinate is already removed.
    // Block receives exactly the remaining pressure areas for its outlet-layout
    // diagonal Schur inverse; K and the removed gauge row are not changed.
    Vec retainedArea=area;if(gauge=="gauge")retainedArea.pop_back();
    Block block(k,retainedArea,"ilu0",method=="aggregation"?"outlet-schur-aggregation":"outlet-schur-diag",nu,0,pk);
    const auto& pressure=*block.pressureSchur;
    const double setup=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    const double initial=linearNorm(rhs);Vec x(rhs.size()),r=rhs;
    std::vector<double> history{initial};std::size_t products=0,restarts=0;bool converged=initial==0;
    while(!converged&&restarts<maximum) {
        const double remaining=target/(linearNorm(r)/initial);
        const auto direction=newtonKrylovDirection2D([&](const Vec& v){++products;return k.apply(block.apply(v));},r,60,remaining);
        if(!direction)break;
        const auto step=block.apply(*direction);
        for(std::size_t i=0;i<x.size();++i)x[i]+=step[i];
        // Same compensated true b-K*x recomputation as the current product API.
        r=k.residual(rhs,x);history.push_back(linearNorm(r));++restarts;
        converged=history.back()/initial<=target;
    }
    const double total=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::ofstream field(output+(converged?".solution":".candidate"),std::ios::binary);
    field.write(reinterpret_cast<const char*>(x.data()),static_cast<std::streamsize>(x.size()*sizeof(double)));
    field.close();linearEnsure(bool(field),"Aggregation field write failed");
    std::size_t positive=0,coefficientCount=0;double maximumPositive=0;
    for(const auto& [ij,value]:pressure.data.off)if(value>0){++positive;maximumPositive=std::max(maximumPositive,value);}
    std::ostringstream out;out<<std::setprecision(17)<<"{\"pressureInverse\":\""<<method<<"\",\"converged\":"<<(converged?"true":"false")
        <<",\"unknowns\":"<<k.n<<",\"pressureUnknowns\":"<<block.np<<",\"target\":"<<target
        <<",\"maximumRestarts\":"<<maximum<<",\"krylovDirections\":60,\"restarts\":"<<restarts<<",\"outerMatrixProducts\":"<<products
        <<",\"trueResidualEvaluations\":"<<restarts<<",\"preconditionerApplications\":"<<block.applications
        <<",\"velocityILUApplications\":"<<block.applications<<",\"pressureInverseApplications\":"<<block.applications
        <<",\"pressureSymmetryError\":"<<pressure.data.symmetryError<<",\"pressurePositiveOffDiagonalCount\":"<<positive
        <<",\"pressurePositiveOffDiagonalMaximum\":"<<maximumPositive<<",\"pressureCoefficientMaximum\":"<<pressure.data.coefficientMaximum
        <<",\"pressureFineCoefficients\":"<<pressure.system.diag.size()+pressure.system.off.size()<<",\"aggregationLevels\":[";
    if(pressure.aggregation)for(std::size_t i=0;i<pressure.aggregation->levels();++i) {
        const auto& level=pressure.aggregation->level(i);coefficientCount+=level.diagonal.size()+level.off.size();
        if(i)out<<',';out<<"{\"unknowns\":"<<level.diagonal.size()<<",\"coefficients\":"<<level.diagonal.size()+level.off.size()<<'}';
    }
    out<<"],\"aggregationTotalCoefficients\":"<<coefficientCount
        <<",\"relativeTrueResidual\":"<<(initial?history.back()/initial:0)<<",\"setupSecondsIncludingRead\":"<<setup
        <<",\"totalSeconds\":"<<total<<",\"compensatedTrueResidual\":true,\"trueResidualHistory\":[";
    for(std::size_t i=0;i<history.size();++i){if(i)out<<',';out<<history[i];}
    out<<"]}\n";std::ofstream log(output+".json");log<<out.str();log.close();linearEnsure(bool(log),"Aggregation record write failed");
    std::cout<<out.str();return converged?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
