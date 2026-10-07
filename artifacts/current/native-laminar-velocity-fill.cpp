// Standard one-level velocity fill, reusing the shared native factorization.
// Original saddle matrix/RHS, pressure inverse, GMRES and gate are unchanged.
#define main archivedBlockProbeMain
#include "native-laminar-block-precondition.cpp"
#undef main
using cartmesh2d::fv::detail::compatible::linear::PressureSchurDiagonal;

void fillWiringControl() {
    using cartmesh2d::fv::detail::compatible::linear::oneLevelFillConnections;
    // Nonsymmetric three-row star becomes an exact LU with one fill level.
    // This detects accidental symmetrization as well as missing numeric fill.
    const Matrix a(3,{{0,0,4},{1,1,6},{2,2,5},{0,1,-1},{1,0,-3},{0,2,-2},{2,0,-1}});
    const SparsePattern2D pattern(3,oneLevelFillConnections(a,3));SparseSystem2D system(pattern);
    for(std::size_t i=0;i<a.n;++i)for(auto j=a.rows[i];j<a.rows[i+1];++j)
        if(a.columns[j]==i)system.diag[i]=a.values[j];else system.off[pattern.slot(i,a.columns[j])]=a.values[j];
    system.factorILU0();const Vec exact{.5,-1,.75},rhs=a.apply(exact);Vec answer(3);
    system.preconditionILU0(rhs,answer);
    linearEnsure(linearNorm(a.residual(rhs,answer))<=256*std::numeric_limits<double>::epsilon()*(1+linearNorm(rhs)),
        "Native one-level ILU failed known nonsymmetric inverse");
    // A generated (1,2) edge must not then create level-2 edge (2,3).
    const Matrix chain(4,{{0,1,1},{0,2,1},{1,3,1}});
    const auto graph=oneLevelFillConnections(chain,4);
    const std::vector<std::pair<std::size_t,std::size_t>> expected{{0,1},{0,2},{1,2},{1,3}};
    linearEnsure(graph==expected,"ILU1 symbolic fill propagated a higher level");
}

int main(int argc,char** argv)try {
    if(argc!=10)throw std::runtime_error("usage: velocity-fill MATRIX CELL_AREAS gauge|outlet nu ilu0|ilu1 ic0|aggregation relative-target maximum-restarts OUTPUT [MATRIX.picard.entries required]");
    fillWiringControl();
    const std::string prefix=argv[1],gauge=argv[3],velocity=argv[5],method=argv[6],output=argv[9];
    const double nu=std::stod(argv[4]),target=std::stod(argv[7]);
    const auto maximum=std::stoull(argv[8]);
    linearEnsure((gauge=="gauge"||gauge=="outlet")&&(method=="ic0"||method=="aggregation")&&(velocity=="ilu0"||velocity=="ilu1")&&
        std::isfinite(nu)&&nu>0&&std::isfinite(target)&&target>0&&target<1&&maximum>0,"Invalid velocity-fill controls");
    linearEnsure(!std::filesystem::exists(output+".json")&&!std::filesystem::exists(output+".solution")&&
        !std::filesystem::exists(output+".candidate"),"Velocity-fill output already exists");
    const auto started=std::chrono::steady_clock::now();
    const auto rhs=binary<double>(prefix+".rhs"),area=areas(argv[2]);
    const Matrix k(rhs.size(),binary<Entry>(prefix+".entries")),pk(rhs.size(),binary<Entry>(prefix+".picard.entries"));
    // For gauge matrices the final pressure coordinate is already removed.
    // Block receives exactly the remaining pressure areas for its outlet-layout
    // diagonal Schur inverse; K and the removed gauge row are not changed.
    Vec retainedArea=area;if(gauge=="gauge")retainedArea.pop_back();
    Block block(k,retainedArea,velocity,method=="aggregation"?"outlet-schur-aggregation":"outlet-schur-diag",nu,0,pk);
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
    field.close();linearEnsure(bool(field),"Velocity-fill field write failed");
    std::size_t positive=0,coefficientCount=0;double maximumPositive=0;
    for(const auto& [ij,value]:pressure.data.off)if(value>0){++positive;maximumPositive=std::max(maximumPositive,value);}
    std::ostringstream out;out<<std::setprecision(17)<<"{\"pressureInverse\":\""<<method<<"\",\"converged\":"<<(converged?"true":"false")
        <<",\"unknowns\":"<<k.n<<",\"pressureUnknowns\":"<<block.np<<",\"target\":"<<target
        <<",\"maximumRestarts\":"<<maximum<<",\"krylovDirections\":60,\"restarts\":"<<restarts<<",\"outerMatrixProducts\":"<<products
        <<",\"trueResidualEvaluations\":"<<restarts<<",\"preconditionerApplications\":"<<block.applications
        <<",\"velocityInverse\":\""<<velocity<<"\",\"velocityOriginalGraphCoefficients\":"<<block.nv+2*cartmesh2d::fv::detail::compatible::linear::connections(pk,block.nv).size()
        <<",\"velocityFactorGraphCoefficients\":"<<block.velocity.diag.size()+block.velocity.off.size()<<",\"velocityILUApplications\":"<<block.applications<<",\"pressureInverseApplications\":"<<block.applications
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
    out<<"]}\n";std::ofstream log(output+".json");log<<out.str();log.close();linearEnsure(bool(log),"Velocity-fill record write failed");
    std::cout<<out.str();return converged?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
