#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/fv/HybridHeat2D.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numbers>
#include <algorithm>
#include <chrono>
#include <limits>
using namespace cartmesh2d;using namespace cartmesh2d::fv;
namespace {
struct Integrals {double temperature=0,inverse=0;};
Integrals analyticTriangle(Point2D a,Point2D b,Point2D c) {
    // Tensor Gauss-Legendre 8 on a Duffy triangle map, signed Jacobian.
    // This integrates the analytic reference only, never numerical fluxes.
    constexpr double x[]={-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498,
                          .1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
    constexpr double w[]={.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620,
                          .3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
    const double cross=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);Integrals out;
    for(unsigned i=0;i<8;++i)for(unsigned j=0;j<8;++j) {
        const double u=.5*(x[i]+1),v=.5*(x[j]+1);
        const Point2D p{a.x+u*(b.x-a.x)+u*v*(c.x-b.x),a.y+u*(b.y-a.y)+u*v*(c.y-b.y)};
        const double t=2+.2/std::log(2.)*std::log(std::hypot(p.x,p.y)/.5);
        const double weight=.25*w[i]*w[j]*u*cross;out.temperature+=weight*t;out.inverse+=weight/t;
    }
    return out;
}
Integrals refinedTriangle(Point2D a,Point2D b,Point2D c) {
    const Point2D ab{.5*(a.x+b.x),.5*(a.y+b.y)},ac{.5*(a.x+c.x),.5*(a.y+c.y)},bc{.5*(b.x+c.x),.5*(b.y+c.y)};
    Integrals result;for(const auto q:{analyticTriangle(a,ab,ac),analyticTriangle(ab,b,bc),analyticTriangle(ac,bc,c),analyticTriangle(ab,bc,ac)}){result.temperature+=q.temperature;result.inverse+=q.inverse;}
    return result;
}
}
int main(int argc,char** argv){try{
    if(argc!=6&&argc!=7)throw std::runtime_error("mesh prefix trace|isothermal|uniform end-time dt [cell-matrix|impulse=N] (0 0 for steady only)");
    const auto start=std::chrono::steady_clock::now();auto read=readCm2dTopology(argv[1]);if(!read.valid())throw std::runtime_error(read.error);
    auto mesh=makeFvMesh2D(read.topology);std::string prefix=argv[2],mode=argv[3];double end=std::stod(argv[4]),dt=std::stod(argv[5]);
    if(mode!="trace"&&mode!="isothermal"&&mode!="uniform")throw std::runtime_error("invalid wall mode");
    if(end<0||!std::isfinite(end)||(end>0&&(!(dt>0)||!std::isfinite(dt))))throw std::runtime_error("invalid physical time");
    const double k=.37,beta=.2/std::log(2.);auto exact=[&](Point2D p){return mode=="uniform"?2.1:2+beta*std::log(std::hypot(p.x,p.y)/.5);};
    std::vector<HeatBoundary2D> bc;for(std::size_t f=0;f<mesh.faces.size();++f)if(!mesh.faces[f].neighbour){const auto& face=mesh.faces[f];const bool inner=dot(Vector2D{face.centre.x,face.centre.y},face.areaVector)<0;bc.push_back({f,HeatBoundaryKind2D::Temperature,mode=="uniform"?2.1:mode=="trace"?exact(face.centre):inner?2:2.2,{}});}
    const auto assemble=std::chrono::steady_clock::now();HybridHeatOperator2D op(mesh,bc,k);const auto built=std::chrono::steady_clock::now();
    auto loose=op.solveSteady(1e-10),steady=op.solveSteady();double sensitivity=0;for(std::size_t i=0;i<mesh.cells.size();++i)sensitivity=std::max(sensitivity,std::abs(loose.temperature[i]-steady.temperature[i]));
    std::ofstream matrix(prefix+".trace-matrix.csv"),local(prefix+".local-matrix.csv");matrix<<std::setprecision(17)<<"row,column,value\n";local<<std::setprecision(17)<<"cell,row,column,value\n";
    op.visitSteadyTraceMatrix([&](auto i,auto j,double v){matrix<<i<<','<<j<<','<<v<<'\n';});op.visitLocalMatrices([&](auto c,auto i,auto j,double v){local<<c<<','<<i<<','<<j<<','<<v<<'\n';});matrix.close();local.close();
    const auto staticStart=std::chrono::steady_clock::now();auto fixedCheck=op.evaluateAtCells(steady.temperature);const double staticSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-staticStart).count();
    double staticCellDifference=0;for(std::size_t i=0;i<steady.temperature.size();++i)staticCellDifference=std::max(staticCellDifference,std::abs(fixedCheck.cellResidual[i]-steady.cellResidual[i]));
    double cellMatrixSeconds=0;std::size_t cellMatrixEvaluations=0;
    if(argc==7&&std::string(argv[6])=="cell-matrix"){
        const auto t0=std::chrono::steady_clock::now();std::ofstream cm(prefix+".cell-matrix.csv");cm<<std::setprecision(17)<<"row,column,value,area\n";
        for(std::size_t j=0;j<mesh.cells.size();++j){auto plus=steady.temperature,minus=steady.temperature;plus[j]+=.5;minus[j]-=.5;auto a=op.evaluateAtCells(plus),b=op.evaluateAtCells(minus);cellMatrixEvaluations+=2;for(std::size_t i=0;i<mesh.cells.size();++i)cm<<i<<','<<j<<','<<a.cellResidual[i]-b.cellResidual[i]<<','<<mesh.cells[i].area<<'\n';}
        cellMatrixSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    }
    double volume=0,error=0,meanError=0,harmonicError=0,near=0,nearArea=0;std::ofstream cells(prefix+".steady.cells.csv");cells<<std::setprecision(17)<<"cell,x,y,area,temperature,analyticCentroidTemperature,analyticMeanTemperature,analyticHarmonicTemperature\n";
    for(std::size_t i=0;i<mesh.cells.size();++i){const auto& c=mesh.cells[i];double t=steady.temperature[i];volume+=c.area;error+=c.area*std::abs(t-exact(c.centre));
        bool wall=std::any_of(c.faces.begin(),c.faces.end(),[&](auto f){return !mesh.faces[f].neighbour;});if(wall){near+=c.area*std::abs(t-exact(c.centre));nearArea+=c.area;}
        Integrals integral;const auto& poly=read.topology.cells[i].vertices;for(std::size_t j=0;j<poly.size();++j){auto a=read.topology.vertices[poly[j]].point,b=read.topology.vertices[poly[(j+1)%poly.size()]].point;auto q=refinedTriangle(c.centre,a,b);integral.inverse+=q.inverse;integral.temperature+=q.temperature;}
        double harmonic=mode=="uniform"?2.1:c.area/integral.inverse,mean=mode=="uniform"?2.1:integral.temperature/c.area;meanError+=c.area*std::abs(t-mean);harmonicError+=c.area*std::abs(t-harmonic);cells<<i<<','<<c.centre.x<<','<<c.centre.y<<','<<c.area<<','<<t<<','<<exact(c.centre)<<','<<mean<<','<<harmonic<<'\n';
    }
    std::ofstream faces(prefix+".steady.faces.csv");faces<<std::setprecision(17)<<"face,owner,neighbour,x,y,sx,sy,trace,heatFlux,analyticChordHeat\n";double fluxError=0,fluxMagnitude=0,wallHeat=0;
    for(std::size_t f=0;f<mesh.faces.size();++f){const auto& x=mesh.faces[f];double q=0;if(!x.neighbour){const Point2D a{x.centre.x+.5*x.areaVector.y,x.centre.y-.5*x.areaVector.x},b{x.centre.x-.5*x.areaVector.y,x.centre.y+.5*x.areaVector.x};q=mode=="uniform"?0:-k*beta*std::atan2(a.x*b.y-a.y*b.x,a.x*b.x+a.y*b.y);fluxError+=std::abs(steady.faceHeatFlux[f]-q);fluxMagnitude+=std::abs(q);wallHeat+=steady.faceHeatFlux[f];}
        faces<<f<<','<<x.owner<<','<<(x.neighbour?std::to_string(*x.neighbour):"-1")<<','<<x.centre.x<<','<<x.centre.y<<','<<x.areaVector.x<<','<<x.areaVector.y<<','<<steady.trace[f]<<','<<steady.faceHeatFlux[f]<<','<<q<<'\n';}
    std::vector<double> state=steady.temperature,cv(state.size(),1);for(std::size_t i=0;i<state.size();++i){auto p=mesh.cells[i].centre;double r=std::hypot(p.x,p.y);if(mode!="uniform")state[i]+=(.02*std::sin(7*std::atan2(p.y,p.x))+.005*std::cos(61*p.x)*std::sin(59*p.y))*std::sin(std::numbers::pi*std::log(r/.5)/std::log(2.));}
    if(argc==7&&std::string(argv[6])!="cell-matrix"){const std::string option=argv[6];if(option.rfind("impulse=",0)!=0)throw std::runtime_error("invalid diagnostic option");const auto colon=option.find(':');const auto cell=std::stoull(option.substr(8,colon==std::string::npos?colon:colon-8));const double amplitude=colon==std::string::npos?.1:std::stod(option.substr(colon+1));if(cell>=state.size()||!std::isfinite(amplitude)||amplitude<=0)throw std::runtime_error("invalid impulse cell/amplitude");state[cell]+=amplitude;}
    double admissibleMinimum=*std::min_element(state.begin(),state.end()),admissibleMaximum=*std::max_element(state.begin(),state.end());for(const auto& b:bc){admissibleMinimum=std::min(admissibleMinimum,b.value);admissibleMaximum=std::max(admissibleMaximum,b.value);}
    auto energy=[&](const std::vector<double>& t){double e=0;for(std::size_t i=0;i<t.size();++i)e+=.5*mesh.cells[i].area*(t[i]-steady.temperature[i])*(t[i]-steady.temperature[i]);return e;};double initialEnergy=energy(state),lastEnergy=initialEnergy,maxEnergyIncrease=0,maxGlobalBalance=0,minTemperature=*std::min_element(state.begin(),state.end()),maxTemperature=*std::max_element(state.begin(),state.end());
    std::ofstream history(prefix+".history.csv");history<<std::setprecision(17)<<"step,time,dt,perturbationEnergy,minimumTemperature,maximumTemperature,faceFluxJump,cellBalance,globalHeatBalance,iterations,localStageCellBalance,sharedStageCellBalance,stageToFinalK,stageToFinalBudgetRatio,acceptedEnergyIdentityDefect\n";
    const auto snapshot=[&](const std::string& name,const std::vector<double>& t,double time){std::ofstream f(prefix+name);f<<std::setprecision(17)<<"time,"<<time<<"\ncell,x,y,area,temperature,volumetricHeatCapacity\n";for(std::size_t i=0;i<t.size();++i)f<<i<<','<<mesh.cells[i].centre.x<<','<<mesh.cells[i].centre.y<<','<<mesh.cells[i].area<<','<<t[i]<<",1\n";};snapshot(".initial.csv",state,0);
    double time=0,maxStageToFinal=0,maxCorrectionRatio=0,maxSharedStageBalance=0,maxLocalStageBalance=0,maxEnergyIdentityDefect=0;
    std::vector<double> finalTrace=steady.trace,finalFlux=steady.faceHeatFlux,finalStage=state;std::size_t steps=0,totalIterations=0;
    while(time<end){double nextTime=std::min(end,(static_cast<double>(steps)+1)*dt);if(end-nextTime<=8*std::numeric_limits<double>::epsilon()*std::max(std::abs(end),std::abs(nextTime)))nextTime=end;const double h=nextTime-time;if(!(h>0)||steps>=100000)throw std::runtime_error("time/step budget exhausted; accepted checkpoint retained");auto candidate=op.backwardEuler(state,cv,h);const double nextEnergy=energy(candidate.temperature);maxEnergyIncrease=std::max(maxEnergyIncrease,nextEnergy-lastEnergy);
        maxStageToFinal=std::max(maxStageToFinal,candidate.maximumStageToFinalDifferenceK);maxCorrectionRatio=std::max(maxCorrectionRatio,candidate.maximumStageToFinalBudgetRatio);
        maxSharedStageBalance=std::max(maxSharedStageBalance,candidate.maximumSharedStageCellBalance);maxLocalStageBalance=std::max(maxLocalStageBalance,candidate.maximumLocalCellBalance);
        double perturbDissipation=0,temporalDefect=0;op.visitLocalMatrices([&](auto c,auto i,auto j,double a){const auto f=mesh.cells[c].faces[i],g=mesh.cells[c].faces[j];const double di=candidate.trace[f]-steady.trace[f]-candidate.stageTemperature[c]+steady.temperature[c],dj=candidate.trace[g]-steady.trace[g]-candidate.stageTemperature[c]+steady.temperature[c];perturbDissipation+=a*di*dj;});
        for(std::size_t i=0;i<state.size();++i)temporalDefect+=.5*mesh.cells[i].area*(candidate.temperature[i]-state[i])*(candidate.temperature[i]-state[i]);
        const double energyIdentityDefect=nextEnergy-lastEnergy+temporalDefect+h*perturbDissipation;maxEnergyIdentityDefect=std::max(maxEnergyIdentityDefect,std::abs(energyIdentityDefect));
        double balance=0;for(std::size_t i=0;i<state.size();++i)balance+=mesh.cells[i].area*(candidate.temperature[i]-state[i])/h;for(auto b:bc)balance+=candidate.faceHeatFlux[b.face];maxGlobalBalance=std::max(maxGlobalBalance,std::abs(balance));
        minTemperature=std::min(minTemperature,*std::min_element(candidate.temperature.begin(),candidate.temperature.end()));maxTemperature=std::max(maxTemperature,*std::max_element(candidate.temperature.begin(),candidate.temperature.end()));
        if(*std::min_element(candidate.temperature.begin(),candidate.temperature.end())<=0||*std::min_element(candidate.stageTemperature.begin(),candidate.stageTemperature.end())<=0){snapshot(".failed-stage.csv",candidate.stageTemperature,time+h);snapshot(".failed-candidate.csv",candidate.temperature,time+h);snapshot(".checkpoint.csv",state,time);throw std::runtime_error("nonpositive Kelvin candidate rejected; accepted checkpoint retained");}
        time=nextTime;++steps;totalIterations+=candidate.iterations;history<<steps<<','<<time<<','<<h<<','<<nextEnergy<<','<<*std::min_element(candidate.temperature.begin(),candidate.temperature.end())<<','<<*std::max_element(candidate.temperature.begin(),candidate.temperature.end())<<','<<candidate.maximumFaceFluxJump<<','<<candidate.maximumCellBalance<<','<<balance<<','<<candidate.iterations<<','<<candidate.maximumLocalCellBalance<<','<<candidate.maximumSharedStageCellBalance<<','<<candidate.maximumStageToFinalDifferenceK<<','<<candidate.maximumStageToFinalBudgetRatio<<','<<energyIdentityDefect<<'\n';finalTrace=candidate.trace;finalFlux=candidate.faceHeatFlux;finalStage=candidate.stageTemperature;state=std::move(candidate.temperature);lastEnergy=nextEnergy;snapshot(".checkpoint.csv",state,time);
    }
    snapshot(".final.csv",state,time);snapshot(".final-stage.csv",finalStage,time);
    std::ofstream transientFaces(prefix+".final.faces.csv");transientFaces<<std::setprecision(17)<<"face,traceAtImplicitStage,sharedHeatFlux\n";for(std::size_t f=0;f<mesh.faces.size();++f)transientFaces<<f<<','<<finalTrace[f]<<','<<finalFlux[f]<<'\n';
    const auto& d=op.diagnostics();std::ofstream report(prefix+".json");report<<std::setprecision(17)<<"{\"scope\":\"isolated native constant-property HMM scalar heat; not Euler coupled\",\"cells\":"<<mesh.cells.size()<<",\"traceUnknowns\":"<<steady.traceUnknowns<<",\"temperatureRelativeL1\":"<<error/volume/.2<<",\"temperatureCellMeanRelativeL1\":"<<meanError/volume/.2<<",\"temperatureHarmonicMeanRelativeL1\":"<<harmonicError/volume/.2<<",\"nearWallRelativeL1\":"<<near/nearArea/.2<<",\"wallFluxRelativeL1\":"<<(fluxMagnitude>0?fluxError/fluxMagnitude:0)<<",\"wallHeatBalance\":"<<wallHeat<<",\"steadyFaceFluxJump\":"<<steady.maximumFaceFluxJump<<",\"steadyCellBalance\":"<<steady.maximumCellBalance<<",\"steadyDissipation\":"<<steady.dissipation<<",\"traceResidualNorm\":"<<steady.traceResidualNorm<<",\"traceRhsNorm\":"<<steady.traceRhsNorm<<",\"steadyIterations\":"<<steady.iterations<<",\"toleranceSensitivityK\":"<<sensitivity<<",\"closureRoundoffRatio\":"<<d.closureRoundoffRatio<<",\"momentRoundoffRatio\":"<<d.momentRoundoffRatio<<",\"minimumNormalDistance\":"<<d.minimumNormalDistance<<",\"affineStabilizationResidual\":"<<d.maximumLocalAffineResidual<<",\"time\":"<<time<<",\"steps\":"<<steps<<",\"temporalOrder\":1,\"initialPerturbationEnergy\":"<<initialEnergy<<",\"finalPerturbationEnergy\":"<<lastEnergy<<",\"maximumPerturbationEnergyIncrease\":"<<maxEnergyIncrease<<",\"maximumGlobalHeatBalance\":"<<maxGlobalBalance<<",\"minimumTemperatureAllStates\":"<<minTemperature<<",\"maximumTemperatureAllStates\":"<<maxTemperature<<",\"maximumPrincipleUndershootK\":"<<std::max(0.,admissibleMinimum-minTemperature)<<",\"maximumPrincipleOvershootK\":"<<std::max(0.,maxTemperature-admissibleMaximum)<<",\"maximumStageToFinalDifferenceK\":"<<maxStageToFinal<<",\"maximumStageToFinalBudgetRatio\":"<<maxCorrectionRatio<<",\"maximumSharedStageCellBalance\":"<<maxSharedStageBalance<<",\"maximumLocalStageCellBalance\":"<<maxLocalStageBalance<<",\"maximumAcceptedEnergyIdentityDefect\":"<<maxEnergyIdentityDefect<<",\"transientIterations\":"<<totalIterations<<",\"staticTraceSolveSeconds\":"<<staticSeconds<<",\"staticCellResidualDifference\":"<<staticCellDifference<<",\"cellMatrixSeconds\":"<<cellMatrixSeconds<<",\"cellMatrixEvaluations\":"<<cellMatrixEvaluations<<",\"geometryBuildSeconds\":"<<std::chrono::duration<double>(built-assemble).count()<<",\"nativeFullSeconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
