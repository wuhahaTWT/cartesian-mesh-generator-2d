#include "cartmesh2d/fv/CompatibleIncompressible2D.hpp"
#include "fixtures/PolygonMesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
template<class F>void rejects(F function){bool rejected=false;try{function();}catch(const std::invalid_argument&){rejected=true;}check(rejected,"Incompatible/corrupt checkpoint accepted");}
FvMesh2D mesh(int n=3,double shear=.3,double L=2.3) {
    std::vector<Polygon2D> polygons;
    const auto point=[&](int i,int j){return Point2D{L*(double(i)/n+shear*double(j)/n),L*double(j)/n};};
    for(int j=0;j<n;++j)for(int i=0;i<n;++i)
        polygons.push_back({{point(i,j),point(i+1,j),point(i+1,j+1),point(i,j+1)}});
    return makeFvMesh2D(cartmesh2d::test::fromPolygons(polygons));
}
CompatibleFlowControls2D controls(const FvMesh2D& mesh) {
    constexpr double U=3.25,L=2.3;
    CompatibleFlowControls2D c;c.viscosity=.1*U*L;c.referenceVelocity=U;c.referenceLength=L;
    c.acceleration=[](Point2D p){return Vector2D{(2-p.x/L)*U*U/L,(-3-p.y/L)*U*U/L};};
    for(std::size_t f=0;f<mesh.faces.size();++f)if(!mesh.faces[f].neighbour) {
        CompatibleBoundary2D b;b.face=f;b.value=[](Point2D p){return Vector2D{-U*p.y/L,U*p.x/L};};c.boundaries.push_back(std::move(b));
    }
    return c;
}
std::string save(const CompatibleFlowCheckpoint2D& checkpoint) {
    std::ostringstream out;
    // The public stream's formatting must not change floating-point bits.
    out<<std::hex<<std::left<<std::showbase<<std::uppercase<<std::setprecision(3);
    writeCompatibleFlowCheckpoint2D(out,checkpoint);return out.str();
}
CompatibleFlowCheckpoint2D read(const std::string& text,const FvMesh2D& mesh) {
    std::istringstream in(text);return readCompatibleFlowCheckpoint2D(in,mesh);
}
bool same(const CompatibleFlowState2D& a,const CompatibleFlowState2D& b){return a.cells==b.cells&&a.faces==b.faces;}
}
int main()try {
    const auto domain=mesh();auto c=controls(domain);
    const auto full=solveCompatibleIncompressible2D(domain,c);
    check(full.converged()&&full.checkpoint,"Uninterrupted reference failed");
    auto limited=c;limited.maximumIterations=2;
    const auto partial=solveCompatibleIncompressible2D(domain,limited);
    check(partial.stop==CompatibleFlowStop2D::NonlinearBudget&&partial.checkpoint&&partial.checkpoint->acceptedIterations()==2,"Budget lost accepted checkpoint");
    const auto serialized=save(*partial.checkpoint);
    const auto checkpoint=read(serialized,domain);
    check(save(checkpoint)==serialized,"Checkpoint read/write changed bits");
    check(checkpoint.nextPseudoStep()==full.iterations[2].pseudoStep,"Next pseudo-step was reset");
    const auto resumed=resumeCompatibleIncompressible2D(domain,c,checkpoint);
    check(resumed.converged()&&resumed.resumed&&!resumed.seed&&resumed.acceptedIterationsBefore==2,"Resume not distinguished from seed");
    check(save(*resumed.checkpoint)==save(*full.checkpoint),"Split/full normalized fields, metrics or pseudo-step differ");
    check(same(*resumed.lastAccepted,*full.lastAccepted),"Split/full physical fields differ");
    check(resumed.iterations.size()+2==full.iterations.size(),"Resume changed iteration trajectory");
    for(std::size_t i=0;i<resumed.iterations.size();++i) {
        const auto& a=resumed.iterations[i];const auto& b=full.iterations[i+2];
        check(a.iteration==b.iteration&&a.pseudoStep==b.pseudoStep&&a.matrixProducts==b.matrixProducts&&
              a.linearInitialRelativeResidual==b.linearInitialRelativeResidual&&a.linearRelativeResidual==b.linearRelativeResidual&&
              a.metrics->residualNorm==b.metrics->residualNorm,"Split/full actual iteration differs");
    }
    auto reordered=c;std::reverse(reordered.boundaries.begin(),reordered.boundaries.end());
    check(save(*resumeCompatibleIncompressible2D(domain,reordered,checkpoint).checkpoint)==save(*full.checkpoint),"Boundary vector order changed physical identity");
    auto changed=c;changed.viscosity*=2;rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    changed=c;changed.referenceVelocity*=2;rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    changed=c;changed.referenceLength*=2;rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    changed=c;changed.quadratureOrder=8;rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    changed=c;changed.equation=CompatibleEquation2D::Stokes;rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    changed=c;changed.globalization=CompatibleGlobalization2D::Backtracking;rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    changed=c;changed.acceleration=[](Point2D){return Vector2D{1,2};};rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    changed=c;changed.boundaries.front().value=[](Point2D){return Vector2D{};};rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    changed=c;changed.boundaries.front().kind=CompatibleBoundaryKind2D::NormalVelocity;changed.boundaries.front().value={};
    rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    changed=c;changed.maximumPseudoStep=c.initialPseudoStep;
    rejects([&]{(void)resumeCompatibleIncompressible2D(domain,changed,checkpoint);});
    const auto distorted=mesh(3,.31);
    // Counts alone match; full geometry/native operator identity must reject.
    const auto wrongMesh=read(serialized,distorted);
    rejects([&]{(void)resumeCompatibleIncompressible2D(distorted,controls(distorted),wrongMesh);});
    rejects([&]{(void)read(serialized,mesh(2));});
    rejects([&]{(void)read(serialized.substr(0,serialized.size()/2),domain);});
    rejects([&]{(void)read(serialized+"EXTRA",domain);});
    rejects([&]{(void)read("{\"kind\":\"seed\"}",domain);});
    auto damaged=serialized;const auto state=damaged.find("STATE\n")+6;
    damaged[state+15]=damaged[state+15]=='0'?'1':'0';rejects([&]{(void)read(damaged,domain);});
    damaged=serialized;const auto accepted=damaged.find("ACCEPTED ")+9;damaged[accepted]='0';
    rejects([&]{(void)read(damaged,domain);});

    // Deliberately sub-roundoff dimensionless equation/change requests:
    // one attempt must re-evaluate convergence, not reuse a saved completion.
    // This negative lifecycle test adds no physical accuracy threshold.
    auto strict=c;strict.equationTolerance=1e-18;strict.stateTolerance=1e-18;strict.maximumIterations=1;
    const auto stricter=resumeCompatibleIncompressible2D(domain,strict,*full.checkpoint);
    check(!stricter.converged()&&stricter.resumed&&stricter.iterations.size()==1,"Saved convergence bypassed new numerical targets");

    auto linear=c;linear.maximumLinearRestarts=1;linear.krylovDirections=1;
    const auto failed=resumeCompatibleIncompressible2D(domain,linear,checkpoint);
    check(failed.stop==CompatibleFlowStop2D::LinearBudget&&same(*failed.lastAccepted,*partial.lastAccepted)&&save(*failed.checkpoint)==serialized,"Failed resume replaced accepted state");
    auto interrupted=c;std::size_t polls=0;interrupted.stopRequested=[&]{return ++polls==60;};
    const auto stopped=resumeCompatibleIncompressible2D(domain,interrupted,checkpoint);
    check(stopped.stop==CompatibleFlowStop2D::Cancelled&&stopped.resumed&&same(*stopped.lastAccepted,*partial.lastAccepted)&&save(*stopped.checkpoint)==serialized,"Cancellation lost resumed state");
    bool cancelAfterSave=false;std::string published;
    auto callback=c;callback.checkpointAccepted=[&](const auto& cp){published=save(cp);cancelAfterSave=true;};
    callback.stopRequested=[&]{return cancelAfterSave;};
    const auto cancelled=solveCompatibleIncompressible2D(domain,callback);
    check(cancelled.stop==CompatibleFlowStop2D::Cancelled&&published==save(*cancelled.checkpoint),"Publication does not correspond to accepted state");
    check(resumeCompatibleIncompressible2D(domain,c,read(published,domain)).converged(),"Published accepted checkpoint not restartable");

    // Supported binary64 subnormals remain exact without decimal underflow.
    auto tiny=c;tiny.equation=CompatibleEquation2D::Stokes;tiny.globalization=CompatibleGlobalization2D::Backtracking;
    tiny.referenceVelocity=1;tiny.referenceLength=1;tiny.acceleration={};
    const double amplitude=1024*std::numeric_limits<double>::denorm_min();
    for(auto& b:tiny.boundaries)b.value=[amplitude](Point2D){return Vector2D{1,amplitude};};
    const auto small=solveCompatibleIncompressible2D(domain,tiny);
    check(small.converged()&&small.checkpoint,"Tiny native field failed");
    bool subnormal=false;for(const auto& face:small.lastAccepted->faces)for(double x:face)subnormal=subnormal||(x!=0&&std::abs(x)<std::numeric_limits<double>::min());
    check(subnormal&&save(read(save(*small.checkpoint),domain))==save(*small.checkpoint),"Subnormal state did not roundtrip");
    std::cout<<"{\"splitFullBitIdentical\":true,\"referenceIterations\":"<<full.iterations.size()
        <<",\"splitAcceptedIterations\":"<<checkpoint.acceptedIterations()<<",\"checkpointBytes\":"<<serialized.size()
        <<",\"mismatchAndCorruptionRejected\":true,\"cancelAndLinearFailurePreserved\":true,\"subnormalRoundtrip\":true}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
