#include "cartmesh2d/fv/ManufacturedFlow2D.hpp"
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/fv/FlowCheckpoint2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/fv/FlowTimeStep2D.hpp"
#include "cartmesh2d/fv/FlowInitialization2D.hpp"
#include "cartmesh2d/fv/TaylorGreen2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cartmesh2d;

namespace {

double number(const std::string& s) {
    std::size_t n = 0;
    const auto v = std::stod(s, &n);
    if (n != s.size() || !std::isfinite(v)) {
        throw std::invalid_argument("expected finite number");
    }
    return v;
}

std::vector<double> viscosityCsv(const std::string& path,std::size_t count) {
    std::ifstream in(path); if(!in)throw std::runtime_error("cannot open face viscosity CSV");
    std::string line; std::getline(in,line);
    if(!line.empty()&&line.back()=='\r')line.pop_back();
    if(line!="face,viscosity")throw std::runtime_error("viscosity CSV header must be face,viscosity");
    std::vector<double> values(count);std::vector<bool> seen(count,false);
    while(std::getline(in,line)) {
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        const auto comma=line.find(',');
        if(comma==std::string::npos || line.find(',',comma+1)!=std::string::npos)
            throw std::runtime_error("viscosity CSV row requires two fields");
        const double index=number(line.substr(0,comma));
        if(index<0 || index>=static_cast<double>(count) || index!=std::floor(index))
            throw std::runtime_error("invalid viscosity face ID");
        const auto id=static_cast<std::size_t>(index);
        if(seen[id])throw std::runtime_error("duplicate viscosity face");
        values[id]=number(line.substr(comma+1));
        if(values[id]<=0)throw std::runtime_error("face viscosity must be positive");
        seen[id]=true;
    }
    if(!in.eof())throw std::runtime_error("cannot read viscosity CSV");
    if(!std::all_of(seen.begin(),seen.end(),[](bool value){return value;}))
        throw std::runtime_error("missing viscosity face");
    return values;
}

std::ofstream out(const std::string& p, const char* ext) {
    std::ofstream s(p + ext);
    s.exceptions(std::ios::badbit | std::ios::failbit);
    s << std::setprecision(17);
    return s;
}

void progress(const fv::FlowIteration2D& h) {
    std::cout << std::setprecision(17)
              << "{\"type\":\"flow-progress\",\"iteration\":" << h.iteration
              << ",\"momentumResidual\":" << h.momentumResidual
              << ",\"continuity\":" << h.continuity
              << ",\"velocityChange\":" << h.velocityChange
              << ",\"pressureChange\":" << h.pressureChange << "}" << std::endl;
}

}

int main(int argc, char** argv) {
    std::string prefix;
    double acceptedTime=0;
    bool transientOutputStarted=false;
    try {
        std::string path,viscosityPath,boundaryPath,boundaryExportPath;
        fv::FlowControls2D controls;
        bool explicitVelocityRelaxation=false;
        double timeStep=0;
        std::size_t requestedSteps=0,completedSteps=0;
        fv::FlowTimeStepControls2D adaptiveControls;
        bool adaptive=false, adaptiveOptions=false, explicitMinimumStep=false;
        double startTime=0;
        std::size_t attemptCount=0,rejectedSteps=0;
        std::string restart;
        fv::FlowInitialVortex2D initialVortex;
        unsigned vortexOptions=0;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--profile") {
                controls.profile = true;
                continue;
            }
            if (a == "--help") {
                std::cout
                    << "Native 2D incompressible laminar SIMPLE (experimental)\n"
            "--mesh FINAL.solver.cm2d --output PREFIX --case external|channel|duct|custom|cavity|manufactured|counterflow\n"
            "--case custom --boundary FILE: named, mesh-bound velocity inlet/pressure outlet/pressure opening/symmetry/wall conditions.\n"
            "--export-boundaries FILE: export channel/duct/cavity/annulus preset without solving; --output optional.\n"
            "--nu 0.01 --speed 1 --max-iterations 1500 --tolerance 1e-6\n"
            "--face-viscosity NU.csv: face,viscosity; all faces, positive kinematic nu; physical cases only.\n"
            "--manufactured-viscosity-slope 0: verification nu(x)=nu*(1+slope*x), steady only.\n"
            "--time-step DT --steps N: backward Euler physical time, converged SIMPLE at each step.\n"
            "--time-step MAX_DT --end-time T: adaptive backward Euler to absolute physical time T.\n"
            "--max-courant 1 --min-time-step MAX_DT/1024 --max-step-retries 10 --max-time-steps 100000: adaptive limits.\n"
            "--velocity-relaxation 0.6: transient inner iterations only; (0,1], larger may be unstable.\n"
            "--restart PREFIX.checkpoint: resume accepted state on identical mesh and physical setup.\n"
            "--case taylor-green: unforced exact slip-box decay; transient verification only.\n"
            "Transient physical cases start at rest; boundary velocities switch on for t>0.\n"
            "--initial-vortex-x X --initial-vortex-y Y --initial-vortex-radius R --initial-vortex-speed V:\n"
            "optional compact interior vortex at t=0, signed peak speed (positive CCW); all four required, no restart.\n"
            "Fixed-step mode reports CFL. Adaptive mode retries unaccepted steps without advancing the saved state.\n"
            "--profile writes extra .performance.json timing/linear iteration diagnostics.\n"
            "--pressure-preconditioner ic0|jacobi|aggregation (default ic0); aggregation experimental; same true-residual tolerance.\n"
            "--viscous-stress symmetric|laplacian (default symmetric); conservative Newtonian stress.\n"
            "--convection upwind|limited-linear|face-limited-linear (default upwind); bounded face reconstruction.\n"
            "--outlet-backflow reject|normal-inlet (default reject).\n"
            "manufactured: unit-square analytic forced vortex; verification only, stationary walls.\n"
            "--manufactured-pressure-slope 0: add Uref^2*slope*(x+y) to the analytic pressure.\n"
            "channel speed=maximum parabolic inlet speed; cavity speed=lid speed.\n"
            "duct: uniform left inlet, right p=0, arbitrary no-slip walls; both openings must be vertical x-extrema.\n"
            "custom: velocity inlets, pressure outlets and walls allow arbitrary orientation; speed is a reference scale only.\n"
            "pressure-opening: axis-aligned static p/rho; normal velocity free, incoming tangential velocity zero.\n"
            "symmetry: axis-aligned impermeable free-slip face, no prescribed velocity or pressure.\n"
            "Ordinary pressure-outlet retains backflow rejection.\n"
            "Other presets require fixed axis-aligned rectangular outer boundaries. Pressure is kinematic.\n"
            "No turbulence/compressibility; outlet backflow policy is explicit.\n";
                return 0;
            }
            if (i + 1 >= argc) {
                throw std::invalid_argument("missing option value");
            }
            std::string v = argv[++i];
            if (a == "--mesh") {
                path = v;
            } else if (a == "--output") {
                prefix = v;
            } else if (a == "--case") {
                controls.scenario = v;
            } else if (a == "--boundary") {
                boundaryPath = v;
            } else if (a == "--export-boundaries") {
                boundaryExportPath = v;
            } else if (a == "--nu") {
                controls.nu = number(v);
            } else if (a == "--face-viscosity") {
                viscosityPath=v;
            } else if (a == "--manufactured-viscosity-slope") {
                controls.manufacturedViscositySlope=number(v);
            } else if (a == "--speed") {
                controls.speed = number(v);
            } else if (a == "--tolerance") {
                controls.tolerance = number(v);
            } else if (a == "--velocity-relaxation") {
                controls.velocityRelaxation=number(v);
                if (!(controls.velocityRelaxation>0 && controls.velocityRelaxation<=1))
                    throw std::invalid_argument("velocity-relaxation must be in (0,1]");
                explicitVelocityRelaxation=true;
            } else if (a == "--time-step") {
                timeStep=number(v);
                if (!(timeStep>0)) throw std::invalid_argument("time-step must be positive");
            } else if (a == "--steps") {
                const double n=number(v);
                if (n<1 || n>1000000 || n!=std::floor(n)) throw std::invalid_argument("bad physical step count");
                requestedSteps=static_cast<std::size_t>(n);
            } else if (a == "--end-time") {
                adaptive=true;adaptiveControls.targetTime=number(v);
            } else if (a == "--max-courant") {
                adaptiveOptions=true;adaptiveControls.maximumCourant=number(v);
            } else if (a == "--min-time-step") {
                adaptiveOptions=true;explicitMinimumStep=true;adaptiveControls.minimumStep=number(v);
            } else if (a == "--max-step-retries" || a == "--max-time-steps") {
                adaptiveOptions=true;
                const double n=number(v);
                if (n<0 || n>1000000 || n!=std::floor(n)) throw std::invalid_argument("bad adaptive iteration budget");
                if (a == "--max-step-retries") adaptiveControls.maximumRetries=static_cast<std::size_t>(n);
                else adaptiveControls.maximumAcceptedSteps=static_cast<std::size_t>(n);
            } else if (a == "--restart") {
                restart=v;
            } else if (a == "--initial-vortex-x" || a == "--initial-vortex-y" ||
                       a == "--initial-vortex-radius" || a == "--initial-vortex-speed") {
                const unsigned bit=a=="--initial-vortex-x"?1:a=="--initial-vortex-y"?2:a=="--initial-vortex-radius"?4:8;
                if (vortexOptions&bit) throw std::invalid_argument("duplicate initial vortex option");
                vortexOptions|=bit;
                const double value=number(v);
                if (bit==1) initialVortex.centre.x=value;
                else if (bit==2) initialVortex.centre.y=value;
                else if (bit==4) initialVortex.radius=value;
                else initialVortex.peakSpeed=value;
            } else if (a == "--manufactured-pressure-slope") {
                controls.manufacturedPressureSlope = number(v);
            } else if (a == "--viscous-stress") {
                if (v != "symmetric" && v != "laplacian") throw std::invalid_argument("viscous-stress must be symmetric or laplacian");
                controls.viscousStress=v=="symmetric"?fv::ViscousStress2D::Symmetric:fv::ViscousStress2D::Laplacian;
            } else if (a == "--convection") {
                if (v != "upwind" && v != "limited-linear" && v != "face-limited-linear")
                    throw std::invalid_argument("convection must be upwind, limited-linear or face-limited-linear");
                controls.convection = v == "face-limited-linear" ? fv::ConvectionScheme2D::FaceLimitedLinearUpwind : v == "limited-linear"
                    ? fv::ConvectionScheme2D::LimitedLinearUpwind : fv::ConvectionScheme2D::Upwind;
            } else if (a == "--outlet-backflow") {
                if (v != "reject" && v != "normal-inlet")
                    throw std::invalid_argument("outlet-backflow must be reject or normal-inlet");
                controls.outletBackflow = v == "normal-inlet"
                    ? fv::OutletBackflow2D::NormalInlet : fv::OutletBackflow2D::Reject;
            } else if (a == "--pressure-preconditioner") {
                if (v != "ic0" && v != "jacobi" && v != "aggregation") {
                    throw std::invalid_argument("pressure preconditioner must be ic0, jacobi or aggregation");
                }
                controls.pressurePreconditioner = v == "aggregation" ? fv::PressurePreconditioner2D::Aggregation : v == "ic0"
                    ? fv::PressurePreconditioner2D::IncompleteCholesky0
                    : fv::PressurePreconditioner2D::Jacobi;
            } else if (a == "--max-iterations") {
                double n = number(v);
                if (n < 1 || n > 100000 || n != std::floor(n)) {
                    throw std::invalid_argument("bad iteration limit");
                }
                controls.maxIterations = static_cast<std::size_t>(n);
            } else {
                throw std::invalid_argument("unknown option " + a);
            }
        }
        if (path.empty() || (prefix.empty() && boundaryExportPath.empty())) {
            throw std::invalid_argument("--mesh and --output required");
        }
        if (!boundaryExportPath.empty() && (timeStep > 0 || requestedSteps > 0 || adaptive || adaptiveOptions || !restart.empty()))
            throw std::invalid_argument("boundary template export does not run time steps or restart");
        if (!path.ends_with(".solver.cm2d") || path.ends_with(".failed.solver.cm2d")) {
            throw std::invalid_argument("requires final *.solver.cm2d");
        }

        if (adaptive) {
            if (timeStep<=0 || requestedSteps) throw std::invalid_argument("--end-time requires --time-step and cannot be combined with --steps");
            adaptiveControls.maximumStep=timeStep;
            if (!explicitMinimumStep) adaptiveControls.minimumStep=timeStep/1024;
            fv::validateFlowTimeStepControls2D(adaptiveControls);
        } else if (adaptiveOptions || (timeStep>0)!=(requestedSteps>0) || (!restart.empty() && timeStep==0))
            throw std::invalid_argument("fixed time mode requires --time-step and --steps; adaptive limits require --end-time");
        if (explicitVelocityRelaxation && timeStep==0)
            throw std::invalid_argument("velocity-relaxation option requires transient flow");
        if (vortexOptions) {
            if (vortexOptions!=15 || timeStep==0 || !restart.empty() || !boundaryExportPath.empty() ||
                (controls.scenario!="external" && controls.scenario!="channel" && controls.scenario!="duct" &&
                 controls.scenario!="custom" && controls.scenario!="cavity"))
                throw std::invalid_argument("initial vortex requires all four options and a fresh physical transient case, without restart/template export");
            fv::validateFlowInitialVortex2D(initialVortex);
        }

        const auto readStart = std::chrono::steady_clock::now();
        const auto read = readCm2dTopology(path);
        if (!read.valid()) {
            throw std::runtime_error(read.error);
        }
        const auto mesh = fv::makeFvMesh2D(read.topology);
        if (controls.scenario == "custom") {
            if (boundaryPath.empty()) throw std::invalid_argument("custom case requires --boundary FILE");
            std::ifstream input(boundaryPath);
            if (!input) throw std::runtime_error("cannot open explicit boundary file");
            controls.boundaryConditions = fv::readFlowBoundaryConditions2D(input, mesh, controls);
        } else if (!boundaryPath.empty()) {
            throw std::invalid_argument("--boundary requires --case custom");
        }
        if (!boundaryExportPath.empty()) {
            controls.boundaryConditions = fv::explicitFlowBoundaryPreset2D(mesh, controls);
            controls.scenario = "custom";
            const auto parent = std::filesystem::path(boundaryExportPath).parent_path();
            if (!parent.empty()) std::filesystem::create_directories(parent);
            std::ofstream output(boundaryExportPath);
            fv::writeFlowBoundaryConditions2D(output, mesh, controls);
            std::cout << "boundary_file=" << boundaryExportPath << '\n';
            return 0;
        }
        if(!viscosityPath.empty()) {
            if(controls.scenario=="manufactured" || controls.scenario=="taylor-green" || controls.scenario=="counterflow")
                throw std::invalid_argument("face-viscosity file requires a physical flow case");
            controls.faceViscosity=viscosityCsv(viscosityPath,mesh.faces.size());
        }
        if(controls.manufacturedViscositySlope!=0) {
            if(controls.scenario!="manufactured" || timeStep!=0 || !viscosityPath.empty() || controls.manufacturedViscositySlope<=-1)
                throw std::invalid_argument("manufactured-viscosity-slope requires steady manufactured case, slope > -1, without viscosity file");
            for(const auto& f:mesh.faces)
                controls.faceViscosity.push_back(controls.nu*(1+controls.manufacturedViscositySlope*f.centre.x));
        }
        const double readSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - readStart).count();
        const auto parent = std::filesystem::path(prefix).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
        if (controls.scenario == "custom") {
            auto boundarySnapshot = out(prefix, ".boundaries");
            fv::writeFlowBoundaryConditions2D(boundarySnapshot, mesh, controls);
        }
        fv::FlowResult2D r;
        std::size_t totalInnerIterations=0;
        fv::FlowPerformance2D totalPerformance;
        if (timeStep==0) r=fv::solveIncompressible2D(mesh,controls,progress);
        else {
            fv::FlowState2D state;
            if (restart.empty()) state=fv::initialIncompressibleState2D(mesh,controls);
            else {
                std::ifstream input(restart);
                if (!input) throw std::runtime_error("cannot open restart checkpoint");
                state=fv::readFlowCheckpoint2D(input,mesh,controls);
            }
            if (vortexOptions) state=fv::withInitialVortex2D(mesh,state,initialVortex);
            startTime=state.time;
            if (adaptive && !(adaptiveControls.targetTime>startTime))
                throw std::invalid_argument("end-time must be after the accepted restart time");
            acceptedTime=state.time;
            auto saveAccepted=[&]() {
                auto checkpoint=out(prefix,".checkpoint.tmp");
                fv::writeFlowCheckpoint2D(checkpoint,mesh,controls,state);
                checkpoint.close();
                std::filesystem::rename(prefix+".checkpoint.tmp",prefix+".checkpoint");
            };
            // Invalidate any old summary before replacing this prefix's files.
            // On interruption/exception, stale fields cannot look like a new pass.
            { auto pending=out(prefix,".json");
              pending << "{\"format\":\"cartmesh2d-flow-summary-v1\",\"status\":\"running\",\"converged\":false}\n"; }
            transientOutputStarted=true;
            if (vortexOptions) {
                auto initial=out(prefix,".initial.checkpoint");
                fv::writeFlowCheckpoint2D(initial,mesh,controls,state);
            } else std::filesystem::remove(prefix+".initial.checkpoint");
            saveAccepted(); // even a failed first step retains the valid initial/restart state
            auto times=out(prefix,".time-history.csv");
            times << "step,time,dt,accepted,innerIterations,momentumResidual,continuity,maxCourant,kineticEnergy,forceX,forceY\n";
            std::ofstream attempts;
            if (adaptive) {
                attempts=out(prefix,".attempt-history.csv");
                attempts << "attempt,step,startTime,time,dt,accepted,reason,innerConverged,innerIterations,momentumResidual,continuity,velocityChange,pressureChange,maxCourant\n";
            }
            for (std::size_t step=1;adaptive ? state.time<adaptiveControls.targetTime : step<=requestedSteps;++step) {
                if (adaptive && completedSteps>=adaptiveControls.maximumAcceptedSteps)
                    throw std::runtime_error("Adaptive accepted-step budget exhausted; last accepted checkpoint retained");
                double trialStep=adaptive?fv::nextAdaptiveFlowTimeStep2D(mesh,state,adaptiveControls):timeStep;
                std::size_t retry=0;
                for (;;) {
                    const double target=state.time+trialStep;
                    auto innerProgress=[&](const fv::FlowIteration2D& h) {
                        std::cout << std::setprecision(17) << "{\"type\":\"flow-progress\",\"time\":" << target
                                  << ",\"timeStep\":" << step << ",\"iteration\":" << h.iteration
                                  << ",\"momentumResidual\":" << h.momentumResidual << ",\"continuity\":" << h.continuity
                                  << ",\"velocityChange\":" << h.velocityChange << ",\"pressureChange\":" << h.pressureChange << "}" << std::endl;
                    };
                    r=fv::advanceIncompressible2D(mesh,controls,state,trialStep,innerProgress);
                    const auto& last=r.history.back();
                    totalInnerIterations+=last.iteration;
                    if (controls.profile) {
                        const auto& p=r.performance;
                        totalPerformance.momentumSolves+=p.momentumSolves;
                        totalPerformance.momentumIterations+=p.momentumIterations;
                        totalPerformance.maxMomentumIterations=std::max(totalPerformance.maxMomentumIterations,p.maxMomentumIterations);
                        totalPerformance.pressureSolves+=p.pressureSolves;
                        totalPerformance.pressureCorrectionPassesSkipped+=p.pressureCorrectionPassesSkipped;
                        totalPerformance.pressureFactorizations+=p.pressureFactorizations;
                        totalPerformance.pressureFactorReuses+=p.pressureFactorReuses;
                        totalPerformance.pressureHierarchyBuilds+=p.pressureHierarchyBuilds;
                        totalPerformance.pressureHierarchyReuses+=p.pressureHierarchyReuses;
                        totalPerformance.pressureHierarchyRefreshes+=p.pressureHierarchyRefreshes;
                        totalPerformance.maxPressureHierarchyLevels=std::max(totalPerformance.maxPressureHierarchyLevels,p.maxPressureHierarchyLevels);
                        totalPerformance.maxPressureCoarseCells=std::max(totalPerformance.maxPressureCoarseCells,p.maxPressureCoarseCells);
                        totalPerformance.pressureIterations+=p.pressureIterations;
                        totalPerformance.maxPressureIterations=std::max(totalPerformance.maxPressureIterations,p.maxPressureIterations);
                        totalPerformance.momentumLinearSolveSeconds+=p.momentumLinearSolveSeconds;
                        totalPerformance.pressureLinearSolveSeconds+=p.pressureLinearSolveSeconds;
                        totalPerformance.solveSeconds+=p.solveSeconds;
                    }
                    if (!adaptive) break;
                    ++attemptCount;
                    const bool accepted=r.converged && r.maxCourant<=adaptiveControls.maximumCourant;
                    const char* reason=accepted?"accepted":r.converged?"courant":"nonconverged";
                    if (accepted && trialStep==adaptiveControls.targetTime-state.time)
                        r.time=adaptiveControls.targetTime; // canonical timestamp of the computed final remainder
                    attempts << attemptCount << ',' << step << ',' << state.time << ',' << r.time << ',' << trialStep
                             << ',' << (accepted?1:0) << ',' << reason << ',' << (r.converged?1:0) << ',' << last.iteration
                             << ',' << last.momentumResidual << ',' << last.continuity << ',' << last.velocityChange
                             << ',' << last.pressureChange << ',' << r.maxCourant << '\n';
                    attempts.flush();
                    if (accepted) break;
                    ++rejectedSteps;
                    if (r.stopped) throw std::runtime_error("Adaptive calculation stopped; last accepted checkpoint retained");
                    const auto smaller=fv::reducedFlowTimeStep2D(trialStep,r.maxCourant,adaptiveControls.targetTime-state.time,adaptiveControls);
                    if (retry>=adaptiveControls.maximumRetries || !smaller)
                        throw std::runtime_error(std::string("Adaptive ")+reason+" rejection exhausted minimum step or retry budget; last accepted checkpoint retained");
                    std::cout << "{\"type\":\"flow-time-retry\",\"step\":" << step << ",\"attempt\":" << attemptCount
                              << ",\"acceptedTime\":" << state.time << ",\"candidateTime\":" << r.time
                              << ",\"dt\":" << trialStep << ",\"nextDt\":" << *smaller
                              << ",\"maxCourant\":" << r.maxCourant << ",\"reason\":\"" << reason << "\"}" << std::endl;
                    trialStep=*smaller;++retry; // state has not been changed by a rejected trial
                }
                const auto& last=r.history.back();
                double energy=0;
                for (std::size_t i=0;i<mesh.cells.size();++i)
                    energy+=.5*mesh.cells[i].area*(r.u[i]*r.u[i]+r.v[i]*r.v[i]);
                times << step << ',' << r.time << ',' << r.timeStep << ',' << (r.converged?1:0) << ',' << last.iteration
                      << ',' << last.momentumResidual << ',' << last.continuity << ',' << r.maxCourant
                      << ',' << energy << ',' << r.forceX << ',' << r.forceY << '\n';
                times.flush();
                if (!r.converged) break; // fixed-step mode preserves its diagnostic failure output
                state={r.time,r.u,r.v,r.p,r.flux};acceptedTime=r.time;++completedSteps;
                saveAccepted();
                std::cout << "{\"type\":\"flow-time-step\",\"time\":" << state.time
                          << ",\"step\":" << step << ",\"maxCourant\":" << r.maxCourant
                          << ",\"kineticEnergy\":" << energy << ",\"forceX\":" << r.forceX << ",\"forceY\":" << r.forceY << "}" << std::endl;
            }
        }
        if (timeStep>0 && controls.profile) r.performance=totalPerformance;
        const auto& last=r.history.back();

        auto cells = out(prefix, ".cells.csv");
        cells << "cell,x,y,area,u,v,p,speed";
        const bool counterflow = controls.scenario == "counterflow";
        if (controls.scenario == "manufactured" || counterflow) cells << ",sourceX,sourceY,exactU,exactV,exactP";
        if (timeStep>0) cells << ",previousU,previousV,temporalX,temporalY";
        cells << '\n';
        auto fields = out(prefix, ".fields.json");
        fields << "{\"format\":\"cartmesh2d-flow-v1\",\"cells\":[\n";
        for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
            const auto& c = mesh.cells[i];
            const double speed = std::hypot(r.u[i], r.v[i]);
            cells << i << ',' << c.centre.x << ',' << c.centre.y << ',' << c.area << ','
                  << r.u[i] << ',' << r.v[i] << ',' << r.p[i] << ',' << speed;
            if (controls.scenario == "manufactured") {
                const auto exact=fv::manufacturedFlow2D(c.centre,controls.speed,controls.nu,controls.manufacturedPressureSlope,controls.manufacturedViscositySlope,controls.viscousStress==fv::ViscousStress2D::Symmetric);
                const double gauge=fv::manufacturedFlow2D(mesh.cells.front().centre,controls.speed,controls.nu,controls.manufacturedPressureSlope,controls.manufacturedViscositySlope,controls.viscousStress==fv::ViscousStress2D::Symmetric).pressure;
                cells << ',' << r.sourceIntegrals[i].x << ',' << r.sourceIntegrals[i].y
                      << ',' << exact.velocity.x << ',' << exact.velocity.y << ',' << exact.pressure-gauge;
            } else if (counterflow) {
                const double pi = std::acos(-1.0);
                const double exactU = controls.speed * (1.0 + 2.0 * std::cos(2.0 * pi * c.centre.y));
                cells << ',' << r.sourceIntegrals[i].x << ',' << r.sourceIntegrals[i].y
                      << ',' << exactU << ',' << 0.0 << ',' << 0.0;
            }
            if (timeStep>0) cells << ',' << r.previousU[i] << ',' << r.previousV[i]
                << ',' << r.temporalIntegrals[i].x << ',' << r.temporalIntegrals[i].y;
            cells << '\n';
            if (i) {
                fields << ",\n";
            }
            fields << "{\"id\":" << i << ",\"u\":" << r.u[i] << ",\"v\":" << r.v[i]
                   << ",\"p\":" << r.p[i] << ",\"speed\":" << speed << '}';
        }
        fields << "\n]}\n";

        auto faces = out(prefix, ".faces.csv");
        faces << "face,owner,neighbour,flux,pressure,advectionX,advectionY,diffusionX,diffusionY,wall";
        if(!controls.faceViscosity.empty())faces<<",viscosity";
        faces<<'\n';
        for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
            const auto& f = mesh.faces[i];
            faces << i << ',' << f.owner << ',';
            if (f.neighbour) {
                faces << *f.neighbour;
            } else {
                faces << -1;
            }
            const auto& fm = r.faceMomentum[i];
            faces << ',' << r.flux[i] << ',' << fm.pressure << ',' << fm.advection.x
                  << ',' << fm.advection.y << ',' << fm.diffusion.x << ',' << fm.diffusion.y << ',' << (fm.wall?1:0);
            if(!controls.faceViscosity.empty())faces<<','<<controls.faceViscosity[i];
            faces<<'\n';
        }

        auto history = out(prefix, ".residuals.csv");
        history << "iteration,momentumResidual,continuity,velocityChange,pressureChange\n";
        for (auto h : r.history) {
            history << h.iteration << ',' << h.momentumResidual << ',' << h.continuity << ','
                    << h.velocityChange << ',' << h.pressureChange << '\n';
        }

        auto summary = out(prefix, ".json");
        const bool custom = controls.scenario == "custom";
        const bool customClosed = custom && std::none_of(controls.boundaryConditions.begin(), controls.boundaryConditions.end(),
            [](const auto& b) { return b.kind == fv::FlowBoundaryKind2D::PressureOutlet || b.kind == fv::FlowBoundaryKind2D::PressureOpening; });
        const char* preconditioner = controls.pressurePreconditioner ==
            fv::PressurePreconditioner2D::IncompleteCholesky0 ? "ic0" :
            (controls.pressurePreconditioner == fv::PressurePreconditioner2D::Aggregation ? "aggregation" : "jacobi");
        const bool symmetric=controls.viscousStress==fv::ViscousStress2D::Symmetric;
        const bool manufactured=controls.scenario == "manufactured";
        const bool counterflowCase=controls.scenario == "counterflow";
        const char* convection = fv::flow_checkpoint_detail::convectionName(controls.convection);
        summary << "{\n";
        if (vortexOptions) summary << "\"initialVortex\":{\"definition\":\"compact-cubic-v1\",\"centre\":["
            << initialVortex.centre.x << ',' << initialVortex.centre.y << "],\"radius\":" << initialVortex.radius
            << ",\"peakSpeed\":" << initialVortex.peakSpeed << ",\"checkpointSuffix\":\".initial.checkpoint\"},\n";
        if (custom) {
            summary << "\"boundaryFileSuffix\":\".boundaries\",\n\"referenceSpeedRole\":\"normalization-only\",\n\"boundaryConditions\":[";
            bool comma = false;
            for (const auto& b : controls.boundaryConditions) {
                if (comma) summary << ',';
                comma = true;
                summary << "{\"face\":" << b.face << ",\"type\":" << std::quoted(fv::flowBoundaryKindName2D(b.kind))
                        << ",\"name\":" << std::quoted(b.name) << ",\"u\":" << b.velocity.x
                        << ",\"v\":" << b.velocity.y << ",\"p\":" << b.pressure << '}';
            }
            summary << "],\n";
            summary << "\"wallLoadReference\":[0,0],\n\"wallLoadDefinition\":\"fluid-on-wall / density / depth; shared-face pressure and selected viscous flux; torque positive counterclockwise\",\n\"namedWallLoads\":[";
            comma = false;
            for (const auto& load : r.namedWallLoads) {
                if (comma) summary << ',';
                comma = true;
                summary << "{\"name\":" << std::quoted(load.name) << ",\"faces\":" << load.faces
                    << ",\"length\":" << load.length
                    << ",\"pressureForceX\":" << load.pressure.x << ",\"pressureForceY\":" << load.pressure.y
                    << ",\"viscousForceX\":" << load.viscous.x << ",\"viscousForceY\":" << load.viscous.y
                    << ",\"forceX\":" << load.pressure.x+load.viscous.x
                    << ",\"forceY\":" << load.pressure.y+load.viscous.y
                    << ",\"pressureTorque\":" << load.pressureTorque << ",\"viscousTorque\":" << load.viscousTorque
                    << ",\"torque\":" << load.pressureTorque+load.viscousTorque << '}';
            }
            summary << "],\n";
        }
        if (timeStep>0) {
            summary << "\"temporalDiscretization\":\"backward-euler\",\n"
                << "\"temporalFaceInterpolation\":\"old-and-iteration-flux-defect-skew-corrected-v2\",\n"
                << "\"time\":" << r.time << ",\n\"dt\":" << r.timeStep
                << ",\n\"acceptedTime\":" << acceptedTime;
            if (adaptive) summary << ",\n\"timeStepControl\":\"adaptive-cfl-retry\",\n\"startTime\":" << startTime
                << ",\n\"targetTime\":" << adaptiveControls.targetTime << ",\n\"maximumTimeStep\":" << timeStep
                << ",\n\"minimumTimeStep\":" << adaptiveControls.minimumStep << ",\n\"targetCourant\":" << adaptiveControls.maximumCourant
                << ",\n\"maximumRetries\":" << adaptiveControls.maximumRetries << ",\n\"maximumAcceptedSteps\":" << adaptiveControls.maximumAcceptedSteps
                << ",\n\"attemptCount\":" << attemptCount << ",\n\"rejectedSteps\":" << rejectedSteps;
            else summary << ",\n\"requestedSteps\":" << requestedSteps;
            summary << ",\n\"completedSteps\":" << completedSteps << ",\n\"maxCourant\":" << r.maxCourant
                << ",\n\"velocityRelaxation\":" << controls.velocityRelaxation << ",\n";
        }
        if (manufactured && controls.manufacturedViscositySlope==0) summary << "\"manufacturedDefinition\":\"psi=(speed/pi)*sin(pi*x)^2*sin(pi*y)^2; p=speed^2*(cos(pi*x)*cos(pi*y)+slope*(x+y)); source=advection+grad(p)-nu*laplacian(U); centroid quadrature\",\n"
                                  << "\"manufacturedPressureSlope\":" << controls.manufacturedPressureSlope << ",\n";
        if(manufactured && controls.manufacturedViscositySlope!=0) summary<<"\"manufacturedDefinition\":\"psi=(speed/pi)*sin(pi*x)^2*sin(pi*y)^2; p=speed^2*(cos(pi*x)*cos(pi*y)+pressureSlope*(x+y)); nu(x)=nu*(1+viscositySlope*x); source=advection+grad(p)-div(nu*stressGradient); stressGradient=grad(U)+grad(U)^T for symmetric, grad(U) for laplacian; centroid quadrature\",\n"
            <<"\"manufacturedPressureSlope\":"<<controls.manufacturedPressureSlope<<",\n";
        if(!controls.faceViscosity.empty()) summary<<"\"viscosityModel\":\"face-values\",\n\"viscosityFile\":"<<std::quoted(viscosityPath)
            <<",\n\"manufacturedViscositySlope\":"<<controls.manufacturedViscositySlope<<",\n";
        if (counterflowCase) summary << "\"counterflowDefinition\":\"u=speed*(1+2*cos(2*pi*y)), v=0, p=0; sourceX=8*pi^2*nu*speed*cos(2*pi*y), sourceY=0\",\n";
        summary << "\"format\":\"cartmesh2d-flow-summary-v1\",\n\"case\":\""
                << controls.scenario << "\",\n\"status\":\""
                << (r.converged ? "converged" : (timeStep>0?"time_step_not_converged":"iteration_limit"))
                << "\",\n\"converged\":" << (r.converged ? "true" : "false")
                << ",\n\"cells\":" << mesh.cells.size()
                << ",\n\"iterations\":" << last.iteration
                << ",\n\"nu\":" << controls.nu
                << ",\n\"speed\":" << controls.speed
                << ",\n\"continuity\":" << last.continuity
                << ",\n\"globalImbalance\":" << r.globalImbalance
                << ",\n\"momentumResidual\":" << last.momentumResidual
                << ",\n\"velocityChange\":" << last.velocityChange
                << ",\n\"pressureChange\":" << last.pressureChange
                << ",\n\"forceX\":" << r.forceX
                << ",\n\"forceY\":" << r.forceY
                << ",\n\"pressureForceX\":" << r.pressureForceX
                << ",\n\"pressureForceY\":" << r.pressureForceY
                << ",\n\"discreteForceX\":" << r.discreteForceX
                << ",\n\"discreteForceY\":" << r.discreteForceY
                << ",\n\"pressureDiscretization\":\"shared-face-gauss\""
                << ",\n\"pressureBoundaryReconstruction\":\"one-sided-linear-adaptive\""
                << ",\n\"convection\":\"" << convection << '"'
                << ",\n\"viscousStress\":\"" << (symmetric?"symmetric":"laplacian") << '\"'
                << ",\n\"outletBackflow\":\"" << (controls.outletBackflow == fv::OutletBackflow2D::NormalInlet ? "normal-inlet" : "reject") << '\"'
                << ",\n\"outletBackflowFaces\":" << r.outletBackflowFaces
                << ",\n\"outletInflow\":" << r.outletInflow
                << ",\n\"forceDefinition\":\"" << (symmetric?"shared-face-newtonian-traction":"reconstructed-newtonian-traction") << '\"'
                << ",\n\"reconstructedForceX\":" << r.reconstructedForceX
                << ",\n\"reconstructedForceY\":" << r.reconstructedForceY
                << ",\n\"wallForceX\":" << r.wallForceX
                << ",\n\"wallForceY\":" << r.wallForceY
                << ",\n\"wallViscousForceX\":" << r.wallViscousForceX
                << ",\n\"wallViscousForceY\":" << r.wallViscousForceY
                << ",\n\"wallForceDefinition\":\"all no-slip walls and moving lid; fluid on boundary; same pressure and viscous flux as momentum\""
                << ",\n\"discreteForceDefinition\":\"pressure plus selected viscous momentum flux; embedded walls\""
                << ",\n\"globalRelativeImbalance\":" << r.globalRelativeImbalance
                << ",\n\"domainHeight\":" << r.domainHeight
                << ",\n\"tolerance\":" << controls.tolerance
                << ",\n\"pressurePreconditioner\":\"" << preconditioner << '"'
                << ",\n\"units\":{\"velocity\":\"m/s\",\"p\":\"m2/s2 (kinematic)\",\"nu\":\"m2/s\",\"faceFlux\":\"m2/s per unit depth\",\"force\":\"m3/s2 (force / density / depth), fluid on stationary embedded walls, positive Cartesian axes\"},\n"
                << "\"pressureReference\":\""
                << ((controls.scenario == "cavity" || controls.scenario=="taylor-green" || manufactured || customClosed)
                        ? "cell 0, kinematic pressure zero"
                       : custom ? (std::any_of(controls.boundaryConditions.begin(), controls.boundaryConditions.end(),
                             [](const auto& b){return b.kind==fv::FlowBoundaryKind2D::PressureOpening;})
                             ? "explicit pressure opening faces, prescribed static kinematic pressure"
                             : "explicit pressure outlet faces, prescribed kinematic pressure")
                                 : "right outlet faces, kinematic pressure zero")
                << "\",\n"
                << "\"method\":\"cell-centred FVM; SIMPLE; Rhie-Chow; shared-face pressure; "
                << convection << " momentum convection; corrected diffusion\",\n"
                << "\"scope\":\"" << (timeStep>0 ? "transient backward-Euler constant-property laminar flow; no turbulence or accuracy certification"
                                                  : "steady constant-property laminar flow; no turbulence or accuracy certification") << "\"\n}\n";

        std::string error;
        if (!writeLegacyVtk2D(read.topology, prefix + ".vtk", &error)) {
            throw std::runtime_error(error);
        }
        std::ofstream vtk(prefix + ".vtk", std::ios::app);
        vtk.exceptions(std::ios::failbit | std::ios::badbit);
        vtk << std::setprecision(17) << "VECTORS velocity double\n";
        for (std::size_t i = 0; i < r.u.size(); ++i) {
            vtk << r.u[i] << ' ' << r.v[i] << " 0\n";
        }
        vtk << "SCALARS pressure_kinematic double 1\nLOOKUP_TABLE default\n";
        for (double p : r.p) {
            vtk << p << '\n';
        }
        vtk << "SCALARS speed double 1\nLOOKUP_TABLE default\n";
        for (std::size_t i = 0; i < r.u.size(); ++i) {
            vtk << std::hypot(r.u[i], r.v[i]) << '\n';
        }

        cells.close();
        faces.close();
        fields.close();
        history.close();
        summary.close();
        vtk.close();
        if (controls.profile) {
            const auto& p = r.performance;
            auto performance = out(prefix, ".performance.json");
            performance << "{\n\"format\":\"cartmesh2d-flow-performance-v1\",\n"
                        << "\"cells\":" << mesh.cells.size()
                        << ",\n\"faces\":" << mesh.faces.size()
                        << ",\n\"simpleIterations\":" << (timeStep>0?totalInnerIterations:last.iteration)
                        << ",\n\"converged\":" << (r.converged ? "true" : "false")
                        << ",\n\"pressurePreconditioner\":\"" << preconditioner << '"'
                        << ",\n\"readAndMeshSeconds\":" << readSeconds
                        << ",\n\"solveSeconds\":" << p.solveSeconds
                        << ",\n\"momentumLinearSolveSeconds\":" << p.momentumLinearSolveSeconds
                        << ",\n\"pressureLinearSolveSeconds\":" << p.pressureLinearSolveSeconds
                        << ",\n\"momentumSolves\":" << p.momentumSolves
                        << ",\n\"momentumIterations\":" << p.momentumIterations
                        << ",\n\"maxMomentumIterations\":" << p.maxMomentumIterations
                        << ",\n\"pressureSolves\":" << p.pressureSolves
                        << ",\n\"pressureCorrectionPassesSkipped\":" << p.pressureCorrectionPassesSkipped
                        << ",\n\"pressureFactorizations\":" << p.pressureFactorizations
                        << ",\n\"pressureFactorReuses\":" << p.pressureFactorReuses
                        << ",\n\"pressureHierarchyBuilds\":" << p.pressureHierarchyBuilds
                        << ",\n\"pressureHierarchyReuses\":" << p.pressureHierarchyReuses
                        << ",\n\"pressureHierarchyRefreshes\":" << p.pressureHierarchyRefreshes
                        << ",\n\"maxPressureHierarchyLevels\":" << p.maxPressureHierarchyLevels
                        << ",\n\"maxPressureCoarseCells\":" << p.maxPressureCoarseCells
                        << ",\n\"pressureIterations\":" << p.pressureIterations
                        << ",\n\"maxPressureIterations\":" << p.maxPressureIterations
                        << ",\n\"scope\":\"steady-clock wall seconds; solve includes validation, assembly, monitoring and callbacks; transient sums all inner solves; linear times include linear setup, exclude assembly; exports/checkpoints excluded; no memory measurement\"\n}\n";
            performance.close();
        }
        if (timeStep==0) progress(last);
        return r.converged ? 0 : 2;
    } catch (const std::exception& e) {
        if (transientOutputStarted) {
            try { auto failed=out(prefix,".json");
                failed << "{\"format\":\"cartmesh2d-flow-summary-v1\",\"status\":\"failed\",\"converged\":false,\"acceptedTime\":" << acceptedTime << "}\n";
            } catch (const std::exception&) { /* preserve the original failure */ }
        }
        std::cerr << "cartmesh2d_flow_cli: " << e.what() << '\n';
        return 1;
    }
}
