#include "cartmesh2d/fv/ManufacturedFlow2D.hpp"
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/fv/FlowCheckpoint2D.hpp"
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
        std::string path;
        fv::FlowControls2D controls;
        bool explicitVelocityRelaxation=false;
        double timeStep=0;
        std::size_t requestedSteps=0,completedSteps=0;
        std::string restart;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--profile") {
                controls.profile = true;
                continue;
            }
            if (a == "--help") {
                std::cout
                    << "Native 2D incompressible laminar SIMPLE (experimental)\n"
            "--mesh FINAL.solver.cm2d --output PREFIX --case external|channel|cavity|manufactured|counterflow\n"
            "--nu 0.01 --speed 1 --max-iterations 1500 --tolerance 1e-6\n"
            "--time-step DT --steps N: backward Euler physical time, converged SIMPLE at each step.\n"
            "--velocity-relaxation 0.6: transient inner iterations only; (0,1], larger may be unstable.\n"
            "--restart PREFIX.checkpoint: resume accepted state on identical mesh and physical setup.\n"
            "--case taylor-green: unforced exact slip-box decay; transient verification only.\n"
            "Transient physical cases start at rest; boundary velocities switch on for t>0.\n"
            "Transient retains inner relaxation flux correction; fixed DT; reports CFL without changing DT.\n"
            "--profile writes extra .performance.json timing/linear iteration diagnostics.\n"
            "--pressure-preconditioner ic0|jacobi|aggregation (default ic0); aggregation experimental; same true-residual tolerance.\n"
            "--viscous-stress symmetric|laplacian (default symmetric); conservative Newtonian stress.\n"
            "--convection upwind|limited-linear (default upwind); bounded face reconstruction.\n"
            "--outlet-backflow reject|normal-inlet (default reject).\n"
            "manufactured: unit-square analytic forced vortex; verification only, stationary walls.\n"
            "--manufactured-pressure-slope 0: add Uref^2*slope*(x+y) to the analytic pressure.\n"
            "channel speed=maximum parabolic inlet speed; cavity speed=lid speed.\n"
            "Only fixed axis-aligned rectangular outer boundaries. Pressure is kinematic.\n"
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
            } else if (a == "--nu") {
                controls.nu = number(v);
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
            } else if (a == "--restart") {
                restart=v;
            } else if (a == "--manufactured-pressure-slope") {
                controls.manufacturedPressureSlope = number(v);
            } else if (a == "--viscous-stress") {
                if (v != "symmetric" && v != "laplacian") throw std::invalid_argument("viscous-stress must be symmetric or laplacian");
                controls.viscousStress=v=="symmetric"?fv::ViscousStress2D::Symmetric:fv::ViscousStress2D::Laplacian;
            } else if (a == "--convection") {
                if (v != "upwind" && v != "limited-linear")
                    throw std::invalid_argument("convection must be upwind or limited-linear");
                controls.convection = v == "limited-linear"
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
        if (path.empty() || prefix.empty()) {
            throw std::invalid_argument("--mesh and --output required");
        }
        if (!path.ends_with(".solver.cm2d") || path.ends_with(".failed.solver.cm2d")) {
            throw std::invalid_argument("requires final *.solver.cm2d");
        }

        if ((timeStep>0)!=(requestedSteps>0) || (!restart.empty() && timeStep==0))
            throw std::invalid_argument("--time-step and --steps must be provided together; restart requires them");
        if (explicitVelocityRelaxation && timeStep==0)
            throw std::invalid_argument("velocity-relaxation option requires transient flow");

        const auto readStart = std::chrono::steady_clock::now();
        const auto read = readCm2dTopology(path);
        if (!read.valid()) {
            throw std::runtime_error(read.error);
        }
        const auto mesh = fv::makeFvMesh2D(read.topology);
        const double readSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - readStart).count();
        const auto parent = std::filesystem::path(prefix).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);
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
            saveAccepted(); // even a failed first step retains the valid initial/restart state
            auto times=out(prefix,".time-history.csv");
            times << "step,time,dt,accepted,innerIterations,momentumResidual,continuity,maxCourant,kineticEnergy,forceX,forceY\n";
            for (std::size_t step=1;step<=requestedSteps;++step) {
                const double target=state.time+timeStep;
                auto innerProgress=[&](const fv::FlowIteration2D& h) {
                    std::cout << std::setprecision(17) << "{\"type\":\"flow-progress\",\"time\":" << target
                              << ",\"timeStep\":" << step << ",\"iteration\":" << h.iteration
                              << ",\"momentumResidual\":" << h.momentumResidual << ",\"continuity\":" << h.continuity
                              << ",\"velocityChange\":" << h.velocityChange << ",\"pressureChange\":" << h.pressureChange << "}" << std::endl;
                };
                r=fv::advanceIncompressible2D(mesh,controls,state,timeStep,innerProgress);
                const auto& last=r.history.back();
                totalInnerIterations+=last.iteration;
                if (controls.profile) {
                    const auto& p=r.performance;
                    totalPerformance.momentumSolves+=p.momentumSolves;
                    totalPerformance.momentumIterations+=p.momentumIterations;
                    totalPerformance.maxMomentumIterations=std::max(totalPerformance.maxMomentumIterations,p.maxMomentumIterations);
                    totalPerformance.pressureSolves+=p.pressureSolves;
                    totalPerformance.pressureFactorizations+=p.pressureFactorizations;
                    totalPerformance.pressureFactorReuses+=p.pressureFactorReuses;
                    totalPerformance.pressureHierarchyBuilds+=p.pressureHierarchyBuilds;
                    totalPerformance.pressureHierarchyReuses+=p.pressureHierarchyReuses;
                    totalPerformance.maxPressureHierarchyLevels=std::max(totalPerformance.maxPressureHierarchyLevels,p.maxPressureHierarchyLevels);
                    totalPerformance.maxPressureCoarseCells=std::max(totalPerformance.maxPressureCoarseCells,p.maxPressureCoarseCells);
                    totalPerformance.pressureIterations+=p.pressureIterations;
                    totalPerformance.maxPressureIterations=std::max(totalPerformance.maxPressureIterations,p.maxPressureIterations);
                    totalPerformance.momentumLinearSolveSeconds+=p.momentumLinearSolveSeconds;
                    totalPerformance.pressureLinearSolveSeconds+=p.pressureLinearSolveSeconds;
                    totalPerformance.solveSeconds+=p.solveSeconds;
                }
                double energy=0;
                for (std::size_t i=0;i<mesh.cells.size();++i)
                    energy+=.5*mesh.cells[i].area*(r.u[i]*r.u[i]+r.v[i]*r.v[i]);
                times << step << ',' << r.time << ',' << timeStep << ',' << (r.converged?1:0) << ',' << last.iteration
                      << ',' << last.momentumResidual << ',' << last.continuity << ',' << r.maxCourant
                      << ',' << energy << ',' << r.forceX << ',' << r.forceY << '\n';
                times.flush();
                if (!r.converged) break; // never advance the physical time with an unconverged candidate
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
                const auto exact=fv::manufacturedFlow2D(c.centre,controls.speed,controls.nu,controls.manufacturedPressureSlope);
                const double gauge=fv::manufacturedFlow2D(mesh.cells.front().centre,controls.speed,controls.nu,controls.manufacturedPressureSlope).pressure;
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
        faces << "face,owner,neighbour,flux,pressure,advectionX,advectionY,diffusionX,diffusionY,wall\n";
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
                  << ',' << fm.advection.y << ',' << fm.diffusion.x << ',' << fm.diffusion.y << ',' << (fm.wall?1:0) << '\n';
        }

        auto history = out(prefix, ".residuals.csv");
        history << "iteration,momentumResidual,continuity,velocityChange,pressureChange\n";
        for (auto h : r.history) {
            history << h.iteration << ',' << h.momentumResidual << ',' << h.continuity << ','
                    << h.velocityChange << ',' << h.pressureChange << '\n';
        }

        auto summary = out(prefix, ".json");
        const char* preconditioner = controls.pressurePreconditioner ==
            fv::PressurePreconditioner2D::IncompleteCholesky0 ? "ic0" :
            (controls.pressurePreconditioner == fv::PressurePreconditioner2D::Aggregation ? "aggregation" : "jacobi");
        const bool symmetric=controls.viscousStress==fv::ViscousStress2D::Symmetric;
        const bool manufactured=controls.scenario == "manufactured";
        const bool counterflowCase=controls.scenario == "counterflow";
        const char* convection = controls.convection == fv::ConvectionScheme2D::LimitedLinearUpwind
            ? "limited-linear" : "upwind";
        summary << "{\n";
        if (timeStep>0) summary << "\"temporalDiscretization\":\"backward-euler\",\n"
            << "\"temporalFaceInterpolation\":\"old-and-iteration-flux-defect-skew-corrected-v2\",\n"
            << "\"time\":" << r.time << ",\n\"dt\":" << timeStep
            << ",\n\"acceptedTime\":" << acceptedTime << ",\n\"requestedSteps\":" << requestedSteps
            << ",\n\"completedSteps\":" << completedSteps << ",\n\"maxCourant\":" << r.maxCourant
            << ",\n\"velocityRelaxation\":" << controls.velocityRelaxation << ",\n";
        if (manufactured) summary << "\"manufacturedDefinition\":\"psi=(speed/pi)*sin(pi*x)^2*sin(pi*y)^2; p=speed^2*(cos(pi*x)*cos(pi*y)+slope*(x+y)); source=advection+grad(p)-nu*laplacian(U); centroid quadrature\",\n"
                                  << "\"manufacturedPressureSlope\":" << controls.manufacturedPressureSlope << ",\n";
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
                << ",\n\"pressureBoundaryReconstruction\":\"one-sided-linear-2ring\""
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
                << ((controls.scenario == "cavity" || controls.scenario=="taylor-green" || manufactured)
                        ? "cell 0, kinematic pressure zero"
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
                        << ",\n\"pressureFactorizations\":" << p.pressureFactorizations
                        << ",\n\"pressureFactorReuses\":" << p.pressureFactorReuses
                        << ",\n\"pressureHierarchyBuilds\":" << p.pressureHierarchyBuilds
                        << ",\n\"pressureHierarchyReuses\":" << p.pressureHierarchyReuses
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
