#include "cartmesh2d/fv/ManufacturedFlow2D.hpp"
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
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
    try {
        std::string path;
        std::string prefix;
        fv::FlowControls2D controls;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--profile") {
                controls.profile = true;
                continue;
            }
            if (a == "--help") {
                std::cout
                    << "Native 2D steady incompressible laminar SIMPLE (experimental)\n"
            "--mesh FINAL.solver.cm2d --output PREFIX --case external|channel|cavity|manufactured\n"
            "--nu 0.01 --speed 1 --max-iterations 1500 --tolerance 1e-6\n"
            "--profile writes extra .performance.json timing/linear iteration diagnostics.\n"
            "--pressure-preconditioner ic0|jacobi (default ic0); same true-residual tolerance.\n"
            "--convection upwind|limited-linear (default upwind); bounded face reconstruction.\n"
            "manufactured: unit-square analytic forced vortex; verification only, stationary walls.\n"
            "--manufactured-pressure-slope 0: add Uref^2*slope*(x+y) to the analytic pressure.\n"
            "channel speed=maximum parabolic inlet speed; cavity speed=lid speed.\n"
            "Only fixed axis-aligned rectangular outer boundaries. Pressure is kinematic.\n"
            "No turbulence/compressibility; outlet backflow explicitly unsupported.\n";
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
            } else if (a == "--manufactured-pressure-slope") {
                controls.manufacturedPressureSlope = number(v);
            } else if (a == "--convection") {
                if (v != "upwind" && v != "limited-linear")
                    throw std::invalid_argument("convection must be upwind or limited-linear");
                controls.convection = v == "limited-linear"
                    ? fv::ConvectionScheme2D::LimitedLinearUpwind : fv::ConvectionScheme2D::Upwind;
            } else if (a == "--pressure-preconditioner") {
                if (v != "ic0" && v != "jacobi") {
                    throw std::invalid_argument("pressure preconditioner must be ic0 or jacobi");
                }
                controls.pressurePreconditioner = v == "ic0"
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

        const auto readStart = std::chrono::steady_clock::now();
        const auto read = readCm2dTopology(path);
        if (!read.valid()) {
            throw std::runtime_error(read.error);
        }
        const auto mesh = fv::makeFvMesh2D(read.topology);
        const double readSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - readStart).count();
        const auto r = fv::solveIncompressible2D(mesh, controls, progress);
        const auto& last = r.history.back();

        const auto parent = std::filesystem::path(prefix).parent_path();
        if (!parent.empty()) {
            std::filesystem::create_directories(parent);
        }

        auto cells = out(prefix, ".cells.csv");
        cells << "cell,x,y,area,u,v,p,speed";
        if (controls.scenario == "manufactured") cells << ",sourceX,sourceY,exactU,exactV,exactP";
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
            }
            cells << '\n';
            if (i) {
                fields << ",\n";
            }
            fields << "{\"id\":" << i << ",\"u\":" << r.u[i] << ",\"v\":" << r.v[i]
                   << ",\"p\":" << r.p[i] << ",\"speed\":" << speed << '}';
        }
        fields << "\n]}\n";

        auto faces = out(prefix, ".faces.csv");
        faces << "face,owner,neighbour,flux,pressure,advectionX,advectionY,diffusionX,diffusionY\n";
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
                  << ',' << fm.advection.y << ',' << fm.diffusion.x << ',' << fm.diffusion.y << '\n';
        }

        auto history = out(prefix, ".residuals.csv");
        history << "iteration,momentumResidual,continuity,velocityChange,pressureChange\n";
        for (auto h : r.history) {
            history << h.iteration << ',' << h.momentumResidual << ',' << h.continuity << ','
                    << h.velocityChange << ',' << h.pressureChange << '\n';
        }

        auto summary = out(prefix, ".json");
        const char* preconditioner = controls.pressurePreconditioner ==
            fv::PressurePreconditioner2D::IncompleteCholesky0 ? "ic0" : "jacobi";
        const bool manufactured=controls.scenario == "manufactured";
        const char* convection = controls.convection == fv::ConvectionScheme2D::LimitedLinearUpwind
            ? "limited-linear" : "upwind";
        summary << "{\n";
        if (manufactured) summary << "\"manufacturedDefinition\":\"psi=(speed/pi)*sin(pi*x)^2*sin(pi*y)^2; p=speed^2*(cos(pi*x)*cos(pi*y)+slope*(x+y)); source=advection+grad(p)-nu*laplacian(U); centroid quadrature\",\n"
                                  << "\"manufacturedPressureSlope\":" << controls.manufacturedPressureSlope << ",\n";
        summary << "\"format\":\"cartmesh2d-flow-summary-v1\",\n\"case\":\""
                << controls.scenario << "\",\n\"status\":\""
                << (r.converged ? "converged" : "iteration_limit")
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
                << ",\n\"pressureBoundaryReconstruction\":\"one-sided-linear\""
                << ",\n\"convection\":\"" << convection << '"'
                << ",\n\"forceDefinition\":\"reconstructed-newtonian-traction\""
                << ",\n\"discreteForceDefinition\":\"pressure plus negative nu grad(U) dot S; embedded walls; Laplacian momentum flux\""
                << ",\n\"globalRelativeImbalance\":" << r.globalRelativeImbalance
                << ",\n\"domainHeight\":" << r.domainHeight
                << ",\n\"tolerance\":" << controls.tolerance
                << ",\n\"pressurePreconditioner\":\"" << preconditioner << '"'
                << ",\n\"units\":{\"velocity\":\"m/s\",\"p\":\"m2/s2 (kinematic)\",\"nu\":\"m2/s\",\"faceFlux\":\"m2/s per unit depth\",\"force\":\"m3/s2 (force / density / depth), fluid on stationary embedded walls, positive Cartesian axes\"},\n"
                << "\"pressureReference\":\""
                << ((controls.scenario == "cavity" || manufactured)
                        ? "cell 0, kinematic pressure zero"
                        : "right outlet faces, kinematic pressure zero")
                << "\",\n"
                << "\"method\":\"cell-centred FVM; SIMPLE; Rhie-Chow; shared-face pressure; "
                << convection << " momentum convection; corrected diffusion\",\n"
                << "\"scope\":\"steady constant-property laminar flow; no turbulence or accuracy certification\"\n}\n";

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
                        << ",\n\"simpleIterations\":" << last.iteration
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
                        << ",\n\"pressureIterations\":" << p.pressureIterations
                        << ",\n\"maxPressureIterations\":" << p.maxPressureIterations
                        << ",\n\"scope\":\"steady-clock wall seconds; solve includes validation, assembly, monitoring and callbacks; linear times include linear setup, exclude assembly; export excluded; no memory measurement\"\n}\n";
            performance.close();
        }
        progress(last);
        return r.converged ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "cartmesh2d_flow_cli: " << e.what() << '\n';
        return 1;
    }
}
