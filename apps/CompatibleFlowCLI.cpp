#include "CompatibleFlowCLI.hpp"
#include "cartmesh2d/fv/CompatibleFlowBoundary2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>

namespace {
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using Clock = std::chrono::steady_clock;
volatile std::sig_atomic_t cancelled = 0;
void cancel(int) { cancelled = 1; }
struct Signals {
    using Handler = void (*)(int);
    Handler interrupt, terminate;
    Signals() : interrupt(std::signal(SIGINT, cancel)), terminate(std::signal(SIGTERM, cancel)) { cancelled = 0; }
    ~Signals() { std::signal(SIGINT, interrupt); std::signal(SIGTERM, terminate); }
};
void jsonString(std::ostream& out, const std::string& text) {
    constexpr char hex[] = "0123456789abcdef";
    out << '"';
    for (const char character : text) {
        const auto c = static_cast<unsigned char>(character);
        if (c == '"' || c == '\\') out << '\\' << static_cast<char>(c);
        else if (c < 32) out << "\\u00" << hex[c >> 4] << hex[c & 15];
        else out << static_cast<char>(c);
    }
    out << '"';
}
void jsonNumber(std::ostream& out, double value) { if (std::isfinite(value)) out << value; else out << "null"; }
std::ofstream output(const std::string& path) {
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out.open(path);
    out << std::setprecision(17);
    return out;
}
double number(const std::string& value) {
    std::size_t used = 0;
    const double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result)) throw std::invalid_argument("Expected finite numeric option");
    return result;
}
std::size_t count(const std::string& value) {
    if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
        throw std::invalid_argument("Expected positive integer option");
    const auto n = std::stoull(value);
    if (!n || n > std::numeric_limits<std::size_t>::max()) throw std::invalid_argument("Integer option out of range");
    return static_cast<std::size_t>(n);
}
const char* stopName(CompatibleFlowStop2D stop) {
    switch (stop) {
    case CompatibleFlowStop2D::Converged: return "converged";
    case CompatibleFlowStop2D::Cancelled: return "cancelled";
    case CompatibleFlowStop2D::NonlinearBudget: return "nonlinear-budget";
    case CompatibleFlowStop2D::LinearBudget: return "linear-budget";
    case CompatibleFlowStop2D::BacktrackingBudget: return "backtracking-budget";
    case CompatibleFlowStop2D::NumericalFailure: return "numerical-failure";
    case CompatibleFlowStop2D::BoundaryFailure: return "boundary-failure";
    }
    return "invalid";
}
void metrics(std::ostream& out, const std::optional<CompatibleFlowMetrics2D>& m) {
    if (!m) { out << "null"; return; }
    out << "{\"cellMomentum\":"; jsonNumber(out,m->cellMomentum);
    out << ",\"faceMomentum\":"; jsonNumber(out,m->faceMomentum);
    out << ",\"divergence\":"; jsonNumber(out,m->divergence);
    out << ",\"stateChange\":"; jsonNumber(out,m->stateChange);
    out << ",\"residualNorm\":"; jsonNumber(out,m->residualNorm);
    out << '}';
}
void progress(const CompatibleFlowIteration2D& it) {
    std::cout << std::setprecision(17) << "{\"type\":\"compatible-flow-progress\",\"iteration\":" << it.iteration
              << ",\"accepted\":" << (it.accepted ? "true" : "false") << ",\"metrics\":";
    metrics(std::cout,it.metrics);
    std::cout << "}" << std::endl;
}
struct Options {
    std::string mesh, prefix, scenario = "external", boundary, restart;
    CompatibleFlowControls2D controls;
};
Options parse(int argc, char** argv) {
    Options opt;
    std::set<std::string> seen;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (!seen.insert(key).second) throw std::invalid_argument("Duplicate option: " + key);
        if (i + 1 == argc) throw std::invalid_argument("Missing value for " + key);
        const std::string value = argv[++i];
        auto& c = opt.controls;
        if (key == "--discretization") { if (value != "compatible") throw std::invalid_argument("Invalid compatible selector"); }
        else if (key == "--mesh") opt.mesh = value;
        else if (key == "--output") opt.prefix = value;
        else if (key == "--case") opt.scenario = value;
        else if (key == "--boundary") opt.boundary = value;
        else if (key == "--restart") opt.restart = value;
        else if (key == "--nu") c.viscosity = number(value);
        else if (key == "--speed") c.referenceVelocity = number(value);
        else if (key == "--reference-length") c.referenceLength = number(value);
        else if (key == "--tolerance") c.equationTolerance = number(value);
        else if (key == "--state-tolerance") c.stateTolerance = number(value);
        else if (key == "--max-iterations") c.maximumIterations = count(value);
        else if (key == "--linear-tolerance") c.linearTolerance = number(value);
        else if (key == "--linear-restarts") c.maximumLinearRestarts = count(value);
        else if (key == "--krylov-directions") c.krylovDirections = count(value);
        else if (key == "--quadrature-order") { const auto n = count(value); if (n > 12) throw std::invalid_argument("Quadrature order must be 4..12"); c.quadratureOrder = static_cast<int>(n); }
        else if (key == "--pseudo-step") {
            const double step = number(value);
            if (step < 0) throw std::invalid_argument("Pseudo-step must be nonnegative");
            c.globalization = step == 0 ? CompatibleGlobalization2D::Backtracking : CompatibleGlobalization2D::PseudoTime;
            if (step > 0) c.initialPseudoStep = step;
        }
        else if (key == "--pseudo-maximum-step") c.maximumPseudoStep = number(value);
        else if (key == "--momentum-inertia") {
            if (value != "0" && value != "1") throw std::invalid_argument("Momentum inertia must be 0 or 1");
            c.equation = value == "0" ? CompatibleEquation2D::Stokes : CompatibleEquation2D::NavierStokes;
        }
        else if (key == "--compatible-pressure-inverse") {
            if (value != "mass" && value != "schur") throw std::invalid_argument("Pressure inverse must be mass or schur");
            c.pressureInverse = value == "mass" ? CompatiblePressureInverse2D::ViscousMass : CompatiblePressureInverse2D::DiagonalSchur;
        }
        else throw std::invalid_argument("Unsupported compatible option: " + key + "; see --discretization compatible --help");
    }
    const auto& c = opt.controls;
    if (opt.mesh.empty() || opt.prefix.empty()) throw std::invalid_argument("Compatible solver requires --mesh and --output");
    if (!opt.mesh.ends_with(".solver.cm2d") || opt.mesh.ends_with(".failed.solver.cm2d"))
        throw std::invalid_argument("Expected a final .solver.cm2d mesh, not a failed mesh");
    if (opt.scenario != "external" && opt.scenario != "cavity" && opt.scenario != "channel" && opt.scenario != "duct" && opt.scenario != "custom")
        throw std::invalid_argument("Compatible CLI supports external, cavity, channel, duct and custom");
    if ((opt.scenario == "custom") != !opt.boundary.empty()) throw std::invalid_argument("--boundary is required exactly for --case custom");
    if (c.viscosity <= 0 || c.referenceVelocity <= 0 || c.referenceLength <= 0 || c.equationTolerance <= 0 || c.stateTolerance <= 0 || c.linearTolerance <= 0 || c.linearTolerance > .01 || c.quadratureOrder < 4 || c.maximumPseudoStep < c.initialPseudoStep)
        throw std::invalid_argument("Invalid positive scale, tolerance, quadrature, Krylov size or pseudo-time limit");
    return opt;
}
void help() {
    std::cout << "Native compatible steady constant-property incompressible flow (development)\n"
        "--discretization compatible --mesh FINAL.solver.cm2d --output FRESH_PREFIX\n"
        "--case external|channel|duct|cavity|custom; custom requires --boundary FILE\n"
        "--nu .01 --speed 1 --reference-length 1 --momentum-inertia 1 (0 selects Stokes)\n"
        "--tolerance 1e-9 --state-tolerance 1e-9 --max-iterations 40\n"
        "--pseudo-step .1 (0 selects residual backtracking) --pseudo-maximum-step 1e6\n"
        "--compatible-pressure-inverse mass|schur --linear-tolerance 1e-13\n"
        "--linear-restarts 50 --krylov-directions 60 --quadrature-order 6\n"
        "Existing explicit boundary values are facewise constant. Symmetry preserves zero normal trace.\n"
        "Pressure outlet uses pseudo-traction -p*n and rejects backflow; it does not separately impose p and all normal velocity derivatives.\n"
        "--restart FILE: restore a matching accepted steady checkpoint and next pseudo-step; iteration budget is additional.\n"
        "No variable material, physical time, pressure opening, smooth moving wall or App integration yet.\n"
        "Outputs retain full P1 cell/face states. Seed and rejected trials are never accepted flow fields.\n"
        "Numerical convergence is not a spatial/physical accuracy certificate. SIGINT/SIGTERM cancel cooperatively.\n";
}
void ensureFresh(const std::string& prefix) {
    for (const auto* suffix : {".checkpoint", ".checkpoint.tmp", ".json", ".summary.json", ".summary.json.tmp", ".cells.csv", ".faces.csv", ".fields.json", ".residuals.csv", ".vtk", ".boundaries", ".seed.json", ".accepted.json", ".rejected.json"})
        if (std::filesystem::exists(prefix + suffix) || std::filesystem::is_symlink(prefix + suffix))
            throw std::invalid_argument("Output already exists: " + prefix + suffix);
    const auto parent = std::filesystem::path(prefix).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
}
void publishCheckpoint(const Options& opt, const CompatibleFlowCheckpoint2D& checkpoint) {
    auto out = output(opt.prefix + ".checkpoint.tmp");
    writeCompatibleFlowCheckpoint2D(out,checkpoint);
    out.close();
    // In the same directory: a process kill before rename leaves the previous
    // complete checkpoint available. Power-loss durability is not claimed.
    std::filesystem::rename(opt.prefix + ".checkpoint.tmp",opt.prefix + ".checkpoint");
}
void stateArrays(std::ostream& out, const CompatibleFlowState2D& state) {
    out << "\"cells\":[";
    for (std::size_t i = 0; i < state.cells.size(); ++i) {
        if (i) out << ',';
        out << '[';
        for (std::size_t j = 0; j < 9; ++j) { if (j) out << ','; out << state.cells[i][j]; }
        out << ']';
    }
    out << "],\"faces\":[";
    for (std::size_t i = 0; i < state.faces.size(); ++i) {
        if (i) out << ',';
        out << '[';
        for (std::size_t j = 0; j < 4; ++j) { if (j) out << ','; out << state.faces[i][j]; }
        out << ']';
    }
    out << ']';
}
void saveState(const Options& opt, const std::optional<CompatibleFlowState2D>& state, const char* suffix, const char* kind) {
    if (!state) return;
    auto out = output(opt.prefix + suffix);
    out << "{\"format\":\"cartmesh2d-compatible-state-v1\",\"kind\":"; jsonString(out,kind);
    out << ",\"meshSource\":"; jsonString(out,opt.mesh);
    out << ",\"contextBoundCheckpoint\":false,\"cellBasis\":\"1,dx/h,dy/h; h=actual cell diameter\","
        "\"cellOrder\":\"u0,uX,uY,v0,vX,vY,p0,pX,pY\",\"faceOrder\":\"u0,us,v0,vs\","
        "\"faceParameter\":\"s in [-1/2,1/2]; point=centre+s*(-Sy,Sx)\","
        "\"units\":\"physical velocity m/s; kinematic pressure m^2/s^2\",";
    stateArrays(out,*state);
    out << "}\n";
    out.close();
}
void history(const Options& opt, const CompatibleFlowResult2D& result) {
    auto out = output(opt.prefix + ".residuals.csv");
    out << "iteration,accepted,linearRestarts,matrixProducts,linearRelativeResidual,trials,alpha,pseudoStep,cellMomentum,faceMomentum,divergence,stateChange,residualNorm\n";
    for (const auto& it : result.iterations) {
        out << it.iteration << ',' << (it.accepted ? 1 : 0) << ',' << it.linearRestarts << ',' << it.matrixProducts << ',' << it.linearRelativeResidual << ',' << it.trials << ',' << it.alpha << ',' << it.pseudoStep;
        if (it.metrics) { const auto& m = *it.metrics; out << ',' << m.cellMomentum << ',' << m.faceMomentum << ',' << m.divergence << ',' << m.stateChange << ',' << m.residualNorm; }
        else out << ",,,,,";
        out << '\n';
    }
    out.close();
}
void fields(const Options& opt, const TopologyMesh2D& topology, const FvMesh2D& mesh, const CompatibleFlowResult2D& result) {
    if (!result.lastAccepted) return;
    const auto& state = *result.lastAccepted;
    auto cells = output(opt.prefix + ".cells.csv");
    cells << "cell,x,y,area,u,v,p,speed,uX,uY,vX,vY,pX,pY\n";
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        const auto& c = mesh.cells[i]; const auto& v = state.cells[i];
        cells << i << ',' << c.centre.x << ',' << c.centre.y << ',' << c.area << ',' << v[0] << ',' << v[3] << ',' << v[6] << ',' << std::hypot(v[0],v[3]) << ',' << v[1] << ',' << v[2] << ',' << v[4] << ',' << v[5] << ',' << v[7] << ',' << v[8] << '\n';
    }
    cells.close();
    auto faces = output(opt.prefix + ".faces.csv");
    faces << "face,owner,neighbour,x,y,Sx,Sy,u0,us,v0,vs,flux,fluxSlope\n";
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
        const auto& f = mesh.faces[i]; const auto& v = state.faces[i];
        faces << i << ',' << f.owner << ','; if (f.neighbour) faces << *f.neighbour; else faces << -1;
        faces << ',' << f.centre.x << ',' << f.centre.y << ',' << f.areaVector.x << ',' << f.areaVector.y;
        for (double x : v) faces << ',' << x;
        faces << ',' << v[0]*f.areaVector.x+v[2]*f.areaVector.y << ',' << v[1]*f.areaVector.x+v[3]*f.areaVector.y << '\n';
    }
    faces.close();
    auto field = output(opt.prefix + ".fields.json");
    field << "{\"format\":\"cartmesh2d-compatible-fields-v1\",\"sourceState\":\"last-accepted\",\"converged\":" << (result.converged() ? "true" : "false")
          << ",\"representation\":\"P1 cell coefficients and independent P1 face traces; see accepted state basis and units\",";
    stateArrays(field,state); field << "}\n"; field.close();
    std::string error;
    if (!writeLegacyVtk2D(topology,opt.prefix + ".vtk",&error)) throw std::runtime_error(error);
    std::ofstream vtk;
    vtk.exceptions(std::ios::failbit | std::ios::badbit);
    vtk.open(opt.prefix + ".vtk",std::ios::app);
    // Mesh writer has already opened CELL_DATA. Centroid P1 values are cell means.
    vtk << std::setprecision(17) << "VECTORS velocity double\n";
    for (const auto& v : state.cells) vtk << v[0] << ' ' << v[3] << " 0\n";
    vtk << "SCALARS pressure_kinematic double 1\nLOOKUP_TABLE default\n";
    for (const auto& v : state.cells) vtk << v[6] << '\n';
    vtk << "SCALARS speed double 1\nLOOKUP_TABLE default\n";
    for (const auto& v : state.cells) vtk << std::hypot(v[0],v[3]) << '\n';
    vtk.close();
}
void summary(const Options& opt, const FvMesh2D& mesh, const CompatibleFlowResult2D* result, const char* status, const std::string& reason, bool complete, double readSeconds, double solveSeconds, double exportSeconds, double totalSeconds) {
    const auto& c = opt.controls;
    const bool savedCheckpoint = std::filesystem::is_regular_file(opt.prefix + ".checkpoint");
    auto out = output(opt.prefix + ".summary.json.tmp");
    out << "{\"format\":\"cartmesh2d-compatible-flow-summary-v1\",\"discretization\":\"compatible\",\"status\":"; jsonString(out,status);
    out << ",\"reason\":"; jsonString(out,reason);
    out << ",\"exportsComplete\":" << (complete ? "true" : "false") << ",\"converged\":" << (complete && result && result->converged() ? "true" : "false")
        << ",\"numericallyConverged\":" << (result && result->converged() ? "true" : "false")
        << ",\"physicalAccuracyQualified\":false,\"contextBoundCheckpoint\":" << (savedCheckpoint ? "true" : "false")
        << ",\"checkpointAvailable\":" << (savedCheckpoint ? "true" : "false") << ",\"checkpointPath\":"; jsonString(out,opt.prefix + ".checkpoint");
    out << ",\"restartInput\":"; if(opt.restart.empty())out << "null"; else jsonString(out,opt.restart);
    out << ",\"resumed\":" << (result && result->resumed ? "true" : "false")
        << ",\"acceptedIterationsBefore\":" << (result ? result->acceptedIterationsBefore : 0)
        << ",\"totalAcceptedIterations\":" << (result && result->checkpoint ? result->checkpoint->acceptedIterations() : 0)
        << ",\"mesh\":"; jsonString(out,opt.mesh);
    out << ",\"case\":"; jsonString(out,opt.scenario);
    out << ",\"cells\":" << mesh.cells.size() << ",\"faces\":" << mesh.faces.size()
        << ",\"lastAcceptedAvailable\":" << (result && result->lastAccepted ? "true" : "false") << ",\"lastRejectedAvailable\":" << (result && result->lastRejected ? "true" : "false")
        << ",\"sourceState\":" << (result && result->lastAccepted ? "\"last-accepted\"" : "null")
        << ",\"attemptedIterations\":" << (result ? result->iterations.size() : 0)
        << ",\"acceptedIterations\":" << (result ? std::count_if(result->iterations.begin(),result->iterations.end(),[](const auto& it){return it.accepted;}) : 0)
        << ",\"lastAcceptedMetrics\":";
    std::optional<CompatibleFlowMetrics2D> acceptedMetrics;
    if(result && result->checkpoint)acceptedMetrics=result->checkpoint->lastAcceptedMetrics();
    if (result) for (const auto& it : result->iterations) if (it.accepted) acceptedMetrics=it.metrics;
    metrics(out,acceptedMetrics);
    out << ",\"controls\":{\"viscosity\":" << c.viscosity << ",\"referenceLength\":" << c.referenceLength << ",\"referenceVelocity\":" << c.referenceVelocity
        << ",\"equation\":\"" << (c.equation == CompatibleEquation2D::Stokes ? "stokes" : "navier-stokes")
        << "\",\"globalization\":\"" << (c.globalization == CompatibleGlobalization2D::PseudoTime ? "pseudo-time" : "backtracking")
        << "\",\"pressureInverse\":\"" << (c.pressureInverse == CompatiblePressureInverse2D::ViscousMass ? "viscous-mass" : "diagonal-schur") << '"'
        << ",\"equationTolerance\":" << c.equationTolerance << ",\"stateTolerance\":" << c.stateTolerance << ",\"linearTolerance\":" << c.linearTolerance
        << ",\"maximumIterations\":" << c.maximumIterations << ",\"maximumLinearRestarts\":" << c.maximumLinearRestarts << ",\"krylovDirections\":" << c.krylovDirections
        << ",\"quadratureOrder\":" << c.quadratureOrder << ",\"initialPseudoStep\":" << c.initialPseudoStep << ",\"maximumPseudoStep\":" << c.maximumPseudoStep
        << ",\"maximumBacktracks\":" << c.maximumBacktracks << ",\"armijo\":" << c.armijo << '}'
        << ",\"timingSeconds\":{\"readAndPrepare\":" << readSeconds << ",\"solve\":" << solveSeconds << ",\"exportBeforeSummary\":" << exportSeconds << ",\"elapsedBeforeSummary\":" << totalSeconds << "}}\n";
    out.close();
    std::filesystem::rename(opt.prefix + ".summary.json.tmp",opt.prefix + ".summary.json");
}
int run(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--help") { help(); return 0; }
    const auto start = Clock::now();
    auto opt = parse(argc,argv);
    const auto read = readCm2dTopology(opt.mesh);
    if (!read.valid()) throw std::invalid_argument("Invalid final mesh: " + read.error);
    const auto mesh = makeFvMesh2D(read.topology);
    FlowControls2D boundaries;
    boundaries.scenario = opt.scenario; boundaries.speed = opt.controls.referenceVelocity;
    if (opt.scenario == "custom") {
        std::ifstream input(opt.boundary);
        if (!input) throw std::invalid_argument("Cannot open boundary file");
        boundaries.boundaryConditions = readFlowBoundaryConditions2D(input,mesh,boundaries);
    } else boundaries.boundaryConditions = explicitFlowBoundaryPreset2D(mesh,boundaries);
    opt.controls.boundaries = compatibleFlowBoundaries2D(mesh,boundaries.boundaryConditions,opt.controls.referenceVelocity);
    boundaries.scenario = "custom";
    opt.controls.stopRequested = [] { return cancelled != 0; };
    opt.controls.iterationAccepted = progress;
    std::optional<CompatibleFlowCheckpoint2D> restart, publishedCheckpoint;
    if(!opt.restart.empty()) {
        std::ifstream input(opt.restart);
        if(!input)throw std::invalid_argument("Cannot open compatible checkpoint");
        restart=readCompatibleFlowCheckpoint2D(input,mesh);
    }
    opt.controls.checkpointAccepted=[&](const CompatibleFlowCheckpoint2D& checkpoint) {
        publishCheckpoint(opt,checkpoint);
        publishedCheckpoint=checkpoint;
    };
    const auto prepared = Clock::now();
    const double readSeconds = std::chrono::duration<double>(prepared-start).count();
    ensureFresh(opt.prefix);
    summary(opt,mesh,nullptr,"running","Preparing native solve",false,readSeconds,0,0,readSeconds);
    CompatibleFlowResult2D result;
    double solveSeconds = 0;
    try {
        Signals signals;
        const auto solveStart = Clock::now();
        result = restart ? resumeCompatibleIncompressible2D(mesh,opt.controls,*restart)
                         : solveCompatibleIncompressible2D(mesh,opt.controls);
        const auto solved = Clock::now();
        solveSeconds = std::chrono::duration<double>(solved-solveStart).count();
        // On a failed resumed attempt, preserve the existing accepted state too.
        if(result.checkpoint)publishCheckpoint(opt,*result.checkpoint);
        // Preserve complete solver evidence before any derived visualization export.
        saveState(opt,result.seed,".seed.json","seed");
        saveState(opt,result.lastAccepted,".accepted.json","accepted-iterate");
        saveState(opt,result.lastRejected,".rejected.json","rejected-trial");
        history(opt,result);
        auto boundary = output(opt.prefix + ".boundaries");
        writeFlowBoundaryConditions2D(boundary,mesh,boundaries); boundary.close();
        fields(opt,read.topology,mesh,result);
        const auto exported = Clock::now();
        summary(opt,mesh,&result,stopName(result.stop),result.reason,true,readSeconds,solveSeconds,std::chrono::duration<double>(exported-solved).count(),std::chrono::duration<double>(exported-start).count());
        std::cout << "{\"type\":\"compatible-flow-result\",\"status\":"; jsonString(std::cout,stopName(result.stop));
        std::cout << ",\"converged\":" << (result.converged() ? "true" : "false") << ",\"reason\":"; jsonString(std::cout,result.reason); std::cout << "}\n";
        return result.converged() ? 0 : result.stop == CompatibleFlowStop2D::Cancelled ? 130 : 2;
    } catch (const std::exception& e) {
        // A callback/export exception must still expose the last successfully
        // published checkpoint; it does not manufacture derived flow fields.
        if(publishedCheckpoint)result.checkpoint=publishedCheckpoint;
        // An incomplete export must never advertise a completed deliverable.
        try { summary(opt,mesh,&result,"failed",e.what(),false,readSeconds,solveSeconds,0,std::chrono::duration<double>(Clock::now()-start).count()); }
        catch (const std::exception& exportError) { std::cerr << "Could not publish failure summary: " << exportError.what() << '\n'; }
        throw;
    }
}
}
std::optional<int> tryCompatibleFlowCLI2D(int argc, char** argv) {
    try {
        std::optional<std::string> selected;
        for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--discretization") {
            if (selected || i + 1 == argc) throw std::invalid_argument("Duplicate or missing --discretization");
            selected = argv[++i];
        }
        if (!selected || *selected == "collocated") return std::nullopt;
        if (*selected != "compatible") throw std::invalid_argument("--discretization must be collocated or compatible");
        return run(argc,argv);
    } catch (const std::exception& e) { std::cerr << "Compatible flow CLI: " << e.what() << '\n'; return 1; }
}
