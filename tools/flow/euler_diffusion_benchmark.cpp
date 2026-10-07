// Research harness: all PDE, characteristic boundaries and SDIRK stages use
// Euler API.
#include "cartmesh2d/fv/Euler2D.hpp"
#include "cartmesh2d/fv/EulerCheckpoint2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
volatile std::sig_atomic_t stopRequested=0;
void requestStop(int){stopRequested=1;}
void require(bool b, const char *s) {
  if (!b)
    throw std::runtime_error(s);
}
std::vector<EulerBoundary2D> read(const std::string &path, const FvMesh2D &m) {
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);
  require(line == "CM2D_EULER_BOUNDARY 3", "research boundary v3 required");
  std::getline(in, line);
  require(line == "MESH " + std::to_string(m.cells.size()) + " " +
                      std::to_string(m.faces.size()),
          "boundary dimensions");
  std::vector<EulerBoundary2D> out;
  bool ended = false;
  while (std::getline(in, line)) {
    if (line == "END") {
      ended = true;
      break;
    }
    EulerBoundary2D b;
    size_t owner;
    double x, y, sx, sy;
    std::string kind, partner, thermal, extra;
    std::istringstream row(line);
    require(bool(row >> b.face >> owner >> x >> y >> sx >> sy >> kind >>
                 b.reference.density >> b.reference.u >> b.reference.v >>
                 b.reference.pressure >> partner >> std::quoted(b.name) >>
                 thermal >> b.thermalValue >> b.wallVelocity.x >>
                 b.wallVelocity.y),
            "boundary row");
    require(!(row >> extra) && b.face < m.faces.size(), "boundary trailing/id");
    const auto &f = m.faces[b.face];
    require(owner == f.owner && x == f.centre.x && y == f.centre.y &&
                sx == f.areaVector.x && sy == f.areaVector.y,
            "boundary geometry binding");
    if (kind == "farfield")
      b.kind = EulerBoundaryKind2D::Farfield;
    else if (kind == "pressure-outlet")
      b.kind = EulerBoundaryKind2D::PressureOutlet;
    else if (kind == "no-slip-wall")
      b.kind = EulerBoundaryKind2D::NoSlipWall;
    else
      throw std::runtime_error("unsupported research boundary kind");
    require(partner == "-", "periodic unsupported in harness");
    if (thermal == "temperature")
      b.thermalKind = HeatBoundaryKind2D::Temperature;
    else if (thermal == "insulated")
      b.thermalKind = HeatBoundaryKind2D::Insulated;
    else if (thermal == "flux")
      b.thermalKind = HeatBoundaryKind2D::OutwardFlux;
    else
      throw std::runtime_error("thermal kind");
    out.push_back(b);
  }
  require(ended && !(in >> line), "boundary terminator");
  return out;
}
} // namespace
int compareAll(int argc, char **argv) {
  auto start = std::chrono::steady_clock::now();
  auto elapsed = [&] {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         start)
        .count();
  };
  auto rd = readCm2dTopology(argv[1]);
  require(rd.valid(), "mesh read");
  auto mesh = makeFvMesh2D(rd.topology);
  auto bc = read(argv[2], mesh);
  IdealGas2D gas;
  EulerTransport2D physics{.025758750694444447, 1.846e-5};
  const std::array<std::string, 3> names{"corrected", "hybrid-heat", "hybrid"};
  const std::array<EulerDiffusionScheme2D, 3> schemes{
      EulerDiffusionScheme2D::Corrected, EulerDiffusionScheme2D::HybridHeat,
      EulerDiffusionScheme2D::Hybrid};
  std::string base = argv[3];
  double end = std::stod(argv[5]), dt = std::stod(argv[6]),
         budget = std::stod(argv[7]);
  require(end > 0 && dt > 0 && budget > 0,
          "positive horizon/step/budget required");
  std::array<EulerState2D, 3> states;
  std::array<EulerStepResult2D, 3> last;
  std::array<std::unique_ptr<EulerStepper2D>, 3> solvers;
  std::array<double, 3> costs{}, construct{};
  std::array<size_t, 3> spatial{}, newton{}, krylov{}, retries{}, discarded{};
  std::array<std::ofstream, 3> hist, bhist;
  std::ofstream attempts(base + ".attempts.csv");
  attempts << std::setprecision(17)
           << "time,scheme,proposedDt,returnedDt,committed,rejected,spatial,"
              "newton,krylov,seconds\n";
  auto save = [&](size_t i) {
    std::ofstream f(base + "." + names[i] + ".checkpoint.tmp");
    writeEulerCheckpoint2D(f, mesh, bc, gas, states[i],
                           "diffusion-research/" + names[i], physics);
    f.close();
    std::filesystem::rename(base + "." + names[i] + ".checkpoint.tmp",
                            base + "." + names[i] + ".checkpoint");
  };
  auto ref = std::find_if(bc.begin(), bc.end(), [](const auto &b) {
    return b.kind == EulerBoundaryKind2D::Farfield;
  });
  require(ref != bc.end(), "farfield reference");
  for (size_t i = 0; i < 3; ++i) {
    if (argc >= 9) {
      std::ifstream f(std::string(argv[8]) + "." + names[i] + ".checkpoint");
      states[i] = readEulerCheckpoint2D(
          f, mesh, bc, gas, "diffusion-research/" + names[i], physics);
    } else
      states[i].cells.assign(mesh.cells.size(),
                             eulerConservative2D(ref->reference, gas));
    save(i);
    hist[i].open(base + "." + names[i] + ".history.csv");
    bhist[i].open(base + "." + names[i] + ".boundary-history.csv");
    hist[i]
        << std::setprecision(17)
        << "step,time,dt,rhoMin,rhoMax,pMin,pMax,tMin,tMax,balanceMass,"
           "balanceMx,balanceMy,balanceEnergy,heat,wallWork,maxCellBalance\n";
    bhist[i] << std::setprecision(17)
             << "step,time,name,mass,mx,my,energy,heat,viscousWork\n";
  }
  require(states[0].time == states[1].time && states[0].time == states[2].time,
          "restart common time mismatch");
  std::string status = "complete", error;
  bool restartExact = true;
  double restartCheckSeconds = 0;
  size_t restartChecks = 0;
  std::string activePhase = "construction";
  size_t activeScheme = 0;
  EulerStepControls2D control;
  control.endTime = end;
  control.maximumStep = dt;
  control.integrator = EulerTimeIntegrator2D::Sdirk2;
  control.order = 2;
  if (argc == 10) {
    require(std::string(argv[9]) == "strict-step", "unknown comparison option");
    control.maximumRetries = 0;
  }
  control.fluxScheme = EulerFluxScheme2D::Hllc;
  control.interrupted = [&] { return stopRequested||elapsed() > budget; };
  try {
    for (size_t i = 0; i < 3; ++i) {
      double t = elapsed();
      activeScheme = i;
      solvers[i] = std::make_unique<EulerStepper2D>(
          mesh, bc, gas, physics, WallGradient2D::Linear, schemes[i]);
      construct[i] = elapsed() - t;
      costs[i] += construct[i];
    }
    while (states[0].time < end) {
      double trial = std::min(dt, end - states[0].time);
      std::array<EulerStepResult2D, 3> candidate;
      bool committed = false;
      for (size_t cycle = 0; !committed; ++cycle) {
        require(cycle < 32, "common-step transaction budget");
        std::array<double, 3> callSeconds{};
        double minimum = trial;
        for (size_t i = 0; i < 3; ++i) {
          control.maximumStep = trial;
          control.diffusionScheme = schemes[i];
          double t = elapsed();
          try {
            activeScheme = i;
            activePhase = "advance";
            candidate[i] = solvers[i]->advance(states[i], control);
          } catch (...) {
            costs[i] += elapsed() - t;
            throw;
          }
          callSeconds[i] = elapsed() - t;
          costs[i] += callSeconds[i];
          minimum = std::min(minimum, candidate[i].step);
          spatial[i] += candidate[i].spatialEvaluations;
          newton[i] += candidate[i].nonlinearIterations;
          krylov[i] += candidate[i].linearIterations;
          retries[i] += candidate[i].rejectedCandidates;
        }
        committed = candidate[0].step == minimum &&
                    candidate[1].step == minimum &&
                    candidate[2].step == minimum;
        for (size_t i = 0; i < 3; ++i) {
          attempts << states[i].time << ',' << names[i] << ',' << trial << ','
                   << candidate[i].step << ',' << committed << ','
                   << candidate[i].rejectedCandidates << ','
                   << candidate[i].spatialEvaluations << ','
                   << candidate[i].nonlinearIterations << ','
                   << candidate[i].linearIterations << ',' << callSeconds[i]
                   << '\n';
          if (!committed)
            ++discarded[i];
        }
        attempts.flush();
        trial = minimum;
      }
      for (size_t i = 0; i < 3; ++i) {
        double t = elapsed();
        const auto &r = candidate[i];
        double rmin = 1e300, rmax = 0, pmin = 1e300, pmax = 0, tmin = 1e300,
               tmax = 0;
        for (const auto &q : r.state.cells) {
          auto p = eulerPrimitive2D(q, gas);
          double temp = p.pressure / (p.density * gas.gasConstant);
          rmin = std::min(rmin, p.density);
          rmax = std::max(rmax, p.density);
          pmin = std::min(pmin, p.pressure);
          pmax = std::max(pmax, p.pressure);
          tmin = std::min(tmin, temp);
          tmax = std::max(tmax, temp);
        }
        hist[i] << r.state.steps << ',' << r.state.time << ',' << r.step << ','
                << rmin << ',' << rmax << ',' << pmin << ',' << pmax << ','
                << tmin << ',' << tmax;
        for (double v : r.balanceError)
          hist[i] << ',' << v;
        hist[i] << ',' << r.boundaryHeat << ',' << r.boundaryViscousWork << ','
                << r.maximumCellBalanceError << '\n';
        for (const auto &b : bc) {
          bhist[i] << r.state.steps << ',' << r.state.time << ',' << b.name;
          for (double v : r.faceFlux[b.face])
            bhist[i] << ',' << v;
          bhist[i] << ',' << r.faceHeatFlux[b.face] << ','
                   << r.faceViscousFlux[b.face][2] << '\n';
        }
        states[i] = r.state;
        last[i] = std::move(candidate[i]);
        save(i);
        hist[i].flush();
        costs[i] += elapsed() - t;
      }
      require(states[0].time == states[1].time &&
                  states[0].time == states[2].time,
              "committed clocks differ");
      // Real mesh checkpoint+solver reconstruction check, outside comparison
      // cost.
      if (states[0].steps == 3) {
        double t = elapsed();
        activePhase = "restart-validation";
        for (size_t i = 0; i < 3; ++i) {
          std::ifstream f(base + "." + names[i] + ".checkpoint");
          auto restored = readEulerCheckpoint2D(
              f, mesh, bc, gas, "diffusion-research/" + names[i], physics);
          control.maximumStep = std::min(dt, end - states[i].time);
          control.diffusionScheme = schemes[i];
          if (control.maximumStep > 0) {
            auto a = solvers[i]->advance(states[i], control);
            auto b = EulerStepper2D(mesh, bc, gas, physics,
                                    WallGradient2D::Linear, schemes[i])
                         .advance(restored, control);
            ++restartChecks;
            restartExact = restartExact && a.state.cells == b.state.cells &&
                           a.state.time == b.state.time;
            require(restartExact, "real mesh restart mismatch");
          }
        }
        restartCheckSeconds = elapsed() - t;
      }
    }
  } catch (const std::exception &e) {
    status = elapsed() > budget ? "budget-exhausted-last-common-accepted"
                                : "failed-last-common-accepted";
    error = names[activeScheme] + "/" + activePhase + ": " + e.what();
    std::ofstream f(base + ".failure.txt");
    f << error << '\n';
  }
  for (size_t i = 0; i < 3; ++i) {
    double t = elapsed();
    save(i);
    std::ofstream cells(base + "." + names[i] + ".cells.csv");
    cells << std::setprecision(17)
          << "cell,x,y,area,rho,rhou,rhov,rhoE,u,v,p,T\n";
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
      auto p = eulerPrimitive2D(states[i].cells[c], gas);
      cells << c << ',' << mesh.cells[c].centre.x << ','
            << mesh.cells[c].centre.y << ',' << mesh.cells[c].area;
      for (double v : states[i].cells[c])
        cells << ',' << v;
      cells << ',' << p.u << ',' << p.v << ',' << p.pressure << ','
            << p.pressure / (p.density * gas.gasConstant) << '\n';
    }
    if (!last[i].faceFlux.empty()) {
      std::ofstream f(base + "." + names[i] + ".faces.csv");
      f << std::setprecision(17)
        << "face,owner,neighbour,x,y,sx,sy,mass,mx,my,energy,heat,viscX,viscY,"
           "work\n";
      for (size_t id = 0; id < mesh.faces.size(); ++id) {
        const auto &g = mesh.faces[id];
        f << id << ',' << g.owner << ','
          << (g.neighbour ? std::to_string(*g.neighbour) : "-1") << ','
          << g.centre.x << ',' << g.centre.y << ',' << g.areaVector.x << ','
          << g.areaVector.y;
        for (double v : last[i].faceFlux[id])
          f << ',' << v;
        f << ',' << last[i].faceHeatFlux[id];
        for (double v : last[i].faceViscousFlux[id])
          f << ',' << v;
        f << '\n';
      }
    }
    costs[i] += elapsed() - t;
    std::ofstream report(base + "." + names[i] + ".json");
    report << std::setprecision(17) << "{\"status\":" << std::quoted(status)
           << ",\"error\":" << std::quoted(error)
           << ",\"scheme\":" << std::quoted(names[i])
           << ",\"time\":" << states[i].time << ",\"steps\":" << states[i].steps
           << ",\"requestedDt\":" << dt
           << ",\"constructionSeconds\":" << construct[i]
           << ",\"spatialEvaluationsIncludingDiscarded\":" << spatial[i]
           << ",\"newtonIncludingDiscarded\":" << newton[i]
           << ",\"krylovIncludingDiscarded\":" << krylov[i]
           << ",\"internalRejectedCandidates\":" << retries[i]
           << ",\"commonTransactionDiscarded\":" << discarded[i]
           << ",\"exclusiveSecondsWithExport\":" << costs[i]
           << ",\"restartChecks\":" << restartChecks
           << ",\"restartExact\":" << (restartExact ? "true" : "false")
           << "}\n";
  }
  std::ofstream overall(base + ".json");
  overall << std::setprecision(17) << "{\"status\":" << std::quoted(status)
          << ",\"error\":" << std::quoted(error)
          << ",\"physicalTime\":" << states[0].time
          << ",\"commonSteps\":" << states[0].steps
          << ",\"fullProcessSeconds\":" << elapsed()
          << ",\"restartCheckSeconds\":" << restartCheckSeconds
          << ",\"restartChecks\":" << restartChecks
          << ",\"restartExact\":" << (restartExact ? "true" : "false") << "}\n";
  std::cout << status << " common time=" << states[0].time
            << " seconds=" << elapsed() << " " << error << '\n';
  return status == "complete" ? 0 : 2;
}
int main(int argc, char **argv) {
  try {
    if ((argc >= 8 && argc <= 10) && std::string(argv[4]) == "all")
      return compareAll(argc, argv);
    require(argc == 8 || argc == 9,
            "mesh boundary prefix corrected|hybrid-heat|hybrid end dt "
            "budgetSeconds [restart]");
    auto start = std::chrono::steady_clock::now();
    auto elapsed = [&] {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                           start)
          .count();
    };
    auto rd = readCm2dTopology(argv[1]);
    require(rd.valid(), "mesh read failed");
    auto m = makeFvMesh2D(rd.topology);
    auto bc = read(argv[2], m);
    IdealGas2D gas;
    EulerTransport2D physics{.025758750694444447, 1.846e-5};
    if(const char* value=std::getenv("CARTMESH_RESEARCH_CONDUCTIVITY"))physics.thermalConductivity=std::stod(value);
    std::string prefix = argv[3], name = argv[4],
                label = "diffusion-research/" + name;
    EulerStepControls2D ctl;
    ctl.integrator = EulerTimeIntegrator2D::Sdirk2;
    if (const char *pc = std::getenv("CARTMESH_RESEARCH_PRECONDITIONER")) {
      require(std::string(pc)=="diagonal" || std::string(pc)=="frozen-flux-ilu0", "unknown research preconditioner");
      if (std::string(pc)=="frozen-flux-ilu0") ctl.implicitPreconditioner=EulerImplicitPreconditioner2D::FrozenFluxIlu0;
    }
    ctl.fluxScheme = EulerFluxScheme2D::Hllc;
    ctl.order = 2;
    if(name=="gauss2")ctl.faceQuadrature=EulerFaceQuadrature2D::Gauss2;
    ctl.wallGradient = WallGradient2D::Linear;
    ctl.endTime = std::stod(argv[5]);
    ctl.maximumStep = std::stod(argv[6]);
    double budget = std::stod(argv[7]);
    std::signal(SIGINT,requestStop);std::signal(SIGTERM,requestStop);
    if(const char* tolerance=std::getenv("CARTMESH_RESEARCH_NONLINEAR_TOLERANCE"))ctl.nonlinearTolerance=std::stod(tolerance);
    ctl.interrupted = [&] { return stopRequested||elapsed() > budget; };
    if (name == "corrected" || name == "gauss2")
      ctl.diffusionScheme = EulerDiffusionScheme2D::Corrected;
    else if (name == "hybrid-heat")
      ctl.diffusionScheme = EulerDiffusionScheme2D::HybridHeat;
    else if (name == "hybrid")
      ctl.diffusionScheme = EulerDiffusionScheme2D::Hybrid;
    else
      throw std::runtime_error("scheme");
    EulerState2D state;
    if (argc == 9) {
      std::ifstream f(argv[8]);
      state = readEulerCheckpoint2D(f, m, bc, gas, std::getenv("CARTMESH_RESEARCH_CHECKPOINT_CASE")?std::getenv("CARTMESH_RESEARCH_CHECKPOINT_CASE"):label, physics);
    } else {
      auto ref = std::find_if(bc.begin(), bc.end(), [](const auto &b) {
        return b.kind == EulerBoundaryKind2D::Farfield;
      });
      require(ref != bc.end(), "farfield reference required");
      state.cells.assign(m.cells.size(),
                         eulerConservative2D(ref->reference, gas));
    }
    auto save = [&] {
      std::ofstream f(prefix + ".checkpoint.tmp");
      writeEulerCheckpoint2D(f, m, bc, gas, state, label, physics);
      f.close();
      std::filesystem::rename(prefix + ".checkpoint.tmp",
                              prefix + ".checkpoint");
    };
    save();
    std::ofstream history(prefix + ".history.csv"),
        bfile(prefix + ".boundary-history.csv");
    history << std::setprecision(17)
            << "step,time,dt,rejected,spatial,newton,krylov,rhoMin,rhoMax,pMin,"
               "pMax,tMin,tMax,balanceMass,balanceMx,balanceMy,balanceEnergy,"
               "heat,wallWork,maxCellBalance,elapsed,stageDefect,outputDefect\n";
    bfile << std::setprecision(17)
          << "step,time,name,mass,mx,my,energy,heat,viscousWork\n";
    double construction = 0;
    std::ofstream rejectionLog(prefix + ".rejections.csv");
    size_t spatial = 0, newton = 0, krylov = 0, rejected = 0;
    EulerStepResult2D last;
    std::ofstream nonlinearTrace,failedWork;
    if(std::getenv("CARTMESH_RESEARCH_IMPLICIT_TRACE")){
      failedWork.open(prefix+".failed-work.csv");failedWork<<"time,spatial,newton,krylov,reason\n"<<std::setprecision(17);
      ctl.implicitFailure=[&](std::size_t spatialCalls,std::size_t newtonIterations,std::size_t linearIterations,const std::string& reason){failedWork<<state.time<<','<<spatialCalls<<','<<newtonIterations<<','<<linearIterations<<','<<std::quoted(reason)<<'\n';};
      nonlinearTrace.open(prefix+".nonlinear.csv");nonlinearTrace<<"time,stage,iteration,h,defect,converged,spatial,krylov\n"<<std::setprecision(17);
      ctl.implicitIteration=[&](const EulerImplicitIteration2D& e){nonlinearTrace<<state.time<<','<<e.stage<<','<<e.iteration<<','<<e.stageStep<<','<<e.maximumScaledDefect<<','<<e.converged<<','<<e.spatialEvaluations<<','<<e.linearIterations<<'\n';};
    }
    const bool readOnly=std::getenv("CARTMESH_RESEARCH_SNAPSHOT")!=nullptr;
    std::string status = readOnly?"read-only-space-snapshot":"complete", error;
    try {
      double begin = elapsed();
      EulerStepper2D solver(m, bc, gas, physics, ctl.wallGradient,
                            ctl.diffusionScheme);
      construction = elapsed() - begin;
      while (!readOnly && state.time < *ctl.endTime) {
        auto old = state;
        auto r = solver.advance(state, ctl);
        require(old.cells == state.cells && old.time == state.time,
                "API mutated input");
        double rmin = 1e300, rmax = 0, pmin = 1e300, pmax = 0, tmin = 1e300,
               tmax = 0;
        for (const auto &q : r.state.cells) {
          auto p = eulerPrimitive2D(q, gas);
          double temp = p.pressure / (p.density * gas.gasConstant);
          rmin = std::min(rmin, p.density);
          rmax = std::max(rmax, p.density);
          pmin = std::min(pmin, p.pressure);
          pmax = std::max(pmax, p.pressure);
          tmin = std::min(tmin, temp);
          tmax = std::max(tmax, temp);
        }
        history << r.state.steps << ',' << r.state.time << ',' << r.step << ','
                << r.rejectedCandidates << ',' << r.spatialEvaluations << ','
                << r.nonlinearIterations << ',' << r.linearIterations << ','
                << rmin << ',' << rmax << ',' << pmin << ',' << pmax << ','
                << tmin << ',' << tmax;
        for (double v : r.balanceError)
          history << ',' << v;
        history << ',' << r.boundaryHeat << ',' << r.boundaryViscousWork << ','
                << r.maximumCellBalanceError << ',' << elapsed() << ',' << r.maximumAcceptedStageDefect << ',' << r.maximumStageOutputDefect << '\n';
        for (const auto &b : bc) {
          const auto f = b.face;
          bfile << r.state.steps << ',' << r.state.time << ',' << b.name;
          for (double v : r.faceFlux[f])
            bfile << ',' << v;
          bfile << ',' << r.faceHeatFlux[f] << ',' << r.faceViscousFlux[f][2]
                << '\n';
        }
        spatial += r.spatialEvaluations;
        newton += r.nonlinearIterations;
        krylov += r.linearIterations;
        rejected += r.rejectedCandidates;
        if (r.rejectedCandidates) rejectionLog << std::setprecision(17) << old.time << "," << r.step << "," << r.rejectedCandidates << "," << std::quoted(r.lastRejectedReason) << "\n";
        state = r.state;
        last = std::move(r);
        save();
        history.flush();
      }
    } catch (const std::exception &e) {
      status = elapsed() > budget ? "budget-exhausted-last-accepted"
                                  : "failed-last-accepted";
      if(stopRequested)status="cancelled-last-accepted";
      error = e.what();
      save();
      std::ofstream failure(prefix + ".failure.txt");
      failure << error << '\n';
    }
    if(std::getenv("CARTMESH_RESEARCH_RECONSTRUCTION")){
      EulerStepper2D diagnosticSolver(m,bc,gas,physics,ctl.wallGradient,ctl.diffusionScheme);
      const auto snapshot=diagnosticSolver.spatialSnapshot(state,ctl);
      std::ofstream faceOutput(prefix+".instantaneous-faces.csv");
      faceOutput<<std::setprecision(17)<<"face,owner,neighbour,x,y,sx,sy,mass,mx,my,energy,heat,viscX,viscY,work\n";
      for(std::size_t f=0;f<m.faces.size();++f){const auto& face=m.faces[f];faceOutput<<f<<','<<face.owner<<','<<(face.neighbour?std::to_string(*face.neighbour):"-1")<<','<<face.centre.x<<','<<face.centre.y<<','<<face.areaVector.x<<','<<face.areaVector.y;
        for(double value:snapshot.faceFlux[f])faceOutput<<','<<value;
        faceOutput<<','<<snapshot.faceHeatFlux[f];for(double value:snapshot.faceViscousFlux[f])faceOutput<<','<<value;faceOutput<<'\n';}
      std::ofstream output(prefix+".reconstruction.csv");
      output<<std::setprecision(17)<<"cell,component,rawX,rawY,limitedX,limitedY,thetaLocalFrame\n";
      for(std::size_t c=0;c<m.cells.size();++c)for(std::size_t k=0;k<4;++k)output<<c<<','<<k<<','<<snapshot.rawGradient[c][k].x<<','<<snapshot.rawGradient[c][k].y<<','<<snapshot.limitedGradient[c][k].x<<','<<snapshot.limitedGradient[c][k].y<<','<<snapshot.limiterTheta[c][k]<<'\n';
    }
    std::ofstream cells(prefix + ".cells.csv");
    cells << std::setprecision(17)
          << "cell,x,y,area,rho,rhou,rhov,rhoE,u,v,p,T\n";
    for (size_t c = 0; c < state.cells.size(); ++c) {
      auto p = eulerPrimitive2D(state.cells[c], gas);
      cells << c << ',' << m.cells[c].centre.x << ',' << m.cells[c].centre.y
            << ',' << m.cells[c].area;
      for (double v : state.cells[c])
        cells << ',' << v;
      cells << ',' << p.u << ',' << p.v << ',' << p.pressure << ','
            << p.pressure / (p.density * gas.gasConstant) << '\n';
    }
    if (!last.faceFlux.empty()) {
      std::ofstream f(prefix + ".faces.csv");
      f << std::setprecision(17)
        << "face,owner,neighbour,x,y,sx,sy,mass,mx,my,energy,heat,viscX,viscY,"
           "work\n";
      for (size_t id = 0; id < m.faces.size(); ++id) {
        const auto &g = m.faces[id];
        f << id << ',' << g.owner << ','
          << (g.neighbour ? std::to_string(*g.neighbour) : "-1") << ','
          << g.centre.x << ',' << g.centre.y << ',' << g.areaVector.x << ','
          << g.areaVector.y;
        for (double v : last.faceFlux[id])
          f << ',' << v;
        f << ',' << last.faceHeatFlux[id];
        for (double v : last.faceViscousFlux[id])
          f << ',' << v;
        f << '\n';
      }
    }
    std::ofstream report(prefix + ".json");
    report << std::setprecision(17) << "{\"status\":" << std::quoted(status)
           << ",\"error\":" << std::quoted(error)
           << ",\"scheme\":" << std::quoted(name) << ",\"time\":" << state.time
           << ",\"steps\":" << state.steps
           << ",\"requestedDt\":" << ctl.maximumStep
           << ",\"constructionSeconds\":" << construction
           << ",\"acceptedSpatialEvaluations\":" << spatial
           << ",\"acceptedNewtonIterations\":" << newton
           << ",\"acceptedKrylovIterations\":" << krylov
           << ",\"rejectedWithinAcceptedCalls\":" << rejected
           << ",\"fullSeconds\":" << elapsed() << "}\n";
    std::cout << status << " time=" << state.time << " steps=" << state.steps
              << " seconds=" << elapsed() << " " << error << '\n';
    return status == "complete" || status == "read-only-space-snapshot" ? 0 : 2;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
