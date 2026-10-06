// Isolated constant-mu viscous manufactured vector problem, NOT curved-wall
// Couette or compressible Navier-Stokes. Velocity is a centroid point estimate.
// Dense cell matrices are sampled from native shared residuals, never rebuilt
// from a separate PDE. Dense export/steady solve are research costs.
#include "cartmesh2d/fv/HybridViscous2D.hpp"
#include "cartmesh2d/fv/detail/FlowLinearSystem2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
using namespace cartmesh2d;
using namespace cartmesh2d::fv;
namespace {
Vector2D exact(Point2D p) {
  double r2 = p.x * p.x + p.y * p.y;
  return {-(.3 + .2 / r2) * p.y, (.3 + .2 / r2) * p.x};
}
const double gx[] = {-.9602898564975363, -.7966664774136267, -.5255324099163290,
                     -.1834346424956498, .1834346424956498,  .5255324099163290,
                     .7966664774136267,  .9602898564975363};
const double gw[] = {.1012285362903763, .2223810344533745, .3137066458778873,
                     .3626837833783620, .3626837833783620, .3137066458778873,
                     .2223810344533745, .1012285362903763};
Vector2D triangle(Point2D a, Point2D b, Point2D c) {
  Vector2D out{};
  double cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  for (size_t i = 0; i < 8; ++i)
    for (size_t j = 0; j < 8; ++j) {
      double u = .5 * (gx[i] + 1), v = .5 * (gx[j] + 1);
      auto q = exact({a.x + u * (b.x - a.x) + u * v * (c.x - b.x),
                      a.y + u * (b.y - a.y) + u * v * (c.y - b.y)});
      double w = .25 * gw[i] * gw[j] * u * cross;
      out.x += w * q.x;
      out.y += w * q.y;
    }
  return out;
}
Vector2D referenceFlux(Point2D p, Vector2D s, double mu) {
  double r2 = p.x * p.x + p.y * p.y;
  return {2 * .2 * mu / r2 *
              (-2 * p.x * p.y / r2 * s.x + (p.x * p.x - p.y * p.y) / r2 * s.y),
          2 * .2 * mu / r2 *
              ((p.x * p.x - p.y * p.y) / r2 * s.x + 2 * p.x * p.y / r2 * s.y)};
}
double seconds(std::chrono::steady_clock::time_point t) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t)
      .count();
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 4) {
      throw std::runtime_error("mesh prefix hybrid|linear");
    }
    auto start = std::chrono::steady_clock::now();
    auto rd = readCm2dTopology(argv[1]);
    if (!rd.valid())
      throw std::runtime_error(rd.error);
    auto m = makeFvMesh2D(rd.topology);
    std::string prefix = argv[2], method = argv[3];
    if (method != "hybrid" && method != "linear")
      throw std::runtime_error("invalid method");
    const double mu = .02;
    size_t n = 2 * m.cells.size();
    std::vector<ViscousBoundary2D> bc, zero;
    for (size_t f = 0; f < m.faces.size(); ++f)
      if (!m.faces[f].neighbour) {
        bc.push_back(
            {f, ViscousBoundaryKind2D::Velocity, exact(m.faces[f].centre), {}});
        zero.push_back({f, ViscousBoundaryKind2D::Velocity, {}, {}});
      }
    auto t = std::chrono::steady_clock::now();
    HybridViscousOperator2D h(m, bc, mu), hz(m, zero, mu);
    ViscousStressOperator2D l(m, bc, mu), lz(m, zero, mu);
    double assembly = seconds(t);
    std::vector<double> rho(m.cells.size(), 1);
    std::vector<Vector2D> u(m.cells.size());
    for (size_t c = 0; c < u.size(); ++c)
      u[c] = exact(m.cells[c].centre);
    auto residual = [&](const std::vector<Vector2D> &v, bool homogeneous) {
      if (method == "hybrid")
        return (homogeneous ? hz : h).evaluate(v, rho, 1e-13).cellResidual;
      return (homogeneous ? lz : l).evaluate(v, rho).cellResidual;
    };
    std::ofstream matrix(prefix + ".cell-matrix.csv");
    matrix << std::setprecision(17) << "row,column,value,area\n";
    std::vector<double> a(n * n);
    t = std::chrono::steady_clock::now();
    for (size_t j = 0; j < n; ++j) {
      std::vector<Vector2D> v(u.size());
      if (j % 2)
        v[j / 2].y = 1;
      else
        v[j / 2].x = 1;
      auto r = residual(v, true);
      for (size_t i = 0; i < n; ++i) {
        a[i * n + j] = r[i / 2][i % 2];
        matrix << i << ',' << j << ',' << a[i * n + j] << ','
               << m.cells[i / 2].area << '\n';
      }
    }
    double matrixSeconds = seconds(t);
    matrix.close();
    t = std::chrono::steady_clock::now();
    std::vector<Vector2D> z(u.size());
    auto forcing = residual(z, false);
    std::vector<std::pair<size_t, size_t>> edges;
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j)
        edges.emplace_back(i, j);
    detail::SparsePattern2D pattern(n, edges);
    detail::SparseSystem2D system(pattern);
    for (size_t i = 0; i < n; ++i) {
      system.diag[i] = a[i * n + i];
      system.rhs[i] = -forcing[i / 2][i % 2];
      for (size_t k = pattern.rows[i]; k < pattern.rows[i + 1]; ++k)
        system.off[k] = a[i * n + pattern.columns[k]];
    }
    detail::LinearVector2D x(n);
    detail::LinearWorkspace2D work(n);
    size_t its =
        method == "hybrid"
            ? system.solvePressure(
                  x, work, detail::LinearPressureMethod2D::Jacobi, 1e-11)
            : system.solve(x, work, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity(),
                           detail::LinearSolveMethod2D::ILU0, 1e-11);
    double steadySeconds = seconds(t);
    std::vector<Vector2D> solved(u.size());
    for (size_t c = 0; c < u.size(); ++c)
      solved[c] = {x[2 * c], x[2 * c + 1]};
    std::ofstream cells(prefix + ".cells.csv");
    cells << std::setprecision(17)
          << "cell,x,y,area,analyticX,analyticY,steadyX,steadyY,meanX,meanY\n";
    double error = 0, meanError = 0, referenceMeanOffset = 0,
           referenceQuadratureGap = 0, norm = 0, balance = 0;
    auto sr = residual(solved, false);
    for (size_t c = 0; c < u.size(); ++c) {
      error += m.cells[c].area *
               std::hypot(solved[c].x - u[c].x, solved[c].y - u[c].y);
      norm += m.cells[c].area * std::hypot(u[c].x, u[c].y);
      Vector2D mean{};
      const auto &poly = rd.topology.cells[c].vertices;
      for (size_t j = 0; j < poly.size(); ++j) {
        auto q =
            triangle(m.cells[c].centre, rd.topology.vertices[poly[j]].point,
                     rd.topology.vertices[poly[(j + 1) % poly.size()]].point);
        mean.x += q.x;
        mean.y += q.y;
      }
      Vector2D refined{};
      for (size_t j = 0; j < poly.size(); ++j) {
        auto aa = m.cells[c].centre, bb = rd.topology.vertices[poly[j]].point,
             cc = rd.topology.vertices[poly[(j + 1) % poly.size()]].point;
        Point2D ab{.5 * (aa.x + bb.x), .5 * (aa.y + bb.y)},
            ac{.5 * (aa.x + cc.x), .5 * (aa.y + cc.y)},
            bcp{.5 * (bb.x + cc.x), .5 * (bb.y + cc.y)};
        for (auto q : {triangle(aa, ab, ac), triangle(ab, bb, bcp),
                       triangle(ac, bcp, cc), triangle(ab, bcp, ac)}) {
          refined.x += q.x;
          refined.y += q.y;
        }
      }
      referenceQuadratureGap = std::max(
          referenceQuadratureGap,
          std::hypot(mean.x - refined.x, mean.y - refined.y) / m.cells[c].area);
      mean = refined;
      mean.x /= m.cells[c].area;
      mean.y /= m.cells[c].area;
      meanError += m.cells[c].area *
                   std::hypot(solved[c].x - mean.x, solved[c].y - mean.y);
      referenceMeanOffset +=
          m.cells[c].area * std::hypot(u[c].x - mean.x, u[c].y - mean.y);
      balance = std::max(balance, std::hypot(sr[c][0], sr[c][1]));
      cells << c << ',' << m.cells[c].centre.x << ',' << m.cells[c].centre.y
            << ',' << m.cells[c].area << ',' << u[c].x << ',' << u[c].y << ','
            << solved[c].x << ',' << solved[c].y << ',' << mean.x << ','
            << mean.y << '\n';
    }
    double fixedTime = 0, traceResidual = 0, jump = 0, diss = 0, heat = 0,
           defect = 0, pertDiss = 0, pertHeat = 0, pertResidual = 0;
    size_t traceIts = 0;
    std::ofstream faces(prefix + ".faces.csv");
    faces << std::setprecision(17)
          << "state,face,neighbour,x,y,sx,sy,fluxX,fluxY,work,referenceX,"
             "referenceY,referenceTorque,referenceWork\n";
    auto diagnostics = [&](const HybridViscousResult2D &q,
                           const std::vector<Vector2D> &v, const char *label) {
      std::ofstream cfile(prefix + "." + label + ".shared-cell.csv"),
          ffile(prefix + "." + label + ".trace-equations.csv");
      cfile << std::setprecision(17)
            << "cell,ux,uy,residualX,residualY,residualWork,dissipation,"
               "sharedHeating\n";
      ffile << std::setprecision(17)
            << "face,traceX,traceY,ownerFluxX,ownerFluxY,neighbourFluxX,"
               "neighbourFluxY,jumpX,jumpY\n";
      std::vector<Vector2D> oq(m.faces.size()), nq(m.faces.size());
      for (size_t c = 0; c < v.size(); ++c) {
        HybridViscousCell2D local(m, c, mu);
        std::vector<Vector2D> trace;
        for (auto f : local.faces())
          trace.push_back(q.trace[f]);
        auto lr = local.evaluate(v[c], trace);
        for (size_t j = 0; j < local.faces().size(); ++j) {
          auto f = local.faces()[j];
          auto &target = m.faces[f].owner == c ? oq[f] : nq[f];
          target = {lr.outwardFlux[j][0], lr.outwardFlux[j][1]};
        }
        cfile << c << ',' << v[c].x << ',' << v[c].y << ','
              << q.cellResidual[c][0] << ',' << q.cellResidual[c][1] << ','
              << q.cellResidual[c][2] << ',' << q.cellDissipation[c] << ','
              << q.cellMechanicalHeating[c] << '\n';
      }
      for (size_t f = 0; f < m.faces.size(); ++f)
        ffile << f << ',' << q.trace[f].x << ',' << q.trace[f].y << ','
              << oq[f].x << ',' << oq[f].y << ',' << nq[f].x << ',' << nq[f].y
              << ',' << (m.faces[f].neighbour ? oq[f].x + nq[f].x : 0) << ','
              << (m.faces[f].neighbour ? oq[f].y + nq[f].y : 0) << '\n';
    };
    auto output = [&](const std::vector<Vector2D> &v, const char *label) {
      std::vector<std::array<double, 3>> flux;
      t = std::chrono::steady_clock::now();
      if (method == "hybrid") {
        auto q = h.evaluate(v, rho, 1e-13);
        flux = q.faceFlux;
        fixedTime = seconds(t);
        traceResidual = q.traceResidualNorm;
        jump = q.maximumTractionJump;
        diss = q.dissipation;
        heat = q.mechanicalHeating;
        defect = q.maximumCellWorkDefect;
        traceIts = q.iterations;
        diagnostics(q, v, label);
      } else {
        flux = l.evaluate(v, rho).faceFlux;
        fixedTime = seconds(t);
      }
      for (size_t f = 0; f < m.faces.size(); ++f) {
        auto p = m.faces[f].centre;
        auto s = m.faces[f].areaVector;
        Vector2D ref{};
        double torque = 0, refWork = 0;
        for (size_t k = 0; k < 8; ++k) {
          Point2D pt{p.x + .5 * gx[k] * s.y, p.y - .5 * gx[k] * s.x};
          auto q = referenceFlux(pt, s, mu);
          double w = .5 * gw[k];
          ref.x += w * q.x;
          ref.y += w * q.y;
          torque += w * (pt.x * q.y - pt.y * q.x);
          auto vel = exact(pt);
          refWork += w * (vel.x * q.x + vel.y * q.y);
        }
        faces << label << ',' << f << ','
              << (m.faces[f].neighbour ? std::to_string(*m.faces[f].neighbour)
                                       : "-1")
              << ',' << p.x << ',' << p.y << ',' << s.x << ',' << s.y << ','
              << flux[f][0] << ',' << flux[f][1] << ',' << flux[f][2] << ','
              << ref.x << ',' << ref.y << ',' << torque << ',' << refWork
              << '\n';
      }
    };
    output(u, "analytic-centroid");
    double analyticTime = fixedTime;
    output(solved, "steady");
    if (method == "hybrid") {
      std::vector<Vector2D> v(u.size());
      for (size_t c = 0; c < v.size(); ++c)
        v[c] = {std::sin(61 * m.cells[c].centre.x) *
                    std::cos(43 * m.cells[c].centre.y),
                std::cos(59 * m.cells[c].centre.x) *
                    std::sin(47 * m.cells[c].centre.y)};
      auto q = hz.evaluate(v, rho, 1e-13);
      pertDiss = q.dissipation;
      pertHeat = q.mechanicalHeating;
      pertResidual = q.traceResidualNorm;
      diagnostics(q, v, "homogeneous-perturbation");
    }
    std::ofstream report(prefix + ".json");
    report << std::setprecision(17) << "{\"cells\":" << u.size()
           << ",\"method\":\"" << method
           << "\",\"velocityCentroidRelativeL1\":" << error / norm
           << ",\"velocityMeanRelativeL1\":" << meanError / norm
           << ",\"centroidMeanReferenceOffsetRelativeL1\":"
           << referenceMeanOffset / norm
           << ",\"referenceQuadratureGapVelocity\":" << referenceQuadratureGap
           << ",\"steadySharedBalanceNPerM\":" << balance
           << ",\"assemblySecondsBothOperators\":" << assembly
           << ",\"matrixSeconds\":" << matrixSeconds
           << ",\"steadySolveSeconds\":" << steadySeconds
           << ",\"steadyIterations\":" << its
           << ",\"analyticTraceSeconds\":" << analyticTime
           << ",\"steadyTraceSeconds\":" << fixedTime
           << ",\"traceIterations\":" << traceIts
           << ",\"traceResidualNPerM\":" << traceResidual
           << ",\"tractionJumpNPerM\":" << jump
           << ",\"dissipationWPerM\":" << diss
           << ",\"sharedMechanicalHeatingWPerM\":" << heat
           << ",\"maximumCellWorkDefectWPerM\":" << defect
           << ",\"perturbationDissipationWPerM\":" << pertDiss
           << ",\"perturbationSharedHeatingWPerM\":" << pertHeat
           << ",\"perturbationTraceResidualNPerM\":" << pertResidual
           << ",\"fullSeconds\":" << seconds(start) << "}\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
