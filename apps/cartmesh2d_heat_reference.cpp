#include "cartmesh2d/fv/ThermalFlow2D.hpp"
#include "cartmesh2d/grid/RectilinearMesh2D.hpp"
#include "cartmesh2d/structured/HeatFlow2D.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

using namespace cartmesh2d;
using namespace cartmesh2d::structured;

namespace {
TopologyMesh2D toTopology(const Grid& g) {
    auto full = makeRectilinearMesh2D(g.x, g.y);
    constexpr auto none = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> cells(full.cells.size(), none);
    std::vector<std::size_t> faces(full.edges.size(), none);
    std::vector<std::size_t> vertices(full.vertices.size(), none);
    TopologyMesh2D result;
    for (auto& c : full.cells) {
        const int i = static_cast<int>(c.id) % g.nx, j = static_cast<int>(c.id) / g.nx;
        if (!g.fluid[g.index(i, j)]) continue;
        cells[c.id] = result.cells.size();
        c.id = result.cells.size();
        for (auto& v : c.vertices) {
            if (vertices[v] == none) {
                vertices[v] = result.vertices.size();
                result.vertices.push_back({vertices[v], full.vertices[v].point});
            }
            v = vertices[v];
        }
        result.cells.push_back(std::move(c));
    }
    for (auto edge : full.edges) {
        const auto owner = cells[edge.owner];
        const auto neighbour = edge.neighbour ? cells[*edge.neighbour] : none;
        if (owner == none && neighbour == none) continue;
        faces[edge.id] = result.edges.size();
        edge.id = result.edges.size();
        if (owner == none) std::swap(edge.v0, edge.v1);
        edge.v0 = vertices[edge.v0];
        edge.v1 = vertices[edge.v1];
        edge.owner = owner == none ? neighbour : owner;
        if (owner != none && neighbour != none) {
            edge.neighbour = neighbour;
        } else {
            if (edge.neighbour) edge.patch = BoundaryPatch2D::EmbeddedBoundary;
            edge.neighbour.reset();
        }
        result.edges.push_back(edge);
    }
    for (auto& c : result.cells)
        for (auto& e : c.edges) e = faces[e];
    return result;
}
}

int main(int argc, char** argv) {
    using Clock = std::chrono::steady_clock;
    if (argc != 3) {
        std::cerr << "Usage: cartmesh2d_heat_reference INPUT.cfg OUTPUT_DIRECTORY\n";
        return 1;
    }
    try {
        const auto p = readParameters(argv[1]);
        if (p.problem != Case::Square) throw std::invalid_argument("Reference expects the square cylinder case");
        const auto g = makeGrid(p);
        const auto topology = toTopology(g);
        const auto mesh = fv::makeFvMesh2D(topology);
        const auto scheme = p.secondOrder ? fv::ConvectionScheme2D::LimitedLinearUpwind : fv::ConvectionScheme2D::Upwind;
        fv::FlowControls2D fc;
        fc.scenario = "custom";
        fc.nu = 1 / p.re;
        fc.velocityRelaxation = p.alphaU;
        fc.pressureRelaxation = p.alphaP;
        fc.tolerance = p.tolerance;
        fc.maxIterations = static_cast<std::size_t>(p.maxIterations);
        fc.convection = scheme;
        fc.pressurePreconditioner = fv::PressurePreconditioner2D::Aggregation;
        fv::ThermalSetup2D thermal;
        thermal.diffusivity = fc.nu / p.pr;
        thermal.sourceDensity.resize(mesh.cells.size());
        thermal.boundary.resize(mesh.faces.size());
        for (std::size_t id = 0; id < mesh.faces.size(); ++id) {
            const auto& face = mesh.faces[id];
            if (face.neighbour) continue;
            const bool cylinder = face.patch == BoundaryPatch2D::EmbeddedBoundary;
            if (face.centre.x == 0) {
                fc.boundaryConditions.push_back({id, fv::FlowBoundaryKind2D::VelocityInlet, {1, 0}, 0, "inlet"});
                thermal.boundary[id] = {fv::ScalarBoundaryKind2D::Value, 0., {}};
            } else if (face.centre.x == g.x.back()) {
                fc.boundaryConditions.push_back({id, fv::FlowBoundaryKind2D::PressureOutlet, {}, 0, "outlet"});
                thermal.boundary[id] = {fv::ScalarBoundaryKind2D::DiffusiveFlux, 0., 0.};
            } else {
                fc.boundaryConditions.push_back({id, fv::FlowBoundaryKind2D::Wall, {}, 0, cylinder ? "cylinder" : "channel_walls"});
                thermal.boundary[id] = cylinder ? fv::ScalarBoundary2D{fv::ScalarBoundaryKind2D::Value, 1., {}}
                    : fv::ScalarBoundary2D{fv::ScalarBoundaryKind2D::DiffusiveFlux, 0., {}};
            }
        }
        fv::ScalarTransportControls2D sc;
        sc.convection = scheme;
        sc.preconditioner = fv::ScalarPreconditioner2D::ILU0;
        const std::filesystem::path output = argv[2];
        std::filesystem::create_directories(output);
        std::ofstream residuals(output / "residuals.csv");
        std::ofstream summary(output / "summary.csv");
        if (!residuals || !summary) throw std::runtime_error("Cannot open reference output files");
        summary << std::setprecision(13) << "quantity,value\nfluid_cells," << mesh.cells.size() << '\n';
        residuals << std::setprecision(13) << "phase,iteration,momentum,continuity,velocity_change,pressure_change,global_mass\n";
        const auto progress = [&](const char* phase) {
            return [&, phase](const fv::FlowIteration2D& r) {
                residuals << phase << ',' << r.iteration << ',' << r.momentumResidual << ',' << r.continuity
                    << ',' << r.velocityChange << ',' << r.pressureChange << ',' << r.globalRelativeImbalance << '\n';
                if (r.iteration == 1 || r.iteration % 100 == 0)
                    std::cout << phase << ' ' << r.iteration << ' ' << r.momentumResidual
                        << ' ' << r.continuity << std::endl;
            };
        };
        const auto start = Clock::now();
        auto flow = fv::solveIncompressible2D(mesh, fc, progress("steady"));
        summary << "flow_seconds," << std::chrono::duration<double>(Clock::now() - start).count() << '\n'
            << "flow_iterations," << flow.history.size() << '\n';
        if (!flow.converged) {
            summary << "status,Steady native flow did not converge\n";
            std::cerr << "Steady native flow did not converge: " << flow.failureReason << '\n';
            return 2;
        }
        auto scalarStart = Clock::now();
        fv::ScalarTransportProblem2D sp;
        sp.diffusivity = thermal.diffusivity;
        sp.sourceDensity = thermal.sourceDensity;
        sp.boundaryData = thermal.boundary;
        sp.volumeFlux = flow.flux;
        auto scalar = fv::solveScalarTransport2D(mesh, sp, sc);
        summary << "thermal_seconds," << std::chrono::duration<double>(Clock::now() - scalarStart).count() << '\n'
            << "thermal_iterations," << scalar.history.size() << '\n';
        std::ofstream scalarLog(output / "thermal-residuals.csv");
        scalarLog << std::setprecision(13) << "phase,iteration,relative_residual,cell_imbalance\n";
        for (const auto& r : scalar.history)
            scalarLog << "steady," << r.iteration << ',' << r.relativeResidual << ',' << r.maxDiagonalScaledImbalance << '\n';
        if (!scalar.converged) {
            summary << "status,Steady native scalar did not converge\n";
            std::cerr << "Steady native scalar did not converge\n";
            return 2;
        }
        // The existing thermal entry point advances physical time. Use the steady
        // native fields as its previous layer and measure the change after one step.
        fv::ThermalFlowState2D previous{{0, flow.u, flow.v, flow.p, flow.flux}, scalar.values};
        const auto stepStart = Clock::now();
        const auto step = fv::advanceThermalFlow2D(mesh, fc, thermal, sc, previous, 1., progress("thermal-step"));
        summary << "step_seconds," << std::chrono::duration<double>(Clock::now() - stepStart).count() << '\n'
            << "step_flow_iterations," << step.flow.history.size() << "\nstep_thermal_iterations," << step.scalar.history.size() << '\n';
        for (const auto& r : step.scalar.history)
            scalarLog << "thermal-step," << r.iteration << ',' << r.relativeResidual << ',' << r.maxDiagonalScaledImbalance << '\n';
        if (!step.accepted) {
            summary << "status,Native thermal step did not converge\n";
            std::cerr << "Native thermal step did not converge: " << step.flow.failureReason << '\n';
            return 2;
        }
        Fields fields;
        fields.u.resize(g.fluid.size());
        fields.v.resize(g.fluid.size());
        fields.p.resize(g.fluid.size());
        fields.temperature.resize(g.fluid.size());
        double du = 0, dt = 0;
        for (const auto& cell : topology.cells) {
            const auto id = cell.id;
            const int i = static_cast<int>(cell.sourceId) % g.nx, j = static_cast<int>(cell.sourceId) / g.nx;
            const int k = g.index(i, j);
            fields.u[k] = step.accepted->flow.u[id];
            fields.v[k] = step.accepted->flow.v[id];
            fields.p[k] = step.accepted->flow.p[id];
            fields.temperature[k] = step.accepted->scalar[id];
            du = std::max({du, std::abs(fields.u[k] - flow.u[id]), std::abs(fields.v[k] - flow.v[id])});
            dt = std::max(dt, std::abs(fields.temperature[k] - scalar.values[id]));
        }
        updateVorticity(g, p.problem, fields);
        writeFields(output, g, fields);
        const auto m = squareMetrics(g, fields, fc.nu);
        summary << "status,Converged\nstep_max_velocity_change," << du << "\nstep_max_theta_change," << dt
            << "\ncd," << m.cd() << "\nnu," << m.meanNu() << "\nlr," << m.reattachment << '\n';
        std::cout << "Native comparison converged: Cd=" << m.cd() << " Nu=" << m.meanNu()
            << " Lr/D=" << m.reattachment << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
