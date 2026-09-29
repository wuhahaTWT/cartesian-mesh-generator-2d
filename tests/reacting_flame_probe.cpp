#include "cartmesh2d/fv/ReactingFlow2D.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::chemistry;
namespace {
void require(bool value, const char* text) { if (!value) throw std::runtime_error(text); }
template<class V> void numbers(std::ostream& out, const V& v) {
    out << '['; bool first = true; for (auto x : v) { if (!first) out << ','; first = false; out << x; } out << ']';
}
template<class V> void matrix(std::ostream& out, const V& v) {
    out << '['; bool first = true; for (const auto& x : v) { if (!first) out << ','; first = false; numbers(out, x); } out << ']';
}
void string(std::ostream& out, const std::string& value) {
    constexpr char hex[] = "0123456789abcdef"; out << '"';
    for (char character : value) {
        const auto c = static_cast<unsigned char>(character);
        if (c == '"' || c == '\\') out << '\\' << static_cast<char>(c);
        else if (c < 32) out << "\\u00" << hex[c / 16] << hex[c % 16];
        else out << static_cast<char>(c);
    }
    out << '"';
}
FvMesh2D strip(const std::vector<double>& x, std::size_t ny, double height) {
    require(x.size() > 2 && ny > 0 && height > 0 && std::isfinite(height), "invalid strip geometry");
    TopologyMesh2D t; const auto nx = x.size() - 1;
    for (std::size_t j = 0; j <= ny; ++j) for (double a : x)
        t.vertices.push_back({t.vertices.size(), {a, height * static_cast<double>(j) / static_cast<double>(ny)}});
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> edges;
    for (std::size_t j = 0; j < ny; ++j) for (std::size_t i = 0; i < nx; ++i) {
        require(std::isfinite(x[i]) && std::isfinite(x[i + 1]) && x[i + 1] > x[i], "non-increasing strip grid");
        TopologyCell2D c; c.id = t.cells.size(); c.geometryArea = (x[i + 1] - x[i]) * height / static_cast<double>(ny);
        const auto a = j * (nx + 1) + i; c.vertices = {a, a + 1, a + nx + 2, a + nx + 1};
        for (std::size_t k = 0; k < 4; ++k) {
            const auto u = c.vertices[k], v = c.vertices[(k + 1) % 4];
            const auto [it, inserted] = edges.emplace(std::minmax(u, v), t.edges.size());
            if (inserted) t.edges.push_back({it->second, u, v, c.id, {}, BoundaryPatch2D::DomainBoundary});
            else { t.edges[it->second].neighbour = c.id; t.edges[it->second].patch = BoundaryPatch2D::None; }
            c.edges.push_back(it->second);
        }
        t.cells.push_back(c);
    }
    return makeFvMesh2D(t);
}
void state(std::ostream& out, DetailedGas& gas, const ReactingState2D& s, const ReactingResidual2D& r) {
    out << "{\"time\":" << s.time << ",\"steps\":" << s.steps << ",\"U\":"; matrix(out, s.cells);
    out << ",\"temperature\":[";
    for (std::size_t i = 0; i < s.cells.size(); ++i) { if (i) out << ','; out << reactingPrimitive2D(gas, s.cells[i]).properties.temperature; }
    out << "],\"pressure\":[";
    for (std::size_t i = 0; i < s.cells.size(); ++i) { if (i) out << ','; out << reactingPrimitive2D(gas, s.cells[i]).properties.pressure; }
    out << "],\"transportDerivative\":"; matrix(out, r.transportDerivative);
    out << ",\"chemistryDerivative\":"; matrix(out, r.chemistryDerivative);
    out << ",\"derivative\":"; matrix(out, r.derivative);
    out << ",\"faceFlux\":"; matrix(out, r.faceFlux);
    out << ",\"transportRate\":"; numbers(out, r.transportRate);
    out << ",\"boundaryFlux\":"; numbers(out, r.boundaryFlux);
    out << ",\"chemistryIntegral\":"; numbers(out, r.chemistryIntegral);
    out << ",\"hlleFallbacks\":" << r.hlleFallbacks << '}';
}
}
int main(int argc, char** argv) {
    try {
        const bool regression = argc == 4 && std::string(argv[3]) == "--regression";
        require(argc == 5 || regression, "expected mechanism, fixture, NEW output directory, physical duration (0 for residual only), or mechanism fixture --regression");
        const auto start = std::chrono::steady_clock::now();
        DetailedGas gas(argv[1]); std::ifstream input(argv[2]);
        std::string token; unsigned version = 0; std::size_t nx = 0, ny = 0, ns = 0; double height = 0;
        require(static_cast<bool>(input >> token >> version) && token == "CM2D_FLAME_FIXTURE" && version == 1, "invalid fixture version");
        require(static_cast<bool>(input >> nx >> ny >> ns >> height) && nx > 1 && nx <= 20000 && ny > 0 && ny <= 8
                && ns == gas.mechanism().species.size(), "invalid fixture sizes");
        for (const auto& species : gas.mechanism().species)
            require(static_cast<bool>(input >> token) && token == species, "fixture species order differs from mechanism");
        std::vector<double> x(nx + 1); for (double& v : x) require(static_cast<bool>(input >> v), "truncated fixture grid");
        std::vector<ReactingConservative2D> columns(nx, ReactingConservative2D(ns + 4));
        for (auto& u : columns) for (double& v : u) require(static_cast<bool>(input >> v), "truncated fixture cells");
        ReactingConservative2D inlet(ns + 4); for (double& v : inlet) require(static_cast<bool>(input >> v), "truncated inlet state");
        require(static_cast<bool>(input >> token) && token == "END" && !(input >> token), "invalid fixture ending");
        const auto mesh = strip(x, ny, height);
        std::vector<ReactingBoundary2D> bc; const auto in = reactingPrimitive2D(gas, inlet);
        for (std::size_t id = 0; id < mesh.faces.size(); ++id) if (!mesh.faces[id].neighbour) {
            ReactingBoundary2D b; b.face = id; const auto& f = mesh.faces[id];
            if (f.areaVector.x < 0) { b.kind = ReactingBoundaryKind2D::Reservoir; b.reservoir = in.gas; b.velocity = in.velocity; }
            else if (f.areaVector.x > 0) b.kind = ReactingBoundaryKind2D::ExtrapolatedOutflow;
            bc.push_back(b);
        }
        ReactingFlowStepper2D solver(gas, mesh, bc);
        std::vector<ReactingConservative2D> cells;
        for (std::size_t j = 0; j < ny; ++j) cells.insert(cells.end(), columns.begin(), columns.end());
        const auto initial = solver.initialState(cells); auto current = initial;
        const auto initialResidual = solver.evaluateResidual(initial);
        if (regression) {
            ReactingStepControls2D controls; controls.endTime = 1e-8;
            while (current.time < *controls.endTime && current.steps < 128) {
                const auto r = solver.advance(current, controls);
                if (!r.accepted) throw std::runtime_error(r.failure);
                require(r.rejectedReasons.empty(), "trace flame fixture requires rejected transport candidates");
                for (std::size_t k = 0; k < ns + 4; ++k) {
                    double scale = std::abs(r.beforeIntegral[k]) + std::abs(r.afterIntegral[k]) + std::abs(r.chemistryChange[k]);
                    for (const auto& f : r.faceFlux) scale += r.step * std::abs(f[k]);
                    require(std::abs(r.balanceError[k]) <= 3e-11 * scale, "trace regression conservation failed");
                }
                for (double v : r.elementalBalanceError) require(std::abs(v) < 2e-8, "trace regression elemental drift");
                current = *r.accepted;
            }
            require(current.time == *controls.endTime, "trace regression did not reach physical endpoint");
            require(initial.cells == cells && initial.steps == 0 && initial.time == 0, "trace regression modified initial state");
            (void)solver.evaluateResidual(current);
            std::cout << "Trace flame regression: " << mesh.cells.size() << " cells, " << current.steps << " accepted steps, no rejections\n";
            return 0;
        }
        std::size_t used = 0; const std::string durationText(argv[4]); const double duration = std::stod(durationText, &used);
        require(used == durationText.size() && duration >= 0 && std::isfinite(duration), "invalid physical duration");
        const std::filesystem::path directory(argv[3]); require(!std::filesystem::exists(directory), "output directory already exists");
        std::filesystem::create_directories(directory);
        { std::ofstream definition(directory / "resolved-mechanism.yaml"); definition << gas.mechanism().resolvedDefinition;
          require(static_cast<bool>(definition), "cannot write resolved mechanism"); }
        std::ofstream log(directory / "steps.jsonl"); log << std::setprecision(17);
        ReactingStepControls2D controls; controls.endTime = duration;
        std::vector<double> boundaryImpulse(ns + 4), chemistryChange(ns + 4), transportClosure(ns + 4), chemistryClosure(ns + 4);
        double maximumClosure = 0, absoluteClosure = 0;
        std::string failure; std::size_t rejected = 0, sourceCalls = 0;
        while (current.time < duration && current.steps < 20000) {
            const auto r = solver.advance(current, controls); sourceCalls += r.sourceCalls; rejected += r.rejectedReasons.size();
            log << "{\"time\":" << current.time << ",\"dt\":" << r.step << ",\"accepted\":" << (r.accepted ? "true" : "false")
                << ",\"courant\":" << r.combinedCourant << ",\"rejectedReasons\":[";
            for (std::size_t j = 0; j < r.rejectedReasons.size(); ++j) { if (j) log << ','; string(log, r.rejectedReasons[j]); }
            log << "],\"failure\":"; string(log, r.failure);
            log << ",\"maximumMassClosureFraction\":" << r.maximumMassClosureFraction
                << ",\"absoluteMassClosureIntegral\":" << r.absoluteMassClosureIntegral
                << ",\"transportMassClosureChange\":"; numbers(log, r.transportMassClosureChange);
            log << ",\"chemistryMassClosureChange\":"; numbers(log, r.chemistryMassClosureChange);
            log << "}\n"; log.flush();
            if (!r.accepted) { failure = r.failure; break; }
            for (std::size_t k = 0; k < ns + 4; ++k) {
                boundaryImpulse[k] += r.step * r.boundaryFlux[k]; chemistryChange[k] += r.chemistryChange[k];
                transportClosure[k] += r.transportMassClosureChange[k]; chemistryClosure[k] += r.chemistryMassClosureChange[k];
            }
            maximumClosure = std::max(maximumClosure, r.maximumMassClosureFraction);
            absoluteClosure += r.absoluteMassClosureIntegral;
            current = *r.accepted;
            if (current.steps % 50 == 0) std::cerr << "accepted " << current.steps << " t=" << current.time << '\n';
        }
        if (current.time < duration && failure.empty()) failure = "step budget exhausted before physical endpoint";
        if (current.steps > 0) { std::ofstream checkpoint(directory / "accepted.checkpoint"); solver.writeCheckpoint(checkpoint, current); }
        // Preserve raw final accepted values before optional residual
        // diagnostics, which can themselves expose an interpolation failure.
        { std::ofstream saved(directory / "accepted-state.json"); saved << std::setprecision(17)
              << "{\"time\":" << current.time << ",\"steps\":" << current.steps << ",\"U\":";
          matrix(saved, current.cells); saved << ",\"failure\":"; string(saved, failure); saved << "}\n";
          require(static_cast<bool>(saved), "cannot save final accepted state"); }
        const auto finalResidual = solver.evaluateResidual(current);
        std::ofstream out(directory / "field.json"); out << std::setprecision(17);
        out << "{\"complete\":" << (current.time == duration ? "true" : "false") << ",\"qualifiedFlame\":false,\"duration\":" << duration
            << ",\"cells\":" << mesh.cells.size() << ",\"nx\":" << nx << ",\"ny\":" << ny << ",\"height\":" << height << ",\"xEdges\":";
        numbers(out, x); out << ",\"areas\":[";
        for (std::size_t i = 0; i < mesh.cells.size(); ++i) { if (i) out << ','; out << mesh.cells[i].area; }
        out << "],\"faces\":[";
        for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
            if (i) out << ','; const auto& f = mesh.faces[i]; out << "{\"owner\":" << f.owner << ",\"neighbour\":";
            if (f.neighbour) out << *f.neighbour; else out << "null";
            out << ",\"centre\":[" << f.centre.x << ',' << f.centre.y << "],\"S\":[" << f.areaVector.x << ',' << f.areaVector.y << "]}";
        }
        out << "],\"boundaryImpulse\":"; numbers(out, boundaryImpulse); out << ",\"chemistryChange\":"; numbers(out, chemistryChange);
        out << ",\"massClosure\":{\"maximumFraction\":" << maximumClosure << ",\"absoluteIntegral\":" << absoluteClosure
            << ",\"transport\":"; numbers(out, transportClosure); out << ",\"chemistry\":"; numbers(out, chemistryClosure); out << '}';
        out << ",\"initial\":"; state(out, gas, initial, initialResidual);
        out << ",\"final\":"; state(out, gas, current, finalResidual);
        out << ",\"rejections\":" << rejected << ",\"sourceCalls\":" << sourceCalls << ",\"failure\":"; string(out, failure);
        out << ",\"elapsedSecondsBeforeFinalFlush\":" << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << "}\n";
        out.close(); require(static_cast<bool>(out) && static_cast<bool>(log), "verification output write failed");
        if (!failure.empty()) std::cerr << failure << '\n';
        return current.time == duration ? 0 : 1;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
