#include "cartmesh2d/fv/ReactingFlow2D.hpp"
#include "FvTestMesh2D.hpp"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;
using namespace cartmesh2d::chemistry;
namespace {
void jsonString(std::ostream& out, const std::string& value) {
    constexpr char hex[] = "0123456789abcdef";
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << static_cast<char>(c);
        else if (c < 32) out << "\\u00" << hex[c / 16] << hex[c % 16];
        else out << static_cast<char>(c);
    }
    out << '"';
}
template<class V> void numbers(std::ostream& out, const V& values) {
    out << '['; bool first = true; for (auto value : values) { if (!first) out << ','; first = false; out << value; } out << ']';
}
void state(std::ostream& out, DetailedGas& gas, const ReactingState2D& s) {
    out << "{\"time\":" << s.time << ",\"steps\":" << s.steps << ",\"U\":[";
    for (std::size_t i = 0; i < s.cells.size(); ++i) { if (i) out << ','; numbers(out, s.cells[i]); }
    out << "],\"temperature\":[";
    for (std::size_t i = 0; i < s.cells.size(); ++i) { if (i) out << ','; out << reactingPrimitive2D(gas, s.cells[i]).properties.temperature; }
    out << "],\"pressure\":[";
    for (std::size_t i = 0; i < s.cells.size(); ++i) { if (i) out << ','; out << reactingPrimitive2D(gas, s.cells[i]).properties.pressure; }
    out << "],\"velocity\":[";
    for (std::size_t i = 0; i < s.cells.size(); ++i) {
        if (i) out << ','; const auto u = reactingPrimitive2D(gas, s.cells[i]).velocity; out << '[' << u.x << ',' << u.y << ']';
    }
    out << "]}";
}
}
int main(int argc, char** argv) {
    try {
        if (argc < 3 || argc > 4) throw std::invalid_argument("expected mechanism path, NEW output directory, optional CFL");
        const std::filesystem::path directory(argv[2]);
        if (std::filesystem::exists(directory)) throw std::invalid_argument("output directory already exists");
        std::filesystem::create_directories(directory);
        const auto started = std::chrono::steady_clock::now();
        DetailedGas gas(argv[1]);
        auto mesh = fv_test::rectangle(8, 6, 1, true);
        const double size = .05;
        for (auto& c : mesh.cells) { c.centre.x *= size; c.centre.y *= size; c.area *= size * size; }
        for (auto& f : mesh.faces) { f.centre.x *= size; f.centre.y *= size; f.areaVector = f.areaVector * size; f.correction = f.correction * size; }
        std::vector<ReactingBoundary2D> boundaries;
        for (std::size_t f = 0; f < mesh.faces.size(); ++f) if (!mesh.faces[f].neighbour)
            boundaries.push_back({f, ReactingBoundaryKind2D::NoSlipWall, {}, {}, 0});
        ReactingFlowStepper2D solver(gas, mesh, boundaries);
        std::vector<double> x(gas.mechanism().species.size());
        for (const auto& [name, value] : {std::pair<std::string,double>{"H2", 2}, {"O2", 1}, {"N2", 3.76}}) {
            const auto it = std::find(gas.mechanism().species.begin(), gas.mechanism().species.end(), name);
            if (it == gas.mechanism().species.end()) throw std::runtime_error("specified mechanism lacks hot-spot fixture species");
            x[static_cast<std::size_t>(it - gas.mechanism().species.begin())] = value;
        }
        std::vector<ReactingConservative2D> cells;
        for (const auto& c : mesh.cells) {
            const double r2 = std::pow((c.centre.x - .022) / .012, 2) + std::pow((c.centre.y - .026) / .012, 2);
            cells.push_back(reactingConservative2D(gas.fromMoleAmounts(1000 + 250 * std::exp(-r2), 101325, x)));
        }
        const auto initial = solver.initialState(cells); auto current = initial;
        ReactingStepControls2D controls; controls.maximumStep = 2e-6; controls.endTime = 4e-5;
        if (argc == 4) controls.courant = std::stod(argv[3]);
        std::vector<double> boundaryImpulse(cells[0].size()), chemicalChange(cells[0].size());
        std::vector<ReactingState2D> frames{initial};
        std::size_t sourceCalls = 0, rejections = 0, fallbacks = 0; double minimumStep = controls.maximumStep;
        std::ofstream log(directory / "steps.jsonl"); log << std::setprecision(17);
        std::string failure;
        while (current.time < *controls.endTime && current.steps < 10000) {
            const auto step = solver.advance(current, controls);
            sourceCalls += step.sourceCalls; rejections += step.rejectedReasons.size(); fallbacks += step.hlleFallbacks;
            log << "{\"previousTime\":" << current.time << ",\"accepted\":" << (step.accepted ? "true" : "false")
                << ",\"dt\":" << step.step << ",\"courant\":" << step.combinedCourant << ",\"rejections\":" << step.rejectedReasons.size()
                << ",\"rejectedReasons\":[";
            for (std::size_t j = 0; j < step.rejectedReasons.size(); ++j) { if (j) log << ','; jsonString(log, step.rejectedReasons[j]); }
            log << "],\"failure\":"; jsonString(log, step.failure); log << "}\n";
            log.flush();
            if (!step.accepted) { failure = step.failure; std::cerr << failure << '\n'; break; }
            minimumStep = std::min(minimumStep, step.step);
            for (std::size_t k = 0; k < boundaryImpulse.size(); ++k) {
                boundaryImpulse[k] += step.step * step.boundaryFlux[k]; chemicalChange[k] += step.chemistryChange[k];
            }
            current = *step.accepted;
            if (current.steps % 10 == 0) {
                frames.push_back(current);
                std::cerr << "accepted " << current.steps << " t=" << current.time << " s\n";
            }
        }
        if (frames.back().time != current.time) frames.push_back(current);
        const bool complete = current.time == *controls.endTime;
        if (!complete && failure.empty()) failure = "accepted step budget exhausted before target time";
        if (current.steps > 0) { std::ofstream checkpoint(directory / "accepted.checkpoint"); solver.writeCheckpoint(checkpoint, current); }
        std::ofstream definition(directory / "resolved-mechanism.yaml"); definition << gas.mechanism().resolvedDefinition;
        std::ofstream output(directory / "field.json"); output << std::setprecision(17);
        output << "{\"complete\":" << (complete ? "true" : "false") << ",\"qualifiedFlame\":false,\"targetTime\":" << *controls.endTime
               << ",\"elapsedSeconds\":" << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()
               << ",\"sourceCalls\":" << sourceCalls << ",\"rejections\":" << rejections << ",\"hlleFallbacks\":" << fallbacks << ",\"minimumStep\":" << minimumStep
               << ",\"controls\":{\"courant\":" << controls.courant << ",\"predictionSafety\":" << controls.predictionSafety
               << ",\"maximumStep\":" << controls.maximumStep << ",\"minimumStep\":" << controls.minimumStep
               << ",\"order\":" << controls.order << ",\"relativeChemistryTolerance\":" << controls.chemistry.relativeTolerance
               << ",\"absoluteChemistryTolerance\":" << controls.chemistry.absoluteTolerance << "}"
               << ",\"failure\":"; jsonString(output, failure); output << ",\"species\":[";
        for (std::size_t k = 0; k < gas.mechanism().species.size(); ++k) { if (k) output << ','; output << std::quoted(gas.mechanism().species[k]); }
        output << "],\"molecularWeights\":"; numbers(output, gas.mechanism().molecularWeights);
        output << ",\"atomicWeights\":"; numbers(output, gas.mechanism().atomicWeights);
        output << ",\"atomCounts\":"; numbers(output, gas.mechanism().atomCounts);
        output << ",\"boundaryImpulse\":"; numbers(output, boundaryImpulse);
        output << ",\"chemistryChange\":"; numbers(output, chemicalChange);
        output << ",\"mesh\":{\"cells\":[";
        for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
            if (i) output << ','; const auto& c = mesh.cells[i];
            output << "{\"centre\":[" << c.centre.x << ',' << c.centre.y << "],\"area\":" << c.area << ",\"faces\":"; numbers(output, c.faces); output << '}';
        }
        output << "],\"faces\":[";
        for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
            if (i) output << ','; const auto& f = mesh.faces[i];
            output << "{\"owner\":" << f.owner << ",\"neighbour\":"; if (f.neighbour) output << *f.neighbour; else output << "null";
            output << ",\"centre\":[" << f.centre.x << ',' << f.centre.y << "],\"S\":[" << f.areaVector.x << ',' << f.areaVector.y << "]}";
        }
        output << "]},\"frames\":[";
        for (std::size_t i = 0; i < frames.size(); ++i) { if (i) output << ','; state(output, gas, frames[i]); }
        output << "]}\n";
        if (!output || !log) throw std::runtime_error("cannot write verification field/log");
        return complete ? 0 : 1;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
