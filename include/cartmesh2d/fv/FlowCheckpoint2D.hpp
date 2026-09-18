#pragma once

#include "cartmesh2d/fv/Incompressible2D.hpp"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cartmesh2d::fv {

namespace flow_checkpoint_detail {

inline void fail(const std::string& message) {
    throw std::runtime_error("Flow checkpoint: " + message);
}

inline void token(std::istream& in, const char* expected) {
    std::string actual;
    if (!(in >> actual) || actual != expected) fail("expected '" + std::string(expected) + "'");
}

inline std::size_t count(std::istream& in, std::size_t expected, const char* what) {
    unsigned long long value = 0;
    if (!(in >> value) || value != static_cast<unsigned long long>(expected))
        fail(std::string("invalid ") + what + " count");
    return expected;
}

inline void finite(double value, const char* what) {
    if (!std::isfinite(value)) fail(std::string("nonfinite ") + what);
}

inline void exact(double actual, double expected, const char* what) {
    finite(actual, what);
    if (actual != expected) fail(std::string("mesh mismatch in ") + what);
}

inline void exact(std::size_t actual, std::size_t expected, const char* what) {
    if (actual != expected) fail(std::string("mesh mismatch in ") + what);
}

inline const char* convectionName(ConvectionScheme2D value) {
    switch (value) {
    case ConvectionScheme2D::Upwind: return "upwind";
    case ConvectionScheme2D::LimitedLinearUpwind: return "limited-linear";
    }
    fail("invalid convection scheme");
    return "";
}

inline const char* stressName(ViscousStress2D value) {
    switch (value) {
    case ViscousStress2D::Laplacian: return "laplacian";
    case ViscousStress2D::Symmetric: return "symmetric";
    }
    fail("invalid viscous stress scheme");
    return "";
}

inline void configuration(std::istream& in, const FlowControls2D& controls) {
    std::string scenario, convection, stress;
    double nu = 0, speed = 0, slope = 0;
    if (!(in >> std::quoted(scenario) >> nu >> speed >> convection >> stress >> slope))
        fail("truncated configuration");
    finite(nu, "configuration nu");
    finite(speed, "configuration speed");
    finite(slope, "configuration pressure slope");
    if (scenario != controls.scenario || nu != controls.nu || speed != controls.speed ||
        convection != convectionName(controls.convection) || stress != stressName(controls.viscousStress) ||
        slope != controls.manufacturedPressureSlope)
        fail("incompatible configuration");
}

inline void writeDouble(std::ostream& out, double value, const char* what) {
    finite(value, what);
    out << std::setprecision(17) << value;
}

} // namespace flow_checkpoint_detail

inline void writeFlowCheckpoint2D(std::ostream& out, const FvMesh2D& mesh,
                                  const FlowControls2D& controls, const FlowState2D& state) {
    struct RestoreFormat {
        std::ostream& stream;
        std::ios_base::fmtflags flags;
        std::streamsize precision;
        ~RestoreFormat() { stream.flags(flags); stream.precision(precision); }
    } restore{out,out.flags(),out.precision()};
    out << std::defaultfloat << std::dec << std::noshowpos << std::noshowbase << std::setprecision(17);
    validateFvMesh2D(mesh);
    flow_checkpoint_detail::finite(controls.nu, "configuration nu");
    flow_checkpoint_detail::finite(controls.speed, "configuration speed");
    flow_checkpoint_detail::finite(controls.manufacturedPressureSlope, "configuration pressure slope");
    flow_checkpoint_detail::finite(state.time, "time");
    if (state.time < 0) flow_checkpoint_detail::fail("time must be nonnegative");
    if (state.u.size() != mesh.cells.size() || state.v.size() != mesh.cells.size() ||
        state.p.size() != mesh.cells.size() || state.flux.size() != mesh.faces.size())
        flow_checkpoint_detail::fail("invalid state vector size");
    for (const auto* values : {&state.u, &state.v, &state.p, &state.flux})
        for (const double value : *values) flow_checkpoint_detail::finite(value, "state value");

    out << "CARTMESH2D_FLOW_CHECKPOINT 1\n"
        << "DISCRETIZATION Euler-RC-v2\n"
        << "CONFIG " << std::quoted(controls.scenario) << ' ';
    flow_checkpoint_detail::writeDouble(out, controls.nu, "configuration nu"); out << ' ';
    flow_checkpoint_detail::writeDouble(out, controls.speed, "configuration speed"); out << ' '
        << flow_checkpoint_detail::convectionName(controls.convection) << ' '
        << flow_checkpoint_detail::stressName(controls.viscousStress) << ' ';
    flow_checkpoint_detail::writeDouble(out, controls.manufacturedPressureSlope, "configuration pressure slope");
    out << "\nCELLS " << mesh.cells.size() << '\n';
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        const auto& cell = mesh.cells[i];
        out << "CELL " << i << ' ';
        flow_checkpoint_detail::writeDouble(out, cell.centre.x, "cell centre"); out << ' ';
        flow_checkpoint_detail::writeDouble(out, cell.centre.y, "cell centre"); out << ' ';
        flow_checkpoint_detail::writeDouble(out, cell.area, "cell area"); out << ' ' << cell.faces.size();
        for (const auto id : cell.faces) out << ' ' << id;
        out << '\n';
    }
    out << "FACES " << mesh.faces.size() << '\n';
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
        const auto& face = mesh.faces[i];
        out << "FACE " << i << ' ' << face.owner << ' ';
        if (face.neighbour) out << *face.neighbour; else out << '-';
        out << ' ' << static_cast<int>(face.patch) << ' ';
        flow_checkpoint_detail::writeDouble(out, face.centre.x, "face centre"); out << ' ';
        flow_checkpoint_detail::writeDouble(out, face.centre.y, "face centre"); out << ' ';
        flow_checkpoint_detail::writeDouble(out, face.areaVector.x, "face area vector"); out << ' ';
        flow_checkpoint_detail::writeDouble(out, face.areaVector.y, "face area vector"); out << ' ';
        flow_checkpoint_detail::writeDouble(out, face.correction.x, "face correction"); out << ' ';
        flow_checkpoint_detail::writeDouble(out, face.correction.y, "face correction"); out << ' ';
        flow_checkpoint_detail::writeDouble(out, face.transmissibility, "face transmissibility"); out << ' ';
        flow_checkpoint_detail::writeDouble(out, face.neighbourWeight, "face neighbour weight"); out << '\n';
    }
    out << "TIME "; flow_checkpoint_detail::writeDouble(out, state.time, "time"); out << '\n';
    const auto writeVector = [&out](const char* name, const std::vector<double>& values) {
        out << name << ' ' << values.size();
        for (const double value : values) { out << ' '; flow_checkpoint_detail::writeDouble(out, value, "state value"); }
        out << '\n';
    };
    writeVector("U", state.u); writeVector("V", state.v); writeVector("P", state.p); writeVector("FLUX", state.flux);
    out << "END\n";
    if (!out) flow_checkpoint_detail::fail("write failed");
}

inline FlowState2D readFlowCheckpoint2D(std::istream& in, const FvMesh2D& mesh,
                                        const FlowControls2D& controls) {
    validateFvMesh2D(mesh);
    flow_checkpoint_detail::token(in, "CARTMESH2D_FLOW_CHECKPOINT");
    flow_checkpoint_detail::token(in, "1");
    flow_checkpoint_detail::token(in, "DISCRETIZATION");
    flow_checkpoint_detail::token(in, "Euler-RC-v2");
    flow_checkpoint_detail::token(in, "CONFIG");
    flow_checkpoint_detail::configuration(in, controls);
    flow_checkpoint_detail::token(in, "CELLS");
    flow_checkpoint_detail::count(in, mesh.cells.size(), "cell");
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        flow_checkpoint_detail::token(in, "CELL");
        std::size_t id = 0, faceCount = 0; double x = 0, y = 0, area = 0;
        if (!(in >> id >> x >> y >> area >> faceCount)) flow_checkpoint_detail::fail("truncated cell");
        flow_checkpoint_detail::exact(id, i, "cell ordering");
        flow_checkpoint_detail::exact(x, mesh.cells[i].centre.x, "cell centre x");
        flow_checkpoint_detail::exact(y, mesh.cells[i].centre.y, "cell centre y");
        flow_checkpoint_detail::exact(area, mesh.cells[i].area, "cell area");
        flow_checkpoint_detail::exact(faceCount, mesh.cells[i].faces.size(), "cell face count");
        for (std::size_t j = 0; j < faceCount; ++j) {
            std::size_t face = 0;
            if (!(in >> face)) flow_checkpoint_detail::fail("truncated cell face list");
            flow_checkpoint_detail::exact(face, mesh.cells[i].faces[j], "cell face ordering");
        }
    }
    flow_checkpoint_detail::token(in, "FACES");
    flow_checkpoint_detail::count(in, mesh.faces.size(), "face");
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
        flow_checkpoint_detail::token(in, "FACE");
        std::size_t id = 0, owner = 0; std::string neighbour; int patch = 0;
        double cx = 0, cy = 0, ax = 0, ay = 0, kx = 0, ky = 0, trans = 0, weight = 0;
        if (!(in >> id >> owner >> neighbour >> patch >> cx >> cy >> ax >> ay >> kx >> ky >> trans >> weight))
            flow_checkpoint_detail::fail("truncated face");
        const auto& expected = mesh.faces[i];
        flow_checkpoint_detail::exact(id, i, "face ordering");
        flow_checkpoint_detail::exact(owner, expected.owner, "face owner");
        if (expected.neighbour) {
            std::size_t value = 0; std::istringstream parsed(neighbour);
            if (!(parsed >> value) || !parsed.eof()) flow_checkpoint_detail::fail("invalid face neighbour");
            flow_checkpoint_detail::exact(value, *expected.neighbour, "face neighbour");
        } else if (neighbour != "-") flow_checkpoint_detail::fail("boundary face has neighbour");
        if (patch != static_cast<int>(expected.patch)) flow_checkpoint_detail::fail("mesh mismatch in face patch");
        flow_checkpoint_detail::exact(cx, expected.centre.x, "face centre x");
        flow_checkpoint_detail::exact(cy, expected.centre.y, "face centre y");
        flow_checkpoint_detail::exact(ax, expected.areaVector.x, "face area vector x");
        flow_checkpoint_detail::exact(ay, expected.areaVector.y, "face area vector y");
        flow_checkpoint_detail::exact(kx, expected.correction.x, "face correction x");
        flow_checkpoint_detail::exact(ky, expected.correction.y, "face correction y");
        flow_checkpoint_detail::exact(trans, expected.transmissibility, "face transmissibility");
        flow_checkpoint_detail::exact(weight, expected.neighbourWeight, "face neighbour weight");
    }
    flow_checkpoint_detail::token(in, "TIME");
    FlowState2D state;
    if (!(in >> state.time)) flow_checkpoint_detail::fail("truncated time");
    flow_checkpoint_detail::finite(state.time, "time");
    if (state.time < 0) flow_checkpoint_detail::fail("time must be nonnegative");
    const auto readVector = [&in](const char* name, std::size_t expected, std::vector<double>& values) {
        flow_checkpoint_detail::token(in, name);
        flow_checkpoint_detail::count(in, expected, name);
        values.resize(expected);
        for (double& value : values) {
            if (!(in >> value)) flow_checkpoint_detail::fail(std::string("truncated ") + name);
            flow_checkpoint_detail::finite(value, name);
        }
    };
    readVector("U", mesh.cells.size(), state.u); readVector("V", mesh.cells.size(), state.v);
    readVector("P", mesh.cells.size(), state.p); readVector("FLUX", mesh.faces.size(), state.flux);
    flow_checkpoint_detail::token(in, "END");
    std::string trailing;
    if (in >> trailing) flow_checkpoint_detail::fail("trailing data");
    if (in.bad()) flow_checkpoint_detail::fail("read failed");
    return state;
}

} // namespace cartmesh2d::fv
