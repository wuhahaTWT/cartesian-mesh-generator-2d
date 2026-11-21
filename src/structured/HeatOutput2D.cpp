#include "cartmesh2d/structured/HeatFlow2D.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <utility>

namespace cartmesh2d::structured {
namespace {
std::ofstream output(const std::filesystem::path& path) {
    std::ofstream file(path);
    if (!file) throw std::runtime_error("Cannot open output: " + path.string());
    file.exceptions(std::ios::badbit | std::ios::failbit);
    file << std::setprecision(13);
    return file;
}

double wallDerivative(double value, double next, double wall, double s, double t) {
    return (value - wall) / s + s / t * ((value - wall) / s - (next - value) / (t - s));
}

double wallPressure(double value, double next, double s, double t) {
    return value - s * (next - value) / (t - s);
}

int wallIndex(const std::vector<double>& axis, double coordinate) {
    return static_cast<int>(std::lower_bound(axis.begin(), axis.end(), coordinate) - axis.begin());
}

enum Side { Front, Rear, Bottom, Top };
struct WallSample {
    double arc, x, y, length, pCell, pWall, nu, shearGradient;
};

std::array<std::vector<WallSample>, 4> squareWalls(const Grid& g, const Fields& f) {
    std::array<std::vector<WallSample>, 4> walls;
    const int left = wallIndex(g.x, 4.5), right = wallIndex(g.x, 5.5);
    const int bottom = wallIndex(g.y, 4.5), top = wallIndex(g.y, 5.5);
    const auto trace = [&](Side side, double arc, double x, double y, double length,
                           Grid::Index p, Grid::Index q, double s, double t) {
        walls[side].push_back({arc, x, y, length, f.p[p], wallPressure(f.p[p], f.p[q], s, t),
            -wallDerivative(f.temperature[p], f.temperature[q], 1., s, t),
            wallDerivative(f.u[p], f.u[q], 0., s, t)});
    };
    const double sf = 4.5 - g.xc(left - 1), tf = 4.5 - g.xc(left - 2);
    const double sr = g.xc(right) - 5.5, tr = g.xc(right + 1) - 5.5;
    for (int j = bottom; j < top; ++j) {
        const double y = g.yc(j), length = g.y[j + 1] - g.y[j];
        trace(Front, y >= 5 ? y - 5 : y - 1, 4.5, y, length,
              g.index(left - 1, j), g.index(left - 2, j), sf, tf);
        trace(Rear, 7 - y, 5.5, y, length, g.index(right, j), g.index(right + 1, j), sr, tr);
    }
    const double sb = 4.5 - g.yc(bottom - 1), tb = 4.5 - g.yc(bottom - 2);
    const double st = g.yc(top) - 5.5, tt = g.yc(top + 1) - 5.5;
    for (int i = left; i < right; ++i) {
        const double x = g.xc(i), length = g.x[i + 1] - g.x[i];
        trace(Bottom, 8 - x, x, 4.5, length,
              g.index(i, bottom - 1), g.index(i, bottom - 2), sb, tb);
        trace(Top, x - 4, x, 5.5, length, g.index(i, top), g.index(i, top + 1), st, tt);
    }
    return walls;
}
}

SquareMetrics squareMetrics(const Grid& g, const Fields& f, double viscosity) {
    SquareMetrics result;
    const auto walls = squareWalls(g, f);
    for (std::size_t i = 0; i < walls[Front].size(); ++i) {
        const auto& a = walls[Front][i];
        const auto& b = walls[Rear][i];
        result.cdPressure += 2 * a.length * (a.pCell - b.pCell);
        result.cdPressureExtrapolated += 2 * a.length * (a.pWall - b.pWall);
    }
    for (std::size_t i = 0; i < walls[Bottom].size(); ++i) {
        const auto& a = walls[Bottom][i];
        const auto& b = walls[Top][i];
        result.cdViscous += 2 * viscosity * a.length * (a.shearGradient + b.shearGradient);
    }
    for (std::size_t side = 0; side < walls.size(); ++side)
        for (const auto& face : walls[side]) result.nu[side] += face.length * face.nu;
    result.leadingCornerPressure = f.p[g.index(wallIndex(g.x, 4.5), wallIndex(g.y, 5.5))];
    const int right = wallIndex(g.x, 5.5), bottom = wallIndex(g.y, 4.5);
    int mid = bottom;
    while (g.yc(mid) < 5.) ++mid;
    const double weight = (5. - g.yc(mid - 1)) / (g.yc(mid) - g.yc(mid - 1));
    double x = 5.5, u = 0;
    for (int i = right; i < g.nx; ++i) {
        const double next = f.u[g.index(i, mid - 1)] * (1 - weight) + f.u[g.index(i, mid)] * weight;
        if (next >= 0) {
            if (u < 0) result.reattachment = x - u * (g.xc(i) - x) / (next - u) - 5.5;
            return result;
        }
        x = g.xc(i);
        u = next;
    }
    result.reattachment = -1; // reverse flow extends to the outlet
    return result;
}

void writeSquareSurface(const std::filesystem::path& path, const Grid& g, const Fields& f) {
    constexpr std::array sideNames{"front", "rear", "bottom", "top"};
    const auto walls = squareWalls(g, f);
    std::vector<std::pair<const WallSample*, const char*>> ordered;
    for (std::size_t side = 0; side < walls.size(); ++side)
        for (const auto& face : walls[side]) ordered.emplace_back(&face, sideNames[side]);
    std::sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) { return a.first->arc < b.first->arc; });
    auto out = output(path);
    out << "s_over_d,x_over_d,y_over_d,face,face_length,nu,cp,cp_cell\n";
    for (const auto& [sample, side] : ordered) {
        const auto& w = *sample;
        out << w.arc << ',' << w.x << ',' << w.y << ',' << side << ',' << w.length << ','
            << w.nu << ',' << 2 * w.pWall << ',' << 2 * w.pCell << '\n';
    }
}

void writeFields(const std::filesystem::path& directory, const Grid& g, const Fields& f) {
    std::filesystem::create_directories(directory);
    auto csv = output(directory / "field.csv");
    auto vtk = output(directory / "field.vtk");
    const std::size_t count = static_cast<std::size_t>(g.nx) * g.ny;
    csv << "i,j,x,y,dx,dy,fluid,u,v,p,theta,vorticity\n";
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            const auto k = g.index(i, j);
            csv << i << ',' << j << ',' << g.xc(i) << ',' << g.yc(j) << ',' << g.x[i + 1] - g.x[i]
                << ',' << g.y[j + 1] - g.y[j] << ',' << static_cast<int>(g.fluid[k]) << ',' << f.u[k]
                << ',' << f.v[k] << ',' << f.p[k] << ',' << f.temperature[k] << ',' << f.vorticity[k] << '\n';
        }
    vtk << "# vtk DataFile Version 3.0\nStructured SIMPLE fields\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS "
        << g.nx + 1 << ' ' << g.ny + 1 << " 1\nX_COORDINATES " << g.x.size() << " double\n";
    for (double x : g.x) vtk << x << '\n';
    vtk << "Y_COORDINATES " << g.y.size() << " double\n";
    for (double y : g.y) vtk << y << '\n';
    vtk << "Z_COORDINATES 1 double\n0\nCELL_DATA " << count << "\nVECTORS velocity double\n";
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            const auto k = g.index(i, j);
            vtk << f.u[k] << ' ' << f.v[k] << " 0\n";
        }
    const auto scalar = [&](const char* name, const auto& values) {
        vtk << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) vtk << static_cast<double>(values[g.index(i, j)]) << '\n';
    };
    scalar("p", f.p);
    scalar("theta", f.temperature);
    scalar("vorticity", f.vorticity);
    scalar("fluid", g.fluid);
}

void writeResult(const std::filesystem::path& directory, const Grid& g, const Parameters& p, const Result& r) {
    writeFields(directory, g, r.fields);
    if (!r.outletU.empty()) {
        auto outlet = output(directory / "outlet.csv");
        outlet << "y,u,v\n";
        for (int j = 0; j < g.ny; ++j)
            outlet << g.yc(j) << ',' << r.outletU[j] << ',' << r.fields.v[g.index(g.nx - 1, j)] << '\n';
    }
    auto history = output(directory / "residuals.csv");
    history << "phase,iteration,continuity,u,v,theta,global_mass,corrected_continuity\n";
    const auto rows = [&](const char* phase, const auto& values) {
        for (const auto& row : values)
            history << phase << ',' << row.iteration << ',' << row.continuity << ',' << row.u << ','
                << row.v << ',' << row.temperature << ',' << row.globalMass << ',' << row.correctedContinuity << '\n';
    };
    rows("flow", r.flowHistory);
    rows("thermal", r.thermalHistory);
    auto summary = output(directory / "summary.csv");
    summary << "quantity,value\nstatus," << r.stop << "\nnx," << g.nx << "\nny," << g.ny
        << "\nfluid_cells," << std::count(g.fluid.begin(), g.fluid.end(), 1) << "\nre," << p.re << "\npr," << p.pr
        << "\nflow_iterations," << r.flowHistory.size() << "\nflow_seconds," << r.flowSeconds
        << "\nthermal_iterations," << r.thermalHistory.size() << "\nthermal_seconds," << r.thermalSeconds
        << "\noutlet_scale," << r.outletScale << "\np_rel_tol," << p.pressureRelativeTolerance
        << "\ninterpolated_initial_field," << !p.initialField.empty()
        << "\nrhie_chow," << p.rhieChow << '\n';
    if (r.flowConverged && r.thermalConverged) {
        const auto m = squareMetrics(g, r.fields, 1. / p.re);
        writeSquareSurface(directory / "surface.csv", g, r.fields);
        summary << "cd," << m.cd() << "\ncd_extrapolated," << m.cdExtrapolated()
            << "\ncd_pressure," << m.cdPressure << "\ncd_pressure_extrapolated," << m.cdPressureExtrapolated
            << "\ncd_viscous," << m.cdViscous << "\nleading_corner_pressure," << m.leadingCornerPressure
            << "\nnu," << m.meanNu() << "\nnu_front," << m.nu[0] << "\nnu_rear," << m.nu[1]
            << "\nnu_bottom," << m.nu[2] << "\nnu_top," << m.nu[3] << "\nlr," << m.reattachment << '\n';
    }
}
}
