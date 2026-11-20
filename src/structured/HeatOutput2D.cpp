#include "cartmesh2d/structured/HeatFlow2D.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>

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
}

SquareMetrics squareMetrics(const Grid& g, const Fields& f, double viscosity) {
    SquareMetrics result;
    const int left = static_cast<int>(std::lower_bound(g.x.begin(), g.x.end(), 4.5) - g.x.begin());
    const int right = static_cast<int>(std::lower_bound(g.x.begin(), g.x.end(), 5.5) - g.x.begin());
    const int bottom = static_cast<int>(std::lower_bound(g.y.begin(), g.y.end(), 4.5) - g.y.begin());
    const int top = static_cast<int>(std::lower_bound(g.y.begin(), g.y.end(), 5.5) - g.y.begin());
    for (int j = bottom; j < top; ++j) {
        const double area = g.y[j + 1] - g.y[j];
        result.cdPressure += 2 * area * (f.p[g.index(left - 1, j)] - f.p[g.index(right, j)]);
        result.nu[0] -= area * wallDerivative(f.temperature[g.index(left - 1, j)],
            f.temperature[g.index(left - 2, j)], 1., 4.5 - g.xc(left - 1), 4.5 - g.xc(left - 2));
        result.nu[1] -= area * wallDerivative(f.temperature[g.index(right, j)],
            f.temperature[g.index(right + 1, j)], 1., g.xc(right) - 5.5, g.xc(right + 1) - 5.5);
    }
    for (int i = left; i < right; ++i) {
        const double area = g.x[i + 1] - g.x[i];
        const double sb = 4.5 - g.yc(bottom - 1), tb = 4.5 - g.yc(bottom - 2);
        const double st = g.yc(top) - 5.5, tt = g.yc(top + 1) - 5.5;
        result.cdViscous += 2 * viscosity * area *
            (wallDerivative(f.u[g.index(i, bottom - 1)], f.u[g.index(i, bottom - 2)], 0., sb, tb) +
             wallDerivative(f.u[g.index(i, top)], f.u[g.index(i, top + 1)], 0., st, tt));
        result.nu[2] -= area * wallDerivative(f.temperature[g.index(i, bottom - 1)],
            f.temperature[g.index(i, bottom - 2)], 1., sb, tb);
        result.nu[3] -= area * wallDerivative(f.temperature[g.index(i, top)],
            f.temperature[g.index(i, top + 1)], 1., st, tt);
    }
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

void writeFields(const std::filesystem::path& directory, const Grid& g, const Fields& f) {
    std::filesystem::create_directories(directory);
    auto csv = output(directory / "field.csv");
    auto vtk = output(directory / "field.vtk");
    const std::size_t count = static_cast<std::size_t>(g.nx) * g.ny;
    csv << "i,j,x,y,dx,dy,fluid,u,v,p,theta,vorticity\n";
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            const int k = g.index(i, j);
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
            const int k = g.index(i, j);
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
    history << "phase,iteration,continuity,u,v,theta,global_mass\n";
    const auto rows = [&](const char* phase, const auto& values) {
        for (const auto& row : values)
            history << phase << ',' << row.iteration << ',' << row.continuity << ',' << row.u << ','
                << row.v << ',' << row.temperature << ',' << row.globalMass << '\n';
    };
    rows("flow", r.flowHistory);
    rows("thermal", r.thermalHistory);
    auto summary = output(directory / "summary.csv");
    summary << "quantity,value\nstatus," << r.stop << "\nnx," << g.nx << "\nny," << g.ny
        << "\nfluid_cells," << std::count(g.fluid.begin(), g.fluid.end(), 1) << "\nre," << p.re << "\npr," << p.pr
        << "\nflow_iterations," << r.flowHistory.size() << "\nflow_seconds," << r.flowSeconds
        << "\nthermal_iterations," << r.thermalHistory.size() << "\nthermal_seconds," << r.thermalSeconds
        << "\noutlet_scale," << r.outletScale << '\n';
    if (r.flowConverged && r.thermalConverged) {
        const auto m = squareMetrics(g, r.fields, 1. / p.re);
        summary << "cd," << m.cd() << "\ncd_pressure," << m.cdPressure << "\ncd_viscous," << m.cdViscous
            << "\nnu," << m.meanNu() << "\nnu_front," << m.nu[0] << "\nnu_rear," << m.nu[1]
            << "\nnu_bottom," << m.nu[2] << "\nnu_top," << m.nu[3] << "\nlr," << m.reattachment << '\n';
    }
}
}
