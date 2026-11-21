#include "cartmesh2d/structured/HeatFlow2D.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cartmesh2d::structured {
Fields readFieldCsv(const std::filesystem::path& path, const Grid& grid, bool interpolate) {
    std::ifstream file(path);
    std::string line;
    std::getline(file, line);
    if (line.ends_with('\r')) line.pop_back();
    if (line != "i,j,x,y,dx,dy,fluid,u,v,p,theta,vorticity") throw std::invalid_argument("Invalid field.csv header");
    using Sample = std::array<double, 4>;
    std::vector<Sample> samples;
    std::vector<double> x, y;
    std::vector<bool> fluid;
    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        int i, j, mask;
        double xc, yc, dx, dy;
        Sample s;
        if (!(row >> i >> j >> xc >> yc >> dx >> dy >> mask >> s[0] >> s[1] >> s[2] >> s[3]))
            throw std::invalid_argument("Invalid field.csv data");
        if (j == 0) x.push_back(xc);
        if (i == 0) y.push_back(yc);
        samples.push_back(s);
        fluid.push_back(mask != 0);
    }
    if (samples.size() != x.size() * y.size() || x.size() < 2 || y.size() < 2 ||
        (!interpolate && (x.size() != static_cast<std::size_t>(grid.nx) || y.size() != static_cast<std::size_t>(grid.ny))))
        throw std::invalid_argument("Field.csv grid dimensions do not match the operation");
    const auto bracket = [](const auto& centres, double value) {
        const auto hi = std::clamp(std::lower_bound(centres.begin(), centres.end(), value), centres.begin() + 1, centres.end() - 1);
        return std::pair{static_cast<std::size_t>(hi - centres.begin() - 1),
                         std::clamp((value - hi[-1]) / (*hi - hi[-1]), 0., 1.)};
    };
    const std::vector<double> zero(grid.fluid.size());
    Fields f{zero, zero, zero, zero, {}};
    const std::array components{&f.u, &f.v, &f.p, &f.temperature};
    for (int j = 0; j < grid.ny; ++j)
        for (int i = 0; i < grid.nx; ++i) {
            const auto k = grid.index(i, j);
            if (!grid.fluid[k]) continue;
            const auto [ix, tx] = interpolate ? bracket(x, grid.xc(i)) : std::pair{static_cast<std::size_t>(i), 0.};
            const auto [iy, ty] = interpolate ? bracket(y, grid.yc(j)) : std::pair{static_cast<std::size_t>(j), 0.};
            double total = 0.;
            for (std::size_t b = 0; b < 2; ++b)
                for (std::size_t a = 0; a < 2; ++a) {
                    const double w = (a ? tx : 1 - tx) * (b ? ty : 1 - ty);
                    if (w == 0.) continue;
                    const auto n = (iy + b) * x.size() + ix + a;
                    if (!fluid[n]) continue;
                    total += w;
                    for (std::size_t c = 0; c < components.size(); ++c) (*components[c])[k] += w * samples[n][c];
                }
            if (total > 0.)
                for (auto* values : components) (*values)[k] /= total;
        }
    return f;
}
}
