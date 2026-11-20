#include "cartmesh2d/structured/HeatFlow2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cartmesh2d::structured {
namespace {
using Sample = std::array<double, 4>;

struct Samples {
    std::vector<double> x, y;
    std::vector<Sample> values;
    std::vector<bool> fluid;
};

Samples readSamples(const std::filesystem::path& path, const Grid& grid) {
    std::ifstream file(path);
    std::string line;
    std::getline(file, line);
    if (line.ends_with('\r')) line.pop_back();
    if (line != "i,j,x,y,dx,dy,fluid,u,v,p,theta,vorticity")
        throw std::invalid_argument("init_from expects a structured field.csv: " + path.string());
    Samples data;
    double minX{}, maxX{}, minY{}, maxY{};
    for (std::size_t row = 0; std::getline(file, line); ++row) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream entry(line);
        std::size_t i, j;
        int fluid;
        double x, y, dx, dy, omega;
        Sample sample;
        if (!(entry >> i >> j >> x >> y >> dx >> dy >> fluid >>
              sample[0] >> sample[1] >> sample[2] >> sample[3] >> omega) ||
            !std::isfinite(x + y + dx + dy + sample[0] + sample[1] + sample[2] + sample[3]) ||
            dx <= 0 || dy <= 0 || (fluid != 0 && fluid != 1))
            throw std::invalid_argument("Invalid row in initial field");
        if (!j) {
            if (i != row || data.y.size() > 1 || (i && x <= data.x.back()))
                throw std::invalid_argument("Invalid initial field x coordinates");
            data.x.push_back(x);
            if (!i) {
                minX = x - .5 * dx;
                minY = y - .5 * dy;
            }
            maxX = x + .5 * dx;
        }
        if (!i) {
            if (j != data.y.size() || (j && y <= data.y.back()))
                throw std::invalid_argument("Invalid initial field y coordinates");
            data.y.push_back(y);
            maxY = y + .5 * dy;
        }
        if (i >= data.x.size() || j >= data.y.size() || row != j * data.x.size() + i ||
            x != data.x[i] || y != data.y[j])
            throw std::invalid_argument("Initial field cells must form a rectilinear grid");
        data.values.push_back(sample);
        data.fluid.push_back(fluid == 1);
    }
    const auto matches = [](double a, double b) {
        // field.csv stores coordinates and widths with 13 significant digits.
        return std::abs(a - b) <= 1e-11 * std::max({1., std::abs(a), std::abs(b)});
    };
    if (data.x.size() < 2 || data.y.size() < 2 || data.values.size() != data.x.size() * data.y.size() ||
        !matches(minX, grid.x.front()) || !matches(maxX, grid.x.back()) ||
        !matches(minY, grid.y.front()) || !matches(maxY, grid.y.back()))
        throw std::invalid_argument("Initial field grid does not cover the computational domain");
    return data;
}

struct Bracket {
    std::size_t lower;
    double weight;
};

Bracket bracket(const std::vector<double>& centres, double x) {
    const auto upper = std::lower_bound(centres.begin(), centres.end(), x);
    const auto k = std::clamp(static_cast<std::size_t>(upper - centres.begin()), std::size_t{1}, centres.size() - 1);
    return {k - 1, std::clamp((x - centres[k - 1]) / (centres[k] - centres[k - 1]), 0., 1.)};
}
}

void interpolateInitialField(const std::filesystem::path& path, const Grid& grid, Fields& f) {
    const auto input = readSamples(path, grid);
    std::vector<Bracket> x, y;
    for (int i = 0; i < grid.nx; ++i) x.push_back(bracket(input.x, grid.xc(i)));
    for (int j = 0; j < grid.ny; ++j) y.push_back(bracket(input.y, grid.yc(j)));
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            const auto k = grid.index(i, j);
            if (!grid.fluid[k]) continue;
            const auto [x0, tx] = x[i];
            const auto [y0, ty] = y[j];
            Sample value{};
            double sum = 0;
            for (std::size_t b = 0; b < 2; ++b)
                for (std::size_t a = 0; a < 2; ++a) {
                    const auto sample = (y0 + b) * input.x.size() + x0 + a;
                    if (!input.fluid[sample]) continue;
                    const double weight = (a ? tx : 1 - tx) * (b ? ty : 1 - ty);
                    sum += weight;
                    for (std::size_t component = 0; component < value.size(); ++component)
                        value[component] += weight * input.values[sample][component];
                }
            if (sum == 0) continue;
            f.u[k] = value[0] / sum;
            f.v[k] = value[1] / sum;
            f.p[k] = value[2] / sum;
            f.temperature[k] = value[3] / sum;
        }
    }
}
}
