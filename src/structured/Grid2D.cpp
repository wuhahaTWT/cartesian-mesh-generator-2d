#include "cartmesh2d/structured/HeatFlow2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace cartmesh2d::structured {
namespace {
void validate(const Parameters& p) {
    const auto positive = [](double v) { return std::isfinite(v) && v > 0; };
    if (!positive(p.re) || !positive(p.pr) || !positive(p.tolerance) || p.tolerance > 1e-6 ||
        !(p.alphaU > 0 && p.alphaU <= 1 && p.alphaP > 0 && p.alphaP <= 1) ||
        !(p.stretch >= 1 && p.stretch <= 1.1) || p.maxIterations < 1 || p.n < 2 ||
        p.baseN < 0 || p.baseN == 1 || (p.baseN && p.n % p.baseN) ||
        p.nx < 4 || p.ny < 4)
        throw std::invalid_argument("Invalid grid, fluid properties or SIMPLE controls");
    const auto cells = p.problem == Case::Square
        ? (20LL * p.n + 2) * static_cast<long double>(10LL * p.n + 2)
        : (p.nx + 2LL) * static_cast<long double>(p.ny + 2LL);
    if (cells > std::numeric_limits<int>::max())
        throw std::invalid_argument("Grid exceeds supported indexing range");
}

std::vector<double> outward(double length, double h, double ratio) {
    int count = 1;
    double width = h, covered = h;
    while (covered < length) {
        width *= ratio;
        covered += width;
        ++count;
    }
    const auto lengthAt = [=](double r) {
        double sum = 0, w = h;
        for (int i = 0; i < count; ++i, w *= r) sum += w;
        return sum;
    };
    // Fit the segment while matching the first cell at the square.
    double low = 1, high = ratio;
    for (int i = 0; i < 48; ++i) {
        const double mid = .5 * (low + high);
        if (lengthAt(mid) < length) low = mid;
        else high = mid;
    }
    const double fitted = .5 * (low + high);
    width = h * length / lengthAt(fitted);
    if (width * ratio < h)
        throw std::invalid_argument("Stretch ratio too small to join the square and outer grid");
    std::vector<double> offsets{0};
    for (int i = 0; i < count; ++i, width *= fitted) offsets.push_back(offsets.back() + width);
    offsets.back() = length;
    return offsets;
}

std::vector<double> squareAxis(double length, int n, int base, double ratio) {
    if (ratio == 1) {
        if (base % 2) throw std::invalid_argument("Unit stretch requires an even base_n to align the half-unit walls");
        std::vector<double> axis;
        for (int k = 0; k <= static_cast<int>(length) * n; ++k)
            axis.push_back(static_cast<double>(k) / n);
        return axis;
    }
    const double h = 1. / base;
    const auto left = outward(4.5, h, ratio);
    const auto right = outward(length - 5.5, h, ratio);
    std::vector<double> coarse;
    for (auto it = left.rbegin(); it != left.rend(); ++it) coarse.push_back(4.5 - *it);
    for (int k = 1; k <= base; ++k) coarse.push_back(4.5 + static_cast<double>(k) / base);
    for (std::size_t k = 1; k < right.size(); ++k) coarse.push_back(5.5 + right[k]);
    std::vector<double> axis{0};
    const int split = n / base;
    for (std::size_t k = 1; k < coarse.size(); ++k)
        for (int sub = 1; sub <= split; ++sub)
            axis.push_back(std::lerp(coarse[k - 1], coarse[k], static_cast<double>(sub) / split));
    return axis;
}

std::vector<double> uniform(double length, int count) {
    std::vector<double> axis;
    for (int k = 0; k <= count; ++k) axis.push_back(length * k / count);
    return axis;
}
}

Parameters readParameters(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot open input file");
    Parameters p;
    std::string line;
    while (std::getline(file, line)) {
        std::replace(line.begin(), line.end(), '=', ' ');
        std::istringstream entry(line.substr(0, line.find('#')));
        std::string key, value, extra;
        if (!(entry >> key)) continue;
        if (!(entry >> value) || (entry >> extra)) throw std::invalid_argument("Invalid entry: " + key);
        std::istringstream number(value);
        if (key == "case") {
            if (value == "square") p.problem = Case::Square;
            else if (value == "cavity") p.problem = Case::Cavity;
            else if (value == "channel") p.problem = Case::Channel;
            else throw std::invalid_argument("Unknown case: " + value);
            continue;
        }
        if (key == "convection") {
            if (value == "central") p.secondOrder = true;
            else if (value == "upwind") p.secondOrder = false;
            else throw std::invalid_argument("Unknown convection scheme: " + value);
            continue;
        }
        if (key == "n") number >> p.n;
        else if (key == "base_n") number >> p.baseN;
        else if (key == "nx") number >> p.nx;
        else if (key == "ny") number >> p.ny;
        else if (key == "stretch") number >> p.stretch;
        else if (key == "re") number >> p.re;
        else if (key == "pr") number >> p.pr;
        else if (key == "alpha_u") number >> p.alphaU;
        else if (key == "alpha_p") number >> p.alphaP;
        else if (key == "tolerance") number >> p.tolerance;
        else if (key == "max_iterations") number >> p.maxIterations;
        else throw std::invalid_argument("Unknown parameter: " + key);
        if (!number || (number >> extra)) throw std::invalid_argument("Invalid value for " + key);
    }
    return p;
}

Grid makeGrid(const Parameters& p) {
    validate(p);
    Grid g;
    if (p.problem == Case::Square) {
        const int base = p.baseN ? p.baseN : p.n;
        g.x = squareAxis(20, p.n, base, p.stretch);
        g.y = squareAxis(10, p.n, base, p.stretch);
    } else {
        g.x = uniform(p.problem == Case::Channel ? 12. : 1., p.nx);
        g.y = uniform(1., p.ny);
    }
    g.nx = static_cast<int>(g.x.size()) - 1;
    g.ny = static_cast<int>(g.y.size()) - 1;
    g.stride = g.nx + 2;
    g.fluid.resize(static_cast<std::size_t>(g.stride) * (g.ny + 2));
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i)
            g.fluid[g.index(i, j)] = !(p.problem == Case::Square &&
                g.xc(i) > 4.5 && g.xc(i) < 5.5 && g.yc(j) > 4.5 && g.yc(j) < 5.5);
    return g;
}
}
