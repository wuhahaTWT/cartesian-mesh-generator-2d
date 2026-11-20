#pragma once

#include <array>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace cartmesh2d::structured {
enum class Case { Square, Cavity, Channel };

struct Parameters {
    Case problem = Case::Square;
    int n = 10, baseN = 0, nx = 240, ny = 40, maxIterations = 12000;
    double stretch = 1.1, re = 20, pr = .71, alphaU = .7, alphaP = .3;
    double tolerance = 1e-6;
    bool secondOrder = true;
};

struct Grid {
    std::vector<double> x, y;
    std::vector<unsigned char> fluid;
    int nx{}, ny{}, stride{};
    int index(int i, int j) const { return (j + 1) * stride + i + 1; }
    double xc(int i) const { return .5 * (x[i] + x[i + 1]); }
    double yc(int j) const { return .5 * (y[j] + y[j + 1]); }
};

struct Fields {
    std::vector<double> u, v, p, temperature, vorticity;
};

struct Residual {
    int iteration{};
    double continuity{}, u{}, v{}, temperature{}, globalMass{};
};

struct Result {
    Fields fields;
    std::vector<Residual> flowHistory, thermalHistory;
    std::vector<double> outletU;
    bool flowConverged{}, thermalConverged{};
    double flowSeconds{}, thermalSeconds{}, outletScale = 1;
    std::string stop;
};

struct SquareMetrics {
    double cdPressure{}, cdViscous{}, reattachment{};
    std::array<double, 4> nu{}; // front, rear, bottom, top
    double cd() const { return cdPressure + cdViscous; }
    double meanNu() const { return .25 * (nu[0] + nu[1] + nu[2] + nu[3]); }
};

Parameters readParameters(const std::filesystem::path&);
Grid makeGrid(const Parameters&);
Result solve(const Grid&, const Parameters&, const std::function<void(bool, const Residual&)>& = {});
SquareMetrics squareMetrics(const Grid&, const Fields&, double viscosity);
void updateVorticity(const Grid&, Case, Fields&);
void writeFields(const std::filesystem::path&, const Grid&, const Fields&);
void writeResult(const std::filesystem::path&, const Grid&, const Parameters&, const Result&);
}
