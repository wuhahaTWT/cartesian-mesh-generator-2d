#include "cartmesh2d/fv/ManufacturedFlow2D.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <string>

using cartmesh2d::Point2D;
using cartmesh2d::Vector2D;
using cartmesh2d::fv::ManufacturedFlowSample2D;
using cartmesh2d::fv::manufacturedFlow2D;

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

void near(double actual, double expected, double tolerance, const std::string& message) {
    check(std::abs(actual - expected) <= tolerance, message);
}

double firstX(const Point2D point, double speed, double nu, double h, bool velocity,
              double pressureSlope = 0.0) {
    const auto plus = manufacturedFlow2D({point.x + h, point.y}, speed, nu, pressureSlope);
    const auto minus = manufacturedFlow2D({point.x - h, point.y}, speed, nu, pressureSlope);
    return velocity ? (plus.velocity.x - minus.velocity.x) / (2.0 * h)
                    : (plus.pressure - minus.pressure) / (2.0 * h);
}

double firstY(const Point2D point, double speed, double nu, double h, bool velocity,
              double pressureSlope = 0.0) {
    const auto plus = manufacturedFlow2D({point.x, point.y + h}, speed, nu, pressureSlope);
    const auto minus = manufacturedFlow2D({point.x, point.y - h}, speed, nu, pressureSlope);
    return velocity ? (plus.velocity.x - minus.velocity.x) / (2.0 * h)
                    : (plus.pressure - minus.pressure) / (2.0 * h);
}

double firstXComponent(const Point2D point, double speed, double nu, double h,
                       bool velocity, bool yComponent, double pressureSlope = 0.0) {
    const auto plus = manufacturedFlow2D({point.x + h, point.y}, speed, nu, pressureSlope);
    const auto minus = manufacturedFlow2D({point.x - h, point.y}, speed, nu, pressureSlope);
    if (!velocity) return (plus.pressure - minus.pressure) / (2.0 * h);
    return (yComponent ? plus.velocity.y - minus.velocity.y
                       : plus.velocity.x - minus.velocity.x) / (2.0 * h);
}

double firstYComponent(const Point2D point, double speed, double nu, double h,
                       bool velocity, bool yComponent, double pressureSlope = 0.0) {
    const auto plus = manufacturedFlow2D({point.x, point.y + h}, speed, nu, pressureSlope);
    const auto minus = manufacturedFlow2D({point.x, point.y - h}, speed, nu, pressureSlope);
    if (!velocity) return (plus.pressure - minus.pressure) / (2.0 * h);
    return (yComponent ? plus.velocity.y - minus.velocity.y
                       : plus.velocity.x - minus.velocity.x) / (2.0 * h);
}

double laplacian(const Point2D point, double speed, double nu, double h, bool yComponent,
                 double pressureSlope = 0.0) {
    const auto centre = manufacturedFlow2D(point, speed, nu, pressureSlope);
    const auto xp = manufacturedFlow2D({point.x + h, point.y}, speed, nu, pressureSlope);
    const auto xm = manufacturedFlow2D({point.x - h, point.y}, speed, nu, pressureSlope);
    const auto xpp = manufacturedFlow2D({point.x + 2.0 * h, point.y}, speed, nu, pressureSlope);
    const auto xmm = manufacturedFlow2D({point.x - 2.0 * h, point.y}, speed, nu, pressureSlope);
    const auto yp = manufacturedFlow2D({point.x, point.y + h}, speed, nu, pressureSlope);
    const auto ym = manufacturedFlow2D({point.x, point.y - h}, speed, nu, pressureSlope);
    const auto ypp = manufacturedFlow2D({point.x, point.y + 2.0 * h}, speed, nu, pressureSlope);
    const auto ymm = manufacturedFlow2D({point.x, point.y - 2.0 * h}, speed, nu, pressureSlope);
    const auto component = [yComponent](const ManufacturedFlowSample2D& sample) {
        return yComponent ? sample.velocity.y : sample.velocity.x;
    };
    const double xPart = (-component(xpp) + 16.0 * component(xp) - 30.0 * component(centre)
                          + 16.0 * component(xm) - component(xmm)) / (12.0 * h * h);
    const double yPart = (-component(ypp) + 16.0 * component(yp) - 30.0 * component(centre)
                          + 16.0 * component(ym) - component(ymm)) / (12.0 * h * h);
    return xPart + yPart;
}

struct Difference {
    double maxRelative = 0.0;
    double rmsRelative = 0.0;
};

Difference finiteDifferenceSource(double h) {
    const std::array<Point2D, 4> points{{{.17, .29}, {.31, .57}, {.63, .78}, {.42, .66}}};
    const std::array<double, 3> speeds{{.3, 1.0, 2.0}};
    const std::array<double, 2> viscosities{{.01, .1}};
    const std::array<double, 3> pressureSlopes{{0.0, 1.0, -.7}};
    double sum = 0.0;
    double maximum = 0.0;
    std::size_t count = 0;
    for (const auto speed : speeds) {
        for (const auto nu : viscosities) {
            for (const auto pressureSlope : pressureSlopes) {
                for (const auto point : points) {
                    const auto sample = manufacturedFlow2D(point, speed, nu, pressureSlope);
                    const double ux = firstXComponent(point, speed, nu, h, true, false, pressureSlope);
                    const double uy = firstYComponent(point, speed, nu, h, true, false, pressureSlope);
                    const double vx = firstXComponent(point, speed, nu, h, true, true, pressureSlope);
                    const double vy = firstYComponent(point, speed, nu, h, true, true, pressureSlope);
                    const double dpdx = firstX(point, speed, nu, h, false, pressureSlope);
                    const double dpdy = firstY(point, speed, nu, h, false, pressureSlope);
                    const double lapU = laplacian(point, speed, nu, h, false, pressureSlope);
                    const double lapV = laplacian(point, speed, nu, h, true, pressureSlope);
                    const double expectedX = sample.velocity.x * ux + sample.velocity.y * uy
                        + dpdx - nu * lapU;
                    const double expectedY = sample.velocity.x * vx + sample.velocity.y * vy
                        + dpdy - nu * lapV;
                    const double scale = std::max({1.0, speed * speed,
                                                   std::abs(expectedX), std::abs(expectedY)});
                    const double error = std::hypot(sample.acceleration.x - expectedX,
                                                    sample.acceleration.y - expectedY) / scale;
                    maximum = std::max(maximum, error);
                    sum += error * error;
                    ++count;
                }
            }
        }
    }
    return {maximum, std::sqrt(sum / static_cast<double>(count))};
}

void knownValuesAndBoundaries() {
    constexpr double pi = std::numbers::pi;
    const auto sample = manufacturedFlow2D({.2, .3}, 1.0, .01);
    near(sample.velocity.x, .328581945074459, 2e-14, "known manufactured u value");
    near(sample.velocity.y, -.622474571220695, 2e-14, "known manufactured v value");
    near(sample.pressure, .475528258147577, 2e-14, "known manufactured pressure value");
    near(manufacturedFlow2D({.2, .3}, 1.0, .01, 1.0).pressure,
         .975528258147577, 2e-14, "positive pressure slope value");
    near(manufacturedFlow2D({.2, .3}, 1.0, .01, -.7).pressure,
         .125528258147577, 2e-14, "negative pressure slope value");
    const std::array<Point2D, 4> walls{{{0.0, .37}, {1.0, .37}, {.37, 0.0}, {.37, 1.0}}};
    for (const auto point : walls) {
        for (const double pressureSlope : {0.0, 1.0, -.7}) {
            const auto wall = manufacturedFlow2D(point, 2.0, .1, pressureSlope);
            near(wall.velocity.x, 0.0, 2e-14, "manufactured boundary u is zero");
            near(wall.velocity.y, 0.0, 2e-14, "manufactured boundary v is zero");
        }
    }
    const double h = 1e-5;
    for (const double speed : {.3, 1.0, 2.0}) {
        for (const double pressureSlope : {0.0, 1.0, -.7}) {
            const double expected = speed * speed * pressureSlope;
            const auto px0 = (manufacturedFlow2D({-h, .37}, speed, .01, pressureSlope).pressure
                             - manufacturedFlow2D({h, .37}, speed, .01, pressureSlope).pressure) / (-2.0 * h);
            const auto px1 = (manufacturedFlow2D({1.0 + h, .37}, speed, .01, pressureSlope).pressure
                             - manufacturedFlow2D({1.0 - h, .37}, speed, .01, pressureSlope).pressure) / (2.0 * h);
            const auto py0 = (manufacturedFlow2D({.37, -h}, speed, .01, pressureSlope).pressure
                             - manufacturedFlow2D({.37, h}, speed, .01, pressureSlope).pressure) / (-2.0 * h);
            const auto py1 = (manufacturedFlow2D({.37, 1.0 + h}, speed, .01, pressureSlope).pressure
                             - manufacturedFlow2D({.37, 1.0 - h}, speed, .01, pressureSlope).pressure) / (2.0 * h);
            near(px0, expected, 2e-8, "x=0 pressure coordinate-normal derivative matches U^2 slope");
            near(px1, expected, 2e-8, "x=1 pressure coordinate-normal derivative matches U^2 slope");
            near(py0, expected, 2e-8, "y=0 pressure coordinate-normal derivative matches U^2 slope");
            near(py1, expected, 2e-8, "y=1 pressure coordinate-normal derivative matches U^2 slope");
        }
    }
    check(std::abs(pi - 3.141592653589793) < 1e-15, "standard pi is available for the fixture");
}

void divergenceAndSpeedChecks() {
    const std::array<Point2D, 4> points{{{.17, .29}, {.31, .57}, {.63, .78}, {.42, .66}}};
    const double h = 5e-5;
    double maximumDivergence = 0.0;
    double maximumSpeed = 0.0;
    for (const auto point : points) {
        const double ux = firstXComponent(point, 1.0, .01, h, true, false);
        const double vy = firstYComponent(point, 1.0, .01, h, true, true);
        maximumDivergence = std::max(maximumDivergence, std::abs(ux + vy));
        const auto sample = manufacturedFlow2D(point, 1.0, .01);
        maximumSpeed = std::max(maximumSpeed, std::hypot(sample.velocity.x, sample.velocity.y));
    }
    std::cout << "manufactured flow: max divergence=" << maximumDivergence
              << ", sampled max speed=" << maximumSpeed << '\n';
    check(maximumDivergence < 2e-8, "manufactured velocity is discretely divergence-free");
    check(maximumSpeed <= 1.0000000001, "speed=1 has unit-bounded manufactured velocity");
}

} // namespace

int main() {
    knownValuesAndBoundaries();
    divergenceAndSpeedChecks();
    const auto coarse = finiteDifferenceSource(1e-4);
    const auto fine = finiteDifferenceSource(5e-5);
    std::cout << "manufactured source: h=1e-4 max=" << coarse.maxRelative
              << ", rms=" << coarse.rmsRelative << "; h=5e-5 max=" << fine.maxRelative
              << ", rms=" << fine.rmsRelative << '\n';
    check(coarse.maxRelative < 2e-5 && fine.maxRelative < 2e-5,
          "manufactured acceleration agrees with independent finite differences");
    check(fine.rmsRelative < coarse.rmsRelative,
          "finite-difference source error improves when h is halved");
    if (failures != 0) std::cerr << "manufactured flow failures=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}
