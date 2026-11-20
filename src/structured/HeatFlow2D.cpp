#include "cartmesh2d/structured/HeatFlow2D.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <utility>

namespace cartmesh2d::structured {
namespace {
using Vector = std::vector<double>;
using Index = Grid::Index;
using Clock = std::chrono::steady_clock;
enum Direction { West, East, South, North };
enum Boundary { Internal, Inlet, Outlet, Wall, Lid, Cylinder };
constexpr std::array opposite{East, West, North, South};

struct Face {
    Index p, q;
    Direction direction;
    Boundary boundary;
    double area, distance, weight, flux{};
};

struct Geometry {
    const Grid& g;
    std::array<Index, 4> offset;
    std::vector<Index> cells;
    Vector volume;
    std::vector<Face> internal, boundary;
    double inflow{};

    Geometry(const Grid& grid, Case problem)
        : g(grid), offset{-1, 1, -g.stride, g.stride}, volume(g.fluid.size()) {
        for (int j = 0; j < g.ny; ++j) {
            for (int i = 0; i < g.nx; ++i) {
                const auto p = g.index(i, j);
                if (!g.fluid[p]) continue;
                cells.push_back(p);
                const double dx = g.x[i + 1] - g.x[i], dy = g.y[j + 1] - g.y[j];
                volume[p] = dx * dy;
                for (Direction d : {West, East, South, North}) {
                    const auto q = p + offset[d];
                    const bool xFace = d <= East;
                    const double area = xFace ? dy : dx;
                    const double half = .5 * (xFace ? dx : dy);
                    if (g.fluid[q]) {
                        if (d != East && d != North) continue;
                        const double distance = d == East ? g.xc(i + 1) - g.xc(i) : g.yc(j + 1) - g.yc(j);
                        internal.push_back({p, q, d, Internal, area, distance, half / distance});
                        continue;
                    }
                    Boundary b = Cylinder;
                    if (i == 0 && d == West) b = problem == Case::Cavity ? Wall : Inlet;
                    if (i == g.nx - 1 && d == East) b = problem == Case::Cavity ? Wall : Outlet;
                    if (j == 0 && d == South) b = Wall;
                    if (j == g.ny - 1 && d == North) b = problem == Case::Cavity ? Lid : Wall;
                    boundary.push_back({p, q, d, b, area, half, 0, b == Inlet ? -area : 0.});
                    if (b == Inlet) inflow += area;
                }
            }
        }
    }
};

struct Coefficients {
    double diagonal{};
    std::array<double, 4> neighbour{};
};
using Matrix = std::vector<Coefficients>;

void clear(Vector& v) { std::fill(v.begin(), v.end(), 0.); }

double product(const Coefficients& a, const Vector& x, Index p, Index stride) {
    return a.diagonal * x[p] - a.neighbour[West] * x[p - 1] - a.neighbour[East] * x[p + 1]
        - a.neighbour[South] * x[p - stride] - a.neighbour[North] * x[p + stride];
}

double wallValue(Boundary b, int component) {
    if (component == 0 && (b == Inlet || b == Lid)) return 1;
    if (component == 2 && b == Cylinder) return 1;
    return 0;
}

bool prescribed(Boundary b, bool thermal) {
    return b != Outlet && (!thermal || b == Inlet || b == Cylinder);
}

double wallWidth(const Face& f, const Geometry& m) {
    const auto inward = f.p - m.offset[f.direction];
    return m.volume[inward] / f.area;
}

void buildTransport(Matrix& matrix, const Geometry& m, double diffusion, bool thermal) {
    std::fill(matrix.begin(), matrix.end(), Coefficients{});
    for (const Face& f : m.internal) {
        const double d = diffusion * f.area / f.distance;
        const double plus = std::max(f.flux, 0.), minus = std::max(-f.flux, 0.);
        matrix[f.p].diagonal += d + plus;
        matrix[f.q].diagonal += d + minus;
        matrix[f.p].neighbour[f.direction] += d + minus;
        matrix[f.q].neighbour[opposite[f.direction]] += d + plus;
    }
    for (const Face& f : m.boundary) {
        auto& row = matrix[f.p];
        if (!prescribed(f.boundary, thermal)) {
            row.diagonal += f.flux;
            continue;
        }
        const double s = f.distance, width = wallWidth(f, m);
        const double d = diffusion * f.area;
        if (width > 0) {
            const double t = 2 * s + .5 * width;
            row.diagonal += d * (1 / s + 1 / (t - s));
            row.neighbour[opposite[f.direction]] += d * s / (t * (t - s));
        } else {
            row.diagonal += d / s;
        }
    }
}

void transportSource(Vector& rhs, const Geometry& m, const Vector& field,
                     double diffusion, bool central, int component) {
    clear(rhs);
    if (central) {
        // Upwind is implicit. Add the difference to the interpolated face value.
        for (const Face& f : m.internal) {
            const double delta = field[f.q] - field[f.p];
            const double correction = f.flux * delta * (f.flux >= 0 ? f.weight : f.weight - 1);
            rhs[f.p] -= correction;
            rhs[f.q] += correction;
        }
    }
    for (const Face& f : m.boundary) {
        if (!prescribed(f.boundary, component == 2)) continue;
        const double width = wallWidth(f, m);
        const double second = width > 0 ? 1 / (2 * f.distance + .5 * width) : 0.;
        const double coefficient = diffusion * f.area * (1 / f.distance + second) - f.flux;
        rhs[f.p] += coefficient * wallValue(f.boundary, component);
    }
}

void pressureGradient(const Geometry& m, const Vector& p, Vector& gx, Vector& gy) {
    clear(gx);
    clear(gy);
    for (const Face& f : m.internal) {
        auto& g = f.direction == East ? gx : gy;
        const double value = std::lerp(p[f.p], p[f.q], f.weight);
        g[f.p] += f.area * (value - p[f.p]);
        g[f.q] -= f.area * (value - p[f.q]);
    }
    for (const Face& f : m.boundary)
        if (f.boundary == Outlet) gx[f.p] -= f.area * p[f.p];
    for (auto k : m.cells) {
        gx[k] /= m.volume[k];
        gy[k] /= m.volume[k];
    }
}

// Alternating line Gauss-Seidel with a tridiagonal solve along each line.
void lines(const Geometry& m, const Matrix& a, const Vector& b, Vector& x,
           const Vector& addedDiagonal, Vector& upper, bool alongX, bool reverse) {
    const Grid& g = m.g;
    const int length = alongX ? g.nx : g.ny;
    const int count = alongX ? g.ny : g.nx;
    const Index step = alongX ? 1 : g.stride;
    const Index lateral = alongX ? g.stride : 1;
    const Direction back = alongX ? West : South, forward = alongX ? East : North;
    const Direction low = alongX ? South : West, high = alongX ? North : East;
    for (int line = 0; line < count; ++line) {
        const int outer = reverse ? count - 1 - line : line;
        const auto first = alongX ? g.index(0, outer) : g.index(outer, 0);
        double previousUpper = 0, previousValue = 0;
        for (Index t = 0, k = first; t < length; ++t, k += step) {
            if (!g.fluid[k]) {
                previousUpper = previousValue = 0;
                upper[t] = 0;
                continue;
            }
            const auto& row = a[k];
            const double lower = row.neighbour[back];
            const double inverse = 1 / (row.diagonal + addedDiagonal[k] - lower * previousUpper);
            upper[t] = previousUpper = row.neighbour[forward] * inverse;
            x[k] = previousValue = (b[k] + lower * previousValue +
                row.neighbour[low] * x[k - lateral] + row.neighbour[high] * x[k + lateral]) * inverse;
        }
        double next = 0;
        for (Index t = length - 1, k = first + t * step; t >= 0; --t, k -= step)
            x[k] = next = x[k] + upper[t] * next;
    }
}

double predictFlux(Geometry& m, const Fields& field, const Vector& uOld, const Vector& vOld,
                   const Vector& gx, const Vector& gy, const Vector& response, double alpha) {
    for (Face& f : m.internal) {
        const bool x = f.direction == East;
        const Vector& velocity = x ? field.u : field.v;
        const Vector& old = x ? uOld : vOld;
        const Vector& grad = x ? gx : gy;
        const auto interpolate = [&](const Vector& z) { return std::lerp(z[f.p], z[f.q], f.weight); };
        const double defect = (field.p[f.q] - field.p[f.p]) / f.distance - interpolate(grad);
        f.flux = f.area * (interpolate(velocity) - interpolate(response) * defect) +
            (1 - alpha) * (f.flux - f.area * interpolate(old));
    }
    double outflow = 0;
    for (const Face& f : m.boundary)
        if (f.boundary == Outlet) outflow += f.area * field.u[f.p];
    const double scale = outflow > 0 ? m.inflow / outflow : 1.;
    for (Face& f : m.boundary)
        if (f.boundary == Outlet) f.flux = scale * f.area * field.u[f.p];
    return scale;
}

double divergence(const Geometry& m, Vector& balance) {
    clear(balance);
    for (const Face& f : m.internal) {
        balance[f.p] += f.flux;
        balance[f.q] -= f.flux;
    }
    double global = 0;
    for (const Face& f : m.boundary) {
        balance[f.p] += f.flux;
        global += f.flux;
    }
    return std::abs(global);
}

double mobility(const Face& f, const Vector& response) {
    const double d = f.boundary == Internal ? std::lerp(response[f.p], response[f.q], f.weight) : response[f.p];
    return d * f.area / f.distance;
}

struct PressureSolver {
    Matrix a;
    Vector inverse, r, z, direction, ad;
    explicit PressureSolver(std::size_t size)
        : a(size), inverse(size), r(size), z(size), direction(size), ad(size) {}

    void assemble(const Geometry& m, const Vector& response) {
        std::fill(a.begin(), a.end(), Coefficients{});
        for (const Face& f : m.internal) {
            const double d = mobility(f, response);
            a[f.p].diagonal += d;
            a[f.q].diagonal += d;
            a[f.p].neighbour[f.direction] = d;
            a[f.q].neighbour[opposite[f.direction]] = d;
        }
        for (const Face& f : m.boundary)
            if (f.boundary == Outlet) a[f.p].diagonal += mobility(f, response);
        if (!m.inflow) a[m.cells.front()].diagonal *= 2; // cavity pressure datum
    }

    void precondition(const Geometry& m) {
        const auto stride = m.g.stride;
        for (auto k : m.cells) {
            const auto& c = a[k].neighbour;
            z[k] = (r[k] + c[West] * z[k - 1] + c[South] * z[k - stride]) * inverse[k];
        }
        for (auto it = m.cells.rbegin(); it != m.cells.rend(); ++it) {
            const auto k = *it;
            const auto& c = a[k].neighbour;
            z[k] += inverse[k] * (c[East] * z[k + 1] + c[North] * z[k + stride]);
        }
    }

    void solve(const Geometry& m, const Vector& balance, Vector& correction, double target) {
        const auto stride = m.g.stride;
        clear(correction);
        if (target == 0.) return;
        for (auto k : m.cells) {
            const auto& c = a[k].neighbour;
            // Modified IC(0): return most of the dropped fill to the diagonal.
            const double west = c[West] * (c[West] + .97 * a[k - 1].neighbour[North]) * inverse[k - 1];
            const double south = c[South] * (c[South] + .97 * a[k - stride].neighbour[East]) * inverse[k - stride];
            inverse[k] = 1 / (a[k].diagonal - west - south);
            r[k] = -balance[k];
        }
        precondition(m);
        direction = z;
        double previous = 0;
        for (auto k : m.cells) previous += r[k] * z[k];
        for (std::size_t iteration = 0; iteration < m.cells.size(); ++iteration) {
            double denominator = 0;
            for (auto k : m.cells) {
                ad[k] = product(a[k], direction, k, stride);
                denominator += direction[k] * ad[k];
            }
            const double step = previous / denominator;
            double norm = 0;
            for (auto k : m.cells) {
                correction[k] += step * direction[k];
                r[k] -= step * ad[k];
                norm += std::abs(r[k]);
            }
            if (!(norm > target)) break; // also return a nonfinite iterate to the residual check
            precondition(m);
            double next = 0;
            for (auto k : m.cells) next += r[k] * z[k];
            const double beta = next / previous;
            for (auto k : m.cells) direction[k] = z[k] + beta * direction[k];
            previous = next;
        }
    }
};

void correct(Geometry& m, Fields& f, const Vector& response, const Vector& pc,
             Vector& gx, Vector& gy, double alphaP) {
    for (Face& e : m.internal) e.flux += mobility(e, response) * (pc[e.p] - pc[e.q]);
    for (Face& e : m.boundary)
        if (e.boundary == Outlet) e.flux += mobility(e, response) * pc[e.p];
    pressureGradient(m, pc, gx, gy);
    for (auto k : m.cells) {
        f.u[k] -= response[k] * gx[k];
        f.v[k] -= response[k] * gy[k];
        f.p[k] += alphaP * pc[k];
    }
}

void vorticity(const Geometry& m, Fields& f) {
    f.vorticity.assign(m.g.fluid.size(), 0.);
    for (const Face& e : m.internal) {
        const bool x = e.direction == East;
        const auto& field = x ? f.v : f.u;
        const double value = (x ? 1. : -1.) * e.area * std::lerp(field[e.p], field[e.q], e.weight);
        f.vorticity[e.p] += value;
        f.vorticity[e.q] -= value;
    }
    for (const Face& e : m.boundary)
        if (e.direction >= South) f.vorticity[e.p] -= (e.direction == North ? 1. : -1.) *
            e.area * wallValue(e.boundary, 0);
        else if (e.boundary == Outlet) f.vorticity[e.p] += e.area * f.v[e.p];
    for (auto k : m.cells) f.vorticity[k] /= m.volume[k];
}
}

void updateVorticity(const Grid& g, Case problem, Fields& f) {
    vorticity(Geometry(g, problem), f);
}

Result solve(const Grid& g, const Parameters& p, const std::function<void(bool, const Residual&)>& progress) {
    const auto start = Clock::now();
    Geometry m(g, p.problem);
    const std::size_t size = g.fluid.size();
    const double diffusion = 1 / p.re, reference = m.inflow ? m.inflow : 1.;
    Result result;
    auto& f = result.fields;
    f = {Vector(size), Vector(size), Vector(size), Vector(size), {}};
    if (!p.initialField.empty()) interpolateInitialField(p.initialField, g, f);
    else if (m.inflow)
        for (auto k : m.cells) f.u[k] = 1.;
    Vector uOld(size), vOld(size), gx(size), gy(size), ru(size), rv(size), relaxation(size);
    Vector response(size), correction(size), balance(size), upper(std::max(g.nx, g.ny));
    Matrix a(size);
    PressureSolver pressure(size);
    result.outletScale = predictFlux(m, f, f.u, f.v, gx, gy, response, 1.);

    const auto assembleMomentum = [&] {
        buildTransport(a, m, diffusion, false);
        pressureGradient(m, f.p, gx, gy);
        transportSource(ru, m, f.u, diffusion, p.secondOrder, 0);
        transportSource(rv, m, f.v, diffusion, p.secondOrder, 1);
        for (auto k : m.cells) {
            ru[k] -= m.volume[k] * gx[k];
            rv[k] -= m.volume[k] * gy[k];
        }
    };
    assembleMomentum();
    for (int it = 1; it <= p.maxIterations; ++it) {
        uOld = f.u;
        vOld = f.v;
        for (auto k : m.cells) {
            relaxation[k] = a[k].diagonal * ((1 - p.alphaU) / p.alphaU);
            response[k] = m.volume[k] / (a[k].diagonal + relaxation[k]);
            ru[k] += relaxation[k] * uOld[k];
            rv[k] += relaxation[k] * vOld[k];
        }
        for (bool alongX : {true, false}) {
            lines(m, a, ru, f.u, relaxation, upper, alongX, it % 2 == 0);
            lines(m, a, rv, f.v, relaxation, upper, alongX, it % 2 == 0);
        }
        result.outletScale = predictFlux(m, f, uOld, vOld, gx, gy, response, p.alphaU);
        divergence(m, balance);
        Residual r{it};
        for (auto k : m.cells) r.continuity += std::abs(balance[k]);
        pressure.assemble(m, response);
        pressure.solve(m, balance, correction, p.pressureRelativeTolerance * r.continuity);
        r.continuity /= reference;
        correct(m, f, response, correction, gx, gy, p.alphaP);
        r.globalMass = divergence(m, balance) / reference;
        assembleMomentum();
        for (auto k : m.cells) {
            r.correctedContinuity += std::abs(balance[k]);
            r.u += std::abs(ru[k] - product(a[k], f.u, k, g.stride));
            r.v += std::abs(rv[k] - product(a[k], f.v, k, g.stride));
        }
        r.correctedContinuity /= reference;
        r.u /= reference;
        r.v /= reference;
        result.flowHistory.push_back(r);
        if (progress) progress(false, r);
        if (!std::isfinite(r.continuity + r.u + r.v + r.globalMass + r.correctedContinuity)) {
            result.stop = "Nonfinite flow residual";
            break;
        }
        if (std::max({r.continuity, r.u, r.v}) < p.tolerance) {
            result.flowConverged = true;
            break;
        }
    }
    result.flowSeconds = std::chrono::duration<double>(Clock::now() - start).count();
    if (result.stop.empty() && !result.flowConverged) result.stop = "Flow iteration limit";
    if (result.flowConverged && p.problem == Case::Square) {
        const auto thermalStart = Clock::now();
        const double diffusivity = diffusion / p.pr;
        buildTransport(a, m, diffusivity, true);
        clear(relaxation);
        transportSource(ru, m, f.temperature, diffusivity, p.secondOrder, 2);
        for (int it = 1; it <= p.maxIterations; ++it) {
            for (bool alongX : {true, false})
                lines(m, a, ru, f.temperature, relaxation, upper, alongX, it % 2 == 0);
            transportSource(ru, m, f.temperature, diffusivity, p.secondOrder, 2);
            Residual r = result.flowHistory.back();
            r.iteration = it;
            for (auto k : m.cells)
                r.temperature += std::abs(ru[k] - product(a[k], f.temperature, k, g.stride));
            r.temperature /= reference;
            result.thermalHistory.push_back(r);
            if (progress) progress(true, r);
            if (!std::isfinite(r.temperature)) {
                result.stop = "Nonfinite thermal residual";
                break;
            }
            if (r.temperature < p.tolerance) {
                result.thermalConverged = true;
                break;
            }
        }
        result.thermalSeconds = std::chrono::duration<double>(Clock::now() - thermalStart).count();
        if (result.stop.empty() && !result.thermalConverged) result.stop = "Thermal iteration limit";
    }
    if (result.stop.empty()) result.stop = "Converged";
    for (const Face& face : m.boundary)
        if (face.boundary == Outlet) result.outletU.push_back(face.flux / face.area);
    vorticity(m, f);
    return result;
}
}
