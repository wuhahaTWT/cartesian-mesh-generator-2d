#include "cartmesh2d/fv/ReactingDiffusion2D.hpp"
#include "cartmesh2d/chemistry/SpeciesMassClosure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <array>
#include <map>
#include <stdexcept>
#include <utility>

namespace cartmesh2d::fv {
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(std::string("reacting diffusion: ") + message);
}
double checked(double value) {
    require(std::isfinite(value), "nonfinite arithmetic");
    return value;
}
}

ReactingDiffusionOperator2D::ReactingDiffusionOperator2D(FvMesh2D mesh,
    std::vector<ReactingDiffusionBoundary2D> boundaries) : mesh_(std::move(mesh)) {
    validateFvMesh2D(mesh_);
    boundaries_.resize(mesh_.faces.size());
    for (auto& b : boundaries) {
        require(b.face < boundaries_.size() && !mesh_.faces[b.face].neighbour && !boundaries_[b.face],
                "invalid or duplicate boundary face");
        require(b.kind == ReactingDiffusionBoundaryKind2D::AdiabaticWall
                || b.kind == ReactingDiffusionBoundaryKind2D::IsothermalWall
                || b.kind == ReactingDiffusionBoundaryKind2D::Reservoir
                || b.kind == ReactingDiffusionBoundaryKind2D::ZeroFluxOutflow, "invalid boundary kind");
        require((b.kind == ReactingDiffusionBoundaryKind2D::Reservoir) == b.reservoir.has_value(),
                "reservoir state missing or inactive");
        require(b.kind == ReactingDiffusionBoundaryKind2D::IsothermalWall
                ? std::isfinite(b.temperature) && b.temperature > 0 : b.temperature == 0,
                "invalid or inactive wall temperature");
        boundaries_[b.face] = std::move(b);
    }
    for (std::size_t f = 0; f < mesh_.faces.size(); ++f)
        require(mesh_.faces[f].neighbour.has_value() || boundaries_[f].has_value(), "missing boundary condition");
    temperatureStencil_ = makeStencil(true);
    compositionStencil_ = makeStencil(false);
}

ReactingDiffusionOperator2D::Stencil ReactingDiffusionOperator2D::makeStencil(bool temperature) const {
    Stencil result(mesh_.cells.size());
    for (std::size_t i = 0; i < mesh_.cells.size(); ++i) {
        double xx = 0, xy = 0, yy = 0;
        std::vector<Sample> samples;
        for (auto id : mesh_.cells[i].faces) {
            const auto& f = mesh_.faces[id];
            const auto other = f.owner == i ? f.neighbour : std::optional<std::size_t>(f.owner);
            const bool prescribed = !other && (boundaries_[id]->kind == ReactingDiffusionBoundaryKind2D::Reservoir
                || (temperature && boundaries_[id]->kind == ReactingDiffusionBoundaryKind2D::IsothermalWall));
            auto d = (other ? mesh_.cells[*other].centre : f.centre) - mesh_.cells[i].centre;
            if (!other && !prescribed) d = f.areaVector; // homogeneous normal derivative in wall stencil
            const double norm = squaredNorm(d);
            require(std::isfinite(norm) && norm > 0, "invalid gradient stencil distance");
            xx += d.x * d.x / norm; xy += d.x * d.y / norm; yy += d.y * d.y / norm;
            if (other || prescribed) samples.push_back({other ? *other : id, !other, {d.x / norm, d.y / norm}});
        }
        const double determinant = xx * yy - xy * xy;
        require(determinant > 64 * std::numeric_limits<double>::epsilon() * (xx + yy) * (xx + yy),
                "rank-deficient gradient stencil");
        for (auto s : samples) {
            s.weight = {checked((yy * s.weight.x - xy * s.weight.y) / determinant),
                        checked((xx * s.weight.y - xy * s.weight.x) / determinant)};
            result[i].push_back(s);
        }
    }
    return result;
}

std::vector<Vector2D> ReactingDiffusionOperator2D::gradients(const Stencil& stencil,
    const std::vector<double>& cells, const std::vector<double>& boundary) const {
    std::vector<Vector2D> result(cells.size());
    for (std::size_t i = 0; i < cells.size(); ++i) {
        for (const auto& s : stencil[i]) {
            const double difference = (s.boundary ? boundary[s.index] : cells[s.index]) - cells[i];
            result[i].x += s.weight.x * difference;
            result[i].y += s.weight.y * difference;
        }
        checked(result[i].x); checked(result[i].y);
    }
    return result;
}

ReactingDiffusionResult2D ReactingDiffusionOperator2D::evaluate(chemistry::DetailedGas& gas,
    const std::vector<chemistry::GasState>& states, bool estimateRate) const {
    require(states.size() == mesh_.cells.size(), "state/mesh cell count mismatch");
    const auto ns = gas.mechanism().species.size(), nc = states.size(), nf = mesh_.faces.size();
    std::vector<chemistry::GasProperties> properties, boundaryProperties(nf);
    std::vector<std::vector<double>> values(ns + 2, std::vector<double>(nc));
    std::vector<std::vector<double>> boundaryValues(ns + 2, std::vector<double>(nf));
    std::vector<long double> abundance(ns);
    const auto fill = [&](const chemistry::GasProperties& p, std::vector<std::vector<double>>& fields, std::size_t i) {
        fields[0][i] = p.temperature;
        fields[1][i] = std::log(p.pressure);
        for (std::size_t k = 0; k < ns; ++k)
            fields[k + 2][i] = p.massFractions[k] * p.meanMolecularWeight / gas.mechanism().molecularWeights[k];
    };
    for (std::size_t i = 0; i < nc; ++i) {
        properties.push_back(gas.properties(states[i]));
        fill(properties.back(), values, i);
        for (std::size_t k = 0; k < ns; ++k) abundance[k] += values[k + 2][i];
    }
    for (std::size_t f = 0; f < nf; ++f) if (boundaries_[f]) {
        const auto& b = *boundaries_[f];
        const auto owner = mesh_.faces[f].owner;
        if (b.reservoir) boundaryProperties[f] = gas.properties(*b.reservoir);
        else {
            boundaryProperties[f] = properties[owner];
            if (b.kind == ReactingDiffusionBoundaryKind2D::IsothermalWall) {
                require(b.temperature >= gas.mechanism().minimumTemperature
                        && b.temperature <= gas.mechanism().maximumTemperature, "wall outside thermodynamic data range");
                boundaryProperties[f].temperature = b.temperature;
            }
        }
        fill(boundaryProperties[f], boundaryValues, f);
    }
    const auto dependent = static_cast<std::size_t>(std::max_element(abundance.begin(), abundance.end()) - abundance.begin());
    std::vector<std::vector<Vector2D>> gradient(ns + 2);
    gradient[0] = gradients(temperatureStencil_, values[0], boundaryValues[0]);
    gradient[1] = gradients(compositionStencil_, values[1], boundaryValues[1]);
    gradient[dependent + 2].resize(nc);
    // Reconstruct N-1 independent mole fractions; the exact identity sum(X)=1
    // defines the remaining gradient. No conserved species density is changed.
    for (std::size_t k = 0; k < ns; ++k) if (k != dependent) {
        gradient[k + 2] = gradients(compositionStencil_, values[k + 2], boundaryValues[k + 2]);
        for (std::size_t i = 0; i < nc; ++i) {
            gradient[dependent + 2][i].x -= gradient[k + 2][i].x;
            gradient[dependent + 2][i].y -= gradient[k + 2][i].y;
        }
    }
    ReactingDiffusionResult2D result;
    result.faces.resize(nf);
    result.speciesResidual.assign(nc, std::vector<double>(ns));
    result.energyResidual.resize(nc);
    result.boundarySpeciesFlux.resize(ns);
    std::vector<std::vector<double>> rateRows;
    if (estimateRate) { result.rate.resize(nc); rateRows.assign(nc, std::vector<double>(ns + 1)); }
    for (std::size_t f = 0; f < nf; ++f) {
        const auto& face = mesh_.faces[f];
        auto& out = result.faces[f];
        out.species.resize(ns);
        const auto* b = boundaries_[f] ? &*boundaries_[f] : nullptr;
        if (b && (b->kind == ReactingDiffusionBoundaryKind2D::AdiabaticWall
                  || b->kind == ReactingDiffusionBoundaryKind2D::ZeroFluxOutflow)) continue;
        const auto owner = face.owner;
        const double length = std::hypot(face.areaVector.x, face.areaVector.y);
        const double w = face.neighbour ? face.neighbourWeight : 1;
        const auto& left = properties[owner];
        const auto& right = face.neighbour ? properties[*face.neighbour] : boundaryProperties[f];
        const double t = (1 - w) * left.temperature + w * right.temperature;
        const double p = std::exp((1 - w) * std::log(left.pressure) + w * std::log(right.pressure));
        std::vector<double> y(ns);
        for (std::size_t k = 0; k < ns; ++k) y[k] = (1 - w) * left.massFractions[k] + w * right.massFractions[k];
        (void)chemistry::closeSpeciesMassRoundoff(1, y);
        const auto faceState = gas.fromMassFractions(t, p, y);
        const auto derivative = [&](std::size_t field) {
            auto g = gradient[field][owner];
            if (face.neighbour) {
                const auto h = gradient[field][*face.neighbour];
                g = {g.x * (1 - w) + w * h.x, g.y * (1 - w) + w * h.y};
            }
            const double delta = (face.neighbour ? values[field][*face.neighbour] : boundaryValues[field][f]) - values[field][owner];
            return checked((face.transmissibility * delta + dot(g, face.correction)) / length);
        };
        if (b && b->kind == ReactingDiffusionBoundaryKind2D::IsothermalWall) {
            // Impermeability imposes j_k.n=0 including thermodiffusion. It is
            // not a zero-concentration-gradient reservoir approximation.
            out.conduction = checked(-gas.transport(faceState).thermalConductivity * derivative(0) * length);
            out.energy = out.conduction;
        } else {
            chemistry::GasDiffusionGradient d;
            d.temperature = derivative(0); d.logPressure = derivative(1); d.moleFractions.resize(ns);
            long double sum = 0;
            for (std::size_t k = 0; k < ns; ++k) if (k != dependent) {
                d.moleFractions[k] = derivative(k + 2); sum += d.moleFractions[k];
            }
            d.moleFractions[dependent] = -static_cast<double>(sum);
            const auto flux = gas.diffusiveFlux(faceState, d);
            for (std::size_t k = 0; k < ns; ++k) out.species[k] = checked(flux.species[k] * length);
            out.conduction = checked(flux.conduction * length);
            out.speciesEnthalpy = checked(flux.speciesEnthalpy * length);
            out.energy = checked(flux.energy * length);
            out.rawMassResidual = checked(flux.rawMassResidual * length);
            out.closureCorrection = checked(flux.closureCorrection * length);
            out.zeroLimitCorrection = checked(flux.zeroLimitCorrection * length);
        }
        if (estimateRate) {
            // Differentiate the SAME frozen constitutive flux and scalar
            // normal-gradient stencils, rather than max(D)/h^2 on a cut cell.
            std::map<std::size_t, std::array<double,2>> weights; // T, composition/ln(p)
            for (std::size_t type = 0; type < 2; ++type) {
                const auto& stencil = type == 0 ? temperatureStencil_ : compositionStencil_;
                weights[owner][type] -= face.transmissibility / length;
                if (face.neighbour) weights[*face.neighbour][type] += face.transmissibility / length;
                const auto add = [&](std::size_t cell, double factor) {
                    for (const auto& s : stencil[cell]) {
                        const double a = factor * dot(s.weight, face.correction) / length;
                        weights[cell][type] -= a;
                        if (!s.boundary) weights[s.index][type] += a;
                    }
                };
                add(owner, face.neighbour ? 1 - w : 1);
                if (face.neighbour) add(*face.neighbour, w);
            }
            const auto fp = gas.properties(faceState);
            const auto ft = gas.transport(faceState);
            const auto& mw = gas.mechanism().molecularWeights;
            std::vector<double> kt(ns + 1), kp(ns + 1);
            std::vector<std::vector<double>> kx(ns + 1, std::vector<double>(ns));
            kt[ns] = -ft.thermalConductivity;
            const bool wall = b && b->kind == ReactingDiffusionBoundaryKind2D::IsothermalWall;
            if (!wall) for (std::size_t k = 0; k < ns; ++k) {
                kt[k] = -ft.thermalDiffusion[k] / fp.temperature;
                for (std::size_t j = 0; j < ns; ++j) {
                    kx[k][j] = faceState.density * mw[k] * mw[j] / (fp.meanMolecularWeight * fp.meanMolecularWeight)
                        * ft.multicomponentDiffusion[k + ns * j];
                    kp[k] += kx[k][j] * fp.massFractions[j] * (fp.meanMolecularWeight / mw[j] - 1);
                    kx[ns][j] += fp.speciesEnthalpies[k] * kx[k][j];
                }
                kt[ns] += fp.speciesEnthalpies[k] * kt[k];
                kp[ns] += fp.speciesEnthalpies[k] * kp[k];
            }
            for (const auto& [cell, weight] : weights) {
                const auto& cp = properties[cell];
                const auto ref = static_cast<std::size_t>(std::max_element(cp.massFractions.begin(), cp.massFractions.end()) - cp.massFractions.begin());
                // Similarity-transform the frozen conserved Jacobian to
                // [delta T/T, independent delta Y] at fixed density. This
                // removes artificial formation-energy scaling stiffness;
                // the conservative operator and its eigenvalues are unchanged.
                for (std::size_t column = 0; column <= ns; ++column) if (column != ref) {
                    const double inverseWeight = column == ns ? 0 : 1 / mw[column] - 1 / mw[ref];
                    const double dt = column == ns ? cp.temperature : 0;
                    const double dp = dt / cp.temperature + cp.meanMolecularWeight * inverseWeight;
                    std::vector<double> fluxDerivative(ns + 1);
                    for (std::size_t row = 0; row <= ns; ++row) {
                        double value = kt[row] * weight[0] * dt + kp[row] * weight[1] * dp;
                        if (column != ns) for (std::size_t j = 0; j < ns; ++j) {
                            const double dx = cp.meanMolecularWeight * ((j == column ? 1 / mw[j] : 0)
                                - (j == ref ? 1 / mw[j] : 0)
                                - cp.massFractions[j] * cp.meanMolecularWeight / mw[j] * inverseWeight);
                            value += kx[row][j] * weight[1] * dx;
                        }
                        fluxDerivative[row] = value * length;
                    }
                    const auto accumulate = [&](std::size_t i) {
                        const auto& pcell = properties[i];
                        const double mass = mesh_.cells[i].area * states[i].density;
                        long double temperatureFlux = fluxDerivative[ns];
                        for (std::size_t k = 0; k < ns; ++k) {
                            rateRows[i][k] += checked(std::abs(fluxDerivative[k]) / mass);
                            const double energy = pcell.speciesEnthalpies[k] - 8314.46261815324 * pcell.temperature / mw[k];
                            temperatureFlux -= static_cast<long double>(energy) * fluxDerivative[k];
                        }
                        rateRows[i][ns] += checked(static_cast<double>(std::abs(temperatureFlux) / (mass * pcell.cv * pcell.temperature)));
                    };
                    accumulate(owner);
                    if (face.neighbour) accumulate(*face.neighbour);
                }
            }
        }
        for (std::size_t k = 0; k < ns; ++k) {
            result.speciesResidual[owner][k] += out.species[k];
            if (face.neighbour) result.speciesResidual[*face.neighbour][k] -= out.species[k];
            else result.boundarySpeciesFlux[k] += out.species[k];
        }
        result.energyResidual[owner] += out.energy;
        if (face.neighbour) result.energyResidual[*face.neighbour] -= out.energy;
        else result.boundaryEnergyFlux += out.energy;
    }
    for (const auto& row : result.speciesResidual) for (double value : row) checked(value);
    for (double value : result.energyResidual) checked(value);
    for (double value : result.boundarySpeciesFlux) checked(value);
    checked(result.boundaryEnergyFlux);
    if (estimateRate) for (std::size_t i = 0; i < nc; ++i)
        result.rate[i] = checked(*std::max_element(rateRows[i].begin(), rateRows[i].end()));
    return result;
}
} // namespace cartmesh2d::fv
