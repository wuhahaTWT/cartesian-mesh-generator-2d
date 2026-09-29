#include "cartmesh2d/chemistry/DetailedGas.hpp"
#include "cartmesh2d/chemistry/SpeciesMassClosure.hpp"
#include "TraceStableMultiTransport.hpp"

#include "cantera/base/Solution.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/global.h"
#include "cantera/base/YamlWriter.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/Transport.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/ReactorNet.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace cartmesh2d::chemistry {
namespace {
void require(bool condition, const std::string& message) {
    if (!condition) throw std::invalid_argument("DetailedGas: " + message);
}
double finite(double value, const char* name) {
    require(std::isfinite(value), std::string("nonfinite ") + name);
    return value;
}
// Only floating-point composition closure is allowed at this API boundary.
// This is a summation/serialization allowance, not a chemistry solver tolerance.
double compositionAllowance(std::size_t count) {
    return 64 * std::numeric_limits<double>::epsilon() * static_cast<double>(count + 1);
}
// Stefan-Maxwell written directly in mass fluxes. No division by trace X/Y,
// and no inverse-matrix differences cancelling O(1) terms to obtain tiny j_k.
std::vector<double> concentrationFlux(const std::vector<double>& x, const std::vector<double>& mw,
    double mean, double density, const std::vector<double>& binary,
    const std::vector<double>& driving, double& backwardError) {
    const auto n = x.size();
    if (n == 1) return {0};
    std::vector<std::vector<double>> full(n, std::vector<double>(n));
    for (std::size_t k = 0; k < n; ++k) for (std::size_t j = 0; j < n; ++j) if (j != k) {
        const double d = binary[k + n * j];
        require(std::isfinite(d) && d > 0, "invalid binary diffusion coefficient");
        full[k][j] = x[k] * mean / (density * mw[j] * d);
        full[k][k] -= x[j] * mean / (density * mw[k] * d);
    }
    const auto constraint = static_cast<std::size_t>(std::max_element(x.begin(), x.end()) - x.begin());
    // Eliminate the abundant flux analytically with sum(j)=0. A separately
    // scaled all-ones constraint row can pollute trace rows during pivoting.
    std::vector<std::size_t> independent, absent;
    std::vector<double> result(n);
    double fluxScale = 0;
    for (std::size_t k = 0; k < n; ++k) if (k != constraint) {
        require(full[k][k] < 0, "invalid Stefan-Maxwell diagonal");
        fluxScale = std::max(fluxScale, std::abs(driving[k] / full[k][k]));
        if (x[k] == 0) {
            // This row has no off-diagonal coefficients, even when grad X_k
            // is nonzero. Solve it exactly before eliminating the constraint.
            absent.push_back(k); result[k] = driving[k] / full[k][k];
        } else independent.push_back(k);
    }
    if (fluxScale == 0) return result;
    const auto reduced = independent.size();
    std::vector<std::vector<double>> a(reduced, std::vector<double>(reduced));
    std::vector<double> rhs(reduced), columnScale(reduced);
    for (std::size_t i = 0; i < reduced; ++i) {
        rhs[i] = driving[independent[i]];
        for (auto k : absent) rhs[i] -= (full[independent[i]][k] - full[independent[i]][constraint]) * result[k];
        // This is the unit of a linear-system unknown, never a physical
        // concentration/flux floor. Keep the unit normal so equilibration
        // does not quantize a subnormal column before the solve. The solved
        // unknown and final stored flux may still be arbitrarily smaller.
        columnScale[i] = std::max({std::numeric_limits<double>::min(), x[independent[i]] * fluxScale,
            std::abs(driving[independent[i]] / full[independent[i]][independent[i]])});
        require(std::isfinite(columnScale[i]) && columnScale[i] > 0, "unrepresentable trace flux scale");
        for (std::size_t j = 0; j < reduced; ++j)
            a[i][j] = full[independent[i]][independent[j]] - full[independent[i]][constraint];
    }
    const auto original = a;
    const auto originalRhs = rhs;
    for (std::size_t i = 0; i < reduced; ++i) {
        // Scale flux unknowns as well as equations. Trace equations then
        // retain their own arithmetic scale through pivoting; no X cutoff.
        for (std::size_t j = 0; j < reduced; ++j) a[i][j] *= columnScale[j];
        double scale = 0; for (double v : a[i]) scale = std::max(scale, std::abs(v));
        require(std::isfinite(scale) && scale > 0, "singular Stefan-Maxwell row");
        for (double& v : a[i]) v /= scale; rhs[i] /= scale;
    }
    for (std::size_t column = 0; column < reduced; ++column) {
        auto pivot = column;
        for (std::size_t i = column + 1; i < reduced; ++i)
            if (std::abs(a[i][column]) > std::abs(a[pivot][column])) pivot = i;
        require(std::isfinite(a[pivot][column]) && a[pivot][column] != 0, "singular Stefan-Maxwell system");
        std::swap(a[column], a[pivot]); std::swap(rhs[column], rhs[pivot]);
        for (std::size_t i = column + 1; i < reduced; ++i) {
            const double factor = a[i][column] / a[column][column]; a[i][column] = 0;
            for (std::size_t j = column + 1; j < reduced; ++j) a[i][j] -= factor * a[column][j];
            rhs[i] -= factor * rhs[column];
        }
    }
    std::vector<double> flux(reduced);
    for (std::size_t i = reduced; i-- > 0;) {
        long double sum = rhs[i];
        for (std::size_t j = i + 1; j < reduced; ++j) sum -= static_cast<long double>(a[i][j]) * flux[j];
        flux[i] = finite(static_cast<double>(sum / a[i][i]), "Stefan-Maxwell flux");
    }
    for (std::size_t i = 0; i < reduced; ++i) flux[i] *= columnScale[i];
    for (std::size_t i = 0; i < reduced; ++i) {
        long double residual = -originalRhs[i], activity = std::abs(originalRhs[i]), coefficientNorm = 0;
        for (std::size_t j = 0; j < reduced; ++j) {
            const long double term = static_cast<long double>(original[i][j]) * flux[j];
            residual += term; activity += std::abs(term);
            coefficientNorm += std::abs(original[i][j]);
        }
        // IEEE subnormal fluxes have an absolute rounding bound, not eps*j.
        // Include exactly the propagated half-ulp of each stored flux. This
        // is neither a concentration floor nor a relaxed relative tolerance.
        const long double quantization = coefficientNorm * std::numeric_limits<double>::denorm_min() / 2;
        if (std::abs(residual) > compositionAllowance(n) * activity + quantization) {
            std::ostringstream error; error.precision(17);
            error << "Stefan-Maxwell backward error at species " << independent[i]
                  << ": residual=" << residual << ", activity=" << activity;
            throw std::runtime_error(error.str());
        }
        if (activity > 0) backwardError = std::max(backwardError, static_cast<double>(std::abs(residual) / activity));
    }
    long double sum = 0; for (auto k : absent) sum += result[k];
    for (std::size_t i = 0; i < reduced; ++i) { result[independent[i]] = flux[i]; sum += flux[i]; }
    result[constraint] = -static_cast<double>(sum);
    return result;
}
class RangeCheckedReactor final : public Cantera::Reactor {
public:
    RangeCheckedReactor(const std::shared_ptr<Cantera::Solution>& solution, double lower, double upper)
        : Cantera::Reactor(solution, false, "chemistry-source"), lower_(lower), upper_(upper) {}
    void eval(double time, double* lhs, double* rhs) override {
        const double temperature = phase()->thermo()->temperature();
        require(std::isfinite(temperature) && temperature >= lower_ && temperature <= upper_,
                "reaction stage outside common thermodynamic data range");
        Cantera::Reactor::eval(time, lhs, rhs);
    }
private:
    double lower_, upper_;
};
}

struct DetailedGas::Impl {
    GasMechanism info;
    std::shared_ptr<Cantera::Solution> solution;

    explicit Impl(const std::string& file, const std::string& phase) {
        require(std::filesystem::is_regular_file(file), "mechanism must be an existing explicit file");
        info.source = std::filesystem::canonical(file).string();
        solution = Cantera::newSolution(info.source, phase, "none");
        auto gas = solution->thermo();
        require(gas->type() == "ideal-gas", "this backend requires an ideal-gas mixture");
        require(solution->kinetics()->nPhases() == 1, "surface/multiphase mechanisms need a separate coupling");
        auto transport = std::make_shared<TraceStableMultiTransport>();
        transport->init(gas.get()); solution->setTransport(transport);
        info.phase = gas->name();
        info.backendVersion = Cantera::version();
        info.species = gas->speciesNames();
        info.elements = gas->elementNames();
        info.reactions = solution->kinetics()->nReactions();
        require(!info.species.empty() && info.reactions > 0, "a reacting mechanism is required");
        info.molecularWeights.resize(info.species.size());
        gas->getMolecularWeights(info.molecularWeights.data());
        info.minimumTemperature = 0;
        info.maximumTemperature = std::numeric_limits<double>::max();
        for (std::size_t k = 0; k < info.species.size(); ++k) {
            require(gas->charge(k) == 0, "charged species require electro-diffusion and are not supported");
            info.minimumTemperature = std::max(info.minimumTemperature, gas->minTemp(k));
            info.maximumTemperature = std::min(info.maximumTemperature, gas->maxTemp(k));
            for (std::size_t m = 0; m < info.elements.size(); ++m)
                info.atomCounts.push_back(gas->nAtoms(k, m));
        }
        for (std::size_t m = 0; m < info.elements.size(); ++m)
            info.atomicWeights.push_back(gas->atomicWeight(m));
        require(info.maximumTemperature > info.minimumTemperature,
                "species thermodynamic data have no common temperature interval");
        Cantera::YamlWriter writer;
        writer.setPrecision(17);
        writer.skipUserDefined();
        writer.addPhase(solution);
        auto definition = Cantera::AnyMap::fromYamlString(writer.toYamlString());
        for (const auto* key : {"generator", "date", "git-commit"}) definition.erase(key);
        for (auto& p : definition["phases"].asVector<Cantera::AnyMap>()) p.erase("state");
        definition.setMetadata("precision", Cantera::AnyValue(17L));
        info.resolvedDefinition = definition.toYamlString();
    }

    void validateFractions(const std::vector<double>& y) const {
        require(y.size() == info.species.size(), "species count mismatch");
        long double sum = 0;
        for (double value : y) {
            require(std::isfinite(value) && value >= 0, "negative or nonfinite mass fraction");
            sum += value;
        }
        if (std::abs(sum - 1) > compositionAllowance(y.size())) {
            std::ostringstream message; message.precision(17);
            message << "mass fractions do not sum to one; normalization is not implicit; sum=" << sum
                    << ", allowance=" << compositionAllowance(y.size());
            require(false, message.str());
        }
    }

    void validateTemperature(double temperature) const {
        require(std::isfinite(temperature) && temperature >= info.minimumTemperature
                && temperature <= info.maximumTemperature, "temperature outside common thermodynamic data range");
    }

    GasState capture() const {
        auto gas = solution->thermo();
        GasState q;
        q.density = gas->density();
        q.internalEnergyDensity = q.density * gas->intEnergy_mass();
        q.speciesDensities.resize(info.species.size());
        gas->getMassFractions(q.speciesDensities.data());
        for (double& value : q.speciesDensities) value *= q.density;
        return q;
    }

    void restore(const GasState& q) {
        require(std::isfinite(q.density) && q.density > 0, "density must be positive");
        finite(q.internalEnergyDensity, "internal energy density");
        require(q.speciesDensities.size() == info.species.size(), "species density count mismatch");
        std::vector<double> y(q.speciesDensities.size());
        for (std::size_t k = 0; k < y.size(); ++k) y[k] = q.speciesDensities[k] / q.density;
        validateFractions(y);
        auto gas = solution->thermo();
        gas->setMassFractions_NoNorm(y.data());
        // Bracket the energy inversion inside the common thermodynamic data
        // interval. Cantera's general UV setter permits extrapolation; here it
        // must fail instead. Formation energy is retained throughout.
        const double target = q.internalEnergyDensity / q.density;
        double lo = info.minimumTemperature, hi = info.maximumTemperature;
        gas->setState_TD(lo, q.density);
        const double lowerEnergy = gas->intEnergy_mass();
        gas->setState_TD(hi, q.density);
        const double upperEnergy = gas->intEnergy_mass();
        const double endpointAllowance = 64 * std::numeric_limits<double>::epsilon()
            * std::max({std::abs(target), std::abs(lowerEnergy), std::abs(upperEnergy), 1.0});
        require(target >= lowerEnergy - endpointAllowance && target <= upperEnergy + endpointAllowance,
                "internal energy outside common thermodynamic data range");
        double temperature = (lo + hi) / 2;
        for (unsigned iteration = 0; iteration < 100; ++iteration) {
            gas->setState_TD(temperature, q.density);
            const double energy = finite(gas->intEnergy_mass(), "internal energy");
            const double cv = finite(gas->cv_mass(), "cv");
            require(cv > 0, "nonpositive mixture cv");
            const double residual = energy - target;
            const double allowance = 64 * std::numeric_limits<double>::epsilon()
                * std::max({std::abs(target), std::abs(energy), cv * temperature, 1.0});
            if (std::abs(residual) <= allowance) return;
            if (residual > 0) hi = temperature;
            else lo = temperature;
            const double next = temperature - residual / cv;
            temperature = (next >= lo && next <= hi && next != temperature) ? next : (lo + hi) / 2;
        }
        throw std::runtime_error("DetailedGas: temperature inversion did not converge");
    }

    GasProperties currentProperties() const {
        auto gas = solution->thermo();
        GasProperties result;
        result.temperature = gas->temperature();
        validateTemperature(result.temperature);
        result.pressure = finite(gas->pressure(), "pressure");
        result.cp = finite(gas->cp_mass(), "cp");
        result.cv = finite(gas->cv_mass(), "cv");
        result.meanMolecularWeight = finite(gas->meanMolecularWeight(), "molecular weight");
        require(result.pressure > 0 && result.cv > 0 && result.cp > result.cv,
                "invalid mixture thermodynamic properties");
        result.massFractions.resize(info.species.size());
        result.speciesEnthalpies.resize(info.species.size());
        result.massProductionRates.resize(info.species.size());
        result.massRateActivities.resize(info.species.size());
        std::vector<double> destruction(info.species.size());
        gas->getMassFractions(result.massFractions.data());
        validateFractions(result.massFractions);
        gas->getPartialMolarEnthalpies(result.speciesEnthalpies.data());
        solution->kinetics()->getNetProductionRates(result.massProductionRates.data());
        solution->kinetics()->getCreationRates(result.massRateActivities.data());
        solution->kinetics()->getDestructionRates(destruction.data());
        long double release = 0;
        for (std::size_t k = 0; k < info.species.size(); ++k) {
            result.speciesEnthalpies[k] /= info.molecularWeights[k];
            result.massProductionRates[k] *= info.molecularWeights[k];
            result.massRateActivities[k] = (result.massRateActivities[k] + destruction[k]) * info.molecularWeights[k];
            finite(result.speciesEnthalpies[k], "species enthalpy");
            finite(result.massProductionRates[k], "mass production rate");
            finite(result.massRateActivities[k], "mass rate activity");
            release -= static_cast<long double>(result.speciesEnthalpies[k]) * result.massProductionRates[k];
        }
        result.enthalpyReleaseRate = finite(static_cast<double>(release), "enthalpy release rate");
        for (std::size_t m = 0; m < info.elements.size(); ++m)
            result.elementalMassFractions.push_back(finite(gas->elementalMassFraction(m), "element fraction"));
        return result;
    }
};

DetailedGas::DetailedGas(const std::string& file, const std::string& phase)
    : impl_(std::make_unique<Impl>(file, phase)) {}
DetailedGas::~DetailedGas() = default;
DetailedGas::DetailedGas(DetailedGas&&) noexcept = default;
DetailedGas& DetailedGas::operator=(DetailedGas&&) noexcept = default;
const GasMechanism& DetailedGas::mechanism() const { return impl_->info; }

GasState DetailedGas::fromMassFractions(double temperature, double pressure, const std::vector<double>& y) {
    impl_->validateTemperature(temperature);
    impl_->validateFractions(y);
    require(std::isfinite(pressure) && pressure > 0, "pressure must be positive");
    auto gas = impl_->solution->thermo();
    gas->setMassFractions_NoNorm(y.data());
    gas->setState_TP(temperature, pressure);
    return impl_->capture();
}

GasState DetailedGas::fromMoleAmounts(double temperature, double pressure, const std::vector<double>& x) {
    require(x.size() == impl_->info.species.size(), "mole amount count mismatch");
    long double mass = 0;
    for (std::size_t k = 0; k < x.size(); ++k) {
        require(std::isfinite(x[k]) && x[k] >= 0, "negative or nonfinite mole amount");
        mass += static_cast<long double>(x[k]) * impl_->info.molecularWeights[k];
    }
    require(std::isfinite(mass) && mass > 0, "empty or overflowing composition");
    std::vector<double> y(x.size());
    for (std::size_t k = 0; k < x.size(); ++k)
        y[k] = static_cast<double>(static_cast<long double>(x[k]) * impl_->info.molecularWeights[k] / mass);
    return fromMassFractions(temperature, pressure, y);
}

GasProperties DetailedGas::properties(const GasState& q) {
    impl_->restore(q);
    return impl_->currentProperties();
}

double DetailedGas::viscosity(const GasState& q) {
    impl_->restore(q);
    const double value = finite(impl_->solution->transport()->viscosity(), "viscosity");
    require(value > 0, "nonpositive viscosity");
    return value;
}

GasTransport DetailedGas::transport(const GasState& q) {
    impl_->restore(q);
    GasTransport result;
    auto transport = impl_->solution->transport();
    result.viscosity = finite(transport->viscosity(), "viscosity");
    const auto count = impl_->info.species.size();
    result.multicomponentDiffusion.resize(count * count);
    result.binaryDiffusion.resize(count * count);
    result.thermalDiffusion.resize(count);
    transport->getMultiDiffCoeffs(count, result.multicomponentDiffusion.data());
    transport->getBinaryDiffCoeffs(count, result.binaryDiffusion.data());
    // getMultiDiffCoeffs overwrites Cantera's L workspace. Query both thermal
    // properties afterwards so their shared full kinetic solve is reused.
    result.thermalConductivity = finite(transport->thermalConductivity(), "thermal conductivity");
    transport->getThermalDiffCoeffs(result.thermalDiffusion.data());
    require(result.viscosity > 0 && result.thermalConductivity > 0, "nonpositive transport property");
    for (double value : result.multicomponentDiffusion) finite(value, "multicomponent diffusion");
    for (double value : result.binaryDiffusion) require(std::isfinite(value) && value > 0, "invalid binary diffusion");
    for (double value : result.thermalDiffusion) finite(value, "thermal diffusion");
    return result;
}

GasDiffusiveFlux DetailedGas::diffusiveFlux(const GasState& state, const GasDiffusionGradient& gradient) {
    impl_->restore(state);
    finite(gradient.temperature, "temperature gradient");
    finite(gradient.logPressure, "log-pressure gradient");
    const auto& info = impl_->info;
    const auto count = info.species.size();
    require(gradient.moleFractions.size() == count, "mole-fraction gradient count mismatch");
    long double gradientSum = 0, gradientMagnitude = 0;
    for (double value : gradient.moleFractions) {
        finite(value, "mole-fraction gradient");
        gradientSum += value;
        gradientMagnitude += std::abs(value);
    }
    require(std::abs(gradientSum) <= compositionAllowance(count) * gradientMagnitude,
            "mole-fraction gradients do not sum to zero");
    auto gas = impl_->solution->thermo();
    std::vector<double> x(count), y(count), enthalpy(count), driving(count);
    gas->getMoleFractions(x.data());
    gas->getMassFractions(y.data());
    gas->getPartialMolarEnthalpies(enthalpy.data());
    const double temperature = gas->temperature(), meanWeight = gas->meanMolecularWeight();
    const auto coefficients = transport(state);
    // Ideal-gas Stefan-Maxwell driving force, including pressure diffusion.
    // The Cantera D_ki matrix follows the positive-sign, W_i-weighted convention;
    // it is not a set of independent species Fick diffusivities.
    for (std::size_t k = 0; k < count; ++k)
        driving[k] = finite(gradient.moleFractions[k] + (x[k] - y[k]) * gradient.logPressure,
                            "diffusion driving force");
    GasDiffusiveFlux result;
    result.species.resize(count);
    const auto concentration = concentrationFlux(x, info.molecularWeights, meanWeight, state.density,
        coefficients.binaryDiffusion, driving, result.concentrationSolveResidual);
    result.pureSpeciesLimit = std::count_if(state.speciesDensities.begin(), state.speciesDensities.end(),
        [](double value) { return value > 0; }) == 1;
    long double massSum = 0, activity = 0;
    for (std::size_t k = 0; k < count; ++k) {
        long double sum = 0;
        const double factor = state.density * info.molecularWeights[k] / (meanWeight * meanWeight);
        for (std::size_t j = 0; j < count; ++j) {
            const long double term = static_cast<long double>(factor) * info.molecularWeights[j]
                * coefficients.multicomponentDiffusion[k + count * j] * driving[j];
            sum += term;
            activity += std::abs(static_cast<long double>(factor) * info.molecularWeights[j]
                * coefficients.multicomponentDiffusion[k + count * j])
                * (std::abs(gradient.moleFractions[j]) + std::abs((x[j] - y[j]) * gradient.logPressure));
        }
        const double backendThermal = finite(coefficients.thermalDiffusion[k] * gradient.temperature / temperature,
                                              "thermal diffusion flux");
        // MultiTransport regularizes X to max(X, 1e-20). At an EXACT zero with
        // zero species gradient the Stefan-Maxwell flux limit is exactly zero;
        // otherwise a strict-positive FV update could create negative traces.
        // A pure gas has no Soret separation. Neither rule is a concentration
        // threshold, clipping, or a reduced-species mechanism.
        // Soret vanishes at an exactly absent face species even if its
        // concentration gradient is nonzero (for example an inlet face).
        // That gradient can still drive a genuine Stefan-Maxwell mass flux.
        const double thermal = result.pureSpeciesLimit || state.speciesDensities[k] == 0 ? 0 : backendThermal;
        const double backendFlux = finite(concentration[k] - backendThermal, "regularized thermal mass flux");
        result.matrixFormDifferenceL1 += std::abs(concentration[k] - static_cast<double>(sum));
        result.species[k] = state.speciesDensities[k] == 0 && gradient.moleFractions[k] == 0
            ? 0 : finite(concentration[k] - thermal, "diffusive mass flux");
        result.zeroLimitCorrection += std::abs(result.species[k] - backendFlux);
        massSum += result.species[k];
        activity += std::abs(thermal);
    }
    result.rawMassResidual = static_cast<double>(massSum);
    // Bound cancellation by the arithmetic activity, not the net flux near
    // equilibrium. No dimensional floor and no user-tunable mass correction.
    require(std::abs(massSum) <= compositionAllowance(count) * activity,
            "multicomponent mass-flux closure failed");
    const auto closureSpecies = static_cast<std::size_t>(std::max_element(y.begin(), y.end()) - y.begin());
    long double others = 0;
    for (std::size_t k = 0; k < count; ++k) if (k != closureSpecies) others += result.species[k];
    const double closed = -static_cast<double>(others);
    result.closureCorrection = closed - result.species[closureSpecies];
    result.species[closureSpecies] = closed;
    finite(result.zeroLimitCorrection, "zero-concentration limit correction");
    long double enthalpyFlux = 0;
    for (std::size_t k = 0; k < count; ++k)
        enthalpyFlux += result.species[k] * static_cast<long double>(enthalpy[k] / info.molecularWeights[k]);
    result.conduction = finite(-coefficients.thermalConductivity * gradient.temperature, "conductive heat flux");
    result.speciesEnthalpy = finite(static_cast<double>(enthalpyFlux), "diffusive enthalpy flux");
    result.energy = finite(result.conduction + result.speciesEnthalpy, "diffusive energy flux");
    return result;
}

ChemistryStep DetailedGas::advanceConstantVolume(const GasState& input, double dt, const ChemistryControls& c) {
    ChemistryStep result;
    result.requestedTime = dt;
    try {
        require(std::isfinite(dt) && dt > 0, "chemistry time step must be positive");
        require(std::isfinite(c.relativeTolerance) && c.relativeTolerance > 0 && c.relativeTolerance < 1,
                "invalid relative integration tolerance");
        require(std::isfinite(c.absoluteTolerance) && c.absoluteTolerance > 0, "invalid absolute integration tolerance");
        require(c.maximumInternalSteps > 0, "invalid internal step budget");
        require(std::isfinite(c.maximumInternalStep) && c.maximumInternalStep >= 0, "invalid maximum internal step");
        require(std::isfinite(c.maximumElementDrift) && c.maximumElementDrift > 0,
                "invalid elemental conservation tolerance");
        require(std::isfinite(c.maximumEnergyDrift) && c.maximumEnergyDrift > 0,
                "invalid energy conservation tolerance");
        impl_->restore(input);
        const auto before = impl_->currentProperties();
        // The generic Reactor evolves total INTERNAL ENERGY, not temperature.
        // No surfaces, walls or flow devices are installed. Its energy is
        // therefore constant without imposing a post-integration correction.
        auto reactor = std::make_shared<RangeCheckedReactor>(impl_->solution,
            impl_->info.minimumTemperature, impl_->info.maximumTemperature);
        Cantera::ReactorNet network(reactor);
        network.setTolerances(c.relativeTolerance, c.absoluteTolerance);
        network.setMaxSteps(c.maximumInternalSteps);
        network.setMaxTimeStep(c.maximumInternalStep);
        try {
            network.advance(dt);
        } catch (...) {
            result.reachedTime = network.time();
            throw;
        }
        result.reachedTime = network.time();
        const auto stats = network.solverStats();
        result.internalSteps = stats.at("steps").asInt();
        result.rhsEvaluations = stats.at("rhs_evals").asInt();
        require(result.reachedTime == dt, "source integrator did not reach requested physical time");
        auto candidate = impl_->capture();
        // Preserve the algebraic mass constraint across repeated source
        // restarts. Validate the raw redundant species first; do not repair
        // a chemistry integration defect outside the summation allowance.
        const auto closure = closeSpeciesMassRoundoff(candidate.density, candidate.speciesDensities);
        result.massClosureChange.assign(candidate.speciesDensities.size(), 0);
        result.massClosureChange[closure.species] = closure.change;
        impl_->restore(candidate);
        const auto after = impl_->currentProperties(); // includes strict positivity, no clipping
        for (std::size_t m = 0; m < before.elementalMassFractions.size(); ++m)
            result.maximumElementDrift = std::max(result.maximumElementDrift,
                std::abs(after.elementalMassFractions[m] - before.elementalMassFractions[m]));
        const double energyScale = std::max(std::abs(input.internalEnergyDensity),
                                            input.density * before.cv * before.temperature);
        result.relativeEnergyDrift = std::abs(candidate.internalEnergyDensity - input.internalEnergyDensity) / energyScale;
        result.relativeDensityDrift = std::abs(candidate.density - input.density) / input.density;
        require(result.relativeDensityDrift <= compositionAllowance(1), "fixed-volume source changed density");
        require(result.maximumElementDrift <= c.maximumElementDrift, "elemental conservation gate failed");
        require(result.relativeEnergyDrift <= c.maximumEnergyDrift, "internal energy conservation gate failed");
        impl_->restore(candidate); // the conserved representation must also be valid
        result.accepted = candidate;
    } catch (const std::exception& error) {
        result.failure = error.what();
    }
    return result;
}

} // namespace cartmesh2d::chemistry
