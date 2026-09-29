#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace cartmesh2d::chemistry {

struct GasMechanism {
    std::string source, phase, backendVersion;
    std::vector<std::string> species, elements;
    std::vector<double> molecularWeights; // kg/kmol
    std::vector<double> atomicWeights; // kg/kmol
    std::vector<double> atomCounts; // species-major: [species * nElements + element]
    std::size_t reactions = 0;
    double minimumTemperature = 0, maximumTemperature = 0; // common data range, K
};

// Volumetric conserved state. Internal energy includes species formation energy
// and may be negative; kinetic energy is owned by the future flow discretization.
struct GasState {
    double density = 0; // kg/m^3
    double internalEnergyDensity = 0; // J/m^3, including chemical energy
    std::vector<double> speciesDensities; // kg/m^3, same ordering as GasMechanism
};

struct GasProperties {
    double temperature = 0, pressure = 0; // K, Pa
    double cp = 0, cv = 0, meanMolecularWeight = 0; // J/(kg K), kg/kmol
    std::vector<double> massFractions, elementalMassFractions;
    std::vector<double> speciesEnthalpies; // J/kg, including formation enthalpy
    std::vector<double> massProductionRates; // kg/(m^3 s)
    // Gross creation + destruction, kg/(m^3 s); roundoff scale for cancellation
    // in net source assembly, including species appearing on both reaction sides.
    std::vector<double> massRateActivities;
    // Diagnostic -sum(h_k * omega_k). Do NOT add it to an energy equation
    // already conserving internal + chemical energy (that would count heat twice).
    double enthalpyReleaseRate = 0; // W/m^3
};

struct GasTransport {
    double viscosity = 0, thermalConductivity = 0; // Pa s, W/(m K)
    // Cantera multicomponent coefficients [m^2/s], column-major D[i + n*j].
    // These require the corresponding multicomponent flux formula; they are
    // not independent Fick coefficients multiplying each species' mass gradient.
    std::vector<double> multicomponentDiffusion;
    std::vector<double> thermalDiffusion; // Soret coefficients, kg/(m s)
};

// Derivatives along one physical unit direction. A caller may evaluate each
// Cartesian component or the outward face normal. sum(dX_k/dn) must be zero.
struct GasDiffusionGradient {
    double temperature = 0; // dT/dn, K/m
    double logPressure = 0; // d(ln p)/dn, 1/m
    std::vector<double> moleFractions; // dX_k/dn, 1/m
};

struct GasDiffusiveFlux {
    std::vector<double> species; // kg/(m^2 s), mass-average reference velocity
    double conduction = 0, speciesEnthalpy = 0, energy = 0; // W/m^2
    // Raw sum before an explicitly bounded roundoff closure correction. This
    // cannot be used to conceal a non-conservative constitutive model.
    double rawMassResidual = 0, closureCorrection = 0; // kg/(m^2 s)
    // L1 difference from the backend's 1e-20 concentration regularization when
    // applying exact absent-species/pure-gas limits; never changes a GasState.
    double zeroLimitCorrection = 0; // kg/(m^2 s)
    bool pureSpeciesLimit = false;
};

struct ChemistryControls {
    // CVODES local error weights, not physical accuracy qualification. The
    // absolute weight applies to its mass, volume, energy and mass-fraction ODEs.
    double relativeTolerance = 1e-8, absoluteTolerance = 1e-18;
    int maximumInternalSteps = 20000;
    double maximumInternalStep = 0; // s; zero leaves the solver's step unbounded
    // Numerical substep gates: absolute elemental mass-fraction change and
    // energy change / max(|rho*e|, rho*cv*T), respectively. Not flame accuracy.
    double maximumElementDrift = 1e-8, maximumEnergyDrift = 1e-8;
};

struct ChemistryStep {
    std::optional<GasState> accepted;
    std::string failure;
    double requestedTime = 0, reachedTime = 0;
    double maximumElementDrift = 0, relativeEnergyDrift = 0;
    double relativeDensityDrift = 0;
    long internalSteps = 0, rhsEvaluations = 0;
};

// A serial, non-copyable native C++ chemistry context, one per worker. Its
// mutable Cantera workspace is never the caller's accepted state. No 3D core,
// Python interpreter, single-step fallback, clipping or hidden normalization.
class DetailedGas {
public:
    explicit DetailedGas(const std::string& mechanismFile, const std::string& phase = "");
    ~DetailedGas();
    DetailedGas(DetailedGas&&) noexcept;
    DetailedGas& operator=(DetailedGas&&) noexcept;
    DetailedGas(const DetailedGas&) = delete;
    DetailedGas& operator=(const DetailedGas&) = delete;

    [[nodiscard]] const GasMechanism& mechanism() const;
    [[nodiscard]] GasState fromMassFractions(double temperature, double pressure,
                                             const std::vector<double>& massFractions);
    // Explicit composition conversion; the inputs are relative mole amounts,
    // not already-normalized fractions. Negative/unknown entries are rejected.
    [[nodiscard]] GasState fromMoleAmounts(double temperature, double pressure,
                                           const std::vector<double>& moleAmounts);
    [[nodiscard]] GasProperties properties(const GasState&);
    [[nodiscard]] GasTransport transport(const GasState&);
    // Full multicomponent concentration/pressure diffusion plus Soret.
    // q = -lambda*dT/dn + sum(h_k*j_k), with formation-inclusive enthalpy.
    // Neutral ideal gases, no differential body forces or Dufour heat term.
    [[nodiscard]] GasDiffusiveFlux diffusiveFlux(const GasState&, const GasDiffusionGradient&);
    // Homogeneous source substep for operator coupling: fixed volume, no mass
    // transfer, no external heat. Uses CVODES with conserved internal energy.
    // Every failure leaves input untouched and returns no accepted candidate.
    [[nodiscard]] ChemistryStep advanceConstantVolume(const GasState&, double timeStep,
                                                       const ChemistryControls& = {});
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace cartmesh2d::chemistry
