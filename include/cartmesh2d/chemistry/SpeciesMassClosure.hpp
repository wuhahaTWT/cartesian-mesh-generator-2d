#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cartmesh2d::chemistry {

struct SpeciesMassClosure {
    std::size_t species = 0; // relative to the first species entry
    double rawResidual = 0, change = 0; // same units as the supplied masses
};

// Explicit N-1 composition representation for integrator stages / face values.
// Check the redundant N-species result BEFORE deriving the most abundant
// component from total mass. Only the existing floating-point summation
// allowance is accepted; negative species and larger defects still fail.
// There is no rescaling, trace threshold, or change to the other N-1 species.
// Stage callers retain `change` in accepted-step numerical error diagnostics.
// Face interpolation uses the closed composition in its shared physical flux.
// This is deliberately not part of the public thermodynamic input conversion.
inline SpeciesMassClosure closeSpeciesMassRoundoff(double total,
    std::vector<double>& masses, std::size_t first = 0) {
    if (!(std::isfinite(total) && total > 0) || first >= masses.size())
        throw std::invalid_argument("species mass closure: invalid total/count");
    const auto count = masses.size() - first;
    auto dependent = first;
    for (std::size_t k = first; k < masses.size(); ++k) {
        if (!(std::isfinite(masses[k]) && masses[k] >= 0))
            throw std::invalid_argument("species mass closure: negative or nonfinite species");
        if (masses[k] > masses[dependent]) dependent = k;
    }
    const auto sum = [&](bool omitDependent) {
        double s = 0, error = 0;
        for (std::size_t k = first; k < masses.size(); ++k) if (!omitDependent || k != dependent) {
            const double next = s + masses[k];
            error += std::abs(s) >= masses[k] ? (s - next) + masses[k] : (masses[k] - next) + s;
            s = next;
        }
        return s + error;
    };
    SpeciesMassClosure result;
    result.species = dependent - first;
    result.rawResidual = sum(false) - total;
    const double allowance = 64 * std::numeric_limits<double>::epsilon() * static_cast<double>(count + 1);
    if (!(std::isfinite(result.rawResidual) && std::abs(result.rawResidual / total) <= allowance))
        throw std::invalid_argument("species mass closure: defect exceeds roundoff allowance");
    const double closed = total - sum(true);
    if (!(std::isfinite(closed) && closed >= 0))
        throw std::invalid_argument("species mass closure: inadmissible dependent species");
    result.change = closed - masses[dependent];
    masses[dependent] = closed;
    return result;
}

} // namespace cartmesh2d::chemistry
