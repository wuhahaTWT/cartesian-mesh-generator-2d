#pragma once

#include "cantera/transport/MultiTransport.h"
#include "cantera/thermo/ThermoPhase.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cartmesh2d::chemistry {

// Same Cantera 3.2 kinetic-theory collision/heat-mode model. Its augmented
// L00 equations subtract O(1) terms to recover O(X_k) trace-species fluxes.
// Write those equations as velocity differences in the mass-average frame,
// and divide each kinetic equation by its species abundance before solving.
class TraceStableMultiTransport final : public Cantera::MultiTransport {
public:
    void getThermalDiffCoeffs(double* dt) override {
        solveLMatrixEquation();
        std::vector<double> x(m_nsp), y(m_nsp);
        m_thermo->getMoleFractions(x.data());
        m_thermo->getMassFractions(y.data());
        long double sum = 0;
        for (std::size_t k = 0; k < m_nsp; ++k) {
            // Back-substitute actual X, including exact zero. The backend
            // regularization remains only in its collision-system assembly,
            // perturbing the coefficient matrix by O(N*Cantera::Tiny).
            dt[k] = 1.6 / Cantera::GasConstant * m_mw[k] * x[k] * m_a[k];
            sum += dt[k];
        }
        const auto dependent = static_cast<std::size_t>(std::max_element(y.begin(), y.end()) - y.begin());
        long double others = 0;
        for (std::size_t k = 0; k < m_nsp; ++k) if (k != dependent) {
            // Transform from the regularized to the actual mass-average frame.
            dt[k] -= y[k] * static_cast<double>(sum);
            others += dt[k];
        }
        dt[dependent] = -static_cast<double>(others);
    }

protected:
    void update_T() override {
        const double previous = m_temp;
        Cantera::MultiTransport::update_T();
        // A viscosity-only visit to B followed by a return to A invalidates
        // binary diffusion even when the old thermal-cache key is still A.
        // updateDiff_T then overwrites the special kinetic self-diffusion
        // diagonal. Invalidate the thermal key on EVERY transport-temperature
        // change so updateThermal_T rebuilds that diagonal and its heat modes.
        // Otherwise transport(A), viscosity(B), transport(A) changes lambda/D_T.
        if (m_temp != previous) m_thermal_tlast = Cantera::Undef;
    }

    void solveLMatrixEquation() override {
        updateThermal_T(); update_C();
        if (m_lmatrix_soln_ok) return;
        const auto n = m_nsp, size = 3 * n;
        m_Lmatrix.resize(size, size, 0.0);
        std::fill(m_Lmatrix.data().begin(), m_Lmatrix.data().end(), 0.0);
        eval_L0010(m_molefracs.data()); eval_L0001(); eval_L1000();
        eval_L1010(m_molefracs.data()); eval_L1001(m_molefracs.data());
        eval_L0100(); eval_L0110(); eval_L0101(m_molefracs.data());
        std::vector<double> rhs(size);
        const double prefactor = 16 * m_temp / 25;
        for (std::size_t i = 0; i < n; ++i) {
            // L00/X_i after using sum(X_j W_j a0_j)=0:
            // (16 T/25) sum_{j!=i} X_j (a0_j-a0_i)/B_ij.
            // Avoid computing the original rank-one term and cancelling it.
            m_Lmatrix(i, i) = 0;
            for (std::size_t j = 0; j < n; ++j) if (i != j) {
                const double v = prefactor * m_molefracs[j] / m_bdiff(i, j);
                m_Lmatrix(i, j) = v; m_Lmatrix(i, i) -= v;
            }
            for (std::size_t j = n; j < size; ++j) m_Lmatrix(i, j) /= m_molefracs[i];
            for (std::size_t j = 0; j < size; ++j) m_Lmatrix(n + i, j) /= m_molefracs[i];
            rhs[n + i] = 1;
            if (hasInternalModes(i)) {
                for (std::size_t j = 0; j < size; ++j) m_Lmatrix(2 * n + i, j) /= m_molefracs[i];
                rhs[2 * n + i] = 1;
            }
            // MultiTransport::thermalConductivity uses the ORIGINAL RHS.
            m_b[i] = 0; m_b[n + i] = m_molefracs[i];
            m_b[2 * n + i] = hasInternalModes(i) ? m_molefracs[i] : 0;
        }
        const auto dependent = static_cast<std::size_t>(std::max_element(m_molefracs.begin(), m_molefracs.end()) - m_molefracs.begin());
        for (std::size_t j = 0; j < size; ++j)
            m_Lmatrix(dependent, j) = j < n ? m_molefracs[j] * m_mw[j] : 0;
        // Row equilibration preserves the kinetic equations. All unknowns
        // remain finite as any individual X approaches zero.
        for (std::size_t i = 0; i < size; ++i) {
            double scale = 0;
            for (std::size_t j = 0; j < size; ++j) scale = std::max(scale, std::abs(m_Lmatrix(i, j)));
            if (!(scale > 0) || !std::isfinite(scale)) throw std::runtime_error("invalid thermal transport row");
            for (std::size_t j = 0; j < size; ++j) m_Lmatrix(i, j) /= scale;
            rhs[i] /= scale;
        }
        const auto original = m_Lmatrix;
        m_a = rhs;
        if (Cantera::solve(m_Lmatrix, m_a.data()) != 0) throw std::runtime_error("thermal transport solve failed");
        const double allowance = 64 * static_cast<double>(size + 1) * std::numeric_limits<double>::epsilon();
        for (std::size_t i = 0; i < size; ++i) {
            long double error = -rhs[i], activity = std::abs(rhs[i]);
            for (std::size_t j = 0; j < size; ++j) {
                const long double term = static_cast<long double>(original(i, j)) * m_a[j];
                error += term; activity += std::abs(term);
            }
            if (!std::isfinite(m_a[i]) || std::abs(error) > allowance * activity)
                throw std::runtime_error("thermal transport linear backward error");
        }
        m_lmatrix_soln_ok = true; m_molefracs_last = m_molefracs; m_l0000_ok = false;
    }
};
} // namespace cartmesh2d::chemistry
