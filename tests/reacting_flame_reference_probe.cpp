// Independent Cantera steady BVP with optional algebraic mass closure and
// conservative species/total-enthalpy flux divergences. The original
// discretization remains the default; the extra selector is explicit.
// This verification executable does not link the native transient FV solver.
#include "cantera/onedim.h"
#include "cantera/oneD/DomainFactory.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <stdexcept>
#include <vector>

using namespace Cantera;

class ClosedFlow : public Flow1D {
public:
    ClosedFlow(std::shared_ptr<Solution> solution, std::size_t points,
               std::size_t dependent, bool enabled, bool conservative)
        : Flow1D(solution, "flame", points), dependent_(dependent), enabled_(enabled), conservative_(conservative) {}
    bool capture = false;
    std::vector<double> originalSpeciesResidual;
    std::vector<double> selectedSpeciesResidual, totalSpeciesFlux, totalEnthalpyFlux;
    std::vector<double> originalEnergyResidual, selectedEnergyResidual;
    // Node j owns a dual volume between interval midpoints. One shared
    // flux is used on each interval; positive mass flow uses its left node.
    double speciesFace(const double* x, std::size_t k, std::size_t j) const {
        const double mdot = 0.5*(rho_u(x,j)+rho_u(x,j+1));
        const std::size_t up = mdot >= 0 ? j : j+1;
        return mdot*Y(x,k,up) + flux(k,j);
    }
    double enthalpyFace(const double* x, std::size_t j) const {
        const double mdot = 0.5*(rho_u(x,j)+rho_u(x,j+1));
        const std::size_t up = mdot >= 0 ? j : j+1;
        double value = -m_tcon[j]*(T(x,j+1)-T(x,j))/m_dz[j];
        for (std::size_t k=0; k<m_nsp; ++k) {
            value += mdot*Y(x,k,up)*m_hk(k,up)/m_wt[k];
            value += flux(k,j)*0.5*(m_hk(k,j)+m_hk(k,j+1))/m_wt[k];
        }
        return value;
    }
    // Read the arrays from the completed steady evaluation. Do not reconstruct
    // or re-evaluate transport when auditing the BVP's discrete element budget.
    double diffusiveMassFlux(std::size_t k, std::size_t j) const { return flux(k,j); }
    double chemicalMassProduction(std::size_t k, std::size_t j) const {
        return m_wt[k] * m_wdot(k,j);
    }
protected:
    void evalEnergy(double* x, double* residual, int* diagonal, double rdt,
                    std::size_t jmin, std::size_t jmax) override {
        Flow1D::evalEnergy(x,residual,diagonal,rdt,jmin,jmax);
        if (capture) {
            originalEnergyResidual.resize(m_points);
            for (std::size_t j=0; j<m_points; ++j) originalEnergyResidual[j]=residual[index(c_offset_T,j)];
        }
        if (conservative_) {
            if (!isFree() || m_do_radiation || m_twoPointControl) throw std::runtime_error("conservative prototype requires free adiabatic flow");
            for (std::size_t j=std::max<std::size_t>(1,jmin); j<=std::min(jmax,m_points-2); ++j) {
                if (!m_do_energy[j]) continue;
                const double width=0.5*(z(j+1)-z(j-1));
                // This is a steady BVP. The diagonal pseudo-time term is
                // nonlinear-solver stabilization, not a physical time model.
                // Formation enthalpy is carried by the complete species
                // mixture. A separate heat source would double count it.
                residual[index(c_offset_T,j)] = -(enthalpyFace(x,j)-enthalpyFace(x,j-1))/(width*m_rho[j]*m_cp[j])
                    - rdt*(T(x,j)-T_prev(j));
                diagonal[index(c_offset_T,j)] = 1;
            }
        }
        if (capture) {
            selectedEnergyResidual.resize(m_points);
            totalEnthalpyFlux.resize(m_points-1);
            for (std::size_t j=0; j<m_points; ++j) selectedEnergyResidual[j]=residual[index(c_offset_T,j)];
            for (std::size_t j=0; j<m_points-1; ++j) totalEnthalpyFlux[j]=enthalpyFace(x,j);
        }
    }
    void evalSpecies(double* x, double* residual, int* diagonal, double rdt,
                     std::size_t jmin, std::size_t jmax) override {
        Flow1D::evalSpecies(x, residual, diagonal, rdt, jmin, jmax);
        if (capture) {
            if (rdt != 0 || jmin != 0 || jmax != m_points - 1) {
                throw std::runtime_error("diagnostic requires complete steady evaluation");
            }
            originalSpeciesResidual.resize(m_points * m_nsp);
            for (std::size_t j = 0; j < m_points; ++j) {
                for (std::size_t k = 0; k < m_nsp; ++k) {
                    originalSpeciesResidual[j*m_nsp+k] = residual[index(c_offset_Y+k,j)];
                }
            }
        }
        if (conservative_) {
            if (jmin==0) {
                for (std::size_t k=0; k<m_nsp; ++k) {
                    if (k!=leftExcessSpecies()) residual[index(c_offset_Y+k,0)]=-speciesFace(x,k,0);
                }
            }
            for (std::size_t j=std::max<std::size_t>(1,jmin); j<=std::min(jmax,m_points-2); ++j) {
                const double width=0.5*(z(j+1)-z(j-1));
                for (std::size_t k=0; k<m_nsp; ++k) {
                    residual[index(c_offset_Y+k,j)] = (m_wt[k]*m_wdot(k,j)
                        -(speciesFace(x,k,j)-speciesFace(x,k,j-1))/width)/m_rho[j]
                        -rdt*(Y(x,k,j)-Y_prev(k,j));
                    diagonal[index(c_offset_Y+k,j)] = 1;
                }
            }
        }
        if (capture) {
            selectedSpeciesResidual.resize(m_points*m_nsp);
            totalSpeciesFlux.resize((m_points-1)*m_nsp);
            for (std::size_t j=0; j<m_points; ++j) for (std::size_t k=0; k<m_nsp; ++k)
                selectedSpeciesResidual[j*m_nsp+k]=residual[index(c_offset_Y+k,j)];
            for (std::size_t j=0; j<m_points-1; ++j) for (std::size_t k=0; k<m_nsp; ++k)
                totalSpeciesFlux[j*m_nsp+k]=speciesFace(x,k,j);
        }
        if (!enabled_) return;
        for (std::size_t j = std::max<std::size_t>(1,jmin); j <= std::min(jmax,m_points-2); ++j) {
            // Mass closure replaces one mathematically redundant species
            // equation. All reaction and diffusion terms are still evaluated.
            double remaining = 1;
            for (std::size_t k = 0; k < m_nsp; ++k) {
                if (k != dependent_) remaining -= Y(x,k,j);
            }
            residual[index(c_offset_Y+dependent_,j)] = remaining - Y(x,dependent_,j);
            diagonal[index(c_offset_Y+dependent_,j)] = 0;
        }
    }
private:
    std::size_t dependent_;
    bool enabled_, conservative_;
};

int main(int argc, char** argv) {
    try {
        if (argc != 5 && argc != 6) throw std::invalid_argument("mechanism input output use_constraint[0/1] [conservative[0/1]]");
        const std::filesystem::path output(argv[3]);
        if (std::filesystem::exists(output)) throw std::invalid_argument("output exists");
        const bool constrained = std::string(argv[4]) == "1";
        if (!constrained && std::string(argv[4]) != "0") throw std::invalid_argument("bad constraint selector");
        const bool conservative = argc == 6 && std::string(argv[5]) == "1";
        if (argc == 6 && !conservative && std::string(argv[5]) != "0") throw std::invalid_argument("bad conservative selector");
        auto solution = newSolution(argv[1], "", "multicomponent");
        auto gas = solution->thermo();
        std::ifstream input(argv[2]);
        input.imbue(std::locale::classic());
        std::string magic;
        int version = 0;
        std::size_t points = 0, species = 0, dependent = 0;
        double pressure = 0, temperature = 0, massFlux = 0, fixedTemperature = 0;
        input >> magic >> version >> points >> species >> dependent >> pressure >> temperature >> massFlux >> fixedTemperature;
        // Even one digit and a separator per value must fit in the input.
        // This protects allocation from a corrupt point count, not grid accuracy.
        if (!input || magic != "CM2D_BVP_GUESS" || version != 1 || species != gas->nSpecies()
            || points < 3 || dependent >= species || points > std::filesystem::file_size(argv[2]) / (2 * (species + 3))
            || !std::isfinite(pressure) || pressure <= 0 || !std::isfinite(temperature) || temperature <= 0
            || !std::isfinite(massFlux) || massFlux <= 0 || !std::isfinite(fixedTemperature) || fixedTemperature <= 0) {
            throw std::invalid_argument("invalid BVP guess header");
        }
        for (std::size_t k = 0; k < species; ++k) {
            std::string name; input >> name;
            if (name != gas->speciesName(k)) throw std::invalid_argument("species order mismatch");
        }
        std::vector<double> inletY(species), inletX(species), grid(points), t(points), u(points);
        std::vector<std::vector<double>> y(species, std::vector<double>(points));
        for (double& v : inletY) input >> v;
        for (std::size_t j = 0; j < points; ++j) {
            input >> grid[j] >> t[j] >> u[j];
            for (std::size_t k = 0; k < species; ++k) input >> y[k][j];
        }
        std::string end; input >> end;
        if (!input || end != "END") throw std::invalid_argument("truncated BVP guess");
        std::string trailing;
        if (input >> trailing) throw std::invalid_argument("trailing BVP guess data");
        for (double value : inletY) {
            if (!std::isfinite(value) || value < 0) throw std::invalid_argument("invalid inlet mass fraction");
        }
        for (std::size_t j = 0; j < points; ++j) {
            if (!std::isfinite(grid[j]) || (j && grid[j] <= grid[j-1])
                || !std::isfinite(t[j]) || t[j] <= 0 || !std::isfinite(u[j]) || u[j] <= 0) {
                throw std::invalid_argument("invalid BVP grid, temperature or velocity");
            }
            for (std::size_t k = 0; k < species; ++k) {
                if (!std::isfinite(y[k][j])) throw std::invalid_argument("nonfinite BVP initial guess");
            }
        }
        gas->setMassFractions_NoNorm(inletY.data()); gas->setState_TP(temperature,pressure);
        gas->getMoleFractions(inletX.data());
        auto flow = std::make_shared<ClosedFlow>(solution,points,dependent,constrained,conservative);
        flow->setFreeFlow(); flow->setupGrid(points,grid.data()); flow->setPressure(pressure);
        flow->setTransportModel("multicomponent"); flow->enableSoret(true);
        flow->setSteadyTolerances(1e-9,1e-15); flow->setTransientTolerances(1e-7,1e-15);
        auto inlet = newDomain<Inlet1D>("inlet",solution);
        inlet->setMoleFractions(inletX.data()); inlet->setTemperature(temperature); inlet->setMdot(massFlux);
        auto outlet = newDomain<Outlet1D>("outlet",solution);
        std::vector<std::shared_ptr<Domain1D>> domains{inlet,flow,outlet};
        Sim1D flame(domains);
        flow->setValues("velocity",u); flow->setValues("T",t);
        for (std::size_t k = 0; k < species; ++k) flow->setValues(gas->speciesName(k),y[k]);
        if (flame.setFixedTemperature(fixedTemperature) != 0) throw std::runtime_error("fixed temperature inserted a grid point");
        flow->setEnergyEnabled(true);
        flame.solve(1,false);
        flow->capture = true; flame.eval(0); flow->capture = false;
        if (flow->grid() != grid) throw std::runtime_error("grid changed");
        t = flow->values("T"); u = flow->values("velocity");
        for (std::size_t k = 0; k < species; ++k) y[k] = flow->values(gas->speciesName(k));
        std::ofstream stream(output);
        stream.imbue(std::locale::classic());
        stream << std::setprecision(17) << "{\"constraint_enabled\":" << (constrained?"true":"false")
               << ",\"steady_relative_tolerance\":1e-9,\"steady_absolute_tolerance\":1e-15"
               << ",\"transient_relative_tolerance\":1e-7,\"transient_absolute_tolerance\":1e-15"
               << ",\"transport\":\"multicomponent\",\"soret\":true,\"energy\":true"
               << ",\"dependent_species\":\"" << gas->speciesName(dependent) << "\",\"p\":" << pressure
               << ",\"species\":[";
        for (std::size_t k=0;k<species;++k) stream << (k?",":"") << '"' << gas->speciesName(k) << '"';
        auto array = [&](const auto& a) { stream << '['; for(std::size_t k=0;k<a.size();++k) stream << (k?",":"") << a[k]; stream << ']'; };
        stream << "],\"grid\":"; array(grid); stream << ",\"T\":"; array(t); stream << ",\"u\":"; array(u);
        stream << ",\"Y\":[";
        std::vector<double> row(species), rho(points);
        for (std::size_t j=0;j<points;++j) {
            for (std::size_t k=0;k<species;++k) row[k]=y[k][j];
            if(j) stream << ','; array(row);
            gas->setMassFractions_NoNorm(row.data()); gas->setState_TP(t[j],pressure); rho[j]=gas->density();
        }
        stream << "],\"rho\":"; array(rho);
        stream << ",\"original_species_residual_per_s_interior\":[";
        for(std::size_t j=1;j<points-1;++j) {
            if(j>1) stream << ',';
            for(std::size_t k=0;k<species;++k) row[k]=flow->originalSpeciesResidual[j*species+k];
            array(row);
        }
        stream << "],\"inlet_prescribed_mass_fractions\":";
        for(std::size_t k=0;k<species;++k) row[k]=inlet->massFraction(k);
        array(row);
        stream << ",\"inlet_mass_flux_kg_per_m2_s\":" << inlet->mdot()
               << ",\"diffusive_flux_convention\":\"positive along increasing grid; row j is interval [j,j+1]\""
               << ",\"diffusive_mass_flux_kg_per_m2_s_intervals\":[";
        for(std::size_t j=0;j<points-1;++j) {
            if(j) stream << ',';
            for(std::size_t k=0;k<species;++k) row[k]=flow->diffusiveMassFlux(k,j);
            array(row);
        }
        stream << "],\"chemical_mass_production_kg_per_m3_s_interior\":[";
        for(std::size_t j=1;j<points-1;++j) {
            if(j>1) stream << ',';
            for(std::size_t k=0;k<species;++k) row[k]=flow->chemicalMassProduction(k,j);
            array(row);
        }
        stream << "],\"conservative_flux_form\":" << (conservative?"true":"false")
               << ",\"total_species_flux_kg_per_m2_s_intervals\":[";
        for (std::size_t j=0; j<points-1; ++j) {
            if(j) stream << ',';
            for (std::size_t k=0; k<species; ++k) row[k]=flow->totalSpeciesFlux[j*species+k];
            array(row);
        }
        stream << "],\"total_enthalpy_flux_W_per_m2_intervals\":"; array(flow->totalEnthalpyFlux);
        stream << ",\"selected_species_residual_per_s_interior\":[";
        for (std::size_t j=1; j<points-1; ++j) {
            if(j>1) stream << ',';
            for (std::size_t k=0; k<species; ++k) row[k]=flow->selectedSpeciesResidual[j*species+k];
            array(row);
        }
        // Boundary temperature equations have different units. Export only
        // interior temperature-equation residuals under the K/s label.
        std::vector<double> energyResidual(points-2);
        for (std::size_t j=1; j<points-1; ++j) energyResidual[j-1]=flow->originalEnergyResidual[j];
        stream << "],\"original_energy_residual_K_per_s_interior\":"; array(energyResidual);
        for (std::size_t j=1; j<points-1; ++j) energyResidual[j-1]=flow->selectedEnergyResidual[j];
        stream << ",\"selected_energy_residual_K_per_s_interior\":"; array(energyResidual);
        stream << "}\n";
        if(!stream) throw std::runtime_error("output write failed");
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
