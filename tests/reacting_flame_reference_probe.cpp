// Independent Cantera steady BVP, with optional algebraic mass closure.
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
               std::size_t dependent, bool enabled)
        : Flow1D(solution, "flame", points), dependent_(dependent), enabled_(enabled) {}
    bool capture = false;
    std::vector<double> originalSpeciesResidual;
protected:
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
    bool enabled_;
};

int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::invalid_argument("mechanism input output use_constraint[0/1]");
        const std::filesystem::path output(argv[3]);
        if (std::filesystem::exists(output)) throw std::invalid_argument("output exists");
        const bool constrained = std::string(argv[4]) == "1";
        if (!constrained && std::string(argv[4]) != "0") throw std::invalid_argument("bad constraint selector");
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
        auto flow = std::make_shared<ClosedFlow>(solution,points,dependent,constrained);
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
        stream << "]}\n";
        if(!stream) throw std::runtime_error("output write failed");
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
