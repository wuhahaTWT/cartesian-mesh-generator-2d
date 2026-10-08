#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

// The archived 3063-cell DFG failure is also a regression for delayed mixing:
// its unstable mode grows before the original first Anderson step (iteration 10).
// This checks equation convergence, not the still-failing DFG physical accuracy.
int main(int argc, char** argv) {
    using namespace cartmesh2d;
    using namespace cartmesh2d::fv;
    try {
        if(argc!=2)throw std::runtime_error("expected extracted fixture directory");
        const std::filesystem::path root=argv[1];
        const auto input=readCm2dTopology(root/"dfg20-l7.solver.cm2d");
        if(!input.valid())throw std::runtime_error(input.error);
        const auto mesh=makeFvMesh2D(input.topology);
        FlowControls2D c;
        c.scenario="custom";c.nu=.001;c.speed=.3;c.tolerance=1e-8;c.maxIterations=2500;
        c.convection=ConvectionScheme2D::LimitedLinearUpwind;
        c.steadyAcceleration=SteadyAcceleration2D::Anderson;
        c.andersonHistory=16;c.andersonStart=1;
        std::ifstream boundaries(root/"dfg20-l7.boundaries");
        c.boundaryConditions=readFlowBoundaryConditions2D(boundaries,mesh,c);
        const auto r=solveIncompressible2D(mesh,c);
        if(!r.converged || r.stopped || r.history.empty() || !r.history.back().strictLinearStep)
            throw std::runtime_error("DFG failed strict final certification");
        const auto& last=r.history.back();
        if(!(last.momentumResidual<c.tolerance && last.velocityChange<c.tolerance &&
             last.pressureChange<c.tolerance && last.continuity<1e-8 && r.globalRelativeImbalance<1e-8))
            throw std::runtime_error("DFG bypassed a convergence gate");
        if(!r.performance.accelerationAccepted)throw std::runtime_error("history was not used");
        std::cout<<"DFG cells="<<mesh.cells.size()<<" iterations="<<r.history.size()
                 <<" momentum="<<last.momentumResidual<<" continuity="<<last.continuity<<'\n';
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';return 1;
    }
}
