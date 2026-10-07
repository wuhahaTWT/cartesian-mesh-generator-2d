#include "FvTestMesh2D.hpp"
#include "cartmesh2d/fv/Euler2D.hpp"
#include <iostream>
#include <iomanip>
using namespace cartmesh2d;using namespace cartmesh2d::fv;
int main(){try{
    std::cout<<std::setprecision(17)<<"[";
    for(bool warped:{false,true}){
        const auto mesh=fv_test::rectangle(8,8,1,warped);const IdealGas2D gas{1.4,1};
        std::vector<EulerBoundary2D> bc;for(std::size_t f=0;f<mesh.faces.size();++f)if(!mesh.faces[f].neighbour)bc.push_back({f,EulerBoundaryKind2D::Transmissive,{}, {},"open"});
        EulerState2D state;for(const auto& c:mesh.cells)state.cells.push_back(eulerConservative2D({1,.4*c.centre.y,0,1},gas));
        EulerStepper2D solver(mesh,bc,gas);double errors[2]{};std::size_t count=0;
        for(int q=0;q<2;++q){EulerStepControls2D ctl;ctl.order=2;ctl.fluxScheme=EulerFluxScheme2D::Hllc;ctl.faceQuadrature=q?EulerFaceQuadrature2D::Gauss2:EulerFaceQuadrature2D::Midpoint;
            const auto op=solver.spatialSnapshot(state,ctl);
            for(std::size_t f=0;f<mesh.faces.size();++f){const auto& face=mesh.faces[f];if(!face.neighbour)continue;const auto a=mesh.cells[face.owner].centre,b=mesh.cells[*face.neighbour].centre;
                if(a.x<.2||a.x>.8||a.y<.2||a.y>.8||b.x<.2||b.x>.8||b.y<.2||b.y>.8)continue;
                const double exact=(1+.16*(face.centre.y*face.centre.y+face.areaVector.x*face.areaVector.x/12))*face.areaVector.x;
                errors[q]=std::max(errors[q],std::abs(op.faceFlux[f][1]-exact));if(q==0)++count;
            }
        }
        // Absolute integrated momentum flux tolerance for this unit-scale algebraic
        // example: roundoff allowance, not a physical accuracy criterion.
        if(count==0||errors[0]<1e-6||errors[1]>1e-13)throw std::runtime_error("affine shear face quadrature test failed");
        std::cout<<(warped?",":"")<<"{\"warped\":"<<warped<<",\"faces\":"<<count<<",\"midpointError\":"<<errors[0]<<",\"gauss2Error\":"<<errors[1]<<"}";
    }
    std::cout<<"]\n";
    const auto mesh=fv_test::rectangle(8,8,1,true);const IdealGas2D gas{1.4,1};std::vector<EulerBoundary2D> bc;
    for(std::size_t f=0;f<mesh.faces.size();++f)if(!mesh.faces[f].neighbour)bc.push_back({f,EulerBoundaryKind2D::SlipWall,{}, {},"wall"});
    EulerState2D state;for(const auto& cell:mesh.cells)state.cells.push_back(eulerConservative2D(cell.centre.x<.5?EulerPrimitive2D{1,.3,0,1}:EulerPrimitive2D{.2,-.2,0,.2},gas));
    EulerStepper2D solver(mesh,bc,gas);EulerStepControls2D ctl;ctl.order=2;ctl.fluxScheme=EulerFluxScheme2D::Hllc;ctl.faceQuadrature=EulerFaceQuadrature2D::Gauss2;ctl.maximumStep=.0005;
    for(int n=0;n<20;++n){const auto before=state;const auto result=solver.advance(state,ctl);
        if(before.cells!=state.cells||result.minimumDensity<=0||result.minimumPressure<=0||std::abs(result.balanceError[0])+std::abs(result.balanceError[3])>1e-13)throw std::runtime_error("Gauss2 stage positivity/shared balance failed");
        state=result.state;}
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
