#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

// Geometry-only consistency spectrum for a final native FVM mesh.  It applies
// the product gradient and corrected scalar-diffusion geometry to the complete
// quadratic basis.  No accepted flow field, residual or case-specific target
// enters the metric; this is a research diagnostic, not a product quality gate.
int main(int argc,char** argv) {
    using namespace cartmesh2d;
    using namespace cartmesh2d::fv;
    if(argc!=2) throw std::runtime_error("usage: topology-spectrum mesh.solver.cm2d");
    const auto read=readCm2dTopology(argv[1]);
    if(!read.valid()) throw std::runtime_error(read.error);
    const auto mesh=makeFvMesh2D(read.topology);
    const std::vector<bool> fixed(mesh.faces.size(),true);
    const auto stencil=detail::buildFlowGradientStencil2D(mesh,fixed);

    const auto laplacian=[&](int basis) {
        std::vector<double> value(mesh.cells.size()),boundary(mesh.faces.size());
        const auto sample=[&](Point2D p) {
            if(basis==0) return p.x*p.x+p.y*p.y;
            if(basis==1) return p.x*p.x-p.y*p.y;
            return 2*p.x*p.y;
        };
        for(std::size_t i=0;i<mesh.cells.size();++i) value[i]=sample(mesh.cells[i].centre);
        for(std::size_t f=0;f<mesh.faces.size();++f) boundary[f]=sample(mesh.faces[f].centre);
        const auto gradient=stencil.apply(value,boundary);
        std::vector<double> integrated(mesh.cells.size());
        for(const auto& face:mesh.faces) {
            const auto i=face.owner;
            double flux=0;
            if(face.neighbour) {
                const auto j=*face.neighbour;
                const Vector2D g{
                    gradient[i].x*(1-face.neighbourWeight)+gradient[j].x*face.neighbourWeight,
                    gradient[i].y*(1-face.neighbourWeight)+gradient[j].y*face.neighbourWeight};
                flux=face.transmissibility*(value[j]-value[i])+dot(g,face.correction);
                integrated[j]-=flux;
            } else {
                flux=face.transmissibility*(boundary[&face-mesh.faces.data()]-value[i])+
                     dot(gradient[i],face.correction);
            }
            integrated[i]+=flux;
        }
        for(std::size_t i=0;i<integrated.size();++i) integrated[i]/=mesh.cells[i].area;
        return integrated;
    };
    const auto radial=laplacian(0);
    const auto normalDifference=laplacian(1);
    const auto shear=laplacian(2);

    std::cout<<std::setprecision(17)
             <<"cell,x,y,area,quadratic_consistency_error,radial_laplacian,normal_difference_laplacian,shear_laplacian,gradient_condition,minimum_neighbour_area_ratio,maximum_nonorthogonal_ratio,minimum_normal_distance_over_sqrt_area,boundary_faces\n";
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        const auto& cell=mesh.cells[i];
        const auto& row=stencil.rows[i];
        const double lambda=.5*(row.xx+row.yy+std::hypot(row.xx-row.yy,2*row.xy));
        const double condition=lambda*lambda/row.det;
        double minimumAreaRatio=1;
        double maximumNonorthogonal=0;
        double minimumNormalDistance=std::numeric_limits<double>::infinity();
        std::size_t boundaryFaces=0;
        for(const auto id:cell.faces) {
            const auto& face=mesh.faces[id];
            const auto other=face.owner==i?face.neighbour:std::optional<std::size_t>(face.owner);
            if(other) minimumAreaRatio=std::min(minimumAreaRatio,
                std::min(cell.area,mesh.cells[*other].area)/std::max(cell.area,mesh.cells[*other].area));
            else ++boundaryFaces;
            const double faceLength=std::hypot(face.areaVector.x,face.areaVector.y);
            maximumNonorthogonal=std::max(maximumNonorthogonal,
                std::hypot(face.correction.x,face.correction.y)/faceLength);
            const Vector2D outward=face.owner==i?face.areaVector:face.areaVector*(-1);
            const auto destination=(other?mesh.cells[*other].centre:face.centre)-cell.centre;
            const double normalDistance=dot(outward,destination)/faceLength;
            minimumNormalDistance=std::min(minimumNormalDistance,normalDistance/std::sqrt(cell.area));
        }
        // These three basis errors form the trace functional error on the
        // symmetric Hessian space.  The Euclidean norm is invariant under a
        // rigid rotation of the mesh coordinates.
        const double consistency=std::hypot(radial[i]-4.,
            std::hypot(normalDifference[i],shear[i]));
        std::cout<<i<<','<<cell.centre.x<<','<<cell.centre.y<<','<<cell.area<<','
                 <<consistency<<','<<radial[i]<<','<<normalDifference[i]<<','<<shear[i]<<','
                 <<condition<<','<<minimumAreaRatio<<','<<maximumNonorthogonal<<','
                 <<minimumNormalDistance<<','<<boundaryFaces<<'\n';
    }
}
