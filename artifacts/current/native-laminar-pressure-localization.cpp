#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Post-process accepted native fields on the product mesh.  This does not
// assemble or solve an independent flow equation.
int main(int argc, char** argv) {
    using namespace cartmesh2d;
    using namespace cartmesh2d::fv;
    if (argc != 4)
        throw std::runtime_error("usage: pressure-localization label mesh.cm2d fields.cells.csv");

    const auto read = readCm2dTopology(argv[2]);
    if (!read.valid()) throw std::runtime_error(read.error);
    const auto mesh = makeFvMesh2D(read.topology);

    std::vector<double> u(mesh.cells.size()), v(u.size()), p(u.size()), exactP(u.size()), error(u.size());
    std::ifstream input(argv[3]);
    if (!input) throw std::runtime_error("cannot read field");
    std::string line;
    std::getline(input, line);
    std::size_t rows = 0;
    while (std::getline(input, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        std::size_t id;
        double x, y, area, speed;
        if (!(row >> id >> x >> y >> area >> u.at(id) >> v.at(id) >> p.at(id) >> speed) || id != rows)
            throw std::runtime_error("invalid field row");
        ++rows;
    }
    if (rows != mesh.cells.size()) throw std::runtime_error("field size mismatch");

    const auto reference = [](Point2D x) {
        const double r2 = x.x*x.x + x.y*x.y;
        return r2/18. - std::log(r2)/9. - 1./(18.*r2);
    };
    double area = 0, gauge = 0;
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        exactP[i] = reference(mesh.cells[i].centre);
        area += mesh.cells[i].area;
        gauge += mesh.cells[i].area * (p[i] - exactP[i]);
    }
    gauge /= area;
    for (std::size_t i = 0; i < error.size(); ++i) error[i] = std::abs(p[i] - exactP[i] - gauge);

    auto weightedQuantile = [&](double fraction) {
        std::vector<std::size_t> order(error.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](auto a, auto b) { return error[a] < error[b]; });
        double cumulative = 0;
        for (auto i : order) {
            cumulative += mesh.cells[i].area;
            if (cumulative >= fraction*area) return error[i];
        }
        return error[order.back()];
    };
    auto concentration = [&](double fraction) {
        std::vector<std::size_t> order(error.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](auto a, auto b) {
            return mesh.cells[a].area*error[a]*error[a] > mesh.cells[b].area*error[b]*error[b];
        });
        double total = 0;
        for (std::size_t i = 0; i < error.size(); ++i)
            total += mesh.cells[i].area*error[i]*error[i];
        double contribution = 0, selectedArea = 0;
        std::size_t count = 0;
        for (auto i : order) {
            contribution += mesh.cells[i].area*error[i]*error[i];
            selectedArea += mesh.cells[i].area;
            ++count;
            if (contribution >= fraction*total) break;
        }
        return std::pair<double,std::size_t>{selectedArea/area,count};
    };

    double l1 = 0, l2 = 0, maximum = 0, wallMaximum = 0, interiorMaximum = 0;
    std::size_t maximumCell = 0, wallCells = 0;
    for (std::size_t i = 0; i < error.size(); ++i) {
        l1 += mesh.cells[i].area*error[i];
        l2 += mesh.cells[i].area*error[i]*error[i];
        if (error[i] > maximum) { maximum = error[i]; maximumCell = i; }
        bool wall = false;
        for (auto face : mesh.cells[i].faces) wall = wall || !mesh.faces[face].neighbour;
        if (wall) { wallMaximum = std::max(wallMaximum,error[i]); ++wallCells; }
        else interiorMaximum = std::max(interiorMaximum,error[i]);
    }

    struct Band { double distance, area=0, square=0, maximum=0; };
    // Distances are metres for the fixed 0.5 m annular gap.  They are
    // diagnostic bands, not acceptance tolerances.
    std::vector<Band> bands{{0.005},{0.01},{0.025},{0.05},{0.10}};
    for (std::size_t i = 0; i < error.size(); ++i) {
        const auto c = mesh.cells[i].centre;
        const double radius = std::hypot(c.x,c.y);
        const double wallDistance = std::min(radius-.5,1.-radius);
        for (auto& band : bands) if (wallDistance >= band.distance) {
            band.area += mesh.cells[i].area;
            band.square += mesh.cells[i].area*error[i]*error[i];
            band.maximum = std::max(band.maximum,error[i]);
        }
    }

    struct Trace { double angle; Vector2D velocity; double length; };
    std::vector<Trace> rotor;
    double normalVelocity = 0;
    for (const auto& bc : rotatingAnnulusBoundaryPreset2D(mesh,.5)) {
        const auto& face = mesh.faces[bc.face];
        const double length = std::hypot(face.areaVector.x,face.areaVector.y);
        normalVelocity = std::max(normalVelocity,
            std::abs(bc.velocity.x*face.areaVector.x+bc.velocity.y*face.areaVector.y)/length);
        if (bc.kind == FlowBoundaryKind2D::SmoothMovingWall)
            rotor.push_back({std::atan2(face.centre.y,face.centre.x),bc.velocity,length});
    }
    std::sort(rotor.begin(),rotor.end(),[](const auto& a,const auto& b){return a.angle<b.angle;});
    double maximumJump = 0, jumpSquare = 0, jumpLength = 0;
    for (std::size_t i = 0; i < rotor.size(); ++i) {
        const auto& a = rotor[i];
        const auto& b = rotor[(i+1)%rotor.size()];
        const double jump = std::hypot(a.velocity.x-b.velocity.x,a.velocity.y-b.velocity.y);
        const double weight = .5*(a.length+b.length);
        maximumJump = std::max(maximumJump,jump);
        jumpSquare += weight*jump*jump;
        jumpLength += weight;
    }

    // Linear patch test on this exact accepted cut-cell mesh.  Rigid-body
    // rotation has identically zero symmetric strain, so the complete product
    // viscous flux should cancel.  The trace is an affine mathematical test,
    // not a proposed annulus boundary condition.
    std::vector<double> patchU(mesh.cells.size()),patchV(mesh.cells.size());
    std::vector<double> boundaryU(mesh.faces.size()),boundaryV(mesh.faces.size());
    std::vector<bool> fixed(mesh.faces.size()),smooth(mesh.faces.size());
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        patchU[i]=-mesh.cells[i].centre.y;
        patchV[i]= mesh.cells[i].centre.x;
    }
    for(std::size_t id=0;id<mesh.faces.size();++id)if(!mesh.faces[id].neighbour) {
        fixed[id]=true;
        boundaryU[id]=-mesh.faces[id].centre.y;
        boundaryV[id]= mesh.faces[id].centre.x;
    }
    const auto patchStencil=detail::buildFlowGradientStencil2D(mesh,fixed);
    const auto patchGu=patchStencil.apply(patchU,boundaryU);
    const auto patchGv=patchStencil.apply(patchV,boundaryV);
    const auto patchTranspose=detail::symmetricViscousCorrection(mesh,patchU,patchV,
        patchGu,patchGv,boundaryU,boundaryV,fixed,fixed,smooth,smooth,.1);
    std::vector<Vector2D> patchResidual(mesh.cells.size());
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& face=mesh.faces[id];
        const auto i=face.owner;
        Vector2D gu=patchGu[i],gv=patchGv[i];
        if(face.neighbour) {
            const auto j=*face.neighbour;
            const double w=face.neighbourWeight;
            gu={gu.x*(1-w)+patchGu[j].x*w,gu.y*(1-w)+patchGu[j].y*w};
            gv={gv.x*(1-w)+patchGv[j].x*w,gv.y*(1-w)+patchGv[j].y*w};
        }
        const double otherU=face.neighbour?patchU[*face.neighbour]:boundaryU[id];
        const double otherV=face.neighbour?patchV[*face.neighbour]:boundaryV[id];
        const Vector2D flux{
            -.1*(face.transmissibility*(otherU-patchU[i])+gu.x*face.correction.x+gu.y*face.correction.y)+patchTranspose[id].x,
            -.1*(face.transmissibility*(otherV-patchV[i])+gv.x*face.correction.x+gv.y*face.correction.y)+patchTranspose[id].y};
        patchResidual[i].x+=flux.x;patchResidual[i].y+=flux.y;
        if(face.neighbour) {
            patchResidual[*face.neighbour].x-=flux.x;
            patchResidual[*face.neighbour].y-=flux.y;
        }
    }
    double patchMaximum=0,patchRms=0;
    for(std::size_t i=0;i<mesh.cells.size();++i) {
        const double acceleration=std::hypot(patchResidual[i].x,patchResidual[i].y)/mesh.cells[i].area;
        patchMaximum=std::max(patchMaximum,acceleration);
        patchRms+=mesh.cells[i].area*acceleration*acceleration;
    }

    const auto c50 = concentration(.5), c90 = concentration(.9);
    const auto& worst = mesh.cells[maximumCell];
    std::cout << std::setprecision(17)
        << "{\"label\":\"" << argv[1] << "\",\"cells\":" << mesh.cells.size()
        << ",\"area_m2\":" << area
        << ",\"pressure_error_m2_s2\":{"
        << "\"l1_area_weighted\":" << l1/area
        << ",\"rms_area_weighted\":" << std::sqrt(l2/area)
        << ",\"maximum\":" << maximum
        << ",\"wall_cell_maximum\":" << wallMaximum
        << ",\"nonwall_cell_maximum\":" << interiorMaximum
        << ",\"q90_area\":" << weightedQuantile(.9)
        << ",\"q99_area\":" << weightedQuantile(.99)
        << ",\"q999_area\":" << weightedQuantile(.999)
        << "},\"worst_cell\":{"
        << "\"id\":" << maximumCell << ",\"x_m\":" << worst.centre.x
        << ",\"y_m\":" << worst.centre.y << ",\"area_m2\":" << worst.area
        << "},\"squared_error_concentration\":{"
        << "\"area_fraction_for_50_percent\":" << c50.first << ",\"cells_for_50_percent\":" << c50.second
        << ",\"area_fraction_for_90_percent\":" << c90.first << ",\"cells_for_90_percent\":" << c90.second
        << "},\"wall_distance_bands\":[";
    for (std::size_t i = 0; i < bands.size(); ++i) {
        const auto& band = bands[i];
        if (i) std::cout << ',';
        std::cout << "{\"minimum_distance_m\":" << band.distance
                  << ",\"retained_area_fraction\":" << band.area/area
                  << ",\"rms_m2_s2\":" << (band.area?std::sqrt(band.square/band.area):0)
                  << ",\"maximum_m2_s2\":" << band.maximum << '}';
    }
    std::cout << "],\"rotor_trace\":{"
              << "\"faces\":" << rotor.size()
              << ",\"maximum_adjacent_velocity_jump_m_s\":" << maximumJump
              << ",\"length_weighted_rms_jump_m_s\":" << std::sqrt(jumpSquare/jumpLength)
              << ",\"maximum_normal_velocity_m_s\":" << normalVelocity
              << "},\"rigid_rotation_symmetric_viscous_patch\":{"
              << "\"rms_acceleration_m_s2\":" << std::sqrt(patchRms/area)
              << ",\"maximum_acceleration_m_s2\":" << patchMaximum
              << "},\"wall_cells\":" << wallCells << "}\n";
}
