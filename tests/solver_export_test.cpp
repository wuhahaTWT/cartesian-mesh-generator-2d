#include "cartmesh2d/io/OpenFoam2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace cartmesh2d;

namespace {
int failures=0;
void check(bool condition,const std::string& message) {
    if (!condition) { ++failures; std::cerr<<"FAIL: "<<message<<'\n'; }
}

CutCell2D fullCell(std::size_t id,const AABB2D& box) {
    CutCell2D cell;
    cell.sourceId=id;
    cell.sourceKey=id;
    cell.backgroundBounds=box;
    cell.kind=CutCellKind::Full;
    cell.fluidPolygon={{{box.min.x,box.min.y},{box.max.x,box.min.y},
                        {box.max.x,box.max.y},{box.min.x,box.max.y}}};
    cell.area=cell.fluidPolygon.area();
    cell.areaFraction=1.0;
    cell.centroid=cell.fluidPolygon.centroid();
    return cell;
}

std::string readText(const std::filesystem::path& path) {
    std::ifstream in(path);
    std::ostringstream text;
    text<<in.rdbuf();
    return text.str();
}
} // namespace

int main() {
    const Domain2D domain{{{0.0,0.0},{2.0,1.0}}};
    const BoundaryLoop embeddedReference({{0.4,0.4},{0.6,0.4},{0.6,0.6},{0.4,0.6}});
    std::vector<CutCell2D> cells{
        fullCell(0,{{0.0,0.0},{1.0,1.0}}),
        fullCell(1,{{1.0,0.0},{2.0,1.0}})
    };
    const auto topology=buildGlobalTopology(cells,domain,embeddedReference);
    check(topology.valid(),"two-prism export fixture topology is valid");

    const auto quality=evaluateSolverQuality2D(topology);
    check(quality.valid(),"orthogonal two-cell fixture passes solver-quality gate");
    check(quality.maxNonOrthogonalityDeg<=1.0e-12 &&
          quality.maxInternalSkewness<=1.0e-12,
          "orthogonal fixture has zero non-orthogonality and skewness");
    check(std::abs(quality.minFaceWeight-0.5)<=1.0e-12 &&
          std::abs(quality.minVolumeRatio-1.0)<=1.0e-12,
          "orthogonal equal-volume fixture matches OpenFOAM face metrics");

    const auto caseDir=std::filesystem::temp_directory_path()/"cartmesh2d-s1-openfoam-fixture";
    std::string error;
    const BoundaryRegion2D boundaryRegion(embeddedReference);
    const auto report=writeExtrudedOpenFoam2D(topology,domain,boundaryRegion,caseDir,0.1,&error);
    check(report.valid(),"OpenFOAM extrusion writer succeeds: "+error);
    check(report.pointCount==12 && report.faceCount==11 &&
          report.internalFaceCount==1 && report.cellCount==2,
          "OpenFOAM extrusion counts match two adjacent prisms");

    const auto polyMesh=caseDir/"constant"/"polyMesh";
    for (const char* name:{"points","faces","owner","neighbour","boundary"}) {
        check(std::filesystem::is_regular_file(polyMesh/name),
              std::string("OpenFOAM file exists: ")+name);
    }
    for (const auto& relative:{std::filesystem::path("system/controlDict"),
                              std::filesystem::path("system/fvSchemes"),
                              std::filesystem::path("system/fvSolution"),
                              std::filesystem::path("constant/transportProperties"),
                              std::filesystem::path("constant/turbulenceProperties"),
                              std::filesystem::path("0/U"),std::filesystem::path("0/p")}) {
        check(std::filesystem::is_regular_file(caseDir/relative),
              "runnable laminar OpenFOAM fixture exists: "+relative.string());
    }
    const std::string boundary=readText(polyMesh/"boundary");
    check(boundary.find("frontAndBack")!=std::string::npos &&
          boundary.find("type empty;")!=std::string::npos,
          "extruded 2D case declares frontAndBack empty patch");
    check(boundary.find("left")!=std::string::npos &&
          boundary.find("right")!=std::string::npos &&
          boundary.find("bottom")!=std::string::npos &&
          boundary.find("top")!=std::string::npos,
          "domain sides retain deterministic boundary patch names");
    check(readText(polyMesh/"neighbour").find("\n1\n(\n1\n)")!=std::string::npos,
          "OpenFOAM neighbour list contains the one internal face");

    SolverQualityPolicy2D strict;
    strict.maxCellAspect=1.5;
    const auto rejected=evaluateSolverQuality2D(topology,strict);
    check(!rejected.valid(),"stricter solver policy rejects measured aspect instead of hiding it");

    // OpenFOAM's allGeometry checks must be represented by the internal gate.
    // A 0.5%-width cell beside a unit cell violates both the 0.05 face-weight
    // and 0.01 neighbouring-volume-ratio policies.
    const Domain2D imbalancedDomain{{{0.0,0.0},{1.005,1.0}}};
    const BoundaryLoop imbalancedReference({{0.2,0.2},{0.3,0.2},{0.3,0.3},{0.2,0.3}});
    const auto imbalanced=buildGlobalTopology({
        fullCell(0,{{0.0,0.0},{0.005,1.0}}),
        fullCell(1,{{0.005,0.0},{1.005,1.0}})
    },imbalancedDomain,imbalancedReference);
    const auto imbalancedQuality=evaluateSolverQuality2D(imbalanced);
    const auto hasIssue=[&](SolverQualityIssueCode2D code) {
        return std::any_of(imbalancedQuality.issues.begin(),imbalancedQuality.issues.end(),
                           [&](const SolverQualityIssue2D& issue) { return issue.code==code; });
    };
    check(hasIssue(SolverQualityIssueCode2D::LowFaceWeight),
          "OpenFOAM-equivalent low interpolation weight is fail-closed");
    check(hasIssue(SolverQualityIssueCode2D::LowVolumeRatio),
          "OpenFOAM-equivalent low neighbouring volume ratio is fail-closed");

    // Minimal cell retained from the 128-segment circle regression that
    // OpenFOAM 2606 reported at skewness 5.7043075454659755.  The very short
    // wall fragment exposes why boundary skewness must use the normal
    // owner-to-face distance, rather than the full owner-to-face distance.
    CutCell2D boundarySkewCell;
    boundarySkewCell.sourceId=0;
    boundarySkewCell.sourceKey=0;
    boundarySkewCell.backgroundBounds={{-0.890625,-0.515625},
                                       {-0.85676465214179398,-0.46875}};
    boundarySkewCell.kind=CutCellKind::Cut;
    boundarySkewCell.fluidPolygon={{{-0.890625,-0.515625},
                                     {-0.85676465214179398,-0.515625},
                                     {-0.85772861000027223,-0.51410274419322155},
                                     {-0.890625,-0.46875}}};
    boundarySkewCell.area=boundarySkewCell.fluidPolygon.area();
    boundarySkewCell.areaFraction=
        boundarySkewCell.area/
        ((boundarySkewCell.backgroundBounds.max.x-boundarySkewCell.backgroundBounds.min.x)*
         (boundarySkewCell.backgroundBounds.max.y-boundarySkewCell.backgroundBounds.min.y));
    boundarySkewCell.centroid=boundarySkewCell.fluidPolygon.centroid();
    const BoundaryLoop skewBoundary(boundarySkewCell.fluidPolygon.vertices);
    const auto boundarySkewTopology=buildGlobalTopology(
        {boundarySkewCell},{boundarySkewCell.backgroundBounds},skewBoundary);
    const auto boundarySkew=evaluateSolverQuality2D(boundarySkewTopology);
    check(std::any_of(boundarySkew.issues.begin(),boundarySkew.issues.end(),
                      [](const SolverQualityIssue2D& issue) {
                          return issue.code==SolverQualityIssueCode2D::ExcessiveBoundarySkewness;
                      }),
          "OpenFOAM-equivalent boundary-face skewness is fail-closed");
    check(std::abs(boundarySkew.maxBoundarySkewness-5.7043075454659755)<1.0e-11,
          "boundary skewness matches the retained OpenFOAM 2606 regression value");

    // Retained NACA0012 leading-edge regression. Convex partitioning the
    // right source Cut-cell alone forces the shared A-F interface onto a tiny
    // triangle (OpenFOAM face weight 0.0352839). Quality-driven agglomeration
    // must remove that source interface before repartitioning.
    const Point2D nacaLeftBottom{-0.038281249999999982,-0.039687523669164407};
    const Point2D nacaA{0.0031250000000000167,-0.039687523669164407};
    const Point2D nacaB{0.044531250000000015,-0.039687523669164407};
    const Point2D nacaC{0.044531250000000015,-0.031574243588545645};
    const Point2D nacaD{0.021529832133895588,-0.024414813345507057};
    const Point2D nacaE{0.0054117450176094928,-0.012689512036677586};
    const Point2D nacaF{0.0031250000000000167,-0.0073275302117124094};
    const Point2D nacaNose{0.0,0.0};
    const Point2D nacaLeftTop{-0.038281249999999982,0.0};
    const auto makeRegressionCell=[](std::size_t id,std::vector<Point2D> vertices) {
        CutCell2D cell;
        cell.sourceId=id;
        cell.sourceKey=id;
        cell.fluidPolygon={std::move(vertices)};
        cell.backgroundBounds=cell.fluidPolygon.bounds();
        cell.kind=CutCellKind::Cut;
        cell.area=cell.fluidPolygon.area();
        cell.areaFraction=0.5;
        cell.centroid=cell.fluidPolygon.centroid();
        return cell;
    };
    const auto nacaRegressionLeft=makeRegressionCell(
        0,{nacaLeftBottom,nacaA,nacaF,nacaNose,nacaLeftTop});
    const auto nacaRegressionRight=makeRegressionCell(
        1,{nacaA,nacaB,nacaC,nacaD,nacaE,nacaF});
    const BoundaryLoop nacaRegressionBoundary(
        {nacaLeftBottom,nacaB,nacaC,nacaD,nacaE,nacaF,nacaNose,nacaLeftTop});
    const BoundaryRegion2D nacaRegressionRegion(nacaRegressionBoundary);
    const Domain2D nacaRegressionDomain{{nacaLeftBottom,{nacaB.x,0.0}}};
    const auto nacaRegressionTopology=buildGlobalTopology(
        {nacaRegressionLeft,nacaRegressionRight},nacaRegressionDomain,
        nacaRegressionBoundary);
    const auto nacaRegressionSolver=buildSolverTopology2D(
        nacaRegressionTopology,nacaRegressionDomain,nacaRegressionRegion);
    const auto nacaRegressionQuality=evaluateSolverQuality2D(
        nacaRegressionSolver.topology);
    check(nacaRegressionTopology.valid() && nacaRegressionSolver.valid() &&
              nacaRegressionSolver.qualityAgglomeratedSourceCellCount==1 &&
              nacaRegressionQuality.valid(),
          "NACA leading-edge low-weight partition is repaired by source agglomeration");
    double minimumSolverTurn=1.0;
    for (const auto& cell:nacaRegressionSolver.topology.cells) {
        for (std::size_t i=0;i<cell.vertices.size();++i) {
            const auto n=cell.vertices.size();
            const Point2D& previous=nacaRegressionSolver.topology.vertices[
                cell.vertices[(i+n-1)%n]].point;
            const Point2D& current=nacaRegressionSolver.topology.vertices[
                cell.vertices[i]].point;
            const Point2D& next=nacaRegressionSolver.topology.vertices[
                cell.vertices[(i+1)%n]].point;
            const Vector2D incoming=current-previous;
            const Vector2D outgoing=next-current;
            const double scale=std::sqrt(squaredNorm(incoming)*squaredNorm(outgoing));
            minimumSolverTurn=std::min(minimumSolverTurn,
                cross(incoming,outgoing)/(scale+std::numeric_limits<double>::min()));
        }
    }
    check(minimumSolverTurn>1.0e-10,
          "near-collinear prism corners are explicitly partitioned for OpenFOAM");

    // This 2:1-style hanging-face fixture has a short shared face but a long
    // cell-centre connector. OpenFOAM normalises internal skewness by at least
    // 0.2*|d|, not by face length alone.
    CutCell2D coarse;
    coarse.sourceId=0;
    coarse.sourceKey=0;
    coarse.backgroundBounds={{0.0,0.0},{1.0,1.0}};
    coarse.kind=CutCellKind::Full;
    coarse.fluidPolygon={{{0.0,0.0},{1.0,0.0},{1.0,0.02},{1.0,1.0},{0.0,1.0}}};
    coarse.area=coarse.fluidPolygon.area();
    coarse.areaFraction=1.0;
    coarse.centroid=coarse.fluidPolygon.centroid();
    const auto hanging=buildGlobalTopology({coarse,
                                            fullCell(1,{{1.0,0.0},{2.0,0.02}}),
                                            fullCell(2,{{1.0,0.02},{2.0,1.0}})},
                                           {{{0.0,0.0},{2.0,1.0}}},embeddedReference);
    const auto hangingQuality=evaluateSolverQuality2D(hanging);
    check(hanging.valid() && hangingQuality.maxInternalSkewness<4.0,
          "short hanging face uses OpenFOAM connector-based skew normalisation");

    CutCell2D transition;
    transition.sourceId=0;
    transition.sourceKey=0;
    transition.backgroundBounds=domain.bounds;
    transition.kind=CutCellKind::Full;
    transition.fluidPolygon={{{0.0,0.0},{1.0,0.0},{2.0,0.0},{2.0,1.0},{0.0,1.0}}};
    transition.area=2.0;
    transition.areaFraction=1.0;
    transition.centroid=transition.fluidPolygon.centroid();
    const auto transitionTopology=buildGlobalTopology({transition},domain,boundaryRegion);
    const auto repaired=buildSolverTopology2D(transitionTopology,domain,boundaryRegion);
    check(repaired.valid() && repaired.partitionedCellCount==1,
          "solver topology partitions a collinear transition cell deterministically");
    for (const auto& cell:repaired.topology.cells) {
        bool strictConvex=true;
        for (std::size_t i=0;i<cell.vertices.size();++i) {
            const auto n=cell.vertices.size();
            const auto& a=repaired.topology.vertices[cell.vertices[(i+n-1)%n]].point;
            const auto& b=repaired.topology.vertices[cell.vertices[i]].point;
            const auto& c=repaired.topology.vertices[cell.vertices[(i+1)%n]].point;
            strictConvex=strictConvex && orientationSign(a,b,c)>0;
        }
        check(strictConvex,"solver partition emits strictly convex cells");
    }

    const Domain2D cornerDomain{{{0.75,-0.0659375},{0.775,-0.053125}}};
    const BoundaryLoop cornerSolid({{0.775,-0.053125},{0.75,-0.0575},{0.75,-0.053125}});
    CutCell2D corner;
    corner.sourceId=0;
    corner.sourceKey=0;
    corner.backgroundBounds=cornerDomain.bounds;
    corner.kind=CutCellKind::Cut;
    corner.fluidPolygon={{{0.75,-0.0659375},{0.775,-0.0659375},
                          {0.775,-0.05328125},{0.775,-0.053125},
                          {0.75,-0.0575}}};
    corner.area=corner.fluidPolygon.area();
    corner.areaFraction=corner.area/
        ((cornerDomain.bounds.max.x-cornerDomain.bounds.min.x)*
         (cornerDomain.bounds.max.y-cornerDomain.bounds.min.y));
    corner.centroid=corner.fluidPolygon.centroid();
    const BoundaryRegion2D cornerBoundary(cornerSolid);
    const auto cornerTopology=buildGlobalTopology({corner},cornerDomain,cornerBoundary);
    const auto cornerRepaired=buildSolverTopology2D(
        cornerTopology,cornerDomain,cornerBoundary);
    check(cornerTopology.valid() && cornerRepaired.valid() &&
          evaluateSolverQuality2D(cornerRepaired.topology).valid(),
          "airfoil collinear-transition minimum regression is solver-valid");

    if (failures!=0) {
        std::cerr<<failures<<" solver export test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout<<"cartmesh2d S1 solver export tests passed\n";
    return EXIT_SUCCESS;
}
