#include "cartmesh2d/io/OpenFoam2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
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

    if (failures!=0) {
        std::cerr<<failures<<" solver export test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout<<"cartmesh2d S1 solver export tests passed\n";
    return EXIT_SUCCESS;
}
