#include "cartmesh2d/io/OpenFoam2D.hpp"
#include "cartmesh2d/quality/PatchLocalQuality2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"
#include "repro/source_halo_circle_patch.hpp"
#include "repro/nozzle_determinant_corner.hpp"
#include "repro/nozzle_nonorth_patch.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
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

CutCell2D polygonCell(std::size_t id,std::vector<Point2D> vertices) {
    CutCell2D cell;
    cell.sourceId=id;cell.sourceKey=id;
    cell.fluidPolygon={std::move(vertices)};
    cell.backgroundBounds=cell.fluidPolygon.bounds();
    cell.kind=CutCellKind::Cut;
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

// Ten-cell one-ring reduction of nozzle r01 at original cell pair 2906/8437.
// The separate rectangle represents the unchanged global compactness minimum.
// It is explicitly part of the measured baseline, not an acceptance override.
void nozzleShortFaceRegression() {
    for (const double scale:{0.001,1.0,1000.0}) {
        for (const double angle:{0.0,0.29670597283903605}) {
            std::vector<CutCell2D> cells;std::vector<bool>immutable,rated;std::vector<double> localH;
            cells.push_back(polygonCell(2885,{{-2.953125,-0.94517625672405037},{-2.9291465205232003,-0.941491481273103},{-2.9297598918682133,-0.9375},{-2.953125,-0.9375}}));
            immutable.push_back(false);rated.push_back(true);localH.push_back(0.0234375);
            cells.push_back(polygonCell(2894,{{-2.953125,-0.9375},{-2.9297598918682133,-0.9375},{-2.9311149786177308,-0.92868184534033438}}));
            immutable.push_back(false);rated.push_back(true);localH.push_back(0.0234375);
            cells.push_back(polygonCell(2896,{{-2.9296875,-0.9284624845517383},{-2.90625,-0.92486084147361325},{-2.90625,-0.9140625}}));
            immutable.push_back(false);rated.push_back(true);localH.push_back(0.0234375);
            cells.push_back(polygonCell(2897,{{-2.90625,-0.9140625},{-2.9296875,-0.9140625},{-2.9311149786177308,-0.92868184534033438},{-2.9296875,-0.9284624845517383}}));
            immutable.push_back(false);rated.push_back(true);localH.push_back(0.0234375);
            cells.push_back(polygonCell(2906,{{-2.90625,-0.92486084147361325},{-2.9061149786177309,-0.92484009272366774},{-2.8828125,-0.92125919839548831},{-2.8828125,-0.9140625},{-2.90625,-0.9140625}}));
            immutable.push_back(false);rated.push_back(true);localH.push_back(0.0234375);
            cells.push_back(polygonCell(2907,{{-2.8828125,-0.92125919839548831},{-2.881114978617731,-0.9209983401070011},{-2.859375,-0.91765755531736337},{-2.859375,-0.9140625},{-2.8828125,-0.9140625}}));
            immutable.push_back(false);rated.push_back(true);localH.push_back(0.0234375);
            cells.push_back(polygonCell(2908,{{-2.90625,-0.9140625},{-2.8828125,-0.9140625},{-2.8828125,-0.890625},{-2.90625,-0.890625}}));
            immutable.push_back(false);rated.push_back(true);localH.push_back(0.0234375);
            cells.push_back(polygonCell(7848,{{-2.9025061387777584,-0.94832442526707694},{-2.9041465205232004,-0.93764972865643637},{-2.9291465205232003,-0.941491481273103},{-2.9275061387777583,-0.95216617788374358}}));
            immutable.push_back(true);rated.push_back(false);localH.push_back(0.025293458900823757);
            cells.push_back(polygonCell(8436,{{-2.8791465205232005,-0.93380797603976973},{-2.881114978617731,-0.9209983401070011},{-2.8828125,-0.92125919839548831},{-2.9061149786177309,-0.92484009272366774},{-2.9041465205232004,-0.93764972865643637}}));
            immutable.push_back(true);rated.push_back(false);localH.push_back(0.025293458900823757);
            cells.push_back(polygonCell(8437,{{-2.9041465205232004,-0.93764972865643637},{-2.9061149786177309,-0.92484009272366774},{-2.90625,-0.92486084147361325},{-2.9296875,-0.9284624845517383},{-2.9311149786177308,-0.92868184534033438},{-2.9297598918682133,-0.9375},{-2.9291465205232003,-0.941491481273103}}));
            immutable.push_back(true);rated.push_back(false);localH.push_back(0.025293458900823757);
            BoundaryLoop outline({{-2.953125,-0.94517625672405037},{-2.9291465205232003,-0.941491481273103},{-2.9275061387777583,-0.95216617788374358},{-2.9025061387777584,-0.94832442526707694},{-2.9041465205232004,-0.93764972865643637},{-2.8791465205232005,-0.93380797603976973},{-2.881114978617731,-0.9209983401070011},{-2.859375,-0.91765755531736337},{-2.859375,-0.9140625},{-2.8828125,-0.9140625},{-2.8828125,-0.890625},{-2.90625,-0.890625},{-2.90625,-0.9140625},{-2.9296875,-0.9140625},{-2.9311149786177308,-0.92868184534033438},{-2.953125,-0.9375}});

            cells.push_back(polygonCell(9000,{{0,0},{1,0},{1,.1},{0,.1}}));
            immutable.push_back(false); rated.push_back(true); localH.push_back(.1);
            BoundaryLoop companion({{0,0},{1,0},{1,.1},{0,.1}});
            const auto transform=[&](Point2D p) {
                return Point2D{scale*(std::cos(angle)*p.x-std::sin(angle)*p.y),
                               scale*(std::sin(angle)*p.x+std::cos(angle)*p.y)};
            };
            for (auto& cell:cells) {
                auto points=cell.fluidPolygon.vertices;
                for (auto& point:points) point=transform(point);
                cell=polygonCell(cell.sourceId,std::move(points));
            }
            const auto transformedLoop=[&](const BoundaryLoop& loop) {
                auto points=loop.vertices();
                for (auto& point:points) point=transform(point);
                return BoundaryLoop(std::move(points));
            };
            BoundaryRegion2D boundary(std::vector<BoundaryLoop>{
                transformedLoop(outline),transformedLoop(companion)});
            for (auto& h:localH) h*=scale;
            const Domain2D domain{{{-6*scale,-6*scale},{6*scale,6*scale}}};
            const auto mesh=buildGlobalTopology(cells,domain,boundary);
            const auto before=evaluateSolverQuality2D(mesh);
            check(mesh.valid() && before.valid(),"nozzle R1 reduced fixture is initially solver-valid");
            const auto repair=repairSolverShortFaces2D(mesh,domain,boundary,immutable,localH,rated,.01);
            check(repair.valid() && repair.accepted && repair.hardFaceCountBefore==1U &&
                  repair.hardFaceCountAfter==0U && repair.localWinnerMatchesGlobalAuthority &&
                  repair.localDeltaMatchesGlobalOracle && repair.patchOutsideStableIdsUnchanged &&
                  repair.candidateGlobalTopologyBuildCount==0U && repair.globalOracleBuildCount==1U,
                  "nozzle R1 rounded support repairs at finite scales and rotation"
                  " scale="+std::to_string(scale)+" angle="+std::to_string(angle)+
                  " detail="+(repair.issues.empty()?std::string("none"):repair.issues.front()));
            if (!repair.accepted) continue;
            double areaBefore=0.0,areaAfter=0.0;
            for (const auto& cell:mesh.cells) areaBefore+=cell.geometryArea;
            for (const auto& cell:repair.topology.cells) areaAfter+=cell.geometryArea;
            check(std::abs(areaBefore-areaAfter)<=1e-10*areaBefore,
                  "nozzle R1 transaction conserves the original fluid area");
            check(std::count(repair.immutableCells.begin(),repair.immutableCells.end(),true)==3,
                  "nozzle R1 keeps every immutable layer cell");
            for (std::size_t i=0;i<cells.size();++i) if (immutable[i]) {
                const auto found=std::find_if(repair.topology.cells.begin(),repair.topology.cells.end(),
                    [&](const auto& cell){return cell.sourceKey==cells[i].sourceKey;});
                check(found!=repair.topology.cells.end(),"immutable source identity survives R1");
                if (found==repair.topology.cells.end()) continue;
                Polygon2D actual;
                for (const auto id:found->vertices) actual.vertices.push_back(repair.topology.vertices[id].point);
                const auto onBoundary=[](Point2D point,const Polygon2D& polygon) {
                    for (std::size_t j=0;j<polygon.vertices.size();++j)
                        if (pointOnSegment(point,{polygon.vertices[j],polygon.vertices[(j+1)%polygon.vertices.size()]})) return true;
                    return false;
                };
                bool same=std::abs(actual.area()-cells[i].area)<=1e-10*cells[i].area;
                for (const auto& point:actual.vertices) same=same && onBoundary(point,cells[i].fluidPolygon);
                for (const auto& point:cells[i].fluidPolygon.vertices) same=same && onBoundary(point,actual);
                check(same,"R1 changes only collinear layer incidences, not layer geometry");
            }
        }
    }
}

void targetNonorthogonalityRegression() {
    auto polygons=nozzleNonorthPolygons();
    const auto outline=nozzleNonorthOutline();
    std::vector<CutCell2D> cells;
    for (const auto& polygon:polygons) cells.push_back(polygonCell(cells.size(),polygon.vertices));
    const BoundaryRegion2D boundary{BoundaryLoop(outline.vertices)};
    const Domain2D domain{boundary.bounds()};
    const auto mesh=buildGlobalTopology(cells,domain,boundary);
    SolverQualityPolicy2D policy;policy.maxNonOrthogonalityDeg=65.0;
    const auto before=evaluateSolverQuality2D(mesh,policy);
    check(!improveSolverForTargetPolicy2D(mesh,domain,boundary,{true},policy).valid(),
          "target repair rejects a truncated immutable mask before indexing");
    const auto repaired=improveSolverForTargetPolicy2D(mesh,domain,boundary,nozzleNonorthLocked(),policy);
    const auto after=evaluateSolverQuality2D(repaired.topology,policy);
    check(mesh.valid() && repaired.valid() && !before.valid() && after.valid(),
          "real nozzle transition face is repaired against target 65 degrees");
    check(repaired.repartitionCount>0U,"target repair changes actual topology");
    double a=0.0,b=0.0;
    for (const auto& cell:mesh.cells) a+=cell.geometryArea;
    for (const auto& cell:repaired.topology.cells) b+=cell.geometryArea;
    check(std::abs(a-b)<1e-9*a,"target repair preserves fluid area");
    const auto frozen=improveSolverForTargetPolicy2D(mesh,domain,boundary,
        std::vector<bool>(cells.size(),true),policy);
    check(frozen.valid() && frozen.repartitionCount==0U &&
          evaluateSolverQuality2D(frozen.topology,policy).issues.size()==before.issues.size(),
          "target repair reports immutable unresolved defects without hiding them");
}

void extrudedDeterminantRegression() {
    check(std::abs(extrudedCellDeterminant2D(Polygon2D{{{0,0},{1,0},{1,1},{0,1}}},1.0)-1.0)<1e-14,
          "all-face determinant is one for a cube");
    for (const double scale:{0.001,1.0,1000.0}) {
        auto polygons=nozzleCornerPolygons();
        auto outline=nozzleCornerOutline();
        for (auto& polygon:polygons) for (auto& point:polygon.vertices) {
            point.x*=scale;point.y*=scale;
        }
        for (auto& point:outline.vertices) {point.x*=scale;point.y*=scale;}
        std::vector<CutCell2D> cells;
        for (const auto& polygon:polygons) cells.push_back(polygonCell(cells.size(),polygon.vertices));
        const BoundaryRegion2D boundary{BoundaryLoop(outline.vertices)};
        const Domain2D domain{boundary.bounds()};
        const auto mesh=buildGlobalTopology(cells,domain,boundary);
        const auto locked=nozzleCornerLocked();
        const double actual=extrudedCellDeterminant2D(polygons[4],.02*scale);
        check(std::abs(actual-0.0009174216495895775)<1e-12,
              "real nozzle corner matches independently measured OpenFOAM all-face determinant across scales");
        SolverQualityPolicy2D policy;policy.maxNonOrthogonalityDeg=65.0;
        const auto repaired=improveSolverExtrudedDeterminant2D(mesh,domain,boundary,locked,.02*scale,policy);
        check(mesh.valid() && repaired.valid() && repaired.repartitionCount>0U &&
              repaired.topology.cells.size()<mesh.cells.size(),
              "actual nozzle corner is repaired through exact union");
        double beforeArea=0.0,afterArea=0.0;
        for (const auto& cell:mesh.cells) beforeArea+=cell.geometryArea;
        for (const auto& cell:repaired.topology.cells) {
            afterArea+=cell.geometryArea;
            Polygon2D polygon;
            for (const auto id:cell.vertices) polygon.vertices.push_back(repaired.topology.vertices[id].point);
            check(extrudedCellDeterminant2D(polygon,.02*scale)>=0.001,
                  "repaired corner keeps all-face determinant above unchanged target");
        }
        check(std::abs(beforeArea-afterArea)<1e-9*beforeArea,
              "corner repair retains fluid area");
        for (std::size_t i=0;i<locked.size();++i) if (locked[i]) {
            const auto found=std::find_if(repaired.topology.cells.begin(),repaired.topology.cells.end(),
                [&](const auto& cell){return cell.sourceLineage==std::vector<std::size_t>{i};});
            bool same=found!=repaired.topology.cells.end();
            if (same) for (const auto vertex:mesh.cells[i].vertices) {
                const auto p=mesh.vertices[vertex].point;
                same=same && std::any_of(found->vertices.begin(),found->vertices.end(),[&](auto id) {
                    const auto q=repaired.topology.vertices[id].point;return p.x==q.x && p.y==q.y;
                });
            }
            check(same,"corner repair retains locked wall geometry and provenance");
        }
        const auto frozen=improveSolverExtrudedDeterminant2D(mesh,domain,boundary,
            std::vector<bool>(cells.size(),true),.02*scale,policy);
        check(frozen.valid() && frozen.repartitionCount==0U,
              "all-locked corner stays unresolved rather than losing wall cells");
    }
}

void directionalRepairRegression() {
    // Eight actual cells around NACA140306 cell108446. Two disconnected
    // unchanged rectangles provide volume ratio .04, representing the full
    // case's still lower .0101983 baseline; no Solver threshold is changed.
    // Without that outside baseline, the eight-cell reduction must reject
    // the local volume-ratio decrease .0718852 -> .0660852.
    for (const double scale:{0.001,1.0,1000.0}) for (const double angle:{0.0,0.29670597283903605}) {
        std::vector<Polygon2D> polygons={{{{0.99974490727337728,-0.00039923769549937834},{1.0000628604999999,0.00084716315351562183},{0.99981489970451842,0.00084716315351562183}}},
            {{{0.99981489970451842,0.00084716315351562183},{1.0000628604999999,0.00084716315351562183},{1.0020160673496092,0.00084716315351562183},{1.0020160673496092,0.0028003700031249967}}},
            {{{1.0000628604999999,0.0012615753090807237},{1.0000838139999999,0.0012572092999999999},{1.0020160673496092,0.0028003700031249967},{1.0000628604999999,0.0028003700031249967}}},
            {{{1.0000838139999999,0.0012572092999999999},{0.99981489970451842,0.00084716315351562183},{1.0020160673496092,0.0028003700031249967}}},
            {{{0.9981096536503905,0.0028003700031249967},{0.99868382814578727,0.0021873431047687948},{1.0000628604999999,0.0012615753090807237},{1.0000628604999999,0.0028003700031249967}}},
            {{{1.0000628604999999,-0.0011060436960937531},{1.0020160673496092,-0.0011060436960937531},{1.0020160673496092,0.00084716315351562183},{1.0000628604999999,0.00084716315351562183}}},
            {{{1.0020160673496092,0.00084716315351562183},{1.0039692741992188,0.00084716315351562183},{1.0039692741992188,0.0028003700031249967},{1.0020160673496092,0.0028003700031249967}}},
            {{{1.0000628604999999,0.0028003700031249967},{1.0020160673496092,0.0028003700031249967},{1.0020160673496092,0.0047535768527343716},{1.0000628604999999,0.0047535768527343716}}}};
        const Polygon2D outline={{{0.9981096536503905,0.0028003700031249967},{0.99868382814578727,0.0021873431047687948},{1.0000628604999999,0.0012615753090807237},{1.0000838139999999,0.0012572092999999999},{0.99981489970451842,0.00084716315351562183},{0.99974490727337728,-0.00039923769549937834},{1.0000628604999999,0.00084716315351562183},{1.0000628604999999,-0.0011060436960937531},{1.0020160673496092,-0.0011060436960937531},{1.0020160673496092,0.00084716315351562183},{1.0039692741992188,0.00084716315351562183},{1.0039692741992188,0.0028003700031249967},{1.0020160673496092,0.0028003700031249967},{1.0020160673496092,0.0047535768527343716},{1.0000628604999999,0.0047535768527343716},{1.0000628604999999,0.0028003700031249967}}};
        polygons.push_back(Polygon2D{{{2,.1},{2.01,.1},{2.01,.11},{2,.11}}});
        polygons.push_back(Polygon2D{{{2.01,.1},{2.02,.1},{2.02,.1004},{2.01,.1004}}});
        const auto transform=[&](Point2D p) {
            return Point2D{scale*(std::cos(angle)*p.x-std::sin(angle)*p.y),
                           scale*(std::sin(angle)*p.x+std::cos(angle)*p.y)};
        };
        std::vector<CutCell2D> cells;
        for (auto polygon:polygons) {
            for (auto& point:polygon.vertices) point=transform(point);
            cells.push_back(polygonCell(cells.size(),polygon.vertices));
        }
        std::vector<Point2D> localLoop;
        for (auto p:outline.vertices) localLoop.push_back(transform(p));
        std::vector<Point2D> outsideLoop;
        for (auto p:std::vector<Point2D>{{2,.1},{2.02,.1},{2.02,.1004},{2.01,.1004},{2.01,.11},{2,.11}})
            outsideLoop.push_back(transform(p));
        const BoundaryRegion2D boundary(std::vector<BoundaryLoop>{BoundaryLoop(localLoop),BoundaryLoop(outsideLoop)});
        const Domain2D domain{boundary.bounds()};
        auto mesh=buildGlobalTopology(cells,domain,boundary);
        std::vector<bool> immutable(cells.size(),true);
        immutable[1]=false;immutable[2]=false;immutable[3]=false;
        const auto repaired=improveSolverDirectionalConnectivity2D(mesh,domain,boundary,immutable);
        check(mesh.valid() && evaluateSolverQuality2D(mesh).valid() && repaired.valid() &&
              repaired.acceptedCount==1U && repaired.topology.cells.size()==9U &&
              repaired.after.failedCells.size()+1U==repaired.before.failedCells.size(),
              "NACA exact union repairs one directional defect under the existing global Solver envelope");
        if (repaired.acceptedCount!=1U) continue;
        check(!repaired.after.valid(),"artificial patch outer cells remain underdetermined and are not reported as a CFD pass");
        const auto merged=std::find_if(repaired.topology.cells.begin(),repaired.topology.cells.end(),
            [](const auto& c){return c.sourceLineage==std::vector<std::size_t>{1,3};});
        check(merged!=repaired.topology.cells.end() &&
              repaired.after.cellDeterminants[merged->id]>=minimumDirectionalDeterminant2D,
              "the repaired cell retains both source identities and passes the actual directional threshold");
        double beforeArea=0,afterArea=0;
        for (const auto& c:mesh.cells) beforeArea+=c.geometryArea;
        for (const auto& c:repaired.topology.cells) afterArea+=c.geometryArea;
        check(std::abs(beforeArea-afterArea)<=1e-10*beforeArea,"directional repair preserves total fluid area");
        for (std::size_t i=0;i<cells.size();++i) if (immutable[i]) {
            const auto found=std::find_if(repaired.topology.cells.begin(),repaired.topology.cells.end(),
                [&](const auto& c){return c.sourceLineage==std::vector<std::size_t>{i};});
            bool unchanged=found!=repaired.topology.cells.end();
            if (unchanged) {
                unchanged=repaired.immutableCells[found->id] && found->vertices.size()==mesh.cells[i].vertices.size();
                for (const auto v:mesh.cells[i].vertices) {
                    const auto p=mesh.vertices[v].point;
                    unchanged=unchanged && std::any_of(found->vertices.begin(),found->vertices.end(),[&](auto w){
                        const auto q=repaired.topology.vertices[w].point;return p.x==q.x && p.y==q.y;});
                }
            }
            check(unchanged,"directional union preserves every immutable polygon and its source lineage");
        }
        const auto frozen=improveSolverDirectionalConnectivity2D(mesh,domain,boundary,std::vector<bool>(cells.size(),true));
        check(frozen.valid() && frozen.acceptedCount==0U && frozen.after.failedCells==frozen.before.failedCells,
              "directional repair does not merge immutable cells or conceal unresolved defects");
        const BoundaryRegion2D localBoundary{BoundaryLoop(localLoop)};
        const auto localMesh=buildGlobalTopology(std::vector<CutCell2D>(cells.begin(),cells.begin()+8),
            Domain2D{localBoundary.bounds()},localBoundary);
        const auto localRepair=improveSolverDirectionalConnectivity2D(localMesh,
            Domain2D{localBoundary.bounds()},localBoundary,std::vector<bool>(immutable.begin(),immutable.begin()+8));
        check(localRepair.valid() && localRepair.acceptedCount==0U,
              "without an outside envelope the local volume-ratio decrease is rejected");
    }
}

} // namespace

int main() {
    targetNonorthogonalityRegression();
    extrudedDeterminantRegression();
    directionalRepairRegression();
    {
        check(directionalDeterminant2D({})==0.0 &&
              directionalDeterminant2D({{1,0}})==0.0 &&
              directionalDeterminant2D({{1,0},{-1,0}})==0.0,
              "absent, single and parallel internal directions are underdetermined");
        check(directionalDeterminant2D({{1,0},{0,1}})==0.125 &&
              directionalDeterminant2D({{1,0},{0,1},{-1,0},{0,-1}})==0.5,
              "directional metric matches analytic orthogonal face tensors");
        check(!std::isfinite(directionalDeterminant2D({{0,0}})),
              "degenerate internal direction is not a connectivity pass");
        for (const double scale:{0.001,1.0,1000.0}) {
            for (const double angle:{0.0,0.29670597283903605}) {
                const auto transform=[&](const Vector2D& v) {
                    return Vector2D{scale*(std::cos(angle)*v.x-std::sin(angle)*v.y),
                                    scale*(std::sin(angle)*v.x+std::cos(angle)*v.y)};
                };
                const Point2D a{1.000083814,.0012572093};
                const Point2D b{.9998148997045184,.0008471631535156218};
                const Point2D c{1.0020160673496092,.0028003700031249967};
                const double measured=directionalDeterminant2D({transform(c-b),transform(c-a)});
                check(std::abs(measured-.0003310395831)<1e-12 &&
                      measured<minimumDirectionalDeterminant2D,
                      "NACA cell 108446 reproduces actual OpenFOAM determinant under finite transforms");
            }
        }
        const Domain2D domain{{{0,0},{2,2}}};
        const BoundaryRegion2D boundary(BoundaryLoop({{0,0},{2,0},{2,2},{0,2}}));
        const auto squareMesh=buildGlobalTopology({
            polygonCell(0,{{0,0},{1,0},{1,1},{0,1}}),
            polygonCell(1,{{1,0},{2,0},{2,1},{1,1}}),
            polygonCell(2,{{0,1},{1,1},{1,2},{0,2}}),
            polygonCell(3,{{1,1},{2,1},{2,2},{1,2}})},domain,boundary);
        const auto connected=evaluateDirectionalConnectivity2D(squareMesh);
        check(connected.valid() && connected.cellDeterminants==std::vector<double>(4,0.125),
              "actual internal face incidence excludes physical boundary directions");
        const auto isolated=buildGlobalTopology({polygonCell(0,{{0,0},{2,0},{2,2},{0,2}})},domain,boundary);
        const auto rankZero=evaluateDirectionalConnectivity2D(isolated);
        check(evaluateSolverQuality2D(isolated).valid() && !rankZero.valid() &&
              rankZero.minimumMeasured==0.0 && rankZero.failedCells==std::vector<std::size_t>{0},
              "directional gate stays separate from legacy Solver policy and reports zero internal faces");
    }
    nozzleShortFaceRegression();
    {
        // Two-cell reduction of the 17-degree rotated NACA failure. The exact
        // physical union is only 0.683 degrees concave; convex-only repair kept
        // a sliver with face weight 0.023898 and volume ratio 0.00855925.
        for (const double scale:{0.001,1.0,1000.0}) {
            const auto point=[&](double x,double y) { return Point2D{x*scale,y*scale}; };
            const auto a=point(.51992455682865402,.078181151383417219);
            const auto b=point(.52341074102393725,.078684924676462265);
            const auto c=point(.5275700288911308,.07954861139755956);
            const auto d=point(.5275700288911308,.085826623445894021);
            const auto e=point(.51992455682865402,.085826623445894021);
            const auto f=point(.51992455682865402,.077917544160927429);
            const BoundaryRegion2D region(BoundaryLoop({a,f,b,c,d,e}));
            const Domain2D localDomain{region.bounds()};
            const auto source=buildGlobalTopology(
                {polygonCell(0,{a,b,c,d,e}),polygonCell(1,{a,f,b})},localDomain,region);
            const auto before=evaluateSolverQuality2D(source);
            check(source.valid() && before.issues.size()==2U &&
                      before.minFaceWeight<before.policy.minFaceWeight &&
                      before.minVolumeRatio<before.policy.minVolumeRatio,
                  "rotated NACA two-cell reduction reproduces weight and volume defects");
            const auto repaired=repartitionSolverTopologyByQuality2D(source,localDomain,region);
            const auto after=evaluateSolverQuality2D(repaired.topology);
            check(repaired.valid() && repaired.topology.cells.size()==1U && after.valid(),
                  "exact concave union repairs sliver within unchanged Solver policy");
            check(after.maxConcavityDeg>.68 && after.maxConcavityDeg<.69,
                  "repair reports real concavity instead of moving physical wall points");
            const auto& output=repaired.topology;
            check(output.cells.size()==1U && output.cells[0].sourceLineage==
                      std::vector<std::size_t>({0U,1U}),
                  "solver agglomeration retains both source identities");
            const double area=source.cells[0].geometryArea+source.cells[1].geometryArea;
            check(output.cells.size()==1U &&
                      std::abs(output.cells[0].geometryArea-area)<=1e-10*area,
                  "sliver area is retained in the exact union");
            const auto same=[](const Point2D& lhs,const Point2D& rhs) {
                return lhs.x==rhs.x && lhs.y==rhs.y;
            };
            std::size_t physicalEdges=0;
            for (const auto& edge:source.edges) {
                if (edge.neighbour) continue;
                ++physicalEdges;
                const auto& x=source.vertices[edge.v0].point;
                const auto& y=source.vertices[edge.v1].point;
                check(std::any_of(output.edges.begin(),output.edges.end(),[&](const auto& target) {
                    const auto& u=output.vertices[target.v0].point;
                    const auto& v=output.vertices[target.v1].point;
                    return !target.neighbour && target.patch==edge.patch &&
                        ((same(x,u) && same(y,v)) || (same(x,v) && same(y,u)));
                }),"every original boundary edge and patch survives agglomeration exactly");
            }
            check(output.edges.size()==physicalEdges,
                  "only the internal sliver interface is removed");
        }
    }
    {
        // One solver face spans two collinear pieces of the same physical wall.
        const Domain2D box{{{-1,-1},{2,2}}};
        const BoundaryLoop wall({{0,0},{0.5,0},{1,0},{1,1},{0,1}});
        auto cell=fullCell(0,{{0,0},{1,1}});
        const auto mesh=buildGlobalTopology({cell},box,wall);
        std::string error;
        const auto output=std::filesystem::temp_directory_path()/"cartmesh2d-collinear-wall-export";
        const auto exported=writeExtrudedOpenFoam2D(mesh,box,BoundaryRegion2D(wall),output,0.1,&error);
        check(exported.valid(),"collinear wall pieces retain a unique export patch: "+error);
    }
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
    const BoundaryRegion2D embeddedRegion(embeddedReference);
    const auto cleanSolverTopology=buildSolverTopology2D(
        topology,domain,embeddedRegion);
    check(cleanSolverTopology.valid() &&
              cleanSolverTopology.profile.candidateTopologyCount==0 &&
              cleanSolverTopology.profile.acceptedSourceRepairs==0 &&
              cleanSolverTopology.profile.acceptedRepartitions==0,
          "quality-clean topology performs no candidate rebuild or repair");
    check(cleanSolverTopology.topology.cells.size()==2U &&
              cleanSolverTopology.topology.cells[0].sourceLineage==
                  std::vector<std::size_t>{0U} &&
              cleanSolverTopology.topology.cells[1].sourceLineage==
                  std::vector<std::size_t>{1U},
          "clean solver topology preserves one-source lineage without geometry lookup");
    // The thin third rectangle is bad, but the first is healthy. Isolated
    // polygon ranking prefers merging 1+2 (near-square), then discovers that
    // its changed centre makes the 0/1 interface fail. Real-neighbour ranking
    // instead merges 2+3 and leaves both healthy sources untouched.
    for (const double scale : {0.001,1.0,1000.0}) {
        std::vector<CutCell2D> strips;
        double x=0.0;
        for (const double width : {0.0529,1.0,0.01,0.3}) {
            strips.push_back(fullCell(strips.size(),{{x,0.0},{x+width*scale,scale}}));
            x+=width*scale;
        }
        const Domain2D stripDomain{{{0.0,0.0},{x,scale}}};
        const BoundaryRegion2D stripBoundary(
            BoundaryLoop({{0.0,0.0},{x,0.0},{x,scale},{0.0,scale}}));
        const auto input=buildGlobalTopology(strips,stripDomain,stripBoundary);
        check(input.valid() && !evaluateSolverQuality2D(input).valid(),
              "source neighbour ranking fixture starts with genuine face-quality failures");
        for (const bool protectFirst : {false,true}) {
            SolverTopologyConstraints2D constraints;
            constraints.immutableInputCells={protectFirst,false,false,false};
            const auto repaired=buildSolverTopology2D(input,stripDomain,stripBoundary,constraints);
            check(repaired.valid() && evaluateSolverQuality2D(repaired.topology).valid(),
                  "source merge with real neighbours passes unchanged full quality gate");
            check(repaired.topology.cells.size()==3 &&
                  repaired.profile.acceptedSourceRepairs==1 &&
                  repaired.profile.candidateTopologyCount==1,
                  "source neighbour ranking avoids a self-induced defect and second global rebuild");
            std::vector<std::vector<std::size_t>> lineage;
            for (const auto& cell:repaired.topology.cells) lineage.push_back(cell.sourceLineage);
            std::sort(lineage.begin(),lineage.end());
            check(lineage==std::vector<std::vector<std::size_t>>{{0},{1},{2,3}},
                  "source neighbour ranking preserves healthy and immutable source polygons");
        }
    }

    for (const auto& fixture:{repro::sourceHaloCirclePatch(),
                              repro::finalRepartitionCirclePatch()}) {
        std::vector<CutCell2D> source;
        for (std::size_t i=0;i<fixture.polygons.size();++i)
            source.push_back(polygonCell(fixture.originalSourceIds[i],fixture.polygons[i].vertices));
        const BoundaryRegion2D patchBoundary(fixture.boundaryLoops);
        const Domain2D patchDomain{patchBoundary.bounds()};
        const auto patch=buildGlobalTopology(source,patchDomain,patchBoundary);
        SolverTopologyConstraints2D constraints;
        constraints.immutableInputCells=fixture.immutable;
        constraints.preserveInputCells=fixture.preserve;
        const auto repaired=buildSolverTopology2D(patch,patchDomain,patchBoundary,constraints);
        check(patch.valid() && repaired.valid() &&
              evaluateSolverQuality2D(repaired.topology).valid(),
              "real circle source-repair patch retains full topology and solver acceptance");
        if (fixture.polygons.size()==21U)
            check(repaired.profile.candidateTopologyCount<30U,
                  "real final-repartition patch avoids exhaustive rebuilding of unchanged convex splits");
        double expectedArea=0.0,actualArea=0.0;
        for (const auto& polygon:fixture.polygons) expectedArea+=polygon.area();
        std::vector<std::size_t> survivingSources;
        for (const auto& cell:repaired.topology.cells) {
            actualArea+=cell.geometryArea;
            survivingSources.insert(survivingSources.end(),cell.sourceLineage.begin(),cell.sourceLineage.end());
        }
        std::sort(survivingSources.begin(),survivingSources.end());
        survivingSources.erase(std::unique(survivingSources.begin(),survivingSources.end()),survivingSources.end());
        auto expectedSources=fixture.originalSourceIds;
        std::sort(expectedSources.begin(),expectedSources.end());
        check(survivingSources==expectedSources &&
              std::abs(actualArea-expectedArea)<=TolerancePolicy{}.relative*expectedArea,
              "real circle patch repair preserves total area and every original source");
    }

    {
        // Independent low-weight pairs require true 2-to-1 unions. A batch
        // may therefore emit fewer cells than it removes; every metadata and
        // source-lineage array must follow the actual replacement count.
        std::vector<CutCell2D> inputCells;
        std::vector<BoundaryLoop> loops;
        for (std::size_t copy=0;copy<3;++copy) {
            const double x=2.0*static_cast<double>(copy);
            inputCells.push_back(fullCell(2*copy,{{x,0},{x+.01,1}}));
            inputCells.push_back(fullCell(2*copy+1,{{x+.01,0},{x+1.01,1}}));
            loops.emplace_back(std::vector<Point2D>{{x,0},{x+1.01,0},{x+1.01,1},{x,1}});
        }
        const BoundaryRegion2D boundary(loops);
        const Domain2D domain{boundary.bounds()};
        const auto input=buildGlobalTopology(inputCells,domain,boundary);
        const auto batch=repartitionSolverTopologyByQuality2D(input,domain,boundary);
        const auto reference=repartitionSolverTopologyByQualitySequentialReference2D(input,domain,boundary);
        check(input.valid() && !evaluateSolverQuality2D(input).valid() &&
              batch.valid() && reference.valid() &&
              evaluateSolverQuality2D(batch.topology).valid() &&
              evaluateSolverQuality2D(reference.topology).valid(),
              "independent union batch and exhaustive reference satisfy the same full solver gate");
        check(batch.topology.cells.size()==3 && batch.immutableCells.size()==3,
              "union batch emits one cell and one protection flag per replacement");
        std::vector<std::vector<std::size_t>> lineage;
        double area=0.0;
        for (const auto& cell:batch.topology.cells) {
            lineage.push_back(cell.sourceLineage);area+=cell.geometryArea;
        }
        std::sort(lineage.begin(),lineage.end());
        check(lineage==std::vector<std::vector<std::size_t>>{{0,1},{2,3},{4,5}} &&
              std::abs(area-3.03)<1e-12,
              "independent union batch preserves all three disjoint regions and their source identities");
    }

    const auto independentPatches=selectIndependentSolverRepairPatches2D(
        {{0,1,2},{2,3},{4,5}});
    check(independentPatches==std::vector<std::size_t>({0,2}),
          "conflicting repair halos cannot be selected in the same batch");

    const double r1Short=0.009;
    const auto r1Topology=buildGlobalTopology({
        polygonCell(0,{{0,0.4},{2,0.4},{2,1},{0,1}}),
        polygonCell(1,{{0,0},{1,0},{r1Short,0.4},{0,0.4}}),
        polygonCell(2,{{1,0},{1,0.4},{r1Short,0.4}}),
        polygonCell(3,{{1,0},{2,0},{2,0.4},{1,0.4}})
    },domain,embeddedReference);
    const std::vector<double> r1LocalH{1,1,1,1};
    const std::vector<bool> r1Rated{false,true,true,true};
    const auto r1Incidence=buildEdgeIncidenceStore2D(r1Topology,0U);
    check(r1Incidence.valid(),"R1 fixture incidence is valid");

    // A scope that covers every cell must reproduce the authoritative global
    // aggregates exactly, otherwise the local evaluator has drifted.
    const std::vector<std::size_t> allCells{0,1,2,3};
    const auto wholeMeshScope=buildPatchLocalScope2D(
        r1Topology,r1Incidence,allCells,r1LocalH,r1Rated);
    const auto wholeMeshLocal=evaluatePatchLocalQuality2D(
        wholeMeshScope.cells,0.01);
    const auto wholeMeshGlobal=evaluateSolverQuality2D(r1Topology);
    check(wholeMeshScope.valid() && wholeMeshLocal.valid() &&
          wholeMeshScope.cells.size()==r1Topology.cells.size() &&
          wholeMeshLocal.ratedCellCount==r1Topology.cells.size() &&
          wholeMeshLocal.ratedFaceCount==r1Topology.edges.size(),
          "whole-mesh patch-local scope covers every cell and face");
    check(wholeMeshLocal.maxNonOrthogonalityDeg==wholeMeshGlobal.maxNonOrthogonalityDeg &&
          wholeMeshLocal.maxInternalSkewness==wholeMeshGlobal.maxInternalSkewness &&
          wholeMeshLocal.maxBoundarySkewness==wholeMeshGlobal.maxBoundarySkewness &&
          wholeMeshLocal.maxConcavityDeg==wholeMeshGlobal.maxConcavityDeg &&
          wholeMeshLocal.maxCellAspect==wholeMeshGlobal.maxCellAspect &&
          wholeMeshLocal.minInteriorAngleDeg==wholeMeshGlobal.minInteriorAngleDeg &&
          wholeMeshLocal.minFaceLength==wholeMeshGlobal.minFaceLength &&
          wholeMeshLocal.minFaceWeight==wholeMeshGlobal.minFaceWeight &&
          wholeMeshLocal.minVolumeRatio==wholeMeshGlobal.minVolumeRatio &&
          wholeMeshLocal.minCompactness==wholeMeshGlobal.minCompactness &&
          wholeMeshLocal.issueCount==wholeMeshGlobal.issues.size(),
          "patch-local quality equals authoritative global quality on the full scope");

    // A patch cell whose neighbour is dropped cannot be evaluated: the halo is
    // what makes the local metrics authoritative, so this must fail closed.
    auto truncated=wholeMeshScope.cells;
    truncated.pop_back();
    check(!evaluatePatchLocalQuality2D(truncated,0.01).valid(),
          "patch-local quality fails closed when a patch face loses its neighbour");

    const auto r1Repair=repairSolverShortFaces2D(
        r1Topology,domain,BoundaryRegion2D(embeddedReference),
        {true,false,false,false},r1LocalH,r1Rated,0.01);
    const auto immovableRepair=repairSolverShortFaces2D(
        r1Topology,domain,BoundaryRegion2D(embeddedReference),
        {true,true,true,true},r1LocalH,r1Rated,0.01);
    check(!immovableRepair.valid() && !immovableRepair.accepted &&
          immovableRepair.applicable && immovableRepair.affectedCells.size()==2U &&
          immovableRepair.candidateCount==0U &&
          immovableRepair.topology.cells.size()==r1Topology.cells.size() &&
          immovableRepair.immutableCells==std::vector<bool>({true,true,true,true}),
          "R1 immutable failure retains rejected topology and affected cell identities");
    check(r1Repair.valid() && r1Repair.accepted &&
          r1Repair.candidateGlobalTopologyBuildCount==0U &&
          r1Repair.candidateFullGlobalQualityEvaluationCount==0U &&
          r1Repair.globalOracleBuildCount==1U &&
          r1Repair.localCandidateCount>0U &&
          r1Repair.localQualityEvaluationCount==2U*r1Repair.localCandidateCount &&
          r1Repair.authoritativeFullQualityEvaluationCount==2U &&
          r1Repair.localWinnerMatchesGlobalAuthority &&
          r1Repair.patchOutsideStableIdsUnchanged &&
          r1Repair.localDeltaMatchesGlobalOracle &&
          r1Repair.hardFaceCountBefore==1U && r1Repair.hardFaceCountAfter==0U &&
          r1Repair.minimumFaceOverLocalHAfter>=0.01,
          "R1 transaction repairs the Q2-B short face with local candidates and one final oracle"
          " valid="+std::to_string(r1Repair.valid())+
          " accepted="+std::to_string(r1Repair.accepted)+
          " candidates="+std::to_string(r1Repair.candidateCount)+
          " local_candidates="+std::to_string(r1Repair.localCandidateCount)+
          " local_quality="+std::to_string(r1Repair.localQualityEvaluationCount)+
          " candidate_global="+std::to_string(r1Repair.candidateGlobalTopologyBuildCount)+
          " candidate_full_quality="+
              std::to_string(r1Repair.candidateFullGlobalQualityEvaluationCount)+
          " oracle="+std::to_string(r1Repair.globalOracleBuildCount)+
          " authoritative_quality="+
              std::to_string(r1Repair.authoritativeFullQualityEvaluationCount)+
          " before="+std::to_string(r1Repair.hardFaceCountBefore)+
          " after="+std::to_string(r1Repair.hardFaceCountAfter)+
          " min_after="+std::to_string(r1Repair.minimumFaceOverLocalHAfter)+
          " issue="+(r1Repair.issues.empty()?std::string("none"):r1Repair.issues.front()));
    const auto repeatedRepair=repairSolverShortFaces2D(
        r1Topology,domain,BoundaryRegion2D(embeddedReference),
        {true,false,false,false},r1LocalH,r1Rated,0.01);
    std::vector<std::size_t> repairedIds(r1Repair.topology.cells.size());
    std::iota(repairedIds.begin(),repairedIds.end(),0U);
    auto repairedRated=r1Repair.immutableCells;
    for (std::size_t i=0;i<repairedRated.size();++i) repairedRated[i]=!repairedRated[i];
    const auto repairedScope=buildPatchLocalScope2D(r1Repair.topology,
        buildEdgeIncidenceStore2D(r1Repair.topology,0U),repairedIds,
        std::vector<double>(repairedIds.size(),1.0),repairedRated);
    const auto repairedLocal=evaluatePatchLocalQuality2D(repairedScope.cells,0.01);
    check(patchLocalQualityWithinGlobalBaseline2D(repairedLocal,wholeMeshLocal,wholeMeshGlobal),
          "actual R1 winner remains inside the measured global quality baseline");
    auto globalRegression=repairedLocal;
    globalRegression.minCompactness=wholeMeshGlobal.minCompactness-0.01;
    check(!patchLocalQualityWithinGlobalBaseline2D(globalRegression,wholeMeshLocal,wholeMeshGlobal),
          "a worse global compactness cannot be hidden by eliminating a short face");
    globalRegression=repairedLocal;
    globalRegression.maxNonOrthogonalityDeg=wholeMeshGlobal.maxNonOrthogonalityDeg+1.0;
    check(!patchLocalQualityWithinGlobalBaseline2D(globalRegression,wholeMeshLocal,wholeMeshGlobal),
          "a worse global nonorthogonality cannot be hidden by eliminating a short face");
    globalRegression=repairedLocal;
    globalRegression.hardShortFaceCount=wholeMeshLocal.hardShortFaceCount+1U;
    check(!patchLocalQualityWithinGlobalBaseline2D(globalRegression,wholeMeshLocal,wholeMeshGlobal),
          "global quality bounds cannot admit a short-face regression");
    globalRegression=repairedLocal;
    globalRegression.issueCount=1U;
    check(!patchLocalQualityWithinGlobalBaseline2D(globalRegression,wholeMeshLocal,wholeMeshGlobal),
          "a candidate with a solver violation cannot use the global quality envelope");
    check(repeatedRepair.accepted &&
          repeatedRepair.candidateCount==r1Repair.candidateCount &&
          repeatedRepair.localCandidateCount==r1Repair.localCandidateCount &&
          repeatedRepair.minimumFaceOverLocalHAfter==r1Repair.minimumFaceOverLocalHAfter,
          "local winner selection is deterministic across repeated runs");

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

    PatchLocalQuality2D q3Base,q3Candidate,q3ShortRegression;
    q3Base.hardVolumeRatioCount=4U;
    q3Base.hardFaceWeightCount=3U;
    q3Base.maximumVolumeRatioSeverity=4.0;
    q3Base.totalVolumeRatioSeverity=10.0;
    q3Base.maximumFaceWeightSeverity=2.0;
    q3Base.totalFaceWeightSeverity=5.0;
    q3Candidate=q3Base;
    q3Candidate.hardVolumeRatioCount=3U;
    q3Candidate.maximumVolumeRatioSeverity=3.0;
    q3Candidate.totalVolumeRatioSeverity=7.0;
    const auto q3Rank=patchLocalTerminationRank2D(q3Base,q3Candidate,7U,9U);
    const auto q3RepeatedRank=patchLocalTerminationRank2D(q3Base,q3Candidate,7U,9U);
    check(patchLocalTerminationQualityNoWorse2D(q3Candidate,q3Base) &&
          q3Rank.hardViolationDelta==-1 &&
          !patchLocalTerminationRankBetter2D(q3Rank,q3RepeatedRank) &&
          !patchLocalTerminationRankBetter2D(q3RepeatedRank,q3Rank),
          "Q3 ranks hard count before volume/weight severity deterministically");
    q3ShortRegression=q3Candidate;
    q3ShortRegression.hardShortFaceCount=1U;
    q3ShortRegression.maximumShortFaceSeverity=2.0;
    q3ShortRegression.totalShortFaceSeverity=2.0;
    check(!patchLocalTerminationQualityNoWorse2D(q3ShortRegression,q3Base),
          "Q3 rejects a target-quality improvement that creates a hard short face");

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
    check(std::all_of(nacaRegressionSolver.topology.cells.begin(),
                      nacaRegressionSolver.topology.cells.end(),
                      [](const TopologyCell2D& cell) {
                          return cell.sourceLineage==std::vector<std::size_t>({0U,1U});
                      }),
          "source agglomeration propagates the union of original source ids");
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

    // Two disconnected copies of the retained NACA defect form independent
    // closed one-ring patches. They must be repaired in one deterministic
    // batch and followed by one authoritative global quality check.
    const auto shifted=[](Point2D point,double dx) {
        point.x+=dx;
        return point;
    };
    const auto shiftedRegressionCell=[&](std::size_t id,const CutCell2D& source,double dx) {
        std::vector<Point2D> vertices;
        for (const auto& point:source.fluidPolygon.vertices) {
            vertices.push_back(shifted(point,dx));
        }
        return makeRegressionCell(id,std::move(vertices));
    };
    const double copyShift=0.2;
    const auto nacaCopyLeft=shiftedRegressionCell(2,nacaRegressionLeft,copyShift);
    const auto nacaCopyRight=shiftedRegressionCell(3,nacaRegressionRight,copyShift);
    std::vector<Point2D> shiftedBoundary;
    for (const auto& point:nacaRegressionBoundary.vertices()) {
        shiftedBoundary.push_back(shifted(point,copyShift));
    }
    BoundaryRegion2D doubleNacaRegion({nacaRegressionBoundary,
                                       BoundaryLoop(std::move(shiftedBoundary))});
    check(doubleNacaRegion.normalizeAlternating(),
          "independent NACA repair fixture has a normalized boundary region");
    const Domain2D doubleNacaDomain{{nacaLeftBottom,{nacaB.x+copyShift,0.0}}};
    const auto doubleNacaTopology=buildGlobalTopology(
        {nacaRegressionLeft,nacaRegressionRight,nacaCopyLeft,nacaCopyRight},
        doubleNacaDomain,doubleNacaRegion);
    const auto doubleNacaSolver=buildSolverTopology2D(
        doubleNacaTopology,doubleNacaDomain,doubleNacaRegion);
    check(doubleNacaTopology.valid() && doubleNacaSolver.valid() &&
              doubleNacaSolver.profile.acceptedSourceRepairs==2 &&
              doubleNacaSolver.profile.sourceRepairIterations==1 &&
              doubleNacaSolver.profile.candidateTopologyCount==1,
          "independent source defects are accepted in one deterministic batch");

    // Retained M11 leading-edge solver partition. The artificial B-C face
    // has weight 0.0333906. Repartitioning the owner with its Cartesian
    // neighbour changes only artificial diagonals and removes that face.
    const Point2D m11A{-0.007226562499999983,-0.0099218809172911017};
    const Point2D m11B{0.0031250000000000167,-0.0099218809172911017};
    const Point2D m11E{0.0036110489721427955,-0.0099218809172911017};
    const Point2D m11C{0.0013547716606548965,-0.0064540311838725203};
    const Point2D m11Nose{0.0,0.0};
    const Point2D m11D{-0.007226562499999983,0.0};
    const Point2D m11A0{-0.007226562499999983,-0.019843761834582203};
    const Point2D m11B0{0.0031250000000000167,-0.019843761834582203};
    const auto m11Owner=makeRegressionCell(0,{m11A,m11B,m11C,m11D});
    const auto m11Sliver=makeRegressionCell(1,{m11B,m11E,m11C});
    const auto m11NoseCell=makeRegressionCell(2,{m11C,m11Nose,m11D});
    const auto m11Below=makeRegressionCell(3,{m11A0,m11B0,m11B,m11A});
    const BoundaryLoop m11Boundary(
        {m11A0,m11B0,m11B,m11E,m11C,m11Nose,m11D,m11A});
    const BoundaryRegion2D m11Region(m11Boundary);
    const Domain2D m11Domain{{m11A0,{m11E.x,0.0}}};
    const auto m11Topology=buildGlobalTopology(
        {m11Owner,m11Sliver,m11NoseCell,m11Below},m11Domain,m11Region);
    const auto m11InitialQuality=evaluateSolverQuality2D(m11Topology);
    std::vector<std::size_t> m11PatchCells(m11Topology.cells.size());
    std::iota(m11PatchCells.begin(),m11PatchCells.end(),0);
    const auto m11PatchQuality=evaluateSolverQualityPatch2D(
        m11Topology,m11PatchCells);
    check(std::abs(m11PatchQuality.maxNonOrthogonalityDeg-
                   m11InitialQuality.maxNonOrthogonalityDeg)<1.0e-12 &&
              std::abs(m11PatchQuality.maxInternalSkewness-
                       m11InitialQuality.maxInternalSkewness)<1.0e-12 &&
              std::abs(m11PatchQuality.minFaceWeight-
                       m11InitialQuality.minFaceWeight)<1.0e-12 &&
              std::abs(m11PatchQuality.minVolumeRatio-
                       m11InitialQuality.minVolumeRatio)<1.0e-12,
          "local closed-patch quality matches full quality on shared metrics");
    // A proper subset has interfaces to cells outside the patch. Its local
    // boundary must still contain those faces, once each, regardless of the
    // unrelated cells and non-contiguous global edge IDs.
    const BoundaryRegion2D m11PairBoundary(BoundaryLoop({m11A,m11B,m11E,m11C,m11D}));
    const Domain2D m11PairDomain{m11PairBoundary.bounds()};
    const auto m11PairOracle=buildGlobalTopology(
        {m11Owner,m11Sliver},m11PairDomain,m11PairBoundary);
    const auto m11PairQuality=evaluateSolverQuality2D(m11PairOracle);
    const auto m11SubsetQuality=evaluateSolverQualityPatch2D(m11Topology,{0,1});
    check(m11PairOracle.valid() &&
              m11SubsetQuality.issues.size()==m11PairQuality.issues.size() &&
              std::abs(m11SubsetQuality.maxNonOrthogonalityDeg-
                       m11PairQuality.maxNonOrthogonalityDeg)<1.0e-12 &&
              std::abs(m11SubsetQuality.maxInternalSkewness-
                       m11PairQuality.maxInternalSkewness)<1.0e-12 &&
              std::abs(m11SubsetQuality.minFaceWeight-
                       m11PairQuality.minFaceWeight)<1.0e-12,
          "proper patch subset agrees with an independently constructed two-cell domain");
    const auto m11Repartitioned=repartitionSolverTopologyByQuality2D(
        m11Topology,m11Domain,m11Region);
    const auto m11Sequential=
        repartitionSolverTopologyByQualitySequentialReference2D(
            m11Topology,m11Domain,m11Region);
    const auto m11Quality=evaluateSolverQuality2D(m11Repartitioned.topology);
    const auto m11SequentialQuality=evaluateSolverQuality2D(m11Sequential.topology);
    const auto topologyArea=[](const TopologyMesh2D& mesh) {
        double area=0.0;
        for (const auto& cell:mesh.cells) area+=cell.geometryArea;
        return area;
    };
    check(m11Topology.valid() && !m11InitialQuality.valid() &&
              m11Repartitioned.valid() && m11Repartitioned.repartitionCount>0 &&
              m11Quality.valid(),
          "M11 low-weight artificial diagonal is repaired by local repartition");
    check(m11Sequential.valid() && m11SequentialQuality.valid() &&
              std::abs(topologyArea(m11Repartitioned.topology)-
                       topologyArea(m11Sequential.topology))<1.0e-14,
          "batch and retained sequential repair both conserve area and pass quality");

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
    transition.sourceLineage={42U};
    transition.backgroundBounds=domain.bounds;
    transition.kind=CutCellKind::Full;
    transition.fluidPolygon={{{0.0,0.0},{1.0,0.0},{2.0,0.0},{2.0,1.0},{0.0,1.0}}};
    transition.area=2.0;
    transition.areaFraction=1.0;
    transition.centroid=transition.fluidPolygon.centroid();
    const auto transitionTopology=buildGlobalTopology({transition},domain,boundaryRegion);
    const auto repaired=buildSolverTopology2D(transitionTopology,domain,boundaryRegion);
    check(repaired.valid() && repaired.partitionedCellCount==0 && repaired.topology.cells.size()==1 && repaired.topology.cells.front().vertices.size()==5,
          "Cartesian rectangle retains its conformal collinear connection without diagonals");
    for (const auto& cell:repaired.topology.cells) {
        check(cell.sourceLineage==std::vector<std::size_t>{42U},
              "every solver partition child retains its original stable lineage");
        bool strictConvex=true;
        for (std::size_t i=0;i<cell.vertices.size();++i) {
            const auto n=cell.vertices.size();
            const auto& a=repaired.topology.vertices[cell.vertices[(i+n-1)%n]].point;
            const auto& b=repaired.topology.vertices[cell.vertices[i]].point;
            const auto& c=repaired.topology.vertices[cell.vertices[(i+1)%n]].point;
            strictConvex=strictConvex && orientationSign(a,b,c)>=0;
        }
        check(strictConvex,"Cartesian transition remains convex with collinear shared vertices");
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
