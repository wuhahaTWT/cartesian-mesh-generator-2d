// Research adapter: fixed original polygonal annulus, independent GEOS cuts,
// and a mixed-side annular layer study using the existing native strip builder.
// This does not extend the product H4 nested-loop capability.
#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/io/OpenFoam2D.hpp"
#include "cartmesh2d/boundary_layer/BoundaryLayer2D.hpp"
#include "cartmesh2d/quality/SolverTopology2D.hpp"
#include "cartmesh2d/stabilization/SmallCell2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <stdexcept>
#include <chrono>

using namespace cartmesh2d;
void require(bool b,const std::string& s){if(!b)throw std::runtime_error(s);}
std::vector<BoundaryLoop> readLoops(const std::string& file) {
  std::ifstream f(file);require(bool(f),"missing XY");std::string line;
  std::vector<BoundaryLoop> loops;std::vector<Point2D> points;
  const auto flush=[&]{if(!points.empty()){loops.emplace_back(points);points.clear();}};
  while(std::getline(f,line)) {std::istringstream row(line);Point2D p;
    if(row>>p.x>>p.y)points.push_back(p);else flush();}
  flush();return loops;
}
int main(int argc,char** argv) {
 try {
  require(argc>=5 && argc<=7,"native|shape|geos|hybrid|crop|ready XY PREFIX N [POLYGONS|TOTAL_ROWS] [RETAINED_ROWS]");
  const std::string mode=argv[1],prefix=argv[3];const int n=std::stoi(argv[4]);
  require(mode=="native" || mode=="shape" || mode=="geos" || mode=="hybrid" || mode=="crop" || mode=="ready","unknown mode");
  require(!std::filesystem::exists(prefix+".metrics.json"),"output exists");
  auto loops=readLoops(argv[2]);BoundaryRegion2D walls(loops);
  require(walls.normalizeAlternating(),"invalid physical boundary");
  const Domain2D domain{{{-1.3,-1.3},{1.3,1.3}}};
  std::vector<CutCell2D> cells;std::size_t layerCount=0;std::set<std::size_t> protectedLayers,preservedSources;
  BoundaryRegion2D remainder=walls;
  if(mode=="hybrid") {
    std::vector<WallChain2D> chains;const auto depths=walls.nestingDepths();
    for(std::size_t i=0;i<walls.loops().size();++i) {
      // Interior of the outer loop, exterior of the inner loop: fluid on right.
      auto c=makeClosedWallChain2D(walls.loops()[i],i,"wall_"+std::to_string(i),
        depths[i]%2==0?WallFluidRegion2D::Interior:WallFluidRegion2D::Exterior);
      require(c.success(),c.message);chains.push_back(*c.chain);
    }
    const auto retainedRows=argc==7?std::stoul(argv[6]):3;
    LayerParameters2D parameters;parameters.nLayers=argc>=6?std::stoul(argv[5]):5;
    require(retainedRows>=1 && parameters.nLayers>=retainedRows,"invalid retained layer count");
    parameters.thickness=.006;parameters.growthRatio=1.2;
    auto layers=buildBoundaryLayerStrips2D(chains,parameters);
    require(layers.success(),layers.failure.message);
    std::vector<BoundaryLoop> fronts;
    for(const auto& strip:layers.strips) {
      fronts.emplace_back(strip.outerEnvelope());
      for(const auto& q:strip.cells) {
        CutCell2D c;c.kind=CutCellKind::Full;c.sourceId=cells.size();c.sourceKey=c.sourceId;
        for(auto id:q.vertices)c.fluidPolygon.vertices.push_back(strip.vertices[id].point);
        if(c.fluidPolygon.signedArea()<0)std::reverse(c.fluidPolygon.vertices.begin(),c.fluidPolygon.vertices.end());
        c.area=c.fluidPolygon.area();c.areaFraction=1;c.centroid=c.fluidPolygon.centroid();
        c.backgroundBounds=c.fluidPolygon.bounds();
        // Keep the requested near-wall rows immutable. Additional rows
        // are an explicitly reported repairable transition buffer, not extra
        // claimed retained wall layers or a replacement for the product H4.
        if(q.layer<retainedRows)protectedLayers.insert(c.sourceId);
        cells.push_back(c);
      }
    }
    layerCount=cells.size();remainder=BoundaryRegion2D(fronts);
    require(remainder.normalizeAlternating(),"invalid layer remainder");
  }
  if(mode=="geos" || mode=="crop" || mode=="ready") {
    require(argc==6,"GEOS polygon file required");std::ifstream f(argv[5]);std::size_t count=0;
    require(bool(f>>count),"invalid GEOS file");
    for(std::size_t i=0;i<count;++i) {
      CutCell2D c;std::size_t nv=0;int kind=0;
      require(bool(f>>c.sourceId>>kind>>c.backgroundBounds.min.x>>c.backgroundBounds.min.y
        >>c.backgroundBounds.max.x>>c.backgroundBounds.max.y>>nv),"truncated GEOS cell");
      c.sourceKey=c.sourceId;c.kind=kind>0?CutCellKind::Cut:CutCellKind::Full;
      if(kind==-1)protectedLayers.insert(c.sourceId);
      if(kind<0)preservedSources.insert(c.sourceId);
      for(std::size_t k=0;k<nv;++k){Point2D p;require(bool(f>>p.x>>p.y),"truncated polygon");c.fluidPolygon.vertices.push_back(p);}
      if(c.fluidPolygon.signedArea()<0)std::reverse(c.fluidPolygon.vertices.begin(),c.fluidPolygon.vertices.end());
      c.area=c.fluidPolygon.area();c.centroid=c.fluidPolygon.centroid();
      const auto b=c.backgroundBounds;c.areaFraction=c.area/((b.max.x-b.min.x)*(b.max.y-b.min.y));
      cells.push_back(c);
    }
    if(mode=="crop")layerCount=protectedLayers.size();
  } else {
    UniformCartesianGrid2D grid(domain,n,n);const auto classification=classifyGrid(grid,remainder);(void)classification;
    const auto offset=layerCount;
    for(const auto& b:grid.cells()) {
      auto parts=buildCutCells(b.bounds,b.classification,remainder,FluidRegion2D::Interior);
      for(auto& c:parts) {
        require(c.valid(),"native cut failed");if(c.kind==CutCellKind::Empty)continue;
        c.sourceId=offset+b.id;c.sourceKey=c.sourceId;cells.push_back(std::move(c));
      }
    }
  }
  auto raw=buildGlobalTopology(cells,domain,walls);
  if(!raw.valid())for(const auto& issue:raw.issues)std::cerr<<"raw issue "<<int(issue.code)<<" object="<<issue.objectId<<" "<<issue.message<<'\n';
  require(raw.valid(),"raw shared topology invalid");
  std::string error;require(writeCm2dTopology(raw,prefix+".raw.cm2d",&error),error);
  const auto rawQuality=evaluateSolverQuality2D(raw);
  std::ofstream(prefix+".raw-quality.json")<<solverQualityReportToJson(rawQuality);
  TopologyMesh2D topology=raw;std::size_t merged=0,repartitions=0;
  SolverTopologyConstraints2D constraints;
  if(layerCount) {
    constraints.immutableInputCells.assign(raw.cells.size(),false);
    constraints.preserveInputCells.assign(raw.cells.size(),false);
    constraints.inputPolygonOverrides.resize(raw.cells.size());
    std::map<std::size_t,Polygon2D> originals;for(const auto& c:cells)originals.emplace(c.sourceId,c.fluidPolygon);
    for(std::size_t i=0;i<raw.cells.size();++i) {
      const auto source=raw.cells[i].sourceId;
      constraints.immutableInputCells[i]=protectedLayers.contains(source);
      constraints.preserveInputCells[i]=preservedSources.contains(source);
      if(constraints.immutableInputCells[i] || constraints.preserveInputCells[i])constraints.inputPolygonOverrides[i]=originals.at(source);
    }
  } else if(mode!="ready" && mode!="shape") {
    const auto analysis=analyzeSmallCells(cells,raw);require(analysis.valid(),"small cell analysis failed");
    const auto agglomerated=agglomerateSmallCells(cells,raw,analysis,domain,walls);
    require(agglomerated.valid(),"small cell agglomeration failed");
    topology=agglomerated.topology;merged=agglomerated.mergedSmallCellCount;
  }
  SolverTopologyResult2D solver;
  if(mode=="shape") {
    const auto started=std::chrono::steady_clock::now();
    const auto shape=agglomerateSolverByShape2D(topology,domain,walls);
    require(shape.valid(),shape.issues.empty()?"shape build failed":shape.issues.front());
    std::ofstream(prefix+".shape-quality.json")<<solverQualityReportToJson(shape.after);
    std::ofstream(prefix+".shape-run.json")<<std::setprecision(17)
      <<"{\"stopReason\":\""<<shape.stopReason<<"\",\"mergedCells\":"<<shape.mergedCellCount
      <<",\"passes\":"<<shape.passes<<",\"candidateCount\":"<<shape.candidateCount
      <<",\"maximumUnionCellCount\":"<<shape.maximumUnionCellCount
      <<",\"remainingIssues\":"<<shape.after.issues.size()
      <<",\"shapeSeconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
    std::cout<<"shape_stop="<<shape.stopReason<<" merges="<<shape.mergedCellCount
      <<" passes="<<shape.passes<<" candidates="<<shape.candidateCount
      <<" remaining="<<shape.after.issues.size()<<'\n';
    solver.topology=shape.topology;solver.inputCellCount=topology.cells.size();
    solver.outputCellCount=shape.topology.cells.size();merged=shape.mergedCellCount;
  }
  else if(mode=="ready") {solver.topology=topology;solver.inputCellCount=topology.cells.size();solver.outputCellCount=topology.cells.size();}
  else solver=buildSolverTopology2D(topology,domain,walls,constraints);
  require(solver.valid(),solver.issues.empty()?"solver build failed":solver.issues.front());
  topology=std::move(solver.topology);repartitions=solver.qualityRepartitionCount;
  const auto quality=evaluateSolverQuality2D(topology);
  const bool shapeTargetMet=mode!="shape" || evaluateSolverQuality2D(
      topology,SolverShapeAgglomerationPolicy2D{}.target).valid();
  require(writeCm2dTopology(topology,prefix+(shapeTargetMet?".solver.cm2d":".failed.solver.cm2d"),&error),error);
  std::ofstream(prefix+".quality.json")<<solverQualityReportToJson(quality);
  double area=0;for(const auto& c:topology.cells)area+=c.geometryArea;
  std::size_t retainedLayers=std::count(solver.immutableOutputCells.begin(),solver.immutableOutputCells.end(),true);
  std::ofstream(prefix+".metrics.json")<<std::setprecision(17)<<"{\"mode\":\""<<mode<<"\",\"n\":"<<n
    <<",\"rawCells\":"<<raw.cells.size()<<",\"cells\":"<<topology.cells.size()<<",\"area\":"<<area
    <<",\"expectedArea\":"<<walls.area()<<",\"requestedLayerCells\":"<<protectedLayers.size()
    <<",\"transitionBufferSourceCells\":"<<layerCount-protectedLayers.size()
    <<",\"retainedLayerCells\":"<<retainedLayers<<",\"mergedSmallCells\":"<<merged
    <<",\"qualityRepartitions\":"<<repartitions<<",\"qualityIssues\":"<<quality.issues.size()<<"}\n";
  if(!quality.valid()){std::cerr<<"quality rejection "<<quality.issues.size()<<'\n';return 3;}
  if(!shapeTargetMet){std::cerr<<"shape target rejection\n";return 3;}
  require(writeExtrudedOpenFoam2D(topology,domain,walls,prefix+".foam",.1,&error).valid(),error);
  // Ready external grids retain their external patch identities. Their adapter
  // writes boundary data from the actual external polyMesh, not a regular-ring
  // preset that would silently impose a different wall geometry.
  if(mode=="ready")return 0;
  const auto mesh=fv::makeFvMesh2D(topology);fv::FlowControls2D c;c.scenario="custom";c.speed=.5;
  c.boundaryConditions=fv::rotatingAnnulusBoundaryPreset2D(mesh,.5);
  std::ofstream boundary(prefix+".boundaries");fv::writeFlowBoundaryConditions2D(boundary,mesh,c);
  std::cout<<"cells="<<topology.cells.size()<<" layerCells="<<retainedLayers<<" area="<<area<<'\n';return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
