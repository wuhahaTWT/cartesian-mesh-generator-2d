#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/topology/EdgeIncidence2D.hpp"
#include <fstream>
#include <iostream>
using namespace cartmesh2d;
int main(int argc,char** argv) {
  if(argc!=2) return 2;
  Domain2D domain{{{0,0},{2,2}}};
  BoundaryLoop wall({{.25,.25},{.5,.25},{.5,.5},{.25,.5}});
  auto cell=buildCutCell(AABB2D{{1,1},{1.5,1.5}},CellClass::Outside,wall);
  cell.sourceId=0; cell.sourceKey=0;
  auto duplicate=cell; duplicate.sourceId=1; duplicate.sourceKey=1;
  auto topology=buildGlobalTopology({cell,duplicate},domain,wall);
  auto quality=evaluateSolverQuality2D(topology);
  std::cout<<"overlapping_cells topology_valid="<<topology.valid()
           <<" cells="<<topology.cells.size()<<" issues="<<topology.issues.size()
           <<" solver_quality_valid="<<quality.valid()
           <<" stronger_incidence_valid="<<buildEdgeIncidenceStore2D(topology,0).valid()<<'\n';
  const auto corrupt=std::filesystem::path(argv[1])/"corrupt.cm2d";
  std::ofstream out(corrupt);
  out<<"CM2D 1\nVERTICES 3\n0 0 0\n1 1 0\n2 0 1\n"
       "EDGES 3\n0 0 1 0 -1 2\n1 1 2 0 -1 2\n2 2 0 0 -1 2\n"
       "CELLS 1\n0 0 0 -42 3 0 0 0 3 0 0 0\n"
       "AUDIT 0 0 0 0 0 0 0\nEND\n";
  out.close();
  auto readback=readCm2dTopology(corrupt);
  std::cout<<"corrupt_readback valid="<<readback.valid()<<" error="<<readback.error<<'\n';
  Domain2D tinyDomain{{{0,0},{1e-6,1e-6}}};
  BoundaryLoop tinyEnvelope({{0,0},{1e-6,0},{1e-6,1e-6},{0,1e-6}});
  auto tiny=buildCutCell(tinyDomain.bounds,CellClass::Inside,tinyEnvelope,FluidRegion2D::Interior);
  tiny.area*=2;
  auto badArea=buildGlobalTopology({tiny},tinyDomain,tinyEnvelope);
  std::cout<<"doubled_area topology_valid="<<badArea.valid()
           <<" polygon_area="<<tiny.fluidPolygon.area()<<" declared_area="<<tiny.area<<'\n';
}
