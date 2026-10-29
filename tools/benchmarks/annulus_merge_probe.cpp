// Research-only controlled reaggregation of the archived annulus. Uses actual
// neighbouring area ratios, NOT an invented reconstruction of background alpha.
#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/io/OpenFoam2D.hpp"
#include "cartmesh2d/stabilization/SmallCell2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>

using namespace cartmesh2d;
void require(bool ok,const std::string& message) {if(!ok)throw std::runtime_error(message);}
int main(int argc,char** argv) {
  try {
    require(argc==4 || argc==5,"MESH OUTPUT_PREFIX NEIGHBOUR_AREA_RATIO [greedy]");
    const auto input=readCm2dTopology(argv[1]);require(input.valid(),input.error);
    const auto& t=input.topology;const std::string prefix=argv[2];const double threshold=std::stod(argv[3]);
    require(threshold>0 && threshold<1,"invalid ratio");
    require(!std::filesystem::exists(prefix+".metrics.json"),"output exists");
    // Reconstruct only the already-classified oriented boundary loops.
    std::map<std::size_t,std::size_t> next;
    for(const auto& e:t.edges)if(!e.neighbour) {
      const auto& owner=t.cells[e.owner];
      const auto k=std::find(owner.edges.begin(),owner.edges.end(),e.id)-owner.edges.begin();
      require(k<owner.vertices.size(),"missing boundary incidence");
      const auto a=owner.vertices[k],b=owner.vertices[(k+1)%owner.vertices.size()];
      require(next.emplace(a,b).second,"branched boundary");
    }
    std::set<std::size_t> used;std::vector<BoundaryLoop> loops;
    for(const auto& [start,to]:next) {
      if(used.contains(start))continue;
      std::vector<Point2D> points;auto v=start;
      do {require(used.insert(v).second,"boundary revisit");points.push_back(t.vertices[v].point);
        require(next.contains(v),"open boundary");v=next.at(v);
      }while(v!=start);
      loops.emplace_back(points);
    }
    BoundaryRegion2D boundary(loops);require(boundary.normalizeAlternating(),"invalid boundary region");
    const Domain2D domain{{{-1.1,-1.1},{1.1,1.1}}};
    std::vector<double> ratio(t.cells.size(),1.);
    for(const auto& e:t.edges)if(e.neighbour) {
      const auto a=e.owner,b=*e.neighbour;
      ratio[a]=std::min(ratio[a],t.cells[a].geometryArea/t.cells[b].geometryArea);
      ratio[b]=std::min(ratio[b],t.cells[b].geometryArea/t.cells[a].geometryArea);
    }
    SmallCellReport2D analysis;analysis.records.resize(t.cells.size());
    std::size_t selected=0;
    for(std::size_t i=0;i<t.cells.size();++i) {
      auto& r=analysis.records[i];r.topologyCellId=i;r.area=t.cells[i].geometryArea;
      if(ratio[i]>=threshold)continue;
      ++selected;double best=0;std::size_t target=i;
      for(const auto eid:t.cells[i].edges) {
        const auto& e=t.edges[eid];if(!e.neighbour)continue;
        const auto other=e.owner==i?*e.neighbour:e.owner;
        if(ratio[other]<threshold)continue;
        const auto d=t.vertices[e.v1].point-t.vertices[e.v0].point;
        const double length=std::hypot(d.x,d.y);
        if(length>best) {best=length;target=other;}
      }
      require(target!=i,"no nonsmall neighbour for experiment");
      r.status=SmallCellStatus2D::CandidateFound;r.targetTopologyCellId=target;r.sharedEdgeLength=best;
    }
    std::size_t candidatesTried=0,candidatesRejected=0;
    auto merged=agglomerateSmallCells({},t,analysis,domain,boundary);
    if(argc==5) {
      require(std::string(argv[4])=="greedy","unknown mode");
      auto current=t;selected=0;
      while(true) {
        std::vector<std::pair<double,std::size_t>> ordered;
        std::vector<double> fractions(current.cells.size(),1.);
        for(const auto& e:current.edges)if(e.neighbour) {
          const auto a=e.owner,b=*e.neighbour;
          fractions[a]=std::min(fractions[a],current.cells[a].geometryArea/current.cells[b].geometryArea);
          fractions[b]=std::min(fractions[b],current.cells[b].geometryArea/current.cells[a].geometryArea);
        }
        for(std::size_t i=0;i<fractions.size();++i)if(fractions[i]<threshold)ordered.emplace_back(fractions[i],i);
        std::sort(ordered.begin(),ordered.end());bool accepted=false;
        for(const auto& [fraction,i]:ordered) {
          std::vector<std::pair<double,std::size_t>> neighbours;
          for(const auto id:current.cells[i].edges) {
            const auto& e=current.edges[id];if(!e.neighbour)continue;
            const auto other=e.owner==i?*e.neighbour:e.owner;
            if(fractions[other]<threshold)continue;
            const auto d=current.vertices[e.v1].point-current.vertices[e.v0].point;
            neighbours.emplace_back(-std::hypot(d.x,d.y),other);
          }
          std::sort(neighbours.begin(),neighbours.end());
          for(const auto& [negativeLength,other]:neighbours) {
            SmallCellReport2D trial;trial.records.resize(current.cells.size());
            for(std::size_t j=0;j<current.cells.size();++j)trial.records[j].topologyCellId=j;
            auto& r=trial.records[i];r.status=SmallCellStatus2D::CandidateFound;
            r.targetTopologyCellId=other;r.sharedEdgeLength=-negativeLength;
            ++candidatesTried;
            auto candidate=agglomerateSmallCells({},current,trial,domain,boundary);
            if(!candidate.valid() || !evaluateSolverQuality2D(candidate.topology).valid()) {++candidatesRejected;continue;}
            current=std::move(candidate.topology);++selected;accepted=true;break;
          }
          if(accepted)break;
        }
        if(!accepted)break;
      }
      // The final no-op uses the same canonical rebuild as every batch case.
      SmallCellReport2D final;final.records.resize(current.cells.size());
      for(std::size_t i=0;i<current.cells.size();++i)final.records[i].topologyCellId=i;
      merged=agglomerateSmallCells({},current,final,domain,boundary);
      merged.inputCellCount=t.cells.size();merged.mergedSmallCellCount=selected;
      merged.totalAreaBefore=0;for(const auto& c:t.cells)merged.totalAreaBefore+=c.geometryArea;
      merged.areaError=merged.totalAreaAfter-merged.totalAreaBefore;
    }
    require(merged.valid(),merged.issues.empty()?"merge invalid":merged.issues.front().message);
    const auto& out=merged.topology;const auto quality=evaluateSolverQuality2D(out);
    std::string error;require(writeCm2dTopology(out,prefix+".solver.cm2d",&error),error);
    std::ofstream(prefix+".quality.json")<<solverQualityReportToJson(quality);
    double minRatio=1;
    for(const auto& e:out.edges)if(e.neighbour)minRatio=std::min(minRatio,
      std::min(out.cells[e.owner].geometryArea,out.cells[*e.neighbour].geometryArea)/
      std::max(out.cells[e.owner].geometryArea,out.cells[*e.neighbour].geometryArea));
    std::ofstream(prefix+".metrics.json")<<std::setprecision(17)
      <<"{\"selectionMetric\":\"cell area / largest face-neighbour area\",\"threshold\":"<<threshold
      <<",\"selected\":"<<selected<<",\"inputCells\":"<<t.cells.size()<<",\"cells\":"<<out.cells.size()
      <<",\"areaBefore\":"<<merged.totalAreaBefore<<",\"areaAfter\":"<<merged.totalAreaAfter
      <<",\"areaError\":"<<merged.areaError<<",\"minNeighbourAreaRatio\":"<<minRatio
      <<",\"qualityGuardedCandidates\":"<<candidatesTried<<",\"rejectedCandidates\":"<<candidatesRejected
      <<",\"qualityIssues\":"<<quality.issues.size()<<"}\n";
    require(writeExtrudedOpenFoam2D(out,domain,boundary,prefix+".foam",.1,&error).valid(),error);
    std::cout<<"selected="<<selected<<" cells="<<out.cells.size()<<" minRatio="<<minRatio<<" qualityIssues="<<quality.issues.size()<<'\n';
    if(!quality.valid())return 3;
    const auto mesh=fv::makeFvMesh2D(out);
    fv::FlowControls2D controls;controls.scenario="custom";controls.speed=.5;
    controls.boundaryConditions=fv::rotatingAnnulusBoundaryPreset2D(mesh,.5);
    std::ofstream bc(prefix+".boundaries");fv::writeFlowBoundaryConditions2D(bc,mesh,controls);
    return 0;
  }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
