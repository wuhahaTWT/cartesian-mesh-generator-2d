// Native 2D EB geometry experiment. AMReX supplies fractions/centroids;
// this adapter only joins the endpoints of its open Cartesian face segments.
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/io/OpenFoam2D.hpp"
#include "cartmesh2d/stabilization/SmallCell2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include "cartmesh2d/fv/FvMesh2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/fv/Incompressible2D.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <tuple>

using namespace cartmesh2d;
constexpr double slope=.5;
void require(bool ok,const std::string& text) {if(!ok)throw std::runtime_error(text);}

int main(int argc,char** argv) {
  // Usage: MODE N INTERCEPT OUTPUT_PREFIX MERGE(0|1); no AMReX input file.
  if(argc!=6) {std::cerr<<"MODE N INTERCEPT OUTPUT_PREFIX MERGE\n";return 2;}
  const std::string mode=argv[1];const int n=std::stoi(argv[2]);
  const double intercept=std::stod(argv[3]);const std::string prefix=argv[4];
  const bool merge=std::stoi(argv[5]);
  int ac=1;char* av[]={argv[0],nullptr};char** ap=av;
  amrex::Initialize(ac,ap,false);
  int status=0;
  try {
    require(n>=8 && intercept>0 && intercept<.5,"invalid geometry");
    require(mode=="native" || mode=="amrex","unknown source");
    require(!std::filesystem::exists(prefix+".metrics.json"),"output already exists");
    const double h=1./n;
    const Domain2D domain{{{0,0},{1,1}}};
    // Explicit interior fluid: above the straight inclined wall.
    const BoundaryLoop wall({{0,intercept},{1,intercept+slope},{1,1},{0,1}});
    std::vector<CutCell2D> cells;
    double momentAreaError=0,momentCentreError=0;
    if(mode=="native") {
      UniformCartesianGrid2D grid(domain,n,n);const auto classification=classifyGrid(grid,wall);
      (void)classification;
      for(const auto& box:grid.cells()) {
        auto c=buildCutCell(box.bounds,box.classification,wall,FluidRegion2D::Interior);
        require(c.valid(),"native cutter rejected cell");
        if(c.kind==CutCellKind::Empty)continue;
        c.sourceId=box.id;c.sourceKey=box.id;cells.push_back(std::move(c));
      }
    } else {
      using namespace amrex;
      // Retain every resolvable EB volume. This is an experiment setting,
      // not a change to AMReX or to the production mesh acceptance policy.
      ParmParse("eb2").add("small_volfrac",0.);
      ParmParse("eb2").add("cover_multiple_cuts",false);
      Box box(IntVect(0),IntVect(n-1));
      RealBox rb({0.,0.},{1.,1.});Geometry geom(box,&rb,0);
      auto shop=EB2::makeShop(EB2::PlaneIF({0.,intercept},{slope,-1.}));
      EB2::Build(shop,geom,0,0,0);
      const auto& level=EB2::IndexSpace::top().getLevel(geom);
      BoxArray ba(box);DistributionMapping dm(ba);
      MultiFab vf(ba,dm,1,0),vc(ba,dm,2,0);
      level.fillVolFrac(vf,geom);level.fillCentroid(vc,geom);
      Array<MultiFab,2> af,fc;
      for(int d=0;d<2;++d) {
        af[d].define(amrex::convert(ba,IntVect::TheDimensionVector(d)),dm,1,0);
        fc[d].define(af[d].boxArray(),dm,1,0);
      }
      level.fillAreaFrac(GetArrOfPtrs(af),geom);
      level.fillFaceCent(GetArrOfPtrs(fc),geom);
      // Symbolic grid-corner / grid-face-endpoint keys give stable identity.
      // Only roundoff-level endpoints (64 eps in grid coordinates) become corners.
      std::map<std::tuple<int,int,int,int>,Point2D> endpointCache;
      const double eps=64*std::numeric_limits<double>::epsilon();
      const auto volume=vf.const_array(0),centroid=vc.const_array(0);
      for(int j=0;j<n;++j)for(int i=0;i<n;++i) {
        const double fraction=volume(i,j,0);
        if(fraction==0)continue;
        std::map<std::tuple<int,int,int,int>,Point2D> points;
        for(int d=0;d<2;++d)for(int side=0;side<2;++side) {
          const int fi=i+(d==0?side:0),fj=j+(d==1?side:0);
          const double a=af[d].const_array(0)(fi,fj,0);
          if(a==0)continue;
          const double c=fc[d].const_array(0)(fi,fj,0);
          for(int end=0;end<2;++end) {
            double t=.5+c+(end?1:-1)*.5*a;
            auto key=std::make_tuple(d+1,fi,fj,end);
            if(std::abs(t)<=eps) {t=0;key={0,fi,fj,0};}
            else if(std::abs(t-1)<=eps) {t=1;key={0,fi+(d==1),fj+(d==0),0};}
            const Point2D p=d==0?Point2D{fi*h,(fj+t)*h}:Point2D{(fi+t)*h,fj*h};
            const auto [it,inserted]=endpointCache.emplace(key,p);
            require(inserted || std::hypot(p.x-it->second.x,p.y-it->second.y)<=eps*h,
                    "inconsistent shared AMReX endpoint");
            points.emplace(key,it->second);
          }
        }
        Point2D centre{};for(const auto& [key,p]:points) {centre.x+=p.x;centre.y+=p.y;}
        centre.x/=points.size();centre.y/=points.size();
        std::vector<Point2D> vertices;for(const auto& [key,p]:points)vertices.push_back(p);
        std::sort(vertices.begin(),vertices.end(),[&](auto a,auto b){
          return std::atan2(a.y-centre.y,a.x-centre.x)<std::atan2(b.y-centre.y,b.x-centre.x);});
        CutCell2D c;c.sourceId=j*n+i;c.sourceKey=c.sourceId;
        c.backgroundBounds={{i*h,j*h},{(i+1)*h,(j+1)*h}};
        c.kind=fraction==1?CutCellKind::Full:CutCellKind::Cut;
        c.fluidPolygon.vertices=vertices;c.area=c.fluidPolygon.area();
        c.areaFraction=fraction;c.centroid=c.fluidPolygon.centroid();
        require(c.centroid.has_value(),"AMReX polygon has no centroid");
        momentAreaError=std::max(momentAreaError,std::abs(c.area/(h*h)-fraction));
        const Point2D ref{(i+.5+centroid(i,j,0,0))*h,(j+.5+centroid(i,j,0,1))*h};
        momentCentreError=std::max(momentCentreError,std::hypot(c.centroid->x-ref.x,c.centroid->y-ref.y)/h);
        for(std::size_t k=0;k<vertices.size();++k) {
          const auto a=vertices[k],b=vertices[(k+1)%vertices.size()];
          if(std::abs(a.y-slope*a.x-intercept)<eps && std::abs(b.y-slope*b.x-intercept)<eps)
            c.embeddedBoundary.push_back({a,b});
        }
        cells.push_back(std::move(c));
      }
      // Adapter consistency budget: 1e-9 of a background area or spacing.
      // Much larger than floating roundoff, far below the smallest sampled cell;
      // used only to reject conversion errors, never as a solver pass tolerance.
      require(momentAreaError<1e-9 && momentCentreError<1e-9,"AMReX moment conversion mismatch");
    }
    auto topology=buildGlobalTopology(cells,domain,wall);
    require(topology.valid(),"imported/native topology invalid");
    const auto raw=topology;
    double minFraction=1;for(const auto& c:cells)minFraction=std::min(minFraction,c.areaFraction);
    std::size_t merged=0;
    if(merge) {
      const auto analysis=analyzeSmallCells(cells,topology);
      require(analysis.valid(),"small cell analysis failed");
      auto result=agglomerateSmallCells(cells,topology,analysis,domain,wall);
      require(result.valid(),"small cell merge failed");
      merged=result.mergedSmallCellCount;topology=std::move(result.topology);
    }
    std::string error;
    require(writeCm2dTopology(raw,prefix+".raw.cm2d",&error),error);
    require(writeCm2dTopology(topology,prefix+".solver.cm2d",&error),error);
    const auto quality=evaluateSolverQuality2D(topology);
    std::ofstream(prefix+".quality.json")<<solverQualityReportToJson(quality);
    double area=0;for(const auto& c:topology.cells)area+=c.geometryArea;
    std::ofstream metrics(prefix+".metrics.json");metrics<<std::setprecision(17)
      <<"{\"source\":\""<<mode<<"\",\"n\":"<<n<<",\"intercept\":"<<intercept
      <<",\"rawCells\":"<<raw.cells.size()<<",\"cells\":"<<topology.cells.size()
      <<",\"rawMinFraction\":"<<minFraction<<",\"mergedSmallCells\":"<<merged
      <<",\"area\":"<<area<<",\"exactArea\":"<<1-intercept-slope/2
      <<",\"momentAreaError\":"<<momentAreaError<<",\"momentCentreError\":"<<momentCentreError
      <<",\"qualityIssues\":"<<quality.issues.size()<<"}\n";
    // Export raw, geometrically valid polygons even when the solver quality gate
    // rejects them, so that checkMesh evidence and the rejection stay distinct.
    const auto foam=writeExtrudedOpenFoam2D(topology,domain,BoundaryRegion2D(wall),prefix+".foam",.1,&error);
    require(foam.valid(),error);
    if(!quality.valid()) {std::cerr<<"Solver quality rejection: "<<quality.issues.size()<<"\n";status=3;}
    else {
      const auto mesh=fv::makeFvMesh2D(topology);
      fv::FlowControls2D controls;controls.scenario="custom";controls.speed=1;
      for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id];if(f.neighbour)continue;
        const double shear=(f.centre.y-slope*f.centre.x-intercept)/(1+slope*slope);
        const Vector2D u{shear,slope*shear};
        fv::FlowBoundaryCondition2D bc;
        bc.face=id;bc.velocity=u;bc.name="outer";
        if(std::abs(f.centre.y-slope*f.centre.x-intercept)<1e-12) {
          bc.kind=fv::FlowBoundaryKind2D::Wall;bc.velocity={};bc.name="inclinedWall";
        } else {
          const bool inlet=dot(u,f.areaVector)<0;
          bc.kind=inlet?fv::FlowBoundaryKind2D::VelocityInlet:fv::FlowBoundaryKind2D::VelocityOutlet;
          bc.name=inlet?"inlet":"velocityOutlet";
        }
        controls.boundaryConditions.push_back(bc);
      }
      std::ofstream bc(prefix+".boundaries");fv::writeFlowBoundaryConditions2D(bc,mesh,controls);
    }
    std::cout<<"cells="<<topology.cells.size()<<" rawMinFraction="<<minFraction<<" merged="<<merged<<" qualityIssues="<<quality.issues.size()<<"\n";
  }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';status=1;}
  amrex::Finalize();return status;
}
