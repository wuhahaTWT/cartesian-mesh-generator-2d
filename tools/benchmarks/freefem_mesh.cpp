// Convert an accepted native annulus mesh to a conforming two-dimensional
// FreeFEM mesh. Original faces remain shared edges; only FV centres are added.
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
#include "cartmesh2d/io/MeshIO2D.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cartmesh2d;
using namespace cartmesh2d::fv;

namespace {
void require(bool ok,const std::string& message) {
    if(!ok)throw std::runtime_error("FreeFEM mesh: "+message);
}
std::ofstream output(const std::filesystem::path& path) {
    std::ofstream out(path);
    require(bool(out),"cannot open "+path.string());
    out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    return out;
}
struct Triangle {
    std::size_t a,b,centre,cell,face;
    long double area;
};
int kindTag(FlowBoundaryKind2D kind) {
    switch(kind) {
    case FlowBoundaryKind2D::Wall:return 1;
    case FlowBoundaryKind2D::MovingWall:return 2;
    case FlowBoundaryKind2D::SmoothMovingWall:return 3;
    default:throw std::runtime_error("FreeFEM mesh: only wall, moving-wall and smooth-moving-wall boundaries are supported");
    }
}
}

int main(int argc,char** argv) {
    if(argc!=4) {
        std::cerr<<"Usage: freefem_mesh MESH BOUNDARIES PREFIX\n";
        return 2;
    }
    try {
        const auto input=readCm2dTopology(argv[1]);
        require(input.valid(),"invalid native topology: "+input.error);
        const auto& topology=input.topology;
        const auto mesh=makeFvMesh2D(topology);
        FlowControls2D controls;controls.scenario="custom";
        std::ifstream boundaryInput(argv[2]);
        require(bool(boundaryInput),"cannot open boundary input");
        controls.boundaryConditions=readFlowBoundaryConditions2D(boundaryInput,mesh,controls);
        const auto nc=mesh.cells.size(),nf=mesh.faces.size(),nv=topology.vertices.size();
        const auto maximumIndex=static_cast<std::size_t>(std::numeric_limits<int>::max());
        require(nv<=maximumIndex&&nc<=maximumIndex-nv&&nf<maximumIndex,
                "indices exceed FreeFEM signed-integer range");
        std::vector<Vector2D> velocity(nf);
        std::vector<int> tags(nf);
        for(const auto& condition:controls.boundaryConditions) {
            tags[condition.face]=kindTag(condition.kind);
            velocity[condition.face]=condition.velocity;
        }
        // The EDP uses this Couette benchmark's analytic solution. Reject other
        // data rather than reporting an unrelated circle error; product gates
        // and the original geometry/trace remain unchanged.
        const auto expected=rotatingAnnulusBoundaryPreset2D(mesh,.5);
        const auto tolerance=TolerancePolicy{};
        const double geometryBudget=16*tolerance.scale(2.);
        const double velocityBudget=16*tolerance.scale(.5);
        double xmin=std::numeric_limits<double>::infinity(),ymin=xmin,xmax=-xmin,ymax=-xmin;
        double innerRadius=0,outerRadius=0;
        for(const auto& condition:expected) {
            const auto face=condition.face;
            const auto& f=mesh.faces[face];
            require(tags[face]==kindTag(condition.kind)&&
                    std::hypot(velocity[face].x-condition.velocity.x,velocity[face].y-condition.velocity.y)<=velocityBudget,
                    "reference requires the original annulus preset with inner speed 0.5");
            for(double sign:{-1.,1.}) {
                const double x=f.centre.x+sign*.5*f.areaVector.y;
                const double y=f.centre.y-sign*.5*f.areaVector.x;
                xmin=std::min(xmin,x);xmax=std::max(xmax,x);
                ymin=std::min(ymin,y);ymax=std::max(ymax,y);
                auto& radius=condition.kind==FlowBoundaryKind2D::SmoothMovingWall?innerRadius:outerRadius;
                radius=std::max(radius,std::hypot(x,y));
            }
        }
        require(std::abs(xmin+xmax)<=2*geometryBudget&&std::abs(ymin+ymax)<=2*geometryBudget&&
                std::abs(innerRadius-.5)<=geometryBudget&&std::abs(outerRadius-1.)<=geometryBudget,
                "reference requires centre (0,0), inner circumradius 0.5 and outer circumradius 1");
        std::size_t nb=0;
        for(std::size_t face=0;face<nf;++face) {
            if(mesh.faces[face].neighbour)require(tags[face]==0,"internal face has a boundary tag");
            else {require(tags[face]!=0,"missing physical wall boundary");++nb;}
        }

        std::vector<Triangle> triangles;
        triangles.reserve(2*nf);
        long double originalArea=0,triangleArea=0;
        long double minimumTriangle=std::numeric_limits<long double>::infinity();
        long double maximumCellAreaError=0;
        for(std::size_t cell=0;cell<nc;++cell) {
            const auto& polygon=topology.cells[cell];
            const auto centre=mesh.cells[cell].centre;
            originalArea+=mesh.cells[cell].area;
            long double cellArea=0;
            AABB2D bounds{{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()},
                          {-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()}};
            for(std::size_t edge=0;edge<polygon.vertices.size();++edge) {
                const auto va=polygon.vertices[edge],vb=polygon.vertices[(edge+1)%polygon.vertices.size()];
                const auto a=topology.vertices[va].point,b=topology.vertices[vb].point;
                require(orientationSign(a,b,centre)>0,
                    "cell "+std::to_string(cell)+" cannot be triangulated from its FV centre at face "+std::to_string(polygon.edges[edge]));
                const long double ax=static_cast<long double>(a.x)-centre.x,ay=static_cast<long double>(a.y)-centre.y;
                const long double bx=static_cast<long double>(b.x)-centre.x,by=static_cast<long double>(b.y)-centre.y;
                const long double area=.5L*(ax*by-ay*bx);
                require(std::isfinite(area)&&area>0,"nonpositive or nonfinite fan triangle at cell "+std::to_string(cell));
                triangles.push_back({va,vb,nv+cell,cell,polygon.edges[edge],area});
                cellArea+=area;minimumTriangle=std::min(minimumTriangle,area);
                bounds.min.x=std::min(bounds.min.x,a.x);bounds.min.y=std::min(bounds.min.y,a.y);
                bounds.max.x=std::max(bounds.max.x,a.x);bounds.max.y=std::max(bounds.max.y,a.y);
            }
            const long double error=std::abs(cellArea-mesh.cells[cell].area);
            maximumCellAreaError=std::max(maximumCellAreaError,error);
            const double extent=std::max(bounds.max.x-bounds.min.x,bounds.max.y-bounds.min.y);
            require(error<=tolerance.areaScale(extent),"fan area differs from native cell area at cell "+std::to_string(cell));
            triangleArea+=cellArea;
        }
        require(triangles.size()<=maximumIndex,"triangle count exceeds FreeFEM signed-integer range");

        const std::filesystem::path prefix(argv[3]);
        if(!prefix.parent_path().empty())std::filesystem::create_directories(prefix.parent_path());
        const auto path=[&](const char* suffix){return std::filesystem::path(prefix.string()+suffix);};
        auto medit=output(path(".mesh"));
        medit<<"MeshVersionFormatted 1\nDimension 2\nVertices\n"<<nv+nc<<'\n';
        for(const auto& vertex:topology.vertices)medit<<vertex.point.x<<' '<<vertex.point.y<<" 0\n";
        for(const auto& cell:mesh.cells)medit<<cell.centre.x<<' '<<cell.centre.y<<" 0\n";
        medit<<"Triangles\n"<<triangles.size()<<'\n';
        for(const auto& triangle:triangles)
            medit<<triangle.a+1<<' '<<triangle.b+1<<' '<<triangle.centre+1<<' '<<triangle.cell+1<<'\n';
        // Medit Edges denotes physical boundary entities. Listing internal
        // faces here would create false boundary labels; they already occur
        // twice with the same original endpoints in adjacent fan triangles.
        medit<<"Edges\n"<<nb<<'\n';
        for(std::size_t face=0;face<nf;++face)if(!mesh.faces[face].neighbour) {
            const auto& owner=topology.cells[mesh.faces[face].owner];
            const auto position=static_cast<std::size_t>(std::find(owner.edges.begin(),owner.edges.end(),face)-owner.edges.begin());
            require(position<owner.edges.size(),"boundary face is absent from owner loop");
            medit<<owner.vertices[position]+1<<' '<<owner.vertices[(position+1)%owner.vertices.size()]+1<<' '<<face+1<<'\n';
        }
        medit<<"End\n";
        auto boundary=output(path(".boundary.dat"));boundary<<nf<<'\n';
        auto cells=output(path(".cells.dat"));cells<<nc<<'\n';
        auto faces=output(path(".faces.dat"));faces<<nf<<'\n';
        auto edges=output(path(".edges.dat"));edges<<nf<<'\n';
        for(std::size_t face=0;face<nf;++face) {
            boundary<<velocity[face].x<<' '<<velocity[face].y<<' '<<tags[face]<<'\n';
            const auto& f=mesh.faces[face];const auto& e=topology.edges[face];
            faces<<face<<' '<<f.centre.x<<' '<<f.centre.y<<' '<<f.areaVector.x<<' '<<f.areaVector.y<<'\n';
            edges<<face<<' '<<e.v0<<' '<<e.v1<<' '<<e.owner<<' '
                 <<(e.neighbour?static_cast<long long>(*e.neighbour):-1)<<' '<<(e.neighbour?0:face+1)<<'\n';
        }
        for(std::size_t cell=0;cell<nc;++cell) {
            const auto& c=mesh.cells[cell];cells<<cell<<' '<<c.centre.x<<' '<<c.centre.y<<' '<<c.area<<'\n';
        }
        auto relation=output(path(".cell-triangles.csv"));
        auto samples=output(path(".triangles.dat"));samples<<triangles.size()<<'\n';
        relation<<"cell_id,triangle_id_zero_based,centre_vertex_mesh_id,face_id,start_vertex_mesh_id,end_vertex_mesh_id,triangle_area\n";
        for(std::size_t id=0;id<triangles.size();++id) {
            const auto& t=triangles[id];relation<<t.cell<<','<<id<<','<<t.centre+1<<','<<t.face<<','<<t.a+1<<','<<t.b+1<<','<<t.area<<'\n';
            const auto a=topology.vertices[t.a].point,b=topology.vertices[t.b].point,c=mesh.cells[t.cell].centre;
            const Point2D sample{c.x+((a.x-c.x)+(b.x-c.x))/3.,c.y+((a.y-c.y)+(b.y-c.y))/3.};
            samples<<id<<' '<<t.cell<<' '<<sample.x<<' '<<sample.y<<' '<<t.area<<'\n';
        }
        auto mapping=output(path(".mapping.json"));
        mapping<<"{\n  \"format\": \"cartmesh2d-freefem-fan-v1\",\n"
            <<"  \"benchmark\": {\"centre\": [0,0], \"innerCircumradius\": 0.5, \"outerCircumradius\": 1, \"innerSpeed\": 0.5},\n"
            <<"  \"nc\": "<<nc<<",\n  \"ntri\": "<<triangles.size()<<",\n"
            <<"  \"nverts\": "<<nv+nc<<",\n  \"originalVertices\": "<<nv<<",\n"
            <<"  \"nf\": "<<nf<<",\n  \"boundaryFaces\": "<<nb<<",\n"
            <<"  \"originalArea\": "<<originalArea<<",\n  \"triangleArea\": "<<triangleArea<<",\n"
            <<"  \"areaDifference\": "<<triangleArea-originalArea<<",\n"
            <<"  \"maximumCellAreaDifference\": "<<maximumCellAreaError<<",\n"
            <<"  \"minimumTriangleArea\": "<<minimumTriangle<<",\n"
            <<"  \"nonpositiveFanTriangles\": 0,\n"
            <<"  \"cellAndFaceIds\": \"zero-based\",\n"
            <<"  \"meshVertexIdsAndBoundaryLabels\": \"one-based\",\n"
            <<"  \"triangleRegion\": \"original cell id + 1\",\n"
            <<"  \"boundaryLabel\": \"original face id + 1\",\n"
            <<"  \"boundaryKindTags\": {\"internal\": 0, \"wall\": 1, \"moving-wall\": 2, \"smooth-moving-wall\": 3},\n"
            <<"  \"meditEdgesContain\": \"physical boundary faces only; all original faces remain shared triangle edges and are also listed in .edges.dat\",\n"
            <<"  \"geometryValidation\": \"native makeFvMesh2D and concentric regular polygon annulus preset\",\n"
            <<"  \"cellTriangleMappingSuffix\": \".cell-triangles.csv\",\n"
            <<"  \"interiorTriangleSamplesSuffix\": \".triangles.dat\",\n"
            <<"  \"triangleSampleMeaning\": \"triangle interior centroid; exact area-mean sample for CR P1nc velocity and P0 pressure, avoiding ambiguous values at the original FV centre vertex\",\n"
            <<"  \"originalEdgesSuffix\": \".edges.dat\"\n}\n";
        for(auto* stream:{&medit,&boundary,&cells,&faces,&edges,&relation,&samples,&mapping}) {
            stream->flush();require(bool(*stream),"output write failed");
        }
        std::cout<<std::setprecision(17)<<"cells="<<nc<<" triangles="<<triangles.size()<<" vertices="<<nv+nc
                 <<" faces="<<nf<<" boundary="<<nb<<" area="<<originalArea<<" triangleArea="<<triangleArea
                 <<" minTriangleArea="<<minimumTriangle<<" nonpositiveFanTriangles=0\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
