#include "cartmesh2d/io/BackgroundGridIO2D.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace cartmesh2d {
namespace {
unsigned classification(CellClass value) {
    switch(value) {
    case CellClass::Outside: return 0;
    case CellClass::Inside: return 1;
    case CellClass::Intersected: return 2;
    }
    return 3;
}
}
bool writeBackgroundGrid2D(const Quadtree2D& tree, const BoundaryRegion2D& boundary,
                            const std::filesystem::path& prefix, bool uniform, std::string* error) {
    const auto fail=[&](const char* message){if(error)*error=message;return false;};
    if(error)error->clear();
    if(tree.leaves().empty() || !tree.deterministicOrderingValid())
        return fail("Invalid background leaf ordering");
    std::array<std::size_t,3> counts{};
    for(const auto& leaf:tree.leaves()) {
        const auto kind=classification(leaf.classification);
        if(kind>2 || !std::isfinite(leaf.area()) || !(leaf.area()>0))
            return fail("Invalid background cell geometry or classification");
        ++counts[kind];
    }
    std::ofstream json(prefix.string()+".background.json");
    std::ofstream vtk(prefix.string()+".background.vtk");
    if(!json || !vtk)return fail("Cannot open background grid output");
    const auto& d=tree.domain().bounds;
    json<<std::setprecision(17)
        <<"{\"format\":\"cartmesh2d-background-v1\",\"mode\":\""<<(uniform?"uniform":"adaptive")
        <<"\",\"solver_ready\":false,\"retains_solid_interior\":true,"
        <<"\"classification_names\":[\"outside\",\"inside\",\"intersected\"],"
        <<"\"domain\":["<<d.min.x<<','<<d.min.y<<','<<d.max.x<<','<<d.max.y<<"],"
        <<"\"counts\":["<<counts[0]<<','<<counts[1]<<','<<counts[2]<<"],\"boundary_loops\":[";
    for(std::size_t l=0;l<boundary.loops().size();++l) {
        if(l)json<<',';
        json<<'[';
        const auto& vertices=boundary.loops()[l].vertices();
        for(std::size_t i=0;i<vertices.size();++i) {
            if(i)json<<',';
            json<<'['<<vertices[i].x<<','<<vertices[i].y<<']';
        }
        json<<']';
    }
    json<<"],\"cells\":[\n";
    const auto n=tree.leaves().size();
    vtk<<std::setprecision(17)<<"# vtk DataFile Version 3.0\nCartMesh2D full background grid; not fluid solver topology\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "<<4*n<<" double\n";
    for(const auto& leaf:tree.leaves()) {
        if(leaf.id)json<<",\n";
        const auto& b=leaf.bounds;
        // Keys are strings because a level-28 key need not fit a JS integer.
        json<<"{\"id\":"<<leaf.id<<",\"key\":\""<<leaf.key<<"\",\"level\":"<<leaf.level
            <<",\"ix\":"<<leaf.ix<<",\"iy\":"<<leaf.iy<<",\"classification\":"<<classification(leaf.classification)
            <<",\"bounds\":["<<b.min.x<<','<<b.min.y<<','<<b.max.x<<','<<b.max.y<<"]}";
        vtk<<b.min.x<<' '<<b.min.y<<" 0\n"<<b.max.x<<' '<<b.min.y<<" 0\n"
           <<b.max.x<<' '<<b.max.y<<" 0\n"<<b.min.x<<' '<<b.max.y<<" 0\n";
    }
    json<<"\n]}\n";
    vtk<<"CELLS "<<n<<' '<<5*n<<'\n';
    for(std::size_t i=0;i<n;++i)vtk<<"4 "<<4*i<<' '<<4*i+1<<' '<<4*i+2<<' '<<4*i+3<<'\n';
    vtk<<"CELL_TYPES "<<n<<'\n';
    for(std::size_t i=0;i<n;++i)vtk<<"9\n";
    vtk<<"CELL_DATA "<<n<<"\nSCALARS classification int 1\nLOOKUP_TABLE default\n";
    for(const auto& leaf:tree.leaves())vtk<<classification(leaf.classification)<<'\n';
    vtk<<"SCALARS level int 1\nLOOKUP_TABLE default\n";
    for(const auto& leaf:tree.leaves())vtk<<leaf.level<<'\n';
    json.close();vtk.close();
    if(!json || !vtk)return fail("Failed while writing background grid");
    return true;
}
}
