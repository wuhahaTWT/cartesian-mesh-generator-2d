#include "cartmesh2d/io/OpenFoam2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <optional>
#include <utility>

namespace cartmesh2d {
namespace {

struct FoamFace2D {
    std::vector<std::size_t> vertices;
    std::size_t owner = 0;
    std::optional<std::size_t> neighbour;
};

struct PatchFaces2D {
    std::string name;
    std::string type;
    std::vector<FoamFace2D> faces;
};

void setError(std::string* error, std::string message) {
    if (error!=nullptr) *error=std::move(message);
}

void writeHeader(std::ofstream& out, const char* className, const char* object) {
    out<<"FoamFile\n{\n"
       <<"    version 2.0;\n"
       <<"    format ascii;\n"
       <<"    arch \"LSB;label=32;scalar=64\";\n"
       <<"    class "<<className<<";\n"
       <<"    location \"constant/polyMesh\";\n"
       <<"    object "<<object<<";\n}\n\n";
}

[[nodiscard]] bool edgeDirectionForCell(const TopologyCell2D& cell,
                                        std::size_t v0,std::size_t v1,
                                        std::size_t& from,std::size_t& to) noexcept {
    for (std::size_t i=0;i<cell.vertices.size();++i) {
        const std::size_t a=cell.vertices[i];
        const std::size_t b=cell.vertices[(i+1)%cell.vertices.size()];
        if ((a==v0 && b==v1) || (a==v1 && b==v0)) {
            from=a;
            to=b;
            return true;
        }
    }
    return false;
}

[[nodiscard]] std::string domainPatchName(const Edge2D& edge,
                                          const TopologyMesh2D& topology,
                                          const Domain2D& domain,
                                          const TolerancePolicy& tol) {
    const auto& a=topology.vertices[edge.v0].point;
    const auto& b=topology.vertices[edge.v1].point;
    const auto near=[&](double x,double y) {
        return std::abs(x-y)<=tol.scale(std::max(std::abs(x),std::abs(y)));
    };
    if (near(a.x,domain.bounds.min.x) && near(b.x,domain.bounds.min.x)) return "left";
    if (near(a.x,domain.bounds.max.x) && near(b.x,domain.bounds.max.x)) return "right";
    if (near(a.y,domain.bounds.min.y) && near(b.y,domain.bounds.min.y)) return "bottom";
    if (near(a.y,domain.bounds.max.y) && near(b.y,domain.bounds.max.y)) return "top";
    return {};
}

[[nodiscard]] PatchFaces2D* findPatch(std::vector<PatchFaces2D>& patches,
                                      const std::string& name) {
    for (auto& patch:patches) if (patch.name==name) return &patch;
    return nullptr;
}

[[nodiscard]] std::optional<std::size_t> embeddedLoopId(
    const Edge2D& edge,const TopologyMesh2D& topology,
    const BoundaryRegion2D& boundary,const TolerancePolicy& tol) {
    const Segment2D face{topology.vertices[edge.v0].point,
                         topology.vertices[edge.v1].point};
    for (std::size_t loopId=0;loopId<boundary.loops().size();++loopId) {
        const auto& vertices=boundary.loops()[loopId].vertices();
        for (std::size_t i=0;i<vertices.size();++i) {
            const Segment2D segment{vertices[i],vertices[(i+1)%vertices.size()]};
            if (pointOnSegment(face.a,segment,tol) && pointOnSegment(face.b,segment,tol)) {
                return loopId;
            }
        }
    }
    return std::nullopt;
}

} // namespace

OpenFoamWriteReport2D writeExtrudedOpenFoam2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const std::filesystem::path& caseDirectory,double thickness,
    std::string* error,const TolerancePolicy& tol) {
    OpenFoamWriteReport2D report;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid() ||
        !(thickness>tol.absolute)) {
        setError(error,"OpenFOAM extrusion requires valid topology, domain, boundary and positive thickness");
        return report;
    }

    const std::size_t layerOffset=topology.vertices.size();
    std::vector<FoamFace2D> internalFaces;
    std::vector<PatchFaces2D> patches;
    for (std::size_t loopId=0;loopId<boundary.loops().size();++loopId) {
        patches.push_back({"wall_"+std::to_string(loopId),"wall",{}});
    }
    const std::vector<PatchFaces2D> domainPatches{
        {"left","patch",{}},
        {"right","patch",{}},
        {"bottom","patch",{}},
        {"top","patch",{}},
        {"frontAndBack","empty",{}}
    };
    patches.insert(patches.end(),domainPatches.begin(),domainPatches.end());

    for (const auto& edge:topology.edges) {
        if (edge.owner>=topology.cells.size()) {
            setError(error,"OpenFOAM edge owner is out of range");
            return {};
        }
        std::size_t from=0,to=0;
        if (!edgeDirectionForCell(topology.cells[edge.owner],edge.v0,edge.v1,from,to)) {
            setError(error,"OpenFOAM edge is absent from owner cell loop");
            return {};
        }
        FoamFace2D face{{from,to,to+layerOffset,from+layerOffset},edge.owner,edge.neighbour};
        if (edge.neighbour) {
            if (*edge.neighbour>=topology.cells.size()) {
                setError(error,"OpenFOAM edge neighbour is out of range");
                return {};
            }
            internalFaces.push_back(std::move(face));
            continue;
        }
        std::string patchName;
        if (edge.patch==BoundaryPatch2D::EmbeddedBoundary) {
            const auto loopId=embeddedLoopId(edge,topology,boundary,tol);
            if (!loopId) {
                setError(error,"OpenFOAM embedded boundary edge does not belong to exactly represented input loop");
                return {};
            }
            patchName="wall_"+std::to_string(*loopId);
        } else if (edge.patch==BoundaryPatch2D::DomainBoundary) {
            patchName=domainPatchName(edge,topology,domain,tol);
        }
        PatchFaces2D* patch=findPatch(patches,patchName);
        if (patch==nullptr) {
            setError(error,"OpenFOAM boundary edge has no deterministic patch");
            return {};
        }
        face.neighbour.reset();
        patch->faces.push_back(std::move(face));
    }

    auto* frontAndBack=findPatch(patches,"frontAndBack");
    for (const auto& cell:topology.cells) {
        if (cell.vertices.size()<3) {
            setError(error,"OpenFOAM cell has fewer than three base vertices");
            return {};
        }
        std::vector<std::size_t> bottomFace(cell.vertices.rbegin(),cell.vertices.rend());
        std::vector<std::size_t> topFace;
        topFace.reserve(cell.vertices.size());
        for (const auto vertex:cell.vertices) topFace.push_back(vertex+layerOffset);
        frontAndBack->faces.push_back({std::move(bottomFace),cell.id,std::nullopt});
        frontAndBack->faces.push_back({std::move(topFace),cell.id,std::nullopt});
    }

    const auto polyMesh=caseDirectory/"constant"/"polyMesh";
    std::error_code filesystemError;
    std::filesystem::create_directories(polyMesh,filesystemError);
    if (filesystemError) {
        setError(error,"failed to create OpenFOAM constant/polyMesh directory");
        return {};
    }

    std::vector<FoamFace2D> faces=internalFaces;
    std::vector<OpenFoamPatchSummary2D> summaries;
    for (const auto& patch:patches) {
        if (patch.faces.empty()) continue;
        summaries.push_back({patch.name,patch.type,patch.faces.size(),faces.size()});
        faces.insert(faces.end(),patch.faces.begin(),patch.faces.end());
    }

    {
        std::ofstream out(polyMesh/"points");
        if (!out) { setError(error,"failed to open OpenFOAM points"); return {}; }
        writeHeader(out,"vectorField","points");
        out<<std::setprecision(17)<<2*topology.vertices.size()<<"\n(\n";
        for (const auto& vertex:topology.vertices) out<<'('<<vertex.point.x<<' '<<vertex.point.y<<" 0)\n";
        for (const auto& vertex:topology.vertices) out<<'('<<vertex.point.x<<' '<<vertex.point.y<<' '<<thickness<<")\n";
        out<<")\n";
        if (!out.good()) { setError(error,"failed while writing OpenFOAM points"); return {}; }
    }
    {
        std::ofstream out(polyMesh/"faces");
        if (!out) { setError(error,"failed to open OpenFOAM faces"); return {}; }
        writeHeader(out,"faceList","faces");
        out<<faces.size()<<"\n(\n";
        for (const auto& face:faces) {
            out<<face.vertices.size()<<'(';
            for (std::size_t i=0;i<face.vertices.size();++i) out<<(i==0?"":" ")<<face.vertices[i];
            out<<")\n";
        }
        out<<")\n";
        if (!out.good()) { setError(error,"failed while writing OpenFOAM faces"); return {}; }
    }
    {
        std::ofstream out(polyMesh/"owner");
        if (!out) { setError(error,"failed to open OpenFOAM owner"); return {}; }
        writeHeader(out,"labelList","owner");
        out<<faces.size()<<"\n(\n";
        for (const auto& face:faces) out<<face.owner<<'\n';
        out<<")\n";
    }
    {
        std::ofstream out(polyMesh/"neighbour");
        if (!out) { setError(error,"failed to open OpenFOAM neighbour"); return {}; }
        writeHeader(out,"labelList","neighbour");
        out<<internalFaces.size()<<"\n(\n";
        for (const auto& face:internalFaces) out<<*face.neighbour<<'\n';
        out<<")\n";
    }
    {
        std::ofstream out(polyMesh/"boundary");
        if (!out) { setError(error,"failed to open OpenFOAM boundary"); return {}; }
        writeHeader(out,"polyBoundaryMesh","boundary");
        out<<summaries.size()<<"\n(\n";
        for (const auto& patch:summaries) {
            out<<patch.name<<"\n{\n"
               <<"    type "<<patch.type<<";\n"
               <<"    nFaces "<<patch.faceCount<<";\n"
               <<"    startFace "<<patch.startFace<<";\n}\n";
        }
        out<<")\n";
    }

    report.pointCount=2*topology.vertices.size();
    report.faceCount=faces.size();
    report.internalFaceCount=internalFaces.size();
    report.cellCount=topology.cells.size();
    report.thickness=thickness;
    report.patches=std::move(summaries);
    return report;
}

} // namespace cartmesh2d
