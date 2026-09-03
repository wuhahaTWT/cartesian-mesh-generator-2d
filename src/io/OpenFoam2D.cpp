#include "cartmesh2d/io/OpenFoam2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
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
    BoundaryConditionRole2D role = BoundaryConditionRole2D::Wall;
    std::string sourceLayer;
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

void writeCaseHeader(std::ofstream& out,const char* className,
                     const char* location,const char* object) {
    out<<"FoamFile\n{\n"
       <<"    version 2.0;\n"
       <<"    format ascii;\n"
       <<"    class "<<className<<";\n"
       <<"    location \""<<location<<"\";\n"
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

[[nodiscard]] double pointSegmentDistance(const Point2D& point,
                                          const Segment2D& segment) noexcept {
    const Vector2D direction=segment.b-segment.a;
    const double lengthSquared=squaredNorm(direction);
    if (!(lengthSquared>0.0)) return std::sqrt(squaredNorm(point-segment.a));
    const double parameter=std::clamp(
        dot(point-segment.a,direction)/lengthSquared,0.0,1.0);
    return std::sqrt(squaredNorm(point-(segment.a+direction*parameter)));
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
    // A shared-construction grid-corner weld is deliberately allowed to move
    // an embedded fragment off the piecewise-linear input by a bounded f*h.
    // The topology has already classified the edge as physical wall; here we
    // only recover its deterministic loop/patch identity within that exact
    // construction budget. Ambiguous loop ownership still fails closed.
    if (!topology.constructionRegistry ||
        topology.canonicalVertexIds.size()!=topology.vertices.size()) {
        return std::nullopt;
    }
    const auto& registry=*topology.constructionRegistry;
    const auto handle0=topology.canonicalVertexIds[edge.v0];
    const auto handle1=topology.canonicalVertexIds[edge.v1];
    if (handle0>=registry.vertices().size() || handle1>=registry.vertices().size()) {
        return std::nullopt;
    }
    const double arithmeticBudget=
        tol.scale(std::max({1.0,std::abs(face.a.x),std::abs(face.a.y),
                           std::abs(face.b.x),std::abs(face.b.y)}));
    double distanceBudgetA=arithmeticBudget;
    double distanceBudgetB=arithmeticBudget;
    // Use the committed event displacement, not mutable local_h metadata from
    // a later convex-partition rebuild. This is a tighter provenance check and
    // remains valid after the same vertex participates in smaller solver cells.
    for (const auto& event:registry.events()) {
        if (event.canonicalVertex==handle0)
            distanceBudgetA=std::max(distanceBudgetA,event.displacement+arithmeticBudget);
        if (event.canonicalVertex==handle1)
            distanceBudgetB=std::max(distanceBudgetB,event.displacement+arithmeticBudget);
    }
    std::optional<std::size_t> matchedLoop;
    for (std::size_t loopId=0;loopId<boundary.loops().size();++loopId) {
        const auto& vertices=boundary.loops()[loopId].vertices();
        bool matches=false;
        for (std::size_t i=0;i<vertices.size();++i) {
            const Segment2D segment{vertices[i],vertices[(i+1)%vertices.size()]};
            if (pointSegmentDistance(face.a,segment)<=distanceBudgetA &&
                pointSegmentDistance(face.b,segment)<=distanceBudgetB) {
                matches=true;
                break;
            }
        }
        if (!matches) continue;
        if (matchedLoop && *matchedLoop!=loopId) return std::nullopt;
        matchedLoop=loopId;
    }
    return matchedLoop;
}

} // namespace

OpenFoamWriteReport2D writeExtrudedOpenFoam2D(
    const TopologyMesh2D& topology,const Domain2D& domain,
    const BoundaryRegion2D& boundary,
    const std::filesystem::path& caseDirectory,double thickness,
    std::string* error,const TolerancePolicy& tol,
    const std::vector<EmbeddedBoundaryPatch2D>& embeddedPatches) {
    OpenFoamWriteReport2D report;
    if (!topology.valid() || !domain.valid(tol) || !boundary.diagnose(tol).valid() ||
        !(thickness>tol.absolute)) {
        setError(error,"OpenFOAM extrusion requires valid topology, domain, boundary and positive thickness");
        return report;
    }

    const std::size_t layerOffset=topology.vertices.size();
    std::vector<FoamFace2D> internalFaces;
    std::vector<PatchFaces2D> patches;
    if (!embeddedPatches.empty() && embeddedPatches.size()!=boundary.loops().size()) {
        setError(error,"embedded boundary patch metadata count does not match loops");
        return report;
    }
    for (std::size_t loopId=0;loopId<boundary.loops().size();++loopId) {
        EmbeddedBoundaryPatch2D spec;
        if (embeddedPatches.empty()) {
            spec={"wall_"+std::to_string(loopId),"wall",
                  BoundaryConditionRole2D::Wall,""};
        } else spec=embeddedPatches[loopId];
        const std::string expectedType=spec.role==BoundaryConditionRole2D::Wall
            ?"wall":(spec.role==BoundaryConditionRole2D::Symmetry?"symmetryPlane":"patch");
        if (!spec.valid() || spec.type!=expectedType ||
            spec.name=="left" || spec.name=="right" || spec.name=="bottom" ||
            spec.name=="top" || spec.name=="frontAndBack") {
            setError(error,"invalid or reserved embedded boundary patch metadata");
            return report;
        }
        if (const auto* existing=findPatch(patches,spec.name);existing!=nullptr) {
            if (existing->type!=spec.type || existing->role!=spec.role ||
                existing->sourceLayer!=spec.sourceLayer) {
                setError(error,"duplicate embedded patch name has conflicting metadata");
                return report;
            }
            continue;
        }
        patches.push_back({spec.name,spec.type,spec.role,spec.sourceLayer,{}});
    }
    const std::vector<PatchFaces2D> domainPatches{
        {"left","patch",BoundaryConditionRole2D::Inlet,"",{}},
        {"right","patch",BoundaryConditionRole2D::Outlet,"",{}},
        {"bottom","patch",BoundaryConditionRole2D::Slip,"",{}},
        {"top","patch",BoundaryConditionRole2D::Slip,"",{}},
        {"frontAndBack","empty",BoundaryConditionRole2D::Slip,"",{}}
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
                const auto& a=topology.vertices[edge.v0].point;
                const auto& b=topology.vertices[edge.v1].point;
                double closestA=std::numeric_limits<double>::infinity();
                double closestB=std::numeric_limits<double>::infinity();
                for (const auto& loop:boundary.loops()) {
                    const auto& vertices=loop.vertices();
                    for (std::size_t i=0;i<vertices.size();++i) {
                        const Segment2D segment{vertices[i],vertices[(i+1)%vertices.size()]};
                        closestA=std::min(closestA,pointSegmentDistance(a,segment));
                        closestB=std::min(closestB,pointSegmentDistance(b,segment));
                    }
                }
                std::ostringstream detail;
                detail<<std::setprecision(17)
                      <<"OpenFOAM embedded boundary edge has no unique input-loop identity within the construction budget"
                      <<" edge="<<edge.id<<" a=("<<a.x<<','<<a.y<<')'
                      <<" b=("<<b.x<<','<<b.y<<')'
                      <<" nearest_a="<<closestA<<" nearest_b="<<closestB;
                setError(error,detail.str());
                return {};
            }
            patchName=embeddedPatches.empty()
                ?"wall_"+std::to_string(*loopId):embeddedPatches[*loopId].name;
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

    // OpenFOAM's upper-triangular addressing requires internal faces to be
    // grouped by owner and then neighbour, in addition to owner < neighbour.
    // Topology edge IDs are vertex-key ordered and therefore cannot be used as
    // the polyMesh face order.
    std::sort(internalFaces.begin(),internalFaces.end(),
              [](const FoamFace2D& lhs,const FoamFace2D& rhs) {
                  if (lhs.owner!=rhs.owner) return lhs.owner<rhs.owner;
                  return *lhs.neighbour<*rhs.neighbour;
              });

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
        summaries.push_back({patch.name,patch.type,patch.role,patch.faces.size(),faces.size()});
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

    // Emit a deliberately low-Reynolds-number laminar exterior-flow fixture.
    // Products with left/right patches are directly consumable by simpleFoam;
    // closed or disconnected regions retain the files for inspection but need
    // problem-specific pressure references before flow validation.
    const auto systemDirectory=caseDirectory/"system";
    const auto zeroDirectory=caseDirectory/"0";
    std::filesystem::create_directories(systemDirectory,filesystemError);
    std::filesystem::create_directories(zeroDirectory,filesystemError);
    if (filesystemError) {
        setError(error,"failed to create OpenFOAM system/0 directories");
        return {};
    }
    {
        std::ofstream out(systemDirectory/"controlDict");
        if (!out) { setError(error,"failed to open OpenFOAM controlDict"); return {}; }
        writeCaseHeader(out,"dictionary","system","controlDict");
        out<<"application simpleFoam;\nstartFrom startTime;\nstartTime 0;\n"
           <<"stopAt endTime;\nendTime 300;\ndeltaT 1;\n"
           <<"writeControl timeStep;\nwriteInterval 300;\npurgeWrite 0;\n"
           <<"writeFormat ascii;\nwritePrecision 10;\nwriteCompression off;\n"
           <<"timeFormat general;\ntimePrecision 6;\nrunTimeModifiable true;\n";
    }
    {
        std::ofstream out(systemDirectory/"fvSchemes");
        if (!out) { setError(error,"failed to open OpenFOAM fvSchemes"); return {}; }
        writeCaseHeader(out,"dictionary","system","fvSchemes");
        out<<"ddtSchemes { default steadyState; }\n"
           <<"gradSchemes { default Gauss linear; }\n"
           <<"divSchemes\n{\n    default none;\n"
           <<"    div(phi,U) bounded Gauss upwind;\n"
           <<"    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n"
           <<"laplacianSchemes { default Gauss linear corrected; }\n"
           <<"interpolationSchemes { default linear; }\n"
           <<"snGradSchemes { default corrected; }\n"
           <<"wallDist { method meshWave; }\n";
    }
    {
        std::ofstream out(systemDirectory/"fvSolution");
        if (!out) { setError(error,"failed to open OpenFOAM fvSolution"); return {}; }
        writeCaseHeader(out,"dictionary","system","fvSolution");
        out<<"solvers\n{\n"
           <<"    p { solver GAMG; tolerance 1e-10; relTol 0.05; smoother GaussSeidel; }\n"
           <<"    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.05; }\n"
           <<"}\nSIMPLE\n{\n    nNonOrthogonalCorrectors 1;\n"
           <<"    residualControl { p 1e-7; U 1e-7; }\n}\n"
           <<"relaxationFactors\n{\n    fields { p 0.2; }\n    equations { U 0.4; }\n}\n";
    }
    {
        std::ofstream out(caseDirectory/"constant"/"transportProperties");
        if (!out) { setError(error,"failed to open OpenFOAM transportProperties"); return {}; }
        writeCaseHeader(out,"dictionary","constant","transportProperties");
        out<<"transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] 0.01;\n";
    }
    {
        std::ofstream out(caseDirectory/"constant"/"turbulenceProperties");
        if (!out) { setError(error,"failed to open OpenFOAM turbulenceProperties"); return {}; }
        writeCaseHeader(out,"dictionary","constant","turbulenceProperties");
        out<<"simulationType laminar;\n";
    }
    {
        std::ofstream out(zeroDirectory/"U");
        if (!out) { setError(error,"failed to open OpenFOAM U field"); return {}; }
        writeCaseHeader(out,"volVectorField","0","U");
        out<<"dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0.1 0 0);\n"
           <<"boundaryField\n{\n";
        for (const auto& patch:summaries) {
            out<<"    "<<patch.name<<" { ";
            if (patch.type=="empty") out<<"type empty;";
            else if (patch.role==BoundaryConditionRole2D::Inlet)
                out<<"type fixedValue; value uniform (0.1 0 0);";
            else if (patch.role==BoundaryConditionRole2D::Outlet) out<<"type zeroGradient;";
            else if (patch.role==BoundaryConditionRole2D::Slip) out<<"type slip;";
            else if (patch.role==BoundaryConditionRole2D::Symmetry) out<<"type symmetryPlane;";
            else out<<"type noSlip;";
            out<<" }\n";
        }
        out<<"}\n";
    }
    {
        std::ofstream out(zeroDirectory/"p");
        if (!out) { setError(error,"failed to open OpenFOAM p field"); return {}; }
        writeCaseHeader(out,"volScalarField","0","p");
        out<<"dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\n"
           <<"boundaryField\n{\n";
        for (const auto& patch:summaries) {
            out<<"    "<<patch.name<<" { ";
            if (patch.type=="empty") out<<"type empty;";
            else if (patch.role==BoundaryConditionRole2D::Outlet)
                out<<"type fixedValue; value uniform 0;";
            else if (patch.role==BoundaryConditionRole2D::Symmetry) out<<"type symmetryPlane;";
            else out<<"type zeroGradient;";
            out<<" }\n";
        }
        out<<"}\n";
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
