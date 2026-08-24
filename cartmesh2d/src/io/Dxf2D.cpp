#include "cartmesh2d/io/Dxf2D.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>

namespace cartmesh2d {
namespace {

struct DxfGroup {
    int code = 0;
    std::string value;
    std::size_t line = 0;
};

struct EntityRecord {
    std::string type;
    std::string layer = "0";
    std::size_t line = 0;
    std::vector<DxfGroup> groups;
};

struct OpenChain {
    std::vector<Point2D> points;
    std::size_t line = 0;
    std::string entity;
    std::string layer;
};

[[nodiscard]] std::string trim(std::string value) {
    const auto first=value.find_first_not_of(" \t\r\n");
    if (first==std::string::npos) return {};
    const auto last=value.find_last_not_of(" \t\r\n");
    return value.substr(first,last-first+1);
}

[[nodiscard]] bool parseInteger(const std::string& text,long long& value) {
    try {
        std::size_t consumed=0;
        value=std::stoll(trim(text),&consumed,10);
        return consumed==trim(text).size();
    } catch (const std::exception&) {
        return false;
    }
}

[[nodiscard]] bool parseDouble(const std::string& text,double& value) {
    try {
        const std::string cleaned=trim(text);
        std::size_t consumed=0;
        value=std::stod(cleaned,&consumed);
        return consumed==cleaned.size() && std::isfinite(value);
    } catch (const std::exception&) {
        return false;
    }
}

[[nodiscard]] std::optional<double> lastDouble(
    const EntityRecord& entity,int code,DxfImportResult2D& result,
    bool required=false) {
    std::optional<double> value;
    for (const auto& group:entity.groups) {
        if (group.code!=code) continue;
        double parsed=0.0;
        if (!parseDouble(group.value,parsed)) {
            result.issues.push_back({DxfIssueCode2D::InvalidNumericValue,group.line,
                                     entity.type,entity.layer,
                                     "invalid finite floating-point group "+
                                     std::to_string(code)});
            return std::nullopt;
        }
        value=parsed;
    }
    if (required && !value) {
        result.issues.push_back({DxfIssueCode2D::MissingRequiredGroup,entity.line,
                                 entity.type,entity.layer,
                                 "missing required group "+std::to_string(code)});
    }
    return value;
}

[[nodiscard]] std::optional<long long> lastInteger(
    const EntityRecord& entity,int code,DxfImportResult2D& result,
    bool required=false) {
    std::optional<long long> value;
    for (const auto& group:entity.groups) {
        if (group.code!=code) continue;
        long long parsed=0;
        if (!parseInteger(group.value,parsed)) {
            result.issues.push_back({DxfIssueCode2D::InvalidNumericValue,group.line,
                                     entity.type,entity.layer,
                                     "invalid integer group "+std::to_string(code)});
            return std::nullopt;
        }
        value=parsed;
    }
    if (required && !value) {
        result.issues.push_back({DxfIssueCode2D::MissingRequiredGroup,entity.line,
                                 entity.type,entity.layer,
                                 "missing required group "+std::to_string(code)});
    }
    return value;
}

[[nodiscard]] bool nearZero(double value,double epsilon) noexcept {
    return std::abs(value)<=epsilon;
}

[[nodiscard]] bool validatePlanarAttributes(
    const EntityRecord& entity,double z,double thickness,
    double extrusionX,double extrusionY,double extrusionZ,
    const DxfImportOptions2D& options,DxfImportResult2D& result) {
    const double epsilon=options.endpointWeldTolerance;
    if (!nearZero(thickness,epsilon)) {
        result.issues.push_back({DxfIssueCode2D::NonZeroWidthOrThickness,entity.line,
                                 entity.type,entity.layer,
                                 "non-zero DXF thickness is not a 2-D boundary"});
        return false;
    }
    if (!nearZero(extrusionX,epsilon) || !nearZero(extrusionY,epsilon) ||
        std::abs(extrusionZ-1.0)>epsilon) {
        result.issues.push_back({DxfIssueCode2D::UnsupportedExtrusion,entity.line,
                                 entity.type,entity.layer,
                                 "only the default +Z object coordinate system is supported"});
        return false;
    }
    if (!result.report.sourcePlaneZDefined) {
        result.report.sourcePlaneZ=z;
        result.report.sourcePlaneZDefined=true;
    } else if (std::abs(result.report.sourcePlaneZ-z)>epsilon) {
        result.issues.push_back({DxfIssueCode2D::NonPlanarEntity,entity.line,
                                 entity.type,entity.layer,
                                 "entities do not share one XY plane"});
        return false;
    }
    return true;
}

[[nodiscard]] std::size_t segmentCountForArc(
    double radius,double sweep,double maximumChordError) {
    const double ratio=std::clamp(maximumChordError/radius,0.0,2.0);
    double maximumAngle=2.0*std::acos(std::clamp(1.0-ratio,-1.0,1.0));
    if (!(maximumAngle>0.0)) maximumAngle=sweep;
    const auto count=static_cast<std::size_t>(std::ceil(sweep/maximumAngle));
    return std::max<std::size_t>(1,count);
}

[[nodiscard]] bool appendCircularArc(
    std::vector<Point2D>& points,const Point2D& centre,double radius,
    double startAngle,double sweep,const DxfImportOptions2D& options,
    DxfImportResult2D& result,const EntityRecord& entity,
    bool includeStart) {
    if (!(radius>0.0) || !std::isfinite(radius) || !(sweep>0.0) ||
        !std::isfinite(sweep)) {
        result.issues.push_back({DxfIssueCode2D::DegenerateEntity,entity.line,
                                 entity.type,entity.layer,
                                 "arc or circle has invalid radius/sweep"});
        return false;
    }
    const std::size_t segments=segmentCountForArc(
        radius,sweep,options.maximumChordError);
    constexpr std::size_t maximumSegments=1000000;
    if (segments>maximumSegments) {
        result.issues.push_back({DxfIssueCode2D::DegenerateEntity,entity.line,
                                 entity.type,entity.layer,
                                 "curve sampling exceeds one million segments"});
        return false;
    }
    const std::size_t first=includeStart?0:1;
    for (std::size_t i=first;i<=segments;++i) {
        const double fraction=static_cast<double>(i)/static_cast<double>(segments);
        const double angle=startAngle+sweep*fraction;
        points.push_back({centre.x+radius*std::cos(angle),
                          centre.y+radius*std::sin(angle)});
    }
    result.report.sampledArcSegmentCount+=segments;
    return true;
}

[[nodiscard]] bool appendBulgeSegment(
    std::vector<Point2D>& points,const Point2D& start,const Point2D& end,
    double bulge,const DxfImportOptions2D& options,DxfImportResult2D& result,
    const EntityRecord& entity) {
    if (nearZero(bulge,options.endpointWeldTolerance)) {
        points.push_back(end);
        return true;
    }
    const Vector2D chord=end-start;
    const double chordLength=std::sqrt(squaredNorm(chord));
    if (!(chordLength>options.endpointWeldTolerance)) {
        result.issues.push_back({DxfIssueCode2D::DegenerateEntity,entity.line,
                                 entity.type,entity.layer,
                                 "bulged polyline segment has coincident endpoints"});
        return false;
    }
    const double sweep=4.0*std::atan(bulge);
    const double radius=chordLength*(1.0+bulge*bulge)/(4.0*std::abs(bulge));
    const Point2D midpoint{0.5*(start.x+end.x),0.5*(start.y+end.y)};
    const Vector2D leftNormal{-chord.y/chordLength,chord.x/chordLength};
    const double centreOffset=chordLength*(1.0-bulge*bulge)/(4.0*bulge);
    const Point2D centre=midpoint+leftNormal*centreOffset;
    const double startAngle=std::atan2(start.y-centre.y,start.x-centre.x);
    const std::size_t segments=segmentCountForArc(
        radius,std::abs(sweep),options.maximumChordError);
    if (segments>1000000) {
        result.issues.push_back({DxfIssueCode2D::DegenerateEntity,entity.line,
                                 entity.type,entity.layer,
                                 "bulge sampling exceeds one million segments"});
        return false;
    }
    for (std::size_t i=1;i<=segments;++i) {
        const double fraction=static_cast<double>(i)/static_cast<double>(segments);
        const double angle=startAngle+sweep*fraction;
        points.push_back({centre.x+radius*std::cos(angle),
                          centre.y+radius*std::sin(angle)});
    }
    // Preserve the exact DXF endpoint so separately authored entities can be
    // welded without depending on trigonometric round-off.
    points.back()=end;
    result.report.sampledArcSegmentCount+=segments;
    return true;
}

[[nodiscard]] bool pointLexLess(const Point2D& lhs,const Point2D& rhs) noexcept {
    return std::tie(lhs.x,lhs.y)<std::tie(rhs.x,rhs.y);
}

void canonicalRotate(std::vector<Point2D>& vertices) {
    if (vertices.empty()) return;
    const auto minimum=std::min_element(vertices.begin(),vertices.end(),pointLexLess);
    std::rotate(vertices.begin(),minimum,vertices.end());
}

[[nodiscard]] std::string jsonEscape(const std::string& value) {
    std::ostringstream out;
    for (const char raw:value) {
        const auto c=static_cast<unsigned char>(raw);
        switch (c) {
            case '\\': out<<"\\\\"; break;
            case '"': out<<"\\\""; break;
            case '\n': out<<"\\n"; break;
            case '\r': out<<"\\r"; break;
            case '\t': out<<"\\t"; break;
            default:
                if (c<0x20U) {
                    out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')
                       <<static_cast<unsigned>(c)<<std::dec;
                } else out<<static_cast<char>(c);
        }
    }
    return out.str();
}

[[nodiscard]] const char* issueCodeName(DxfIssueCode2D code) noexcept {
    switch (code) {
        case DxfIssueCode2D::CannotOpen: return "cannot_open";
        case DxfIssueCode2D::BinaryDxfUnsupported: return "binary_dxf_unsupported";
        case DxfIssueCode2D::MalformedGroupPair: return "malformed_group_pair";
        case DxfIssueCode2D::MissingEntitiesSection: return "missing_entities_section";
        case DxfIssueCode2D::UnsupportedEntity: return "unsupported_entity";
        case DxfIssueCode2D::MissingRequiredGroup: return "missing_required_group";
        case DxfIssueCode2D::InvalidNumericValue: return "invalid_numeric_value";
        case DxfIssueCode2D::NonPlanarEntity: return "non_planar_entity";
        case DxfIssueCode2D::UnsupportedExtrusion: return "unsupported_extrusion";
        case DxfIssueCode2D::NonZeroWidthOrThickness: return "non_zero_width_or_thickness";
        case DxfIssueCode2D::OpenOrBranchedBoundary: return "open_or_branched_boundary";
        case DxfIssueCode2D::DegenerateEntity: return "degenerate_entity";
        case DxfIssueCode2D::InvalidBoundaryRegion: return "invalid_boundary_region";
        case DxfIssueCode2D::OutputFailure: return "output_failure";
    }
    return "unknown";
}

class DisjointSet {
public:
    explicit DisjointSet(std::size_t size):parent_(size),rank_(size,0) {
        std::iota(parent_.begin(),parent_.end(),0);
    }
    [[nodiscard]] std::size_t find(std::size_t value) {
        if (parent_[value]!=value) parent_[value]=find(parent_[value]);
        return parent_[value];
    }
    void unite(std::size_t lhs,std::size_t rhs) {
        lhs=find(lhs); rhs=find(rhs);
        if (lhs==rhs) return;
        if (rank_[lhs]<rank_[rhs]) std::swap(lhs,rhs);
        parent_[rhs]=lhs;
        if (rank_[lhs]==rank_[rhs]) ++rank_[lhs];
    }
private:
    std::vector<std::size_t> parent_;
    std::vector<unsigned> rank_;
};

[[nodiscard]] bool assembleOpenChains(
    const std::vector<OpenChain>& chains,double weldTolerance,
    std::vector<std::vector<Point2D>>& loops,DxfImportResult2D& result) {
    if (chains.empty()) return true;
    struct Endpoint { Point2D point; std::size_t chain=0; bool end=false; };
    std::vector<Endpoint> endpoints;
    endpoints.reserve(chains.size()*2);
    for (std::size_t i=0;i<chains.size();++i) {
        if (chains[i].points.size()<2) {
            result.issues.push_back({DxfIssueCode2D::DegenerateEntity,chains[i].line,
                                     chains[i].entity,chains[i].layer,
                                     "entity emitted fewer than two points"});
            return false;
        }
        endpoints.push_back({chains[i].points.front(),i,false});
        endpoints.push_back({chains[i].points.back(),i,true});
    }
    std::vector<std::size_t> order(endpoints.size());
    std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](std::size_t lhs,std::size_t rhs) {
        return std::tie(endpoints[lhs].point.x,endpoints[lhs].point.y,lhs)<
               std::tie(endpoints[rhs].point.x,endpoints[rhs].point.y,rhs);
    });
    DisjointSet sets(endpoints.size());
    for (std::size_t oi=0;oi<order.size();++oi) {
        const auto lhs=order[oi];
        for (std::size_t oj=oi+1;oj<order.size();++oj) {
            const auto rhs=order[oj];
            if (endpoints[rhs].point.x-endpoints[lhs].point.x>weldTolerance) break;
            if (std::abs(endpoints[rhs].point.y-endpoints[lhs].point.y)<=weldTolerance &&
                std::hypot(endpoints[rhs].point.x-endpoints[lhs].point.x,
                           endpoints[rhs].point.y-endpoints[lhs].point.y)<=weldTolerance) {
                sets.unite(lhs,rhs);
            }
        }
    }
    std::map<std::size_t,std::vector<std::size_t>> clusterMembers;
    for (std::size_t endpoint=0;endpoint<endpoints.size();++endpoint)
        clusterMembers[sets.find(endpoint)].push_back(endpoint);
    for (const auto& [root,members]:clusterMembers) {
        (void)root;
        for (std::size_t i=0;i<members.size();++i) {
            for (std::size_t j=i+1;j<members.size();++j) {
                const auto& a=endpoints[members[i]].point;
                const auto& b=endpoints[members[j]].point;
                if (std::hypot(a.x-b.x,a.y-b.y)>weldTolerance) {
                    result.issues.push_back({DxfIssueCode2D::OpenOrBranchedBoundary,0,
                                             "","",
                                             "endpoint weld cluster exceeds the configured "
                                             "diameter; chaining would over-weld geometry"});
                    return false;
                }
            }
        }
    }
    std::map<std::size_t,std::size_t> rootToNode;
    std::vector<Point2D> nodePoints;
    std::vector<std::size_t> endpointNode(endpoints.size());
    for (const auto endpointIndex:order) {
        const auto root=sets.find(endpointIndex);
        auto [it,inserted]=rootToNode.emplace(root,nodePoints.size());
        if (inserted) nodePoints.push_back(endpoints[endpointIndex].point);
        else if (pointLexLess(endpoints[endpointIndex].point,nodePoints[it->second]))
            nodePoints[it->second]=endpoints[endpointIndex].point;
        endpointNode[endpointIndex]=it->second;
    }
    std::vector<std::array<std::size_t,2>> chainNodes(chains.size());
    std::vector<std::vector<std::size_t>> adjacency(nodePoints.size());
    for (std::size_t chain=0;chain<chains.size();++chain) {
        chainNodes[chain]={endpointNode[2*chain],endpointNode[2*chain+1]};
        adjacency[chainNodes[chain][0]].push_back(chain);
        adjacency[chainNodes[chain][1]].push_back(chain);
    }
    for (std::size_t node=0;node<adjacency.size();++node) {
        std::sort(adjacency[node].begin(),adjacency[node].end());
        if (adjacency[node].size()!=2) {
            std::ostringstream detail;
            detail<<"boundary endpoint ("<<std::setprecision(17)<<nodePoints[node].x<<','
                  <<nodePoints[node].y<<") has degree "<<adjacency[node].size()
                  <<"; expected 2";
            result.issues.push_back({DxfIssueCode2D::OpenOrBranchedBoundary,0,"","",
                                     detail.str()});
            return false;
        }
    }
    std::vector<bool> used(chains.size(),false);
    for (std::size_t seed=0;seed<chains.size();++seed) {
        if (used[seed]) continue;
        std::size_t startNode=chainNodes[seed][0];
        if (pointLexLess(nodePoints[chainNodes[seed][1]],nodePoints[startNode]))
            startNode=chainNodes[seed][1];
        std::size_t currentNode=startNode;
        std::vector<Point2D> loop{nodePoints[startNode]};
        while (true) {
            std::optional<std::size_t> nextChain;
            for (const auto candidate:adjacency[currentNode]) {
                if (!used[candidate]) { nextChain=candidate; break; }
            }
            if (!nextChain) {
                if (currentNode==startNode) break;
                result.issues.push_back({DxfIssueCode2D::OpenOrBranchedBoundary,0,"","",
                                         "chain traversal terminated before closing"});
                return false;
            }
            const auto chain=*nextChain;
            used[chain]=true;
            const bool forward=chainNodes[chain][0]==currentNode;
            const auto& points=chains[chain].points;
            if (forward) {
                for (std::size_t i=1;i<points.size();++i) loop.push_back(points[i]);
                currentNode=chainNodes[chain][1];
            } else {
                for (std::size_t i=points.size()-1;i-->0;) loop.push_back(points[i]);
                currentNode=chainNodes[chain][0];
            }
            loop.back()=nodePoints[currentNode];
            if (currentNode==startNode) break;
        }
        if (loop.size()>1 && nearlyEqual(loop.front(),loop.back(),
                TolerancePolicy{weldTolerance,0.0})) loop.pop_back();
        loops.push_back(std::move(loop));
    }
    return true;
}

} // namespace

DxfImportResult2D readAsciiDxfBoundary2D(
    const std::filesystem::path& path,const DxfImportOptions2D& options,
    const TolerancePolicy& tol) {
    DxfImportResult2D result;
    result.report.maximumChordError=options.maximumChordError;
    result.report.endpointWeldTolerance=options.endpointWeldTolerance;
    if (!(options.maximumChordError>0.0) || !std::isfinite(options.maximumChordError) ||
        !(options.endpointWeldTolerance>0.0) ||
        !std::isfinite(options.endpointWeldTolerance)) {
        result.issues.push_back({DxfIssueCode2D::InvalidNumericValue,0,"","",
                                 "DXF tolerances must be finite and positive"});
        return result;
    }

    std::ifstream in(path,std::ios::binary);
    if (!in) {
        result.issues.push_back({DxfIssueCode2D::CannotOpen,0,"","",
                                 "cannot open DXF file: "+path.string()});
        return result;
    }
    std::string signature(22,'\0');
    in.read(signature.data(),static_cast<std::streamsize>(signature.size()));
    if (signature.rfind("AutoCAD Binary DXF",0)==0) {
        result.issues.push_back({DxfIssueCode2D::BinaryDxfUnsupported,1,"","",
                                 "DXF-1 accepts ASCII DXF only"});
        return result;
    }
    in.clear();
    in.seekg(0);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in,line)) lines.push_back(line);
    if (lines.size()%2!=0) {
        result.issues.push_back({DxfIssueCode2D::MalformedGroupPair,lines.size(),"","",
                                 "ASCII DXF must contain code/value line pairs"});
        return result;
    }
    std::vector<DxfGroup> groups;
    groups.reserve(lines.size()/2);
    for (std::size_t i=0;i<lines.size();i+=2) {
        long long code=0;
        if (!parseInteger(lines[i],code) || code<std::numeric_limits<int>::min() ||
            code>std::numeric_limits<int>::max()) {
            result.issues.push_back({DxfIssueCode2D::MalformedGroupPair,i+1,"","",
                                     "invalid DXF group code"});
            return result;
        }
        groups.push_back({static_cast<int>(code),trim(lines[i+1]),i+1});
    }

    bool inEntities=false;
    bool inHeader=false;
    bool sawEntities=false;
    std::vector<EntityRecord> entities;
    for (std::size_t i=0;i<groups.size();) {
        if (groups[i].code==0 && groups[i].value=="SECTION" && i+1<groups.size() &&
            groups[i+1].code==2) {
            inEntities=groups[i+1].value=="ENTITIES";
            inHeader=groups[i+1].value=="HEADER";
            sawEntities=sawEntities || inEntities;
            i+=2;
            continue;
        }
        if (groups[i].code==0 && groups[i].value=="ENDSEC") {
            inEntities=false;
            inHeader=false;
            ++i;
            continue;
        }
        if (inHeader && groups[i].code==9 && groups[i].value=="$INSUNITS") {
            if (i+1>=groups.size() || groups[i+1].code!=70) {
                result.issues.push_back({DxfIssueCode2D::MalformedGroupPair,
                                         groups[i].line,"","",
                                         "$INSUNITS is not followed by group 70"});
                return result;
            }
            long long units=0;
            if (!parseInteger(groups[i+1].value,units)) {
                result.issues.push_back({DxfIssueCode2D::InvalidNumericValue,
                                         groups[i+1].line,"","",
                                         "invalid $INSUNITS value"});
                return result;
            }
            result.report.insertionUnitsCode=units;
            result.report.insertionUnitsCodeDefined=true;
            i+=2;
            continue;
        }
        if (!inEntities) { ++i; continue; }
        if (groups[i].code==999) { ++i; continue; }
        if (groups[i].code!=0) {
            result.issues.push_back({DxfIssueCode2D::MalformedGroupPair,groups[i].line,"","",
                                     "entity does not begin with group code 0"});
            return result;
        }
        EntityRecord entity;
        entity.type=groups[i].value;
        entity.line=groups[i].line;
        ++i;
        while (i<groups.size() && groups[i].code!=0) {
            entity.groups.push_back(groups[i]);
            if (groups[i].code==8) entity.layer=groups[i].value;
            ++i;
        }
        entities.push_back(std::move(entity));
    }
    if (!sawEntities) {
        result.issues.push_back({DxfIssueCode2D::MissingEntitiesSection,0,"","",
                                 "DXF has no ENTITIES section"});
        return result;
    }

    std::vector<OpenChain> openChains;
    std::vector<std::vector<Point2D>> closedLoops;
    std::set<std::string> layers;
    for (const auto& entity:entities) {
        layers.insert(entity.layer);
        const auto extrusionX=lastDouble(entity,210,result).value_or(0.0);
        const auto extrusionY=lastDouble(entity,220,result).value_or(0.0);
        const auto extrusionZ=lastDouble(entity,230,result).value_or(1.0);
        const auto thickness=lastDouble(entity,39,result).value_or(0.0);
        if (!result.issues.empty()) return result;

        if (entity.type=="LINE") {
            const auto x0=lastDouble(entity,10,result,true);
            const auto y0=lastDouble(entity,20,result,true);
            const auto x1=lastDouble(entity,11,result,true);
            const auto y1=lastDouble(entity,21,result,true);
            const double z0=lastDouble(entity,30,result).value_or(0.0);
            const double z1=lastDouble(entity,31,result).value_or(0.0);
            if (!result.issues.empty()) return result;
            if (std::abs(z0-z1)>options.endpointWeldTolerance) {
                result.issues.push_back({DxfIssueCode2D::NonPlanarEntity,entity.line,
                                         entity.type,entity.layer,
                                         "LINE endpoints do not lie in one XY plane"});
                return result;
            }
            if (!validatePlanarAttributes(entity,z0,thickness,extrusionX,extrusionY,
                                           extrusionZ,options,result)) return result;
            const Point2D a{*x0,*y0},b{*x1,*y1};
            if (std::hypot(a.x-b.x,a.y-b.y)<=options.endpointWeldTolerance) {
                result.issues.push_back({DxfIssueCode2D::DegenerateEntity,entity.line,
                                         entity.type,entity.layer,"zero-length LINE"});
                return result;
            }
            openChains.push_back({{a,b},entity.line,entity.type,entity.layer});
            ++result.report.lineEntityCount;
        } else if (entity.type=="ARC") {
            const auto cx=lastDouble(entity,10,result,true);
            const auto cy=lastDouble(entity,20,result,true);
            const double z=lastDouble(entity,30,result).value_or(0.0);
            const auto radius=lastDouble(entity,40,result,true);
            const auto start=lastDouble(entity,50,result,true);
            const auto end=lastDouble(entity,51,result,true);
            if (!result.issues.empty()) return result;
            if (!validatePlanarAttributes(entity,z,thickness,extrusionX,extrusionY,
                                           extrusionZ,options,result)) return result;
            double sweep=std::fmod(*end-*start,360.0);
            if (sweep<=0.0) sweep+=360.0;
            if (std::abs(sweep-360.0)<1.0e-12) {
                result.issues.push_back({DxfIssueCode2D::DegenerateEntity,entity.line,
                                         entity.type,entity.layer,
                                         "ARC start and end angles coincide; use CIRCLE"});
                return result;
            }
            std::vector<Point2D> points;
            if (!appendCircularArc(points,{*cx,*cy},*radius,
                    *start*std::numbers::pi/180.0,sweep*std::numbers::pi/180.0,
                    options,result,entity,true)) return result;
            openChains.push_back({std::move(points),entity.line,entity.type,entity.layer});
            ++result.report.arcEntityCount;
        } else if (entity.type=="CIRCLE") {
            const auto cx=lastDouble(entity,10,result,true);
            const auto cy=lastDouble(entity,20,result,true);
            const double z=lastDouble(entity,30,result).value_or(0.0);
            const auto radius=lastDouble(entity,40,result,true);
            if (!result.issues.empty()) return result;
            if (!validatePlanarAttributes(entity,z,thickness,extrusionX,extrusionY,
                                           extrusionZ,options,result)) return result;
            std::vector<Point2D> points;
            if (!appendCircularArc(points,{*cx,*cy},*radius,0.0,
                    2.0*std::numbers::pi,options,result,entity,true)) return result;
            if (!points.empty()) points.pop_back();
            closedLoops.push_back(std::move(points));
            ++result.report.circleEntityCount;
        } else if (entity.type=="LWPOLYLINE") {
            const auto expected=lastInteger(entity,90,result,true);
            const auto flags=lastInteger(entity,70,result).value_or(0);
            const double z=lastDouble(entity,38,result).value_or(0.0);
            const double constantWidth=lastDouble(entity,43,result).value_or(0.0);
            if (!result.issues.empty()) return result;
            if (!validatePlanarAttributes(entity,z,thickness,extrusionX,extrusionY,
                                           extrusionZ,options,result)) return result;
            if (!nearZero(constantWidth,options.endpointWeldTolerance)) {
                result.issues.push_back({DxfIssueCode2D::NonZeroWidthOrThickness,entity.line,
                                         entity.type,entity.layer,
                                         "wide LWPOLYLINE is not a boundary centreline"});
                return result;
            }
            struct Vertex { Point2D point; double bulge=0.0; bool hasY=false; };
            std::vector<Vertex> vertices;
            for (const auto& group:entity.groups) {
                if (group.code==10) {
                    double x=0.0;
                    if (!parseDouble(group.value,x)) {
                        result.issues.push_back({DxfIssueCode2D::InvalidNumericValue,
                                                 group.line,entity.type,entity.layer,
                                                 "invalid LWPOLYLINE x coordinate"});
                        return result;
                    }
                    vertices.push_back({{x,0.0},0.0,false});
                } else if (group.code==20) {
                    double y=0.0;
                    if (vertices.empty() || vertices.back().hasY ||
                        !parseDouble(group.value,y)) {
                        result.issues.push_back({DxfIssueCode2D::MalformedGroupPair,
                                                 group.line,entity.type,entity.layer,
                                                 "LWPOLYLINE y coordinate is out of order"});
                        return result;
                    }
                    vertices.back().point.y=y;
                    vertices.back().hasY=true;
                } else if (group.code==42) {
                    double bulge=0.0;
                    if (vertices.empty() || !parseDouble(group.value,bulge)) {
                        result.issues.push_back({DxfIssueCode2D::InvalidNumericValue,
                                                 group.line,entity.type,entity.layer,
                                                 "invalid LWPOLYLINE bulge"});
                        return result;
                    }
                    vertices.back().bulge=bulge;
                } else if (group.code==40 || group.code==41) {
                    double width=0.0;
                    if (!parseDouble(group.value,width) ||
                        !nearZero(width,options.endpointWeldTolerance)) {
                        result.issues.push_back({DxfIssueCode2D::NonZeroWidthOrThickness,
                                                 group.line,entity.type,entity.layer,
                                                 "variable-width LWPOLYLINE is unsupported"});
                        return result;
                    }
                }
            }
            if (*expected<0 || static_cast<unsigned long long>(*expected)!=vertices.size() ||
                vertices.size()<2 ||
                std::any_of(vertices.begin(),vertices.end(),[](const auto& vertex) {
                    return !vertex.hasY;
                })) {
                result.issues.push_back({DxfIssueCode2D::MissingRequiredGroup,entity.line,
                                         entity.type,entity.layer,
                                         "LWPOLYLINE vertex count/group structure is invalid"});
                return result;
            }
            const bool closed=(flags&1LL)!=0;
            std::vector<Point2D> points{vertices.front().point};
            const std::size_t segmentCount=closed?vertices.size():vertices.size()-1;
            for (std::size_t i=0;i<segmentCount;++i) {
                const std::size_t next=(i+1)%vertices.size();
                if (!appendBulgeSegment(points,vertices[i].point,vertices[next].point,
                        vertices[i].bulge,options,result,entity)) return result;
            }
            if (closed) {
                if (points.size()>1 && nearlyEqual(points.front(),points.back(),
                        TolerancePolicy{options.endpointWeldTolerance,0.0})) points.pop_back();
                closedLoops.push_back(std::move(points));
            } else {
                openChains.push_back({std::move(points),entity.line,entity.type,entity.layer});
            }
            ++result.report.lightweightPolylineCount;
        } else if (options.rejectUnsupportedEntities) {
            result.issues.push_back({DxfIssueCode2D::UnsupportedEntity,entity.line,
                                     entity.type,entity.layer,
                                     "DXF-1 supports LINE, LWPOLYLINE, ARC and CIRCLE only"});
            return result;
        }
    }
    result.report.sourceEntityCount=result.report.lineEntityCount+
        result.report.lightweightPolylineCount+result.report.arcEntityCount+
        result.report.circleEntityCount;
    result.report.layers.assign(layers.begin(),layers.end());
    if (result.report.sourceEntityCount==0) {
        result.issues.push_back({DxfIssueCode2D::InvalidBoundaryRegion,0,"","",
                                 "DXF ENTITIES section contains no supported boundary entities"});
        return result;
    }
    if (!assembleOpenChains(openChains,options.endpointWeldTolerance,closedLoops,result))
        return result;

    std::vector<BoundaryLoop> loops;
    loops.reserve(closedLoops.size());
    for (auto& vertices:closedLoops) loops.emplace_back(std::move(vertices));
    BoundaryRegion2D boundary(std::move(loops));
    const auto diagnostics=boundary.diagnose(tol);
    if (!diagnostics.valid() || !boundary.normalizeAlternating(tol)) {
        for (const auto& issue:diagnostics.issues) {
            result.issues.push_back({DxfIssueCode2D::InvalidBoundaryRegion,0,"","",
                                     issue.message});
        }
        if (diagnostics.issues.empty()) {
            result.issues.push_back({DxfIssueCode2D::InvalidBoundaryRegion,0,"","",
                                     "failed to normalize DXF loop orientation/nesting"});
        }
        return result;
    }
    std::vector<std::vector<Point2D>> canonical;
    canonical.reserve(boundary.loops().size());
    for (const auto& loop:boundary.loops()) {
        auto vertices=loop.vertices();
        canonicalRotate(vertices);
        canonical.push_back(std::move(vertices));
    }
    std::sort(canonical.begin(),canonical.end(),[](const auto& lhs,const auto& rhs) {
        return std::lexicographical_compare(lhs.begin(),lhs.end(),rhs.begin(),rhs.end(),
                                            pointLexLess);
    });
    std::vector<BoundaryLoop> canonicalLoops;
    canonicalLoops.reserve(canonical.size());
    for (auto& vertices:canonical) canonicalLoops.emplace_back(std::move(vertices));
    result.boundary.emplace(std::move(canonicalLoops));
    if (!result.boundary->diagnose(tol).valid()) {
        result.boundary.reset();
        result.issues.push_back({DxfIssueCode2D::InvalidBoundaryRegion,0,"","",
                                 "canonical DXF boundary failed final diagnostics"});
        return result;
    }
    result.report.outputLoopCount=result.boundary->loops().size();
    for (const auto& loop:result.boundary->loops())
        result.report.outputVertexCount+=loop.vertices().size();
    return result;
}

bool writeBoundaryXy2D(const BoundaryRegion2D& boundary,
                       const std::filesystem::path& path,std::string* error) {
    if (!boundary.diagnose().valid()) {
        if (error) *error="cannot write invalid boundary region";
        return false;
    }
    const auto parent=path.parent_path();
    std::error_code ec;
    if (!parent.empty()) std::filesystem::create_directories(parent,ec);
    if (ec) {
        if (error) *error="cannot create output directory: "+ec.message();
        return false;
    }
    std::ofstream out(path);
    if (!out) {
        if (error) *error="cannot open boundary output: "+path.string();
        return false;
    }
    out<<"# normalized from ASCII DXF by cartmesh2d\n"<<std::setprecision(17);
    for (std::size_t loop=0;loop<boundary.loops().size();++loop) {
        if (loop!=0) out<<'\n';
        for (const auto& point:boundary.loops()[loop].vertices())
            out<<point.x<<' '<<point.y<<'\n';
    }
    if (!out.good()) {
        if (error) *error="failed while writing boundary output";
        return false;
    }
    return true;
}

std::string dxfImportReportToJson(const DxfImportResult2D& result,int indentSpaces) {
    const std::string indent(static_cast<std::size_t>(std::max(0,indentSpaces)),' ');
    std::ostringstream out;
    out<<std::setprecision(17)<<"{\n"
       <<indent<<"\"valid\": "<<(result.valid()?"true":"false")<<",\n"
       <<indent<<"\"source_entity_count\": "<<result.report.sourceEntityCount<<",\n"
       <<indent<<"\"line_entity_count\": "<<result.report.lineEntityCount<<",\n"
       <<indent<<"\"lwpolyline_entity_count\": "
       <<result.report.lightweightPolylineCount<<",\n"
       <<indent<<"\"arc_entity_count\": "<<result.report.arcEntityCount<<",\n"
       <<indent<<"\"circle_entity_count\": "<<result.report.circleEntityCount<<",\n"
       <<indent<<"\"output_loop_count\": "<<result.report.outputLoopCount<<",\n"
       <<indent<<"\"output_vertex_count\": "<<result.report.outputVertexCount<<",\n"
       <<indent<<"\"sampled_arc_segment_count\": "
       <<result.report.sampledArcSegmentCount<<",\n"
       <<indent<<"\"maximum_chord_error\": "<<result.report.maximumChordError<<",\n"
       <<indent<<"\"endpoint_weld_tolerance\": "
       <<result.report.endpointWeldTolerance<<",\n"
       <<indent<<"\"source_plane_z\": ";
    if (result.report.sourcePlaneZDefined) out<<result.report.sourcePlaneZ;
    else out<<"null";
    out<<",\n"<<indent<<"\"insertion_units_code\": ";
    if (result.report.insertionUnitsCodeDefined) out<<result.report.insertionUnitsCode;
    else out<<"null";
    out<<",\n"<<indent<<"\"layers\": [";
    for (std::size_t i=0;i<result.report.layers.size();++i) {
        if (i!=0) out<<", ";
        out<<'"'<<jsonEscape(result.report.layers[i])<<'"';
    }
    out<<"],\n"<<indent<<"\"issues\": [";
    for (std::size_t i=0;i<result.issues.size();++i) {
        if (i!=0) out<<',';
        const auto& issue=result.issues[i];
        out<<"\n"<<indent<<indent<<"{\"code\": \""<<issueCodeName(issue.code)
           <<"\", \"line\": "<<issue.line<<", \"entity\": \""
           <<jsonEscape(issue.entity)<<"\", \"layer\": \""
           <<jsonEscape(issue.layer)<<"\", \"message\": \""
           <<jsonEscape(issue.message)<<"\"}";
    }
    if (!result.issues.empty()) out<<'\n'<<indent;
    out<<"]\n}\n";
    return out.str();
}

} // namespace cartmesh2d
