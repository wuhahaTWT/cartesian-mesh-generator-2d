#include "FlowSolverDetail2D.hpp"
#include <algorithm>
#include <map>
#include <limits>
#include "cartmesh2d/fv/ManufacturedFlow2D.hpp"
#include "cartmesh2d/fv/TaylorGreen2D.hpp"
#include "cartmesh2d/fv/FlowBoundaryIO2D.hpp"
namespace cartmesh2d::fv::solver_detail {
double counterflowSpeed(double y,double speed){return speed*(1+2*std::cos(2*std::acos(-1.)*y));}
Boundary boundaries(const FvMesh2D& m, const FlowControls2D& c) {
    const auto nf = m.faces.size();
    Boundary b;
    b.role.assign(nf, Role::Wall);
    b.u.resize(nf);
    b.v.resize(nf);
    b.fixedU.assign(nf, false);
    b.fixedV.assign(nf, false);
    b.fixedP.assign(nf, false);
    b.constantU.assign(nf, false);
    b.constantV.assign(nf, false);
    b.gradientGroups.resize(nf);
    b.xmin = b.ymin = std::numeric_limits<double>::infinity();
    b.xmax = b.ymax = -b.xmin;
    for (const auto& f : m.faces) {
        if (!f.neighbour) {
            // Curved duct extrema may occur at vertices, not face centres.
            // For a 2-D edge, rotating its area vector gives its tangent.
            const bool vertexBounds = c.scenario == "duct" || c.scenario == "custom";
            const double dx = vertexBounds ? .5 * std::abs(f.areaVector.y) : 0;
            const double dy = vertexBounds ? .5 * std::abs(f.areaVector.x) : 0;
            b.xmin = std::min(b.xmin, f.centre.x - dx);
            b.xmax = std::max(b.xmax, f.centre.x + dx);
            b.ymin = std::min(b.ymin, f.centre.y - dy);
            b.ymax = std::max(b.ymax, f.centre.y + dy);
        }
    }
    const double width = b.xmax - b.xmin;
    const double height = b.ymax - b.ymin;
    ensure(width > 0 && height > 0, "Flow empty domain");
    const double eps = TolerancePolicy{}.scale(std::max(width, height));
    const auto equal = [&](double a, double d) { return std::abs(a - d) <= eps; };
    if (c.scenario == "custom") {
        ensure(std::isfinite(c.speed) && c.speed > 0, "Custom boundary reference speed must be positive");
        ensure(c.outletBackflow == OutletBackflow2D::Reject,
               "Custom pressure outlets currently require explicit backflow rejection");
        std::vector<bool> seen(nf, false);
        std::map<std::string, FlowBoundaryKind2D> namedKinds;
        std::map<std::string, std::size_t> gradientGroups;
        b.p.resize(nf);
        std::size_t inletCount = 0, outletCount = 0, openingCount = 0;
        std::size_t velocityOutletCount = 0;
        double prescribedFlux = 0, prescribedFluxMagnitude = 0;
        double inletLength = 0, outletLength = 0;
        for (const auto& condition : c.boundaryConditions) {
            const auto id = condition.face;
            ensure(id < nf, "Custom boundary face ID is out of range");
            const auto& face = m.faces[id];
            ensure(!face.neighbour, "Custom boundary condition refers to an internal face");
            ensure(!seen[id], "Duplicate custom boundary face");
            seen[id] = true;
            ensure(!condition.name.empty() && condition.name.size() <= 128 &&
                   std::none_of(condition.name.begin(), condition.name.end(), [](unsigned char ch) {
                       return ch < 32 || ch == 127 || ch == ',' || ch == '"';
                   }), "Invalid custom boundary name");
            const auto [it, inserted] = namedKinds.emplace(condition.name, condition.kind);
            ensure(inserted || it->second == condition.kind,
                   "One named boundary cannot mix physical condition types");
            b.gradientGroups[id]=gradientGroups.emplace(condition.name,gradientGroups.size()).first->second;
            ensure(std::isfinite(condition.velocity.x) && std::isfinite(condition.velocity.y) &&
                   std::isfinite(condition.pressure), "Non-finite custom boundary value");
            const double length = std::hypot(face.areaVector.x, face.areaVector.y);
            const double q = finite(dot(condition.velocity, face.areaVector));
            switch (condition.kind) {
            case FlowBoundaryKind2D::VelocityOutlet:
                ensure(condition.pressure == 0 && q > 0,
                       "Velocity outlet must point out of the fluid and cannot prescribe pressure");
                // Role::Inlet is the internal prescribed-vector/flux stencil;
                // external boundary kind/name retain the outlet semantics.
                b.role[id] = Role::Inlet;
                b.fixedU[id] = b.fixedV[id] = true;
                b.u[id] = condition.velocity.x; b.v[id] = condition.velocity.y;
                prescribedFlux += q; prescribedFluxMagnitude += std::abs(q);
                ++velocityOutletCount;
                break;
            case FlowBoundaryKind2D::VelocityInlet:
                ensure(condition.pressure == 0 && q < 0,
                       "Velocity inlet must point into the fluid and cannot prescribe pressure");
                b.role[id] = Role::Inlet;
                b.fixedU[id] = b.fixedV[id] = true;
                b.u[id] = condition.velocity.x; b.v[id] = condition.velocity.y;
                b.initialU += condition.velocity.x * length;
                b.initialV += condition.velocity.y * length;
                inletLength += length;
                ++inletCount;
                prescribedFlux += q; prescribedFluxMagnitude += std::abs(q);
                break;
            case FlowBoundaryKind2D::PressureOutlet:
            case FlowBoundaryKind2D::PressureOpening:
                ensure(condition.velocity.x == 0 && condition.velocity.y == 0,
                       "Pressure boundary cannot also prescribe velocity");
                if (condition.kind == FlowBoundaryKind2D::PressureOpening) {
                    ensure(std::min(std::abs(face.areaVector.x), std::abs(face.areaVector.y)) <=
                               TolerancePolicy{}.scale(1.) * length,
                           "Pressure opening must be axis aligned");
                    ++openingCount;
                    b.pressureOpenings = true;
                }
                b.role[id] = condition.kind == FlowBoundaryKind2D::PressureOpening ? Role::Opening : Role::Outlet;
                b.fixedP[id] = true; b.p[id] = condition.pressure;
                b.initialP += condition.pressure * length;
                outletLength += length;
                ++outletCount;
                break;
            case FlowBoundaryKind2D::Symmetry:
                ensure(condition.velocity.x == 0 && condition.velocity.y == 0 && condition.pressure == 0,
                       "Symmetry cannot prescribe velocity or pressure");
                ensure(std::min(std::abs(face.areaVector.x), std::abs(face.areaVector.y)) <=
                           TolerancePolicy{}.scale(1.) * length,
                       "Symmetry face must be axis aligned");
                b.role[id] = Role::Slip;
                if (std::abs(face.areaVector.x) > std::abs(face.areaVector.y))
                    b.fixedU[id] = b.constantU[id] = true;
                else b.fixedV[id] = b.constantV[id] = true;
                break;
            case FlowBoundaryKind2D::Wall:
            case FlowBoundaryKind2D::MovingWall:
            case FlowBoundaryKind2D::SmoothMovingWall:
                ensure(condition.pressure == 0, "Wall cannot prescribe pressure");
                if (condition.kind == FlowBoundaryKind2D::Wall)
                    ensure(condition.velocity.x == 0 && condition.velocity.y == 0,
                           "Stationary wall velocity must be zero");
                ensure(std::abs(q) <= TolerancePolicy{}.scale(std::max(c.speed,
                           std::hypot(condition.velocity.x, condition.velocity.y))) * length,
                       "Moving wall velocity must be tangential to its face");
                b.role[id] = condition.kind == FlowBoundaryKind2D::Wall ? Role::Wall : Role::Lid;
                b.fixedU[id] = b.fixedV[id] = true;
                b.constantU[id] = b.constantV[id] = condition.kind != FlowBoundaryKind2D::SmoothMovingWall;
                b.u[id] = condition.velocity.x; b.v[id] = condition.velocity.y;
                break;
            default:
                throw std::runtime_error("Unknown custom boundary condition type");
            }
        }
        for (std::size_t id = 0; id < nf; ++id)
            ensure(m.faces[id].neighbour || seen[id], "Missing custom boundary face");
        ensure(openingCount > 0 || (inletCount > 0 && (outletCount > 0 || velocityOutletCount > 0)) ||
               (inletCount == 0 && outletCount == 0 && velocityOutletCount == 0),
               "Custom open flow needs a pressure opening or velocity inlet and pressure outlet; closed flow needs only walls");
        b.closed = outletCount == 0;
        if (b.closed && velocityOutletCount > 0)
            ensure(std::abs(prescribedFlux) <= TolerancePolicy{}.scale(prescribedFluxMagnitude),
                   "Fully prescribed port fluxes must balance for a pressure-gauge solve");
        if (!b.closed) {
            if (inletLength > 0) {
                b.initialU = finite(b.initialU / inletLength);
                b.initialV = finite(b.initialV / inletLength);
            }
            b.initialP = finite(b.initialP / outletLength);
        }
        return b;
    }
    ensure(c.boundaryConditions.empty(), "Explicit boundary conditions require the custom scenario");
    if(c.scenario=="flatplate")
        ensure(c.flatPlateLeadingEdge>=b.xmin && c.flatPlateLeadingEdge<b.xmax,
               "Flat plate leading edge must lie within the bottom boundary");

    std::size_t inlets = 0;
    std::size_t outlets = 0;
    std::size_t walls = 0;
    for (std::size_t id = 0; id < nf; ++id) {
        const auto& f = m.faces[id];
        if (f.neighbour) {
            continue;
        }
        const bool left = equal(f.centre.x, b.xmin);
        const bool right = equal(f.centre.x, b.xmax);
        const bool top = equal(f.centre.y, b.ymax);
        const bool bottom = equal(f.centre.y, b.ymin);
        const bool embedded =
            c.scenario == "external" && f.patch == BoundaryPatch2D::EmbeddedBoundary;
        if (embedded) {
            b.role[id] = Role::Wall;
            ++walls;
        } else if (c.scenario == "duct") {
            // A user-selected duct has planar x-end openings and arbitrary
            // stationary walls between them. Unlike the channel preset, its
            // fluid region need not fill the bounding rectangle.
            const bool vertical = std::abs(f.areaVector.y) <= eps;
            if (left && vertical && f.areaVector.x < 0) {
                b.role[id] = Role::Inlet;
                ++inlets;
            } else if (right && vertical && f.areaVector.x > 0) {
                b.role[id] = Role::Outlet;
                ++outlets;
            } else {
                b.role[id] = Role::Wall;
                ++walls;
            }
        } else {
            ensure(left || right || top || bottom,
                   "Selected flow case requires rectangular outer boundary");
            // An axis-position match alone is insufficient for a sloping edge.
            ensure((left || right) ? std::abs(f.areaVector.y) <= eps
                                   : std::abs(f.areaVector.x) <= eps,
                   "Flow outer edge is not axis aligned");
            if (c.scenario == "taylor-green") {
                b.role[id] = Role::Slip;
            } else if (c.scenario == "manufactured") {
                b.role[id] = Role::Wall;
            } else if (c.scenario == "cavity") {
                b.role[id] = top ? Role::Lid : Role::Wall;
            } else if (left) {
                b.role[id] = Role::Inlet;
                ++inlets;
            } else if (right) {
                b.role[id] = Role::Outlet;
                ++outlets;
            } else if(c.scenario=="flatplate") {
                b.role[id]=top && c.flatPlateTop==FlatPlateTop2D::PressureFarfield ? Role::Farfield : Role::Slip;
                if(bottom) {
                    const double halfLength=.5*std::abs(f.areaVector.y);
                    const double lo=f.centre.x-halfLength,hi=f.centre.x+halfLength;
                    ensure(!(lo<c.flatPlateLeadingEdge-eps && hi>c.flatPlateLeadingEdge+eps),
                           "Flat plate leading edge crosses a boundary face; split the mesh edge");
                    if(f.centre.x>=c.flatPlateLeadingEdge) {b.role[id]=Role::Wall;++walls;}
                }
            } else {
                b.role[id] = (c.scenario == "external" || c.scenario == "counterflow") ? Role::Slip : Role::Wall;
            }
        }

        b.gradientGroups[id]=static_cast<std::size_t>(b.role[id]);
        b.constantU[id]=b.role[id]==Role::Wall || b.role[id]==Role::Lid;
        b.constantV[id]=b.constantU[id];
        switch (b.role[id]) {
        case Role::Wall:
            b.fixedU[id] = b.fixedV[id] = true;
            break;
        case Role::Lid:
            b.fixedU[id] = b.fixedV[id] = true;
            b.u[id] = c.speed;
            break;
        case Role::Inlet:
            b.fixedU[id] = b.fixedV[id] = true;
            b.u[id] = c.scenario == "channel"
                          ? 4 * c.speed * (f.centre.y - b.ymin) * (b.ymax - f.centre.y) /
                                (height * height)
                          : c.speed;
            if (c.scenario == "counterflow") b.u[id] = counterflowSpeed(f.centre.y, c.speed);
            break;
        case Role::Outlet:
        case Role::Opening:
            b.fixedP[id] = true;
            break;
        case Role::Farfield:
            b.fixedP[id] = true;
            b.u[id] = c.speed;
            break;
        case Role::Slip:
            // Outer boundaries were checked to be axis aligned above.
            if (left || right) b.fixedU[id] = b.constantU[id] = true;
            else b.fixedV[id] = b.constantV[id] = true;
            break;
        }
    }
    b.closed = c.scenario == "cavity" || c.scenario == "manufactured" || c.scenario == "taylor-green";
    if (c.scenario == "manufactured" || c.scenario == "taylor-green" || c.scenario == "counterflow")
        ensure(equal(b.xmin,0) && equal(b.ymin,0) && equal(b.xmax,1) && equal(b.ymax,1),
               "Verification flow requires the unit square [0,1]^2");
    if (!b.closed) {
        ensure(inlets && outlets, "Flow needs left inlet and right pressure outlet");
    }
    if (c.scenario == "external") {
        ensure(walls > 0, "External case requires embedded solid wall");
    } else if (c.scenario == "duct") {
        ensure(walls > 0, "Duct requires stationary wall faces between its openings");
    } else {
        if(c.scenario=="flatplate")ensure(walls>0,"Flat plate requires resolved no-slip faces");
        double area = 0;
        for (const auto& cell : m.cells) {
            area += cell.area;
        }
        ensure(std::abs(area - width * height) <
                   TolerancePolicy{}.areaScale(std::max(width, height)) * 10,
               "Channel/cavity requires a filled rectangle without holes");
    }
    return b;
}

void updateOutletBoundary(Boundary& b, const FvMesh2D& m,
                          const FlowControls2D& c, const Vec& flux) {
    if(c.scenario!="flatplate" && c.scenario!="custom" && c.outletBackflow!=OutletBackflow2D::NormalInlet)return;
    for (std::size_t id = 0; id < m.faces.size(); ++id) {
        if (!m.faces[id].neighbour && b.role[id] == Role::Opening) {
            const bool normalX = std::abs(m.faces[id].areaVector.x) > std::abs(m.faces[id].areaVector.y);
            b.fixedU[id] = b.constantU[id] = !normalX && flux[id] < 0;
            b.fixedV[id] = b.constantV[id] = normalX && flux[id] < 0;
            continue;
        }
        if(!m.faces[id].neighbour && b.role[id]==Role::Farfield) {
            // Horizontal open top: pressure controls normal flux; only the
            // entering tangential component is prescribed from the freestream.
            b.fixedU[id]=b.constantU[id]=flux[id]<0;
            continue;
        }
        if (c.outletBackflow != OutletBackflow2D::NormalInlet) continue;
        if (m.faces[id].neighbour || b.role[id] != Role::Outlet) continue;
        // Current cases have a right, axis-aligned pressure outlet. Keep its
        // normal component free; constrain only the incoming tangential one.
        b.fixedV[id] = b.constantV[id] = flux[id] < 0;
        b.v[id] = 0;
    }
}

void initializeCaseSources(const FvMesh2D& m,const FlowControls2D& c,FlowResult2D& r){
    const auto n=m.cells.size();
    if (c.scenario == "manufactured") {
        r.sourceIntegrals.reserve(n);
        for (const auto& cell : m.cells) {
            const auto acceleration=manufacturedFlow2D(cell.centre,c.speed,c.nu,c.manufacturedPressureSlope,c.manufacturedViscositySlope,c.viscousStress==ViscousStress2D::Symmetric).acceleration;
            r.sourceIntegrals.push_back({finite(cell.area*acceleration.x),finite(cell.area*acceleration.y)});
        }
    }
    if (c.scenario == "counterflow") {
        const double pi=std::acos(-1.);
        for (const auto& cell : m.cells)
            r.sourceIntegrals.push_back({finite(cell.area*8*pi*pi*c.nu*c.speed*std::cos(2*pi*cell.centre.y)),0});
    }
}
void initializeCaseVelocity(const FvMesh2D& m,const FlowControls2D& c,const Boundary& b,FlowResult2D& r){
    const auto n=m.cells.size();const double h=b.ymax-b.ymin;
    for (std::size_t i = 0; i < n; ++i) {
        const double y = m.cells[i].centre.y;
        r.u[i] = b.closed ? 0
                          : (c.scenario == "channel"
                                 ? 4 * c.speed * (y - b.ymin) * (b.ymax - y) / (h * h)
                                 : c.speed);
        if (c.scenario == "counterflow") r.u[i] = counterflowSpeed(y,c.speed);
        if (c.scenario == "custom") {
            r.u[i] = b.initialU; r.v[i] = b.initialV; r.p[i] = b.initialP;
        }
    }
}
} // namespace cartmesh2d::fv::solver_detail
namespace cartmesh2d::fv {
using namespace solver_detail;
FlowState2D initialIncompressibleState2D(const FvMesh2D& m, const FlowControls2D& c) {
    validateFvMesh2D(m);
    ensure(c.flatPlateLeadingEdge==0 && c.flatPlateTop==FlatPlateTop2D::PressureFarfield,
           "Flat plate controls are not supported by transient initialization");
    ensure((c.scenario=="external" || c.scenario=="channel" || c.scenario=="duct" || c.scenario=="custom" || c.scenario=="cavity" || c.scenario=="taylor-green") &&
           std::isfinite(c.nu) && c.nu>0 && std::isfinite(c.speed) && c.speed>0,
           "Invalid transient initial-state controls");
    validateViscosity(m,c);
    const auto b=boundaries(m,c); (void)b;
    FlowState2D s;
    s.u.resize(m.cells.size());s.v.resize(m.cells.size());s.p.resize(m.cells.size());s.flux.resize(m.faces.size());
    if (c.scenario == "custom") std::fill(s.p.begin(), s.p.end(), b.initialP);
    if (c.scenario=="taylor-green") {
        ensure(std::isfinite(c.nu) && c.nu>0 && std::isfinite(c.speed) && c.speed>0,"Invalid vortex controls");
        const double gauge=taylorGreen2D(m.cells.front().centre,0,c.speed,c.nu).pressure;
        for (std::size_t i=0;i<m.cells.size();++i) {
            const auto q=taylorGreen2D(m.cells[i].centre,0,c.speed,c.nu);
            s.u[i]=q.velocity.x;s.v[i]=q.velocity.y;s.p[i]=q.pressure-gauge;
        }
        // Streamfunction differences give an exactly conservative integral flux.
        for (std::size_t id=0;id<m.faces.size();++id) {
            const auto& f=m.faces[id];
            if (!f.neighbour) continue;
            const Point2D a{f.centre.x+.5*f.areaVector.y,f.centre.y-.5*f.areaVector.x};
            const Point2D z{f.centre.x-.5*f.areaVector.y,f.centre.y+.5*f.areaVector.x};
            s.flux[id]=taylorGreen2D(z,0,c.speed,c.nu).streamfunction-taylorGreen2D(a,0,c.speed,c.nu).streamfunction;
        }
    }
    return s;
}

void validateFlowBoundaryConditions2D(const FvMesh2D& mesh, const FlowControls2D& controls) {
    (void)boundaries(mesh, controls);
}

std::vector<FlowBoundaryCondition2D> explicitFlowBoundaryPreset2D(
    const FvMesh2D& mesh, const FlowControls2D& controls) {
    validateFvMesh2D(mesh);
    if (controls.scenario == "custom") {
        validateFlowBoundaryConditions2D(mesh, controls);
        return controls.boundaryConditions;
    }
    if (controls.scenario == "annulus")
        return rotatingAnnulusBoundaryPreset2D(mesh, controls.speed);
    ensure(controls.scenario == "duct" || controls.scenario == "channel" || controls.scenario == "cavity",
           "Explicit boundary template currently supports duct, channel and cavity only");
    const auto b=boundaries(mesh, controls);
    std::vector<FlowBoundaryCondition2D> result;
    for (std::size_t id=0;id<mesh.faces.size();++id) {
        if (mesh.faces[id].neighbour) continue;
        FlowBoundaryCondition2D value{id,FlowBoundaryKind2D::Wall,{b.u[id],b.v[id]},0,"wall"};
        if (b.role[id]==Role::Inlet) {value.kind=FlowBoundaryKind2D::VelocityInlet;value.name="inlet";}
        else if (b.role[id]==Role::Outlet) {value.kind=FlowBoundaryKind2D::PressureOutlet;value.name="outlet";}
        else if (b.role[id]==Role::Lid) {value.kind=FlowBoundaryKind2D::MovingWall;value.name="lid";}
        else ensure(b.role[id]==Role::Wall,"Unsupported condition in explicit boundary template");
        result.push_back(std::move(value));
    }
    return result;
}
}
