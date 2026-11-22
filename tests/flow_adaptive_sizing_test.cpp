#include "cartmesh2d/fv/FlowAdaptiveSizing2D.hpp"
#include "cartmesh2d/sizing/SizeField2D.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

using namespace cartmesh2d;

namespace {

int failures = 0;

void check(bool pass, const char* message) {
    if (!pass) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

CutCell2D makeQuadCell(std::size_t id, std::vector<Point2D> pts) {
    CutCell2D c;
    c.sourceId = id;
    c.sourceKey = id;
    c.sourceLineage = {id};
    c.kind = CutCellKind::Cut;
    c.fluidPolygon = {std::move(pts)};
    c.area = c.fluidPolygon.area();
    c.areaFraction = 1.0;
    c.centroid = c.fluidPolygon.centroid();
    c.backgroundBounds = c.fluidPolygon.bounds();
    return c;
}

TopologyMesh2D buildTransformedGrid(std::size_t nx, std::size_t ny, double scale, double angleRad) {
    const double c = std::cos(angleRad);
    const double s = std::sin(angleRad);
    const auto xform = [&](double x, double y) -> Point2D {
        return {scale * (c * x - s * y), scale * (s * x + c * y)};
    };
    const std::vector<Point2D> outer{
        xform(0.0, 0.0), xform(1.0, 0.0), xform(1.0, 1.0), xform(0.0, 1.0)
    };
    const BoundaryRegion2D boundary{BoundaryLoop(outer)};
    const Domain2D domain{{{-3.0 * scale, -3.0 * scale}, {3.0 * scale, 3.0 * scale}}};
    std::vector<CutCell2D> cells;
    cells.reserve(nx * ny);
    const double dx = 1.0 / static_cast<double>(nx);
    const double dy = 1.0 / static_cast<double>(ny);
    std::size_t id = 0;
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const double x0 = static_cast<double>(i) * dx;
            const double x1 = static_cast<double>(i + 1U) * dx;
            const double y0 = static_cast<double>(j) * dy;
            const double y1 = static_cast<double>(j + 1U) * dy;
            cells.push_back(makeQuadCell(id++, {
                xform(x0, y0), xform(x1, y0), xform(x1, y1), xform(x0, y1)
            }));
        }
    }
    return buildGlobalTopology(cells, domain, boundary);
}

} // namespace

int main() {
    // 1. Uniform flow produces zero refinement boxes.
    const auto topo = buildTransformedGrid(10U, 10U, 1.0, 0.0);
    const auto fv = fv::makeFvMesh2D(topo);
    const std::size_t n = fv.cells.size();
    {
        std::vector<double> u(n, 2.5), v(n, -0.5), p(n, 1.0);
        fv::FlowAdaptivePolicy2D policy;
        policy.embeddedWallsNoSlip = false;
        const auto res = fv::suggestFlowAdaptiveRefinement2D(topo, fv, u, v, p, policy);
        check(res.valid(), "uniform flow suggestion should be valid");
        check(res.selectedCellCount == 0U && res.boxes.empty(),
              "uniform flow must not produce spurious refinement boxes");
    }

    // 2. Localized shear layer + pressure feature produces tight boxes that round-trip via JSON.
    {
        std::vector<double> u(n, 0.0), v(n, 0.0), p(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const auto pt = fv.cells[i].centre;
            const double r2 = (pt.x - 0.72) * (pt.x - 0.72) + (pt.y - 0.50) * (pt.y - 0.50);
            u[i] = std::tanh((pt.y - 0.50) / 0.05) * std::exp(-((pt.x - 0.55) * (pt.x - 0.55)) / 0.04);
            p[i] = 0.5 * std::exp(-r2 / 0.01);
        }
        fv::FlowAdaptivePolicy2D policy;
        policy.embeddedWallsNoSlip = false;
        policy.refineFraction = 0.15;
        policy.maxBoxes = 4U;
        policy.maxTargetLevel = 6U;
        policy.wallLevel = 7U;
        const auto res = fv::suggestFlowAdaptiveRefinement2D(topo, fv, u, v, p, policy);
        check(res.valid() && !res.boxes.empty() && res.boxes.size() <= 4U,
              "localized shear/pressure feature should emit 1..4 boxes");
        check(res.selectedCellCount > 0U && res.selectedAreaFraction < 0.25,
              "selected cells should stay concentrated on the feature");
        bool coversPeak = false;
        double unionBoxArea = 0.0;
        for (const auto& box : res.boxes) {
            const auto& b = box.region.bounds;
            unionBoxArea += (b.max.x - b.min.x) * (b.max.y - b.min.y);
            if (b.min.x <= 0.55 && b.max.x >= 0.55 && b.min.y <= 0.50 && b.max.y >= 0.50)
                coversPeak = true;
            check(box.region.targetLevel == 4U,
                  "level-3 uniform cells should upgrade by +1 to target level 4");
            check(box.levelsBelowWall == 3U,
                  "levelsBelowWall should match wallLevel(7) - targetLevel(4)");
        }
        check(coversPeak, "emitted boxes must cover the peak shear region (0.55, 0.50)");
        check(unionBoxArea < 0.65, "emitted boxes should stay localized rather than filling the domain");

        const std::string json1 = fv::flowAdaptiveRefinementToJson2D(res);
        const std::string json2 = fv::flowAdaptiveRefinementToJson2D(
            fv::suggestFlowAdaptiveRefinement2D(topo, fv, u, v, p, policy));
        check(json1 == json2, "adaptive refinement JSON must be deterministic");

        std::vector<BoxRefinementRegion2D> parsed;
        std::string parseError;
        check(parseBoxRefinementRegionsJson2D(json1, parsed, &parseError),
              "parseBoxRefinementRegionsJson2D must parse emitted JSON");
        check(parsed.size() == res.boxes.size(), "parsed box count must match emitted boxes");
        if (parsed.size() == res.boxes.size()) {
            for (std::size_t k = 0; k < parsed.size(); ++k) {
                const auto& a = parsed[k];
                const auto& b = res.boxes[k].region;
                check(a.targetLevel == b.targetLevel &&
                      std::abs(a.bounds.min.x - b.bounds.min.x) < 1e-12 &&
                      std::abs(a.bounds.max.y - b.bounds.max.y) < 1e-12,
                      "round-tripped box coordinates and targetLevel must match");
            }
        }
    }

    // 3. Scale and rotation invariance of the dimensionless indicator.
    {
        const double scale = 100.0;
        const double uScale = 25.0;
        const double angle = 0.63;
        const auto topoScaled = buildTransformedGrid(8U, 8U, scale, angle);
        const auto fvScaled = fv::makeFvMesh2D(topoScaled);
        const auto topoUnit = buildTransformedGrid(8U, 8U, 1.0, 0.0);
        const auto fvUnit = fv::makeFvMesh2D(topoUnit);

        const double c = std::cos(angle);
        const double s = std::sin(angle);
        const auto evalField = [&](const fv::FvMesh2D& m, double L, double U, bool rotated) {
            std::vector<double> uu(m.cells.size()), vv(m.cells.size()), pp(m.cells.size());
            for (std::size_t i = 0; i < m.cells.size(); ++i) {
                const auto pt = m.cells[i].centre;
                const double lx = rotated ? (c * pt.x + s * pt.y) / L : pt.x / L;
                const double ly = rotated ? (-s * pt.x + c * pt.y) / L : pt.y / L;
                const double dx = lx - 0.5;
                const double dy = ly - 0.5;
                const double r2 = dx * dx + dy * dy;
                const double e = std::exp(-r2 / 0.02);
                const double ul = -U * dy * e;
                const double vl = U * dx * e;
                uu[i] = rotated ? (c * ul - s * vl) : ul;
                vv[i] = rotated ? (s * ul + c * vl) : vl;
                pp[i] = U * U * e;
            }
            return std::tuple{std::move(uu), std::move(vv), std::move(pp)};
        };

        auto [u1, v1, p1] = evalField(fvUnit, 1.0, 1.0, false);
        auto [u2, v2, p2] = evalField(fvScaled, scale, uScale, true);
        fv::FlowAdaptivePolicy2D policy;
        policy.embeddedWallsNoSlip = false;
        const auto r1 = fv::suggestFlowAdaptiveRefinement2D(topoUnit, fvUnit, u1, v1, p1, policy);
        const auto r2 = fv::suggestFlowAdaptiveRefinement2D(topoScaled, fvScaled, u2, v2, p2, policy);
        check(r1.valid() && r2.valid(), "scaled/rotated suggestions must be valid");
        double maxDiff = 0.0;
        for (std::size_t i = 0; i < r1.cellIndicators.size(); ++i) {
            maxDiff = std::max(maxDiff, std::abs(r1.cellIndicators[i] - r2.cellIndicators[i]));
        }
        check(maxDiff < 1e-11, "dimensionless indicator must be scale- and rotation-invariant");
    }

    // 3b. Multi-physics scalar/thermal gradient drives localized refinement even
    // when velocity and pressure are completely uniform, and VTK / fields JSON export work.
    {
        std::vector<double> u(n, 1.0), v(n, 0.0), p(n, 0.0), theta(n, 300.0);
        for (std::size_t i = 0; i < n; ++i) {
            const auto pt = fv.cells[i].centre;
            const double r2 = (pt.x - 0.30) * (pt.x - 0.30) + (pt.y - 0.70) * (pt.y - 0.70);
            theta[i] = 300.0 + 50.0 * std::exp(-r2 / 0.01);
        }
        fv::FlowAdaptivePolicy2D policy;
        policy.embeddedWallsNoSlip = false;
        policy.refineFraction = 0.15;
        policy.maxBoxes = 4U;
        policy.maxTargetLevel = 6U;
        policy.wallLevel = 7U;
        const auto res = fv::suggestFlowAdaptiveRefinement2D(
            topo, fv, u, v, p, policy, {}, {}, {}, theta);
        check(res.valid() && !res.boxes.empty(), "localized thermal plume should emit adaptive boxes");
        check(res.referenceScalarSpan > 25.0 && res.referenceScalarSpan < 50.0,
              "referenceScalarSpan should reflect thermal range");
        check(res.cellSelected.size() == n && res.selectedCellCount > 0U,
              "cellSelected bitmask must match cell count");
        bool coversPlume = false;
        for (const auto& box : res.boxes) {
            const auto& b = box.region.bounds;
            if (b.min.x <= 0.30 && b.max.x >= 0.30 && b.min.y <= 0.70 && b.max.y >= 0.70)
                coversPlume = true;
        }
        check(coversPlume, "emitted thermal adaptive boxes must cover the plume peak (0.30, 0.70)");

        const std::string fieldsJson = fv::flowAdaptiveFieldsToJson2D(res);
        check(fieldsJson.find("\"cartmesh2d-adaptive-fields-v1\"") != std::string::npos &&
              fieldsJson.find("\"indicator\":") != std::string::npos &&
              fieldsJson.find("\"selected\":") != std::string::npos,
              "flowAdaptiveFieldsToJson2D must emit v1 indicator/selected cell records");

        const auto tmpVtk = std::filesystem::temp_directory_path() / "cartmesh2d_adaptive_test.vtk";
        std::string vtkErr;
        check(fv::writeFlowAdaptiveVtk2D(topo, res, tmpVtk, &vtkErr),
              "writeFlowAdaptiveVtk2D must succeed");
        {
            std::ifstream in(tmpVtk);
            std::ostringstream ss;
            ss << in.rdbuf();
            const std::string content = ss.str();
            check(content.find("SCALARS adaptive_indicator double 1") != std::string::npos &&
                  content.find("SCALARS adaptive_selected int 1") != std::string::npos,
                  "adaptive VTK must contain adaptive_indicator and adaptive_selected scalars");
        }
        std::error_code ec;
        std::filesystem::remove(tmpVtk, ec);
    }

    // 4. Fail-closed checks for non-finite fields, invalid policy, and malformed JSON.
    {
        std::vector<double> u(n, 0.0), v(n, 0.0), p(n, 0.0);
        u[0] = std::numeric_limits<double>::quiet_NaN();
        check(!fv::suggestFlowAdaptiveRefinement2D(topo, fv, u, v, p).valid(),
              "NaN velocity must fail closed");
        u[0] = 0.0;
        fv::FlowAdaptivePolicy2D badPolicy;
        badPolicy.refineFraction = 0.9;
        check(!fv::suggestFlowAdaptiveRefinement2D(topo, fv, u, v, p, badPolicy).valid(),
              "refineFraction > 0.5 must fail closed");

        std::vector<BoxRefinementRegion2D> dummy;
        std::string err;
        check(!parseBoxRefinementRegionsJson2D("{\"boxes\":[{\"xmin\":1,\"ymin\":0,\"xmax\":0,\"ymax\":1,\"target_level\":5}]}", dummy, &err),
              "inverted box bounds must be rejected by JSON parser");
        check(!parseBoxRefinementRegionsJson2D("{\"boxes\":[{\"xmin\":0,\"ymin\":0,\"xmax\":1,\"ymax\":1}]}", dummy, &err),
              "missing target_level must be rejected by JSON parser");
        check(!parseBoxRefinementRegionsJson2D("{}", dummy, &err),
              "missing boxes key must be rejected by JSON parser");
    }

    return failures ? 1 : 0;
}
