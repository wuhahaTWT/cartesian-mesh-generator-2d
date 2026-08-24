#include "cartmesh2d/io/Dxf2D.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <sstream>
#include <string>

using namespace cartmesh2d;

namespace {

int failures=0;

void check(bool condition,const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr<<"FAIL: "<<message<<'\n';
    }
}

void near(double actual,double expected,double epsilon,const std::string& message) {
    check(std::abs(actual-expected)<=epsilon,message);
}

std::filesystem::path writeText(const std::filesystem::path& directory,
                                const std::string& name,const std::string& text) {
    const auto path=directory/name;
    std::ofstream out(path);
    out<<text;
    return path;
}

std::string readText(const std::filesystem::path& path) {
    std::ifstream in(path);
    std::ostringstream out;
    out<<in.rdbuf();
    return out.str();
}

const std::string prefix="0\nSECTION\n2\nHEADER\n9\n$INSUNITS\n70\n6\n"
                         "0\nENDSEC\n0\nSECTION\n2\nENTITIES\n999\nDXF test\n";
const std::string prefixMm="0\nSECTION\n2\nHEADER\n9\n$INSUNITS\n70\n4\n"
                           "0\nENDSEC\n0\nSECTION\n2\nENTITIES\n";
const std::string prefixUnitless="0\nSECTION\n2\nHEADER\n0\nENDSEC\n"
                                 "0\nSECTION\n2\nENTITIES\n";
const std::string suffix="0\nENDSEC\n0\nEOF\n";

} // namespace

int main() {
    const auto directory=std::filesystem::temp_directory_path()/"cartmesh2d_dxf_tests";
    std::error_code ec;
    std::filesystem::remove_all(directory,ec);
    std::filesystem::create_directories(directory,ec);
    check(!ec,"create DXF test directory");

    const auto rectangle=writeText(directory,"rectangle.dxf",prefix+
        "0\nLWPOLYLINE\n8\nwall\n90\n4\n70\n1\n"
        "10\n0\n20\n0\n10\n2\n20\n0\n10\n2\n20\n1\n10\n0\n20\n1\n"+suffix);
    DxfImportOptions2D options;
    options.maximumChordError=0.01;
    const auto rectangleResult=readAsciiDxfBoundary2D(rectangle,options);
    check(rectangleResult.valid(),"closed LWPOLYLINE imports");
    if (rectangleResult.valid()) {
        check(rectangleResult.boundary->loops().size()==1,"rectangle has one loop");
        check(rectangleResult.boundary->loops()[0].vertices().size()==4,
              "rectangle retains four vertices");
        near(rectangleResult.boundary->area(),2.0,1e-12,"rectangle area exact");
        check(rectangleResult.report.layers==std::vector<std::string>{"wall"},
              "DXF layer reported deterministically");
        check(rectangleResult.report.insertionUnitsCodeDefined &&
              rectangleResult.report.insertionUnitsCode==6,
              "DXF insertion units code is preserved in the report");
        check(rectangleResult.report.effectiveUnitsCode==6 &&
              rectangleResult.report.coordinateScaleToMetres==1.0,
              "metre DXF units are applied");
        check(rectangleResult.embeddedPatches.size()==1 &&
              rectangleResult.embeddedPatches[0].name=="wall" &&
              rectangleResult.embeddedPatches[0].role==BoundaryConditionRole2D::Wall,
              "wall layer maps to OpenFOAM wall metadata");
        const auto first=directory/"rectangle-a.xy";
        const auto second=directory/"rectangle-b.xy";
        check(writeBoundaryXy2D(*rectangleResult.boundary,first,nullptr,
                                rectangleResult.embeddedPatches),"write first XY");
        const auto repeat=readAsciiDxfBoundary2D(rectangle,options);
        check(repeat.valid() && writeBoundaryXy2D(*repeat.boundary,second,nullptr,
                                                  repeat.embeddedPatches),
              "repeat import/write succeeds");
        check(readText(first)==readText(second),"repeat DXF conversion is byte-identical");
    }

    const auto circle=writeText(directory,"circle.dxf",prefix+
        "0\nCIRCLE\n8\nsolid\n10\n0\n20\n0\n30\n0\n40\n1\n"+suffix);
    const auto circleResult=readAsciiDxfBoundary2D(circle,options);
    check(circleResult.valid(),"CIRCLE imports");
    if (circleResult.valid()) {
        const auto& vertices=circleResult.boundary->loops()[0].vertices();
        check(vertices.size()>=16,"circle uses a resolved polygon");
        near(circleResult.boundary->area(),std::numbers::pi,0.05,
             "circle polygon area is within chord-controlled error");
        double maximumSagitta=0.0;
        for (std::size_t i=0;i<vertices.size();++i) {
            const auto& a=vertices[i];
            const auto& b=vertices[(i+1)%vertices.size()];
            const double midpointRadius=std::hypot(0.5*(a.x+b.x),0.5*(a.y+b.y));
            maximumSagitta=std::max(maximumSagitta,1.0-midpointRadius);
        }
        check(maximumSagitta<=options.maximumChordError+1e-12,
              "CIRCLE tessellation respects maximum chord error");
    }

    const auto lineLoop=writeText(directory,"line_loop.dxf",prefix+
        "0\nLINE\n8\nwall\n10\n1\n20\n0\n11\n1\n21\n1\n"
        "0\nLINE\n8\nwall\n10\n0\n20\n1\n11\n0\n21\n0\n"
        "0\nLINE\n8\nwall\n10\n1\n20\n1\n11\n0\n21\n1\n"
        "0\nLINE\n8\nwall\n10\n0\n20\n0\n11\n1\n21\n0\n"+suffix);
    const auto lineResult=readAsciiDxfBoundary2D(lineLoop,options);
    check(lineResult.valid(),"unordered/reversed LINE entities assemble into a loop");
    if (lineResult.valid()) near(lineResult.boundary->area(),1.0,1e-12,"LINE loop area exact");

    const auto bulge=writeText(directory,"bulge.dxf",prefix+
        "0\nLWPOLYLINE\n8\nwall\n90\n2\n70\n1\n"
        "10\n-1\n20\n0\n42\n1\n10\n1\n20\n0\n42\n0\n"+suffix);
    const auto bulgeResult=readAsciiDxfBoundary2D(bulge,options);
    check(bulgeResult.valid(),"LWPOLYLINE bulge imports");
    if (bulgeResult.valid()) near(bulgeResult.boundary->area(),0.5*std::numbers::pi,0.03,
                                  "bulge semicircle area is correct");

    const auto open=writeText(directory,"open.dxf",prefix+
        "0\nLINE\n8\nwall\n10\n0\n20\n0\n11\n1\n21\n0\n"+suffix);
    const auto openResult=readAsciiDxfBoundary2D(open,options);
    check(!openResult.valid(),"open LINE boundary fails closed");
    check(!openResult.issues.empty() &&
          openResult.issues.front().code==DxfIssueCode2D::OpenOrBranchedBoundary,
          "open LINE reports endpoint degree issue");

    const auto selfIntersecting=writeText(directory,"bowtie.dxf",prefix+
        "0\nLWPOLYLINE\n8\nwall\n90\n4\n70\n1\n"
        "10\n0\n20\n0\n10\n1\n20\n1\n10\n0\n20\n1\n10\n1\n20\n0\n"+suffix);
    const auto bowtieResult=readAsciiDxfBoundary2D(selfIntersecting,options);
    check(!bowtieResult.valid(),"self-intersecting DXF loop fails closed");
    check(!bowtieResult.issues.empty() &&
          bowtieResult.issues.front().code==DxfIssueCode2D::InvalidBoundaryRegion,
          "self-intersection is reported as invalid region");

    const auto unsupported=writeText(directory,"hatch.dxf",prefix+
        "0\nHATCH\n8\nwall\n70\n0\n"+suffix);
    const auto unsupportedResult=readAsciiDxfBoundary2D(unsupported,options);
    check(!unsupportedResult.valid(),"unsupported HATCH fails closed");
    check(!unsupportedResult.issues.empty() &&
          unsupportedResult.issues.front().code==DxfIssueCode2D::UnsupportedEntity,
          "unsupported entity issue is explicit");

    const auto ellipse=writeText(directory,"ellipse_mm.dxf",prefixMm+
        "0\nELLIPSE\n8\nwall_ellipse\n10\n500\n20\n0\n30\n0\n"
        "11\n500\n21\n0\n31\n0\n40\n0.2\n41\n0\n42\n6.283185307179586\n"+suffix);
    auto curveOptions=options;
    curveOptions.maximumChordError=5.0e-4;
    const auto ellipseResult=readAsciiDxfBoundary2D(ellipse,curveOptions);
    check(ellipseResult.valid(),"millimetre ELLIPSE imports");
    if (ellipseResult.valid()) {
        near(ellipseResult.boundary->area(),0.05*std::numbers::pi,0.003,
             "ELLIPSE area and millimetre-to-metre scaling are correct");
        near(ellipseResult.boundary->bounds().max.x,1.0,1e-12,
             "ELLIPSE maximum x is converted to metres");
        check(ellipseResult.report.effectiveUnitsCode==4 &&
              ellipseResult.report.coordinateScaleToMetres==1.0e-3,
              "millimetre conversion is reported");
        check(ellipseResult.report.ellipseEntityCount==1 &&
              ellipseResult.report.sampledEllipseSegmentCount>4,
              "ELLIPSE is chord-error sampled");
    }

    const std::string splineKnots=
        "40\n0\n40\n0\n40\n0\n40\n0\n40\n1\n40\n1\n40\n1\n40\n1\n";
    const auto splineAirfoil=writeText(directory,"spline_airfoil_mm.dxf",prefixMm+
        "0\nSPLINE\n8\nwall_airfoil\n70\n8\n71\n3\n72\n8\n73\n4\n74\n0\n"+
        splineKnots+
        "10\n0\n20\n0\n30\n0\n10\n250\n20\n180\n30\n0\n"
        "10\n750\n20\n100\n30\n0\n10\n1000\n20\n0\n30\n0\n"
        "0\nSPLINE\n8\nwall_airfoil\n70\n8\n71\n3\n72\n8\n73\n4\n74\n0\n"+
        splineKnots+
        "10\n1000\n20\n0\n30\n0\n10\n750\n20\n-80\n30\n0\n"
        "10\n250\n20\n-120\n30\n0\n10\n0\n20\n0\n30\n0\n"+suffix);
    const auto splineResult=readAsciiDxfBoundary2D(splineAirfoil,curveOptions);
    check(splineResult.valid(),"two open cubic SPLINE entities form a closed airfoil");
    if (splineResult.valid()) {
        check(splineResult.report.splineEntityCount==2 &&
              splineResult.report.sampledSplineSegmentCount>4,
              "SPLINE entities are adaptively sampled");
        check(splineResult.embeddedPatches.size()==1 &&
              splineResult.embeddedPatches[0].name=="wall_airfoil",
              "SPLINE loop layer reaches patch metadata");
        near(splineResult.boundary->bounds().max.x,1.0,1e-12,
             "SPLINE control points are converted from millimetres to metres");
    }

    const auto mixedLayers=writeText(directory,"mixed_layers.dxf",prefix+
        "0\nLINE\n8\nwall_a\n10\n0\n20\n0\n11\n1\n21\n0\n"
        "0\nLINE\n8\nwall_b\n10\n1\n20\n0\n11\n1\n21\n1\n"
        "0\nLINE\n8\nwall_a\n10\n1\n20\n1\n11\n0\n21\n1\n"
        "0\nLINE\n8\nwall_a\n10\n0\n20\n1\n11\n0\n21\n0\n"+suffix);
    const auto mixedResult=readAsciiDxfBoundary2D(mixedLayers,options);
    check(!mixedResult.valid(),"mixed boundary layers on one loop fail closed");
    check(!mixedResult.issues.empty() &&
          mixedResult.issues.front().code==DxfIssueCode2D::BoundaryMetadataConflict,
          "mixed boundary roles report metadata conflict");

    const auto collidingLayers=writeText(directory,"colliding_layers.dxf",prefix+
        "0\nCIRCLE\n8\nwall-a\n10\n-2\n20\n0\n40\n0.5\n"
        "0\nCIRCLE\n8\nwall_a\n10\n2\n20\n0\n40\n0.5\n"+suffix);
    const auto collidingResult=readAsciiDxfBoundary2D(collidingLayers,options);
    check(!collidingResult.valid(),
          "distinct DXF layers that sanitize to one patch name fail closed");
    check(!collidingResult.issues.empty() &&
          collidingResult.issues.front().code==DxfIssueCode2D::BoundaryMetadataConflict,
          "sanitized patch-name collision reports metadata conflict");

    const auto unitless=writeText(directory,"unitless.dxf",prefixUnitless+
        "0\nCIRCLE\n8\nwall\n10\n0\n20\n0\n40\n1\n"+suffix);
    const auto unitlessResult=readAsciiDxfBoundary2D(unitless,options);
    check(!unitlessResult.valid(),"unitless DXF fails without an explicit override");
    check(!unitlessResult.issues.empty() &&
          unitlessResult.issues.front().code==DxfIssueCode2D::UnknownOrUnitlessUnits,
          "unitless failure is explicit");
    auto overrideOptions=curveOptions;
    overrideOptions.sourceUnitsOverrideCode=4;
    const auto overrideResult=readAsciiDxfBoundary2D(unitless,overrideOptions);
    check(overrideResult.valid() && overrideResult.report.unitsOverrideApplied,
          "explicit millimetre override imports a unitless DXF");

    const auto inletLayer=writeText(directory,"inlet_layer.dxf",prefix+
        "0\nCIRCLE\n8\ninlet_feed\n10\n0\n20\n0\n40\n1\n"+suffix);
    const auto inletResult=readAsciiDxfBoundary2D(inletLayer,options);
    check(inletResult.valid() && inletResult.embeddedPatches.size()==1 &&
          inletResult.embeddedPatches[0].name=="inlet_feed" &&
          inletResult.embeddedPatches[0].type=="patch" &&
          inletResult.embeddedPatches[0].role==BoundaryConditionRole2D::Inlet,
          "inlet_* layer maps to inlet patch metadata");

    const auto symmetryLayer=writeText(directory,"symmetry_layer.dxf",prefix+
        "0\nCIRCLE\n8\nsymmetry_axis\n10\n0\n20\n0\n40\n1\n"+suffix);
    const auto symmetryResult=readAsciiDxfBoundary2D(symmetryLayer,options);
    check(symmetryResult.valid() && symmetryResult.embeddedPatches.size()==1 &&
          symmetryResult.embeddedPatches[0].type=="symmetryPlane" &&
          symmetryResult.embeddedPatches[0].role==BoundaryConditionRole2D::Symmetry,
          "symmetry_* layer maps to symmetryPlane metadata");

    const auto extrusion=writeText(directory,"extrusion.dxf",prefix+
        "0\nCIRCLE\n8\nwall\n10\n0\n20\n0\n40\n1\n210\n1\n220\n0\n230\n0\n"+suffix);
    const auto extrusionResult=readAsciiDxfBoundary2D(extrusion,options);
    check(!extrusionResult.valid(),"non-default OCS extrusion fails closed");
    check(!extrusionResult.issues.empty() &&
          extrusionResult.issues.front().code==DxfIssueCode2D::UnsupportedExtrusion,
          "unsupported extrusion issue is explicit");

    const auto nonPlanar=writeText(directory,"nonplanar.dxf",prefix+
        "0\nLINE\n8\nwall\n10\n0\n20\n0\n30\n0\n11\n1\n21\n0\n31\n1\n"+suffix);
    const auto nonPlanarResult=readAsciiDxfBoundary2D(nonPlanar,options);
    check(!nonPlanarResult.valid(),"non-planar LINE fails closed");
    check(!nonPlanarResult.issues.empty() &&
          nonPlanarResult.issues.front().code==DxfIssueCode2D::NonPlanarEntity,
          "non-planar LINE issue is explicit");

    const auto binary=writeText(directory,"binary.dxf",
        "AutoCAD Binary DXF\r\n\x1a\0ignored");
    const auto binaryResult=readAsciiDxfBoundary2D(binary,options);
    check(!binaryResult.valid(),"binary DXF fails closed");
    check(!binaryResult.issues.empty() &&
          binaryResult.issues.front().code==DxfIssueCode2D::BinaryDxfUnsupported,
          "binary DXF issue is explicit");

    const auto malformed=writeText(directory,"malformed.dxf","0\nSECTION\n2\n");
    const auto malformedResult=readAsciiDxfBoundary2D(malformed,options);
    check(!malformedResult.valid(),"odd ASCII group line fails closed");
    check(!malformedResult.issues.empty() &&
          malformedResult.issues.front().code==DxfIssueCode2D::MalformedGroupPair,
          "malformed group-pair issue is explicit");

    const auto missingEntities=writeText(directory,"missing_entities.dxf",
        "0\nSECTION\n2\nHEADER\n0\nENDSEC\n0\nEOF\n");
    const auto missingEntitiesResult=readAsciiDxfBoundary2D(missingEntities,options);
    check(!missingEntitiesResult.valid(),"missing ENTITIES section fails closed");
    check(!missingEntitiesResult.issues.empty() &&
          missingEntitiesResult.issues.front().code==DxfIssueCode2D::MissingEntitiesSection,
          "missing ENTITIES issue is explicit");

    if (failures!=0) {
        std::cerr<<failures<<" DXF test(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout<<"cartmesh2d DXF-2 tests passed\n";
    return EXIT_SUCCESS;
}
