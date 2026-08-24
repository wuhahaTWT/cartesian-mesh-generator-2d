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
        const auto first=directory/"rectangle-a.xy";
        const auto second=directory/"rectangle-b.xy";
        check(writeBoundaryXy2D(*rectangleResult.boundary,first),"write first XY");
        const auto repeat=readAsciiDxfBoundary2D(rectangle,options);
        check(repeat.valid() && writeBoundaryXy2D(*repeat.boundary,second),
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

    const auto unsupported=writeText(directory,"spline.dxf",prefix+
        "0\nSPLINE\n8\nwall\n70\n8\n"+suffix);
    const auto unsupportedResult=readAsciiDxfBoundary2D(unsupported,options);
    check(!unsupportedResult.valid(),"unsupported SPLINE fails closed");
    check(!unsupportedResult.issues.empty() &&
          unsupportedResult.issues.front().code==DxfIssueCode2D::UnsupportedEntity,
          "unsupported entity issue is explicit");

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
    std::cout<<"cartmesh2d DXF-1 tests passed\n";
    return EXIT_SUCCESS;
}
