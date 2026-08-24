#include "cartmesh2d/io/Dxf2D.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>

using namespace cartmesh2d;

namespace {

void usage() {
    std::cerr<<"usage: cartmesh2d_dxf_cli <input-ascii.dxf> <output.xy> "
               "<maximum-chord-error-metres> [report.json] "
               "[endpoint-weld-tolerance-metres=1e-10] "
               "[source-unit=auto|in|ft|mm|cm|m|km|1..24]\n";
}

std::optional<long long> unitCode(const std::string& value,bool& valid) {
    valid=true;
    if (value=="auto") return std::nullopt;
    if (value=="in" || value=="inch") return 1;
    if (value=="ft" || value=="foot") return 2;
    if (value=="mile") return 3;
    if (value=="mm") return 4;
    if (value=="cm") return 5;
    if (value=="m") return 6;
    if (value=="km") return 7;
    if (value=="yard") return 10;
    try {
        std::size_t consumed=0;
        const long long code=std::stoll(value,&consumed);
        if (consumed==value.size() && code>=1 && code<=24) return code;
    } catch (const std::exception&) {}
    valid=false;
    return std::nullopt;
}

} // namespace

int main(int argc,char** argv) {
    if (argc<4 || argc>7) {
        usage();
        return EXIT_FAILURE;
    }
    const std::filesystem::path input=argv[1];
    const std::filesystem::path output=argv[2];
    const std::filesystem::path report=argc>=5
        ?std::filesystem::path(argv[4])
        :std::filesystem::path(output.string()+".dxf.json");
    DxfImportOptions2D options;
    try {
        std::size_t consumed=0;
        const std::string chordText=argv[3];
        options.maximumChordError=std::stod(chordText,&consumed);
        if (consumed!=chordText.size()) throw std::invalid_argument("trailing chord token");
        if (argc>=6) {
            consumed=0;
            const std::string weldText=argv[5];
            options.endpointWeldTolerance=std::stod(weldText,&consumed);
            if (consumed!=weldText.size()) throw std::invalid_argument("trailing weld token");
        }
        if (argc>=7) {
            bool valid=false;
            options.sourceUnitsOverrideCode=unitCode(argv[6],valid);
            if (!valid) throw std::invalid_argument("invalid source unit");
        }
    } catch (const std::exception&) {
        std::cerr<<"invalid numeric DXF tolerance\n";
        return EXIT_FAILURE;
    }
    if (!(options.maximumChordError>0.0) ||
        !std::isfinite(options.maximumChordError) ||
        !(options.endpointWeldTolerance>0.0) ||
        !std::isfinite(options.endpointWeldTolerance)) {
        std::cerr<<"DXF tolerances must be finite and positive\n";
        return EXIT_FAILURE;
    }

    auto imported=readAsciiDxfBoundary2D(input,options);
    std::string outputError;
    if (imported.valid() && !writeBoundaryXy2D(*imported.boundary,output,&outputError,
                                               imported.embeddedPatches)) {
        imported.issues.push_back({DxfIssueCode2D::OutputFailure,0,"","",outputError});
        imported.boundary.reset();
    }
    std::error_code ec;
    if (!report.parent_path().empty())
        std::filesystem::create_directories(report.parent_path(),ec);
    std::ofstream reportOut(report);
    if (!reportOut) {
        std::cerr<<"cannot write DXF report: "<<report.string()<<'\n';
        return EXIT_FAILURE;
    }
    reportOut<<dxfImportReportToJson(imported);
    if (!reportOut.good()) {
        std::cerr<<"failed while writing DXF report\n";
        return EXIT_FAILURE;
    }
    if (!imported.valid()) {
        std::cerr<<"DXF import failed with "<<imported.issues.size()<<" issue(s)\n";
        for (std::size_t i=0;i<imported.issues.size();++i) {
            const auto& issue=imported.issues[i];
            std::cerr<<"dxf_issue["<<i<<"] code="<<static_cast<int>(issue.code)
                     <<" line="<<issue.line<<" entity="<<issue.entity
                     <<" layer="<<issue.layer<<" message="<<issue.message<<'\n';
        }
        std::cerr<<"report="<<report.string()<<'\n';
        return EXIT_FAILURE;
    }
    std::cout<<std::setprecision(17)
             <<"cartmesh2d DXF import PASS\n"
             <<"source_entities="<<imported.report.sourceEntityCount<<'\n'
             <<"output_loops="<<imported.report.outputLoopCount<<'\n'
             <<"output_vertices="<<imported.report.outputVertexCount<<'\n'
             <<"sampled_arc_segments="<<imported.report.sampledArcSegmentCount<<'\n'
             <<"ellipse_entities="<<imported.report.ellipseEntityCount<<'\n'
             <<"spline_entities="<<imported.report.splineEntityCount<<'\n'
             <<"sampled_ellipse_segments="
             <<imported.report.sampledEllipseSegmentCount<<'\n'
             <<"sampled_spline_segments="
             <<imported.report.sampledSplineSegmentCount<<'\n'
             <<"maximum_chord_error="<<imported.report.maximumChordError<<'\n'
             <<"endpoint_weld_tolerance="<<imported.report.endpointWeldTolerance<<'\n'
             <<"output_units=metre\n"
             <<"effective_source_units="<<imported.report.effectiveUnitsName<<'\n'
             <<"coordinate_scale_to_metres="
             <<imported.report.coordinateScaleToMetres<<'\n'
             <<"insertion_units_code=";
    if (imported.report.insertionUnitsCodeDefined)
        std::cout<<imported.report.insertionUnitsCode;
    else std::cout<<"unspecified";
    std::cout<<'\n'
             <<"boundary_xy="<<output.string()<<'\n'
             <<"report="<<report.string()<<'\n';
    return EXIT_SUCCESS;
}
