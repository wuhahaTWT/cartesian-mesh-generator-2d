#include "cartmesh2d/io/Dxf2D.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace cartmesh2d;

namespace {

void usage() {
    std::cerr<<"usage: cartmesh2d_dxf_cli <input-ascii.dxf> <output.xy> "
               "<maximum-chord-error> [report.json] [endpoint-weld-tolerance=1e-10]\n";
}

} // namespace

int main(int argc,char** argv) {
    if (argc<4 || argc>6) {
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
    if (imported.valid() && !writeBoundaryXy2D(*imported.boundary,output,&outputError)) {
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
             <<"maximum_chord_error="<<imported.report.maximumChordError<<'\n'
             <<"endpoint_weld_tolerance="<<imported.report.endpointWeldTolerance<<'\n'
             <<"insertion_units_code=";
    if (imported.report.insertionUnitsCodeDefined)
        std::cout<<imported.report.insertionUnitsCode;
    else std::cout<<"unspecified";
    std::cout<<'\n'
             <<"boundary_xy="<<output.string()<<'\n'
             <<"report="<<report.string()<<'\n';
    return EXIT_SUCCESS;
}
