// Read-only shape screening using the existing native quality implementation.
// Candidate limits are research targets, never installed as product defaults.
#include "cartmesh2d/io/MeshIO2D.hpp"
#include "cartmesh2d/quality/SolverQuality2D.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
using namespace cartmesh2d;
const char* name(SolverQualityIssueCode2D code) {
    switch(code) {
    case SolverQualityIssueCode2D::ExcessiveNonOrthogonality:return "nonOrthogonality";
    case SolverQualityIssueCode2D::ExcessiveAspect:return "hydraulicAspect";
    case SolverQualityIssueCode2D::SmallInteriorAngle:return "interiorAngle";
    default:return "otherExistingChecks";
    }
}
int main(int argc,char** argv) {
    try {
        if(argc!=3)throw std::runtime_error("MESH OUTPUT_JSON");
        if(std::filesystem::exists(argv[2]))throw std::runtime_error("output exists");
        const auto input=readCm2dTopology(argv[1]);
        if(!input.valid())throw std::runtime_error(input.error);
        const auto& mesh=input.topology;
        std::ofstream out(argv[2]);if(!out)throw std::runtime_error("cannot create output");
        out<<std::setprecision(17)<<"{\"cells\":"<<mesh.cells.size()<<",\"currentNativeReport\":"
           <<solverQualityReportToJson(evaluateSolverQuality2D(mesh))<<",\"screens\":[";
        for(int target=0;target<2;++target) {
            SolverQualityPolicy2D p;
            p.maxNonOrthogonalityDeg=target?45.:60.;
            p.maxCellAspect=target?10.:20.;
            p.minInteriorAngleDeg=target?20.:10.;
            const auto q=evaluateSolverQuality2D(mesh,p);
            std::map<std::string,std::size_t> counts;
            for(const auto& issue:q.issues)++counts[name(issue.code)];
            if(target)out<<',';
            out<<"{\"name\":\""<<(target?"generation-target-candidate":"shape-rejection-candidate")
               <<"\",\"installedAsDefault\":false,\"maxNonOrthogonalityDeg\":"<<p.maxNonOrthogonalityDeg
               <<",\"maxHydraulicAspect\":"<<p.maxCellAspect<<",\"minInteriorAngleDeg\":"<<p.minInteriorAngleDeg
               <<",\"passes\":"<<(q.valid()?"true":"false")<<",\"issueEvents\":"<<q.issues.size()
               <<",\"issueCounts\":{";
            bool first=true;for(const auto& [key,value]:counts){if(!first)out<<',';first=false;out<<'"'<<key<<"\":"<<value;}
            out<<"}}";
        }
        const Polygon2D square{{{0,0},{1,0},{1,1},{0,1}}};
        const Polygon2D rectangle{{{0,0},{10,0},{10,1},{0,1}}};
        out<<"],\"nativeMetricCalibration\":{\"squareHydraulicAspect\":"<<evaluateSolverCellMetrics2D(square).hydraulicAspect
           <<",\"rectangle10By1HydraulicAspect\":"<<evaluateSolverCellMetrics2D(rectangle).hydraulicAspect
           <<"},\"scope\":\"Read-only shape screening; existing native kernels; no CFD run or convergence qualification. Issue events are not unique cell counts. Other policy values unchanged.\"}\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
