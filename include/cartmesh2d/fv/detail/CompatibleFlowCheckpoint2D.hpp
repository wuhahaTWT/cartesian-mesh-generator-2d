#pragma once
#include "cartmesh2d/fv/CompatibleIncompressible2D.hpp"
#include <cstdint>
#include <stdexcept>
namespace cartmesh2d::fv::detail {
struct CompatibleCheckpointData2D {
    std::shared_ptr<const std::vector<std::uint64_t>> context;
    std::vector<double> normalizedState;
    std::size_t cells=0,faces=0,acceptedIterations=0;
    double nextPseudoStep=0;
    CompatibleFlowMetrics2D metrics;
};
struct CompatibleCheckpointAccess2D {
    static const CompatibleCheckpointData2D& get(const CompatibleFlowCheckpoint2D& checkpoint) {
        if(!checkpoint.data_)throw std::invalid_argument("Empty compatible checkpoint");
        return *checkpoint.data_;
    }
    static CompatibleFlowCheckpoint2D make(CompatibleCheckpointData2D data) {
        return CompatibleFlowCheckpoint2D(std::make_shared<const CompatibleCheckpointData2D>(std::move(data)));
    }
};
// Exact storage extent from mesh identity, boundary kinds, projected traces
// and the already assembled native load. No PDE is reimplemented here.
inline std::size_t compatibleCheckpointContextWords2D(const FvMesh2D& mesh) {
    std::size_t incidence=0;
    for(const auto& cell:mesh.cells)incidence+=cell.faces.size();
    return 9+14*mesh.cells.size()+17*mesh.faces.size()+5*incidence;
}
}
