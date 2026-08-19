#pragma once
#include "cartmesh2d/grid/CartesianGrid2D.hpp"
#include <cstddef>
#include <cstdint>
#include <vector>
namespace cartmesh2d {
struct DistanceRefinementBand2D { double distance = 0.0; std::size_t targetLevel = 0; };
struct QuadtreeRefinementPolicy2D { std::size_t boundaryLevel = 0; std::vector<DistanceRefinementBand2D> distanceBands; };
struct QuadtreeLeaf2D {
    std::size_t id = 0; std::uint64_t key = 0; std::size_t level = 0; std::uint64_t ix = 0; std::uint64_t iy = 0;
    AABB2D bounds; CellClass classification = CellClass::Outside;
    [[nodiscard]] Point2D center() const noexcept;
    [[nodiscard]] double area() const noexcept;
};
struct FaceNeighborPair2D { std::size_t first = 0; std::size_t second = 0; };
struct QuadtreeBalanceReport2D { std::size_t iterations = 0; std::size_t refinedLeaves = 0; std::size_t violationsBefore = 0; std::size_t violationsAfter = 0; };
class Quadtree2D {
public:
    Quadtree2D(Domain2D domain, std::size_t maxLevel, const BoundaryLoop& boundary, const TolerancePolicy& tol = {});
    [[nodiscard]] const Domain2D& domain() const noexcept { return domain_; }
    [[nodiscard]] std::size_t maxLevel() const noexcept { return maxLevel_; }
    [[nodiscard]] const std::vector<QuadtreeLeaf2D>& leaves() const noexcept { return leaves_; }
    void refine(const BoundaryLoop& boundary, const QuadtreeRefinementPolicy2D& policy, const TolerancePolicy& tol = {});
    [[nodiscard]] bool refineLeafByKey(std::uint64_t key, const BoundaryLoop& boundary, const TolerancePolicy& tol = {});
    [[nodiscard]] std::vector<FaceNeighborPair2D> faceNeighbors() const;
    [[nodiscard]] std::size_t countBalanceViolations() const;
    [[nodiscard]] QuadtreeBalanceReport2D enforceTwoToOneBalance(const BoundaryLoop& boundary, const TolerancePolicy& tol = {});
    [[nodiscard]] double totalLeafArea() const noexcept;
    [[nodiscard]] bool deterministicOrderingValid() const noexcept;
private:
    Domain2D domain_; std::size_t maxLevel_ = 0; std::vector<QuadtreeLeaf2D> leaves_;
    [[nodiscard]] static std::uint64_t mortonPath(std::size_t level, std::uint64_t ix, std::uint64_t iy) noexcept;
    [[nodiscard]] static std::uint64_t makeKey(std::size_t level, std::uint64_t ix, std::uint64_t iy) noexcept;
    [[nodiscard]] AABB2D boundsFor(std::size_t level, std::uint64_t ix, std::uint64_t iy) const noexcept;
    [[nodiscard]] QuadtreeLeaf2D makeLeaf(std::size_t level, std::uint64_t ix, std::uint64_t iy, const BoundaryLoop& boundary, const TolerancePolicy& tol) const;
    [[nodiscard]] bool splitLeafAt(std::size_t index, const BoundaryLoop& boundary, const TolerancePolicy& tol);
    void sortAndAssignIds();
};
[[nodiscard]] double distancePointToSegment(const Point2D& point, const Segment2D& segment) noexcept;
[[nodiscard]] double distanceAABBToBoundary(const AABB2D& box, const BoundaryLoop& boundary, const TolerancePolicy& tol = {});
} // namespace cartmesh2d
