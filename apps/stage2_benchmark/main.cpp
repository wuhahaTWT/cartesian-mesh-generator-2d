#include "cartmesh/grid/LinearOctree.hpp"
#include "cartmesh/util/Performance.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

using Clock = std::chrono::steady_clock;

[[nodiscard]] std::uint64_t parse_target(std::string_view text) {
    std::size_t parsed = 0;
    const auto value = std::stoull(std::string(text), &parsed);
    if (parsed != text.size() || value < 8 || value > 16'777'216ULL) {
        throw std::invalid_argument("目标叶数必须位于 8..16777216");
    }
    return value;
}

[[nodiscard]] double elapsed(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
}

[[nodiscard]] constexpr std::string_view build_type() noexcept {
#if defined(NDEBUG)
    return "Release";
#else
    return "Debug";
#endif
}

} // 匿名命名空间

int main(int argc, char** argv) {
    try {
        const std::uint64_t requested = argc > 1 ? parse_target(argv[1]) : 10'000'000ULL;
        std::uint8_t base_level = 0;
        while (base_level < 8U &&
               (std::uint64_t{1} << (3U * static_cast<std::uint32_t>(base_level + 1U))) <=
                   requested) {
            ++base_level;
        }
        if (base_level == 8U) {
            base_level = 7U;
        }
        const auto initial_leaf_count =
            std::uint64_t{1} << (3U * static_cast<std::uint32_t>(base_level));
        const auto split_count =
            requested > initial_leaf_count ? (requested - initial_leaf_count + 6U) / 7U : 0U;
        if (split_count > initial_leaf_count) {
            throw std::invalid_argument("请求叶数无法用相邻两层基准生成");
        }
        const auto total_start = Clock::now();
        const auto construction_start = Clock::now();
        cartmesh::LinearOctree tree(cartmesh::AABB({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}),
                                    base_level,
                                    static_cast<std::uint8_t>(base_level + 1U));
        const auto construction_end = Clock::now();
        const auto refinement_start = Clock::now();
        const std::uint64_t anchor_threshold = split_count * 8U;
        const auto refinement = tree.refine_to_desired_levels(
            [&](cartmesh::OctreeNodeCode code, const cartmesh::AABB&) {
                return cartmesh::octree_anchor_morton(code, tree.maximum_level()) <
                               anchor_threshold
                           ? tree.maximum_level()
                           : tree.base_level();
            });
        const auto refinement_end = Clock::now();
        const auto balance_start = Clock::now();
        const auto balance = tree.balance_faces_2_to_1();
        const auto balance_end = Clock::now();
        const auto report = tree.check_face_balance();
        if (!tree.validate_partition() || !report.balanced ||
            tree.leaf_count() < requested || tree.leaf_count() > requested + 6U) {
            throw std::runtime_error("千万叶八叉树基准不变量失败");
        }
        std::cout << std::setprecision(17)
                  << "{\n  \"schemaVersion\": 1,\n  \"projectStage\": 2,\n"
                  << "  \"benchmark\": \"compact_linear_octree_construction\",\n"
                  << "  \"solverReadyCutCellMesh\": false,\n"
                  << "  \"requestedLeafCount\": " << requested
                  << ",\n  \"initialLeafCount\": " << initial_leaf_count
                  << ",\n  \"finalLeafCount\": " << tree.leaf_count()
                  << ",\n  \"compactLeafStorageBytes\": " << tree.compact_storage_bytes()
                  << ",\n  \"bytesPerLeafCode\": 8,\n"
                  << "  \"ruleSplitCount\": " << refinement.split_count
                  << ",\n  \"balanceSplitCount\": " << balance.split_count
                  << ",\n  \"partitionValid\": true,\n  \"faceBalanced2To1\": true,\n"
                  << "  \"resultHashFnv1a64Decimal\": " << tree.result_hash_fnv1a64()
                  << ",\n  \"threads\": 1,\n  \"buildType\": \"" << build_type()
                  << "\",\n  \"peakRssBytes\": " << cartmesh::peak_rss_bytes()
                  << ",\n  \"timingsSeconds\": {\"initialConstruction\": "
                  << elapsed(construction_start, construction_end) << ", \"refinement\": "
                  << elapsed(refinement_start, refinement_end) << ", \"balanceCheck\": "
                  << elapsed(balance_start, balance_end) << ", \"total\": "
                  << elapsed(total_start, Clock::now()) << "}\n}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "cartmesh_stage2_benchmark 错误：" << error.what() << '\n';
        return 1;
    }
}
