#include "cartmesh/geometry/SurfaceMesh.hpp"
#include "cartmesh/cutcell/TriangulatedSurfaceCutter.hpp"
#include "cartmesh/cutcell/CutCellMeshFingerprint.hpp"
#include "cartmesh/grid/LinearOctree.hpp"
#include "cartmesh/grid/OctreeRefinement.hpp"
#include "cartmesh/incremental/IncrementalRemesher.hpp"
#include "cartmesh/incremental/MeshMapping.hpp"
#include "cartmesh/incremental/SurfaceChangeSet.hpp"
#include "cartmesh/io/StlReader.hpp"
#include "cartmesh/spatial/TriangleBvh.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct TestFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

void expect(bool condition, const std::string& message) {
    if (!condition) throw TestFailure(message);
}

cartmesh::SurfaceMesh make_box_surface(const cartmesh::Vec3& minimum,
                                       const cartmesh::Vec3& maximum,
                                       bool reorder = false) {
    using cartmesh::Triangle;
    const cartmesh::Vec3 p000{minimum.x, minimum.y, minimum.z};
    const cartmesh::Vec3 p100{maximum.x, minimum.y, minimum.z};
    const cartmesh::Vec3 p010{minimum.x, maximum.y, minimum.z};
    const cartmesh::Vec3 p110{maximum.x, maximum.y, minimum.z};
    const cartmesh::Vec3 p001{minimum.x, minimum.y, maximum.z};
    const cartmesh::Vec3 p101{maximum.x, minimum.y, maximum.z};
    const cartmesh::Vec3 p011{minimum.x, maximum.y, maximum.z};
    const cartmesh::Vec3 p111{maximum.x, maximum.y, maximum.z};
    std::vector<Triangle> triangles{
        Triangle(p000, p010, p110), Triangle(p000, p110, p100),
        Triangle(p001, p101, p111), Triangle(p001, p111, p011),
        Triangle(p000, p100, p101), Triangle(p000, p101, p001),
        Triangle(p010, p011, p111), Triangle(p010, p111, p110),
        Triangle(p000, p001, p011), Triangle(p000, p011, p010),
        Triangle(p100, p110, p111), Triangle(p100, p111, p101)};
    if (reorder) {
        std::reverse(triangles.begin(), triangles.end());
        for (auto& triangle : triangles) {
            const auto vertices = triangle.vertices();
            triangle = Triangle(vertices[1], vertices[0], vertices[2]);
        }
    }
    return cartmesh::SurfaceMesh(std::move(triangles),
                                 cartmesh::SurfaceFormat::ascii_stl,
                                 "box");
}

void test_surface_change_is_order_independent() {
    const auto first = make_box_surface({-0.4, -0.4, -0.4},
                                        {0.4, 0.4, 0.4});
    const auto reordered = make_box_surface({-0.4, -0.4, -0.4},
                                            {0.4, 0.4, 0.4}, true);
    const auto unchanged = cartmesh::detect_surface_changes(first, reordered);
    expect(unchanged.empty() && !unchanged.bounds,
           "三角片顺序、顶点起点和绕序变化不得伪造几何变化");

    const auto moved = make_box_surface({-0.4, -0.4, -0.4},
                                        {0.55, 0.4, 0.4});
    const auto changed = cartmesh::detect_surface_changes(first, moved);
    expect(!changed.empty() && changed.bounds,
           "局部轮廓坐标变化必须形成显式变化集");
    expect(changed.bounds->minimum().x <= 0.4 &&
               changed.bounds->maximum().x >= 0.55,
           "变化包围盒必须同时覆盖旧、新表面位置");
}

void test_incremental_octree_matches_full_rebuild() {
    const cartmesh::AABB domain({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
    const auto old_surface = make_box_surface({-0.45, -0.35, -0.35},
                                              {0.35, 0.35, 0.35});
    const auto new_surface = make_box_surface({-0.45, -0.35, -0.35},
                                              {0.55, 0.35, 0.35});
    cartmesh::OctreeRefinementConfiguration configuration;
    configuration.surface_target_level = 5;
    configuration.distance_bands.push_back({0.08, 4});

    const cartmesh::TriangleBvh old_bvh(old_surface);
    cartmesh::LinearOctree old_tree(domain, 2, 5);
    static_cast<void>(cartmesh::OctreeRefinementEngine(configuration, &old_bvh)
                          .apply(old_tree));

    const cartmesh::TriangleBvh new_bvh(new_surface);
    const cartmesh::OctreeRefinementEngine new_engine(configuration, &new_bvh);
    cartmesh::LinearOctree full_tree(domain, 2, 5);
    static_cast<void>(new_engine.apply(full_tree));

    const auto changes = cartmesh::detect_surface_changes(old_surface, new_surface);
    const auto affected = cartmesh::conservative_affected_bounds(
        changes, domain, 0.12);
    const auto incremental = cartmesh::update_octree_incrementally(
        old_tree, new_engine, affected);

    expect(incremental.tree.validate_partition(),
           "增量更新后必须保持完整无重叠分区");
    expect(incremental.tree.check_face_balance().balanced,
           "增量更新后必须保持面 2:1 平衡");
    expect(incremental.tree.leaf_codes().size() == full_tree.leaf_codes().size() &&
               std::equal(incremental.tree.leaf_codes().begin(),
                          incremental.tree.leaf_codes().end(),
                          full_tree.leaf_codes().begin()),
           "增量八叉树必须与新几何全量重构逐节点码一致");
    expect(incremental.statistics.preserved_leaf_count > 0 &&
               incremental.statistics.created_leaf_count > 0 &&
               incremental.statistics.removed_leaf_count > 0,
           "局部轮廓变化必须同时保留未影响叶并记录创建/删除叶");
}

void expect_mesh_equivalent(const cartmesh::ConvexCutCellMesh& incremental,
                            const cartmesh::ConvexCutCellMesh& full) {
    const double tolerance = 1.0e-12;
    expect(incremental.fluid_cells.size() == full.fluid_cells.size(),
           "增量与全量流体单元数必须一致");
    expect(incremental.internal_faces.size() == full.internal_faces.size(),
           "增量与全量内部面数必须一致");
    expect(incremental.component_internal_faces.size() ==
               full.component_internal_faces.size(),
           "增量与全量分量内部面数必须一致");
    expect(incremental.full_fluid_cell_count == full.full_fluid_cell_count &&
               incremental.full_solid_cell_count == full.full_solid_cell_count &&
               incremental.cut_cell_count == full.cut_cell_count &&
               incremental.nonclosed_cell_count == full.nonclosed_cell_count &&
               incremental.negative_volume_cell_count ==
                   full.negative_volume_cell_count &&
               incremental.shared_face_mismatch_count ==
                   full.shared_face_mismatch_count &&
               incremental.classification_conflict_count ==
                   full.classification_conflict_count,
           "增量与全量 Cut-cell 分类和拓扑计数必须一致");
    expect(std::abs(incremental.total_fluid_volume - full.total_fluid_volume) <
               tolerance &&
               std::abs(incremental.total_embedded_boundary_area -
                        full.total_embedded_boundary_area) < tolerance,
           "增量与全量总体积和嵌入边界面积必须一致");
    for (std::size_t index = 0; index < full.fluid_cells.size(); ++index) {
        const auto& first = incremental.fluid_cells[index];
        const auto& second = full.fluid_cells[index];
        expect(first.background_cell_id == second.background_cell_id &&
                   first.cut == second.cut &&
                   first.fluid_component_count == second.fluid_component_count &&
                   first.fluid_piece_count == second.fluid_piece_count &&
                   std::abs(first.volume - second.volume) < tolerance &&
                   std::abs(first.volume_fraction - second.volume_fraction) <
                       tolerance &&
                   cartmesh::norm(first.centroid - second.centroid) < tolerance,
               "增量与全量逐流体单元几何必须一致");
        expect(first.embedded_boundary_faces.size() ==
                   second.embedded_boundary_faces.size() &&
                   first.fluid_component_region_ids ==
                       second.fluid_component_region_ids,
               "增量与全量边界面和全局 region 必须一致");
        for (std::size_t face = 0; face < 6; ++face) {
            expect(std::abs(first.cartesian_faces[face].area -
                            second.cartesian_faces[face].area) < tolerance &&
                       cartmesh::norm(first.cartesian_faces[face].centroid -
                                      second.cartesian_faces[face].centroid) <
                           tolerance,
                   "增量与全量 Cartesian 面开口必须一致");
        }
    }
    for (std::size_t index = 0; index < full.internal_faces.size(); ++index) {
        const auto& first = incremental.internal_faces[index];
        const auto& second = full.internal_faces[index];
        expect(first.first_background_cell_id ==
                       second.first_background_cell_id &&
                   first.second_background_cell_id ==
                       second.second_background_cell_id &&
                   std::abs(first.area - second.area) < tolerance &&
                   cartmesh::norm(first.centroid - second.centroid) < tolerance,
               "增量与全量内部面邻接必须一致");
    }
}

void test_incremental_cutcell_reuses_unchanged_geometry() {
    const cartmesh::AABB domain({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
    const auto old_surface = make_box_surface({-0.4, -0.3, -0.3},
                                              {0.25, 0.3, 0.3});
    const auto new_surface = make_box_surface({-0.4, -0.3, -0.3},
                                              {0.48, 0.3, 0.3});
    cartmesh::OctreeRefinementConfiguration configuration;
    configuration.surface_target_level = 3;

    const cartmesh::TriangulatedSurfaceCutter old_cutter(old_surface, 0);
    cartmesh::LinearOctree old_tree(domain, 1, 3);
    static_cast<void>(cartmesh::OctreeRefinementEngine(
                          configuration, &old_cutter.bvh())
                          .apply(old_tree));
    const auto old_mesh = cartmesh::build_triangulated_cut_cell_mesh(
        old_tree, old_cutter);

    const cartmesh::TriangulatedSurfaceCutter new_cutter(new_surface, 0);
    const cartmesh::OctreeRefinementEngine new_engine(configuration,
                                                       &new_cutter.bvh());
    cartmesh::LinearOctree full_tree(domain, 1, 3);
    static_cast<void>(new_engine.apply(full_tree));
    const auto full_mesh = cartmesh::build_triangulated_cut_cell_mesh(
        full_tree, new_cutter);

    const auto changes = cartmesh::detect_surface_changes(old_surface, new_surface);
    const auto affected = cartmesh::conservative_affected_bounds(
        changes, domain, 0.02);
    auto updated_tree = cartmesh::update_octree_incrementally(
        old_tree, new_engine, affected);
    expect(std::equal(updated_tree.tree.leaf_codes().begin(),
                      updated_tree.tree.leaf_codes().end(),
                      full_tree.leaf_codes().begin(),
                      full_tree.leaf_codes().end()),
           "Cut-cell 增量案例的叶码必须与全量树一致");

    const auto incremental =
        cartmesh::build_incremental_triangulated_cut_cell_mesh(
            old_tree, old_mesh, updated_tree.tree, new_cutter, affected);
    expect(incremental.statistics.reused_leaf_count > 0 &&
               incremental.statistics.rebuilt_leaf_count > 0 &&
               incremental.statistics.geometry_reuse_fraction > 0.0 &&
               incremental.statistics.geometry_reuse_fraction < 1.0,
           "局部变化必须真实复用部分旧单元并重构受影响单元");
    expect_mesh_equivalent(incremental.mesh, full_mesh);
    expect(cartmesh::cut_cell_mesh_fingerprint_fnv1a64(
               updated_tree.tree, incremental.mesh) ==
               cartmesh::cut_cell_mesh_fingerprint_fnv1a64(full_tree,
                                                            full_mesh),
           "增量与全量 Cut-cell 结果指纹必须逐字节一致");
    expect(incremental.mesh.nonclosed_cell_count == 0 &&
               incremental.mesh.negative_volume_cell_count == 0 &&
               incremental.mesh.shared_face_mismatch_count == 0 &&
               incremental.mesh.classification_conflict_count == 0,
           "增量 Cut-cell 结果必须保持阶段三/四硬不变量");

    const auto mapping = cartmesh::build_incremental_mesh_mapping(
        old_tree, old_mesh, updated_tree.tree, incremental.mesh, affected);
    expect(!mapping.entries.empty() && mapping.preserved_pair_count > 0 &&
               mapping.rebuilt_pair_count > 0,
           "增量映射必须同时记录保留和受影响重构单元");
    expect(std::abs(mapping.total_background_overlap_volume - domain.volume()) <
               1.0e-12,
           "旧、新八叉树映射必须无缝覆盖固定计算域");
    expect(mapping.total_shared_fluid_overlap_volume >= 0.0 &&
               mapping.total_shared_fluid_overlap_volume <=
                   std::min(mapping.old_fluid_volume,
                            mapping.new_fluid_volume) + 1.0e-12,
           "真实流体重叠体积必须位于旧、新流体总体积范围内");
    expect(mapping.removed_fluid_volume > 0.0,
           "固体右侧外移必须显式报告被新固体占据的旧流体体积");
    for (const auto& entry : mapping.entries) {
        expect(entry.exact_fluid_overlap &&
                   entry.old_volume_preserved_fraction >= 0.0 &&
                   entry.old_volume_preserved_fraction <= 1.0 + 1.0e-10 &&
                   entry.new_volume_from_old_fraction >= 0.0 &&
                   entry.new_volume_from_old_fraction <= 1.0 + 1.0e-10,
               "每条增量映射必须提供有界的真实几何重叠权重");
    }
}

void test_coplanar_tunnel_corner_regression() {
    const std::string path =
        std::string(CARTMESH_SOURCE_DIR) +
        "/tests/data/stage5_coplanar_tunnel_corner_ascii.stl";
    const auto surface = cartmesh::read_stl(path);
    const cartmesh::TriangulatedSurfaceCutter cutter(surface, 0);
    const cartmesh::AABB domain({-0.2, -0.2, -0.2}, {1.2, 1.2, 1.2});
    cartmesh::LinearOctree tree(domain, 2, 4);
    cartmesh::OctreeRefinementConfiguration configuration;
    configuration.surface_target_level = 4;
    static_cast<void>(cartmesh::OctreeRefinementEngine(
                          configuration, &cutter.bvh())
                          .apply(tree));
    const auto mesh = cartmesh::build_triangulated_cut_cell_mesh(tree, cutter);
    expect(mesh.cut_cell_count > 0 && mesh.nonclosed_cell_count == 0 &&
               mesh.shared_face_mismatch_count == 0 &&
               mesh.negative_volume_cell_count == 0,
           "共面通道壁与开口凹角必须保持闭合且共享面一致");
}

} // namespace

int main() {
    try {
        test_surface_change_is_order_independent();
        test_incremental_octree_matches_full_rebuild();
        test_incremental_cutcell_reuses_unchanged_geometry();
        test_coplanar_tunnel_corner_regression();
        std::cout << "阶段5增量测试通过\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "阶段5增量测试失败：" << error.what() << '\n';
        return 1;
    }
}
