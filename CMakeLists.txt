cmake_minimum_required(VERSION 3.20)

project(cartesian_mesher VERSION 0.1.0 LANGUAGES CXX)

option(CARTMESH_BUILD_TESTS "构建项目测试套件" ON)
option(CARTMESH_BUILD_BENCHMARKS "构建基准程序" ON)

function(cartmesh_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    endif()
endfunction()

add_library(cartmesh
    src/classify/SurfaceClassifier.cpp
    src/cutcell/ConvexCutCellMesh.cpp
    src/cutcell/CutCellMeshFingerprint.cpp
    src/cutcell/CutCellTopology.cpp
    src/cutcell/CutCellRegions.cpp
    src/cutcell/ConvexPolyhedron.cpp
    src/cutcell/ConvexSurfaceCutter.cpp
    src/cutcell/TriangulatedSurfaceCutter.cpp
    src/cutcell/TriangulatedSurfaceCutCellMesh.cpp
    src/geometry/SurfaceDiagnostics.cpp
    src/geometry/SurfaceMesh.cpp
    src/geometry/TriangleBoxIntersection.cpp
    src/geometry/TriangleTriangleIntersection.cpp
    src/grid/LinearOctree.cpp
    src/grid/OctreeRefinement.cpp
    src/grid/UniformCartesianGrid.cpp
    src/incremental/IncrementalRemesher.cpp
    src/incremental/MeshMapping.cpp
    src/incremental/SurfaceChangeSet.cpp
    src/io/CutCellJsonWriter.cpp
    src/io/CutCellVtkWriter.cpp
    src/io/StlReader.cpp
    src/io/DiagnosticVtkWriter.cpp
    src/io/OpenFoamWriter.cpp
    src/io/ScalableOpenFoamWriter.cpp
    src/io/VtkWriter.cpp
    src/scalable/CompactUniformCutCellMesh.cpp
    src/spatial/TriangleBvh.cpp
)
add_library(cartmesh::cartmesh ALIAS cartmesh)

target_compile_features(cartmesh PUBLIC cxx_std_20)
target_include_directories(cartmesh
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

cartmesh_enable_warnings(cartmesh)

add_executable(cartmesh_cli apps/cli/main.cpp)
target_link_libraries(cartmesh_cli PRIVATE cartmesh::cartmesh)
target_compile_features(cartmesh_cli PRIVATE cxx_std_20)
cartmesh_enable_warnings(cartmesh_cli)

add_executable(cartmesh_octree_cli apps/octree_cli/main.cpp)
target_link_libraries(cartmesh_octree_cli PRIVATE cartmesh::cartmesh)
target_compile_features(cartmesh_octree_cli PRIVATE cxx_std_20)
cartmesh_enable_warnings(cartmesh_octree_cli)

add_executable(cartmesh_cutcell_cli apps/cutcell_cli/main.cpp)
target_link_libraries(cartmesh_cutcell_cli PRIVATE cartmesh::cartmesh)
target_compile_features(cartmesh_cutcell_cli PRIVATE cxx_std_20)
cartmesh_enable_warnings(cartmesh_cutcell_cli)

add_executable(cartmesh_incremental_cli apps/incremental_cli/main.cpp)
target_link_libraries(cartmesh_incremental_cli PRIVATE cartmesh::cartmesh)
target_compile_features(cartmesh_incremental_cli PRIVATE cxx_std_20)
cartmesh_enable_warnings(cartmesh_incremental_cli)

add_executable(cartmesh_stage6_cli apps/stage6_cli/main.cpp)
target_link_libraries(cartmesh_stage6_cli PRIVATE cartmesh::cartmesh)
target_compile_features(cartmesh_stage6_cli PRIVATE cxx_std_20)
cartmesh_enable_warnings(cartmesh_stage6_cli)

if(CARTMESH_BUILD_BENCHMARKS)
    add_executable(cartmesh_benchmark apps/benchmark/main.cpp)
    target_link_libraries(cartmesh_benchmark PRIVATE cartmesh::cartmesh)
    target_compile_features(cartmesh_benchmark PRIVATE cxx_std_20)
    cartmesh_enable_warnings(cartmesh_benchmark)
    add_executable(cartmesh_stage2_benchmark apps/stage2_benchmark/main.cpp)
    target_link_libraries(cartmesh_stage2_benchmark PRIVATE cartmesh::cartmesh)
    target_compile_features(cartmesh_stage2_benchmark PRIVATE cxx_std_20)
    cartmesh_enable_warnings(cartmesh_stage2_benchmark)
endif()

include(CTest)
if(BUILD_TESTING AND CARTMESH_BUILD_TESTS)
    add_executable(cartmesh_tests tests/test_main.cpp)
    target_link_libraries(cartmesh_tests PRIVATE cartmesh::cartmesh)
    target_compile_features(cartmesh_tests PRIVATE cxx_std_20)
    cartmesh_enable_warnings(cartmesh_tests)
    add_test(NAME cartmesh_stage0_tests COMMAND cartmesh_tests)
    add_executable(cartmesh_stage1_tests tests/stage1_test.cpp)
    target_link_libraries(cartmesh_stage1_tests PRIVATE cartmesh::cartmesh)
    target_compile_features(cartmesh_stage1_tests PRIVATE cxx_std_20)
    cartmesh_enable_warnings(cartmesh_stage1_tests)
    add_test(NAME cartmesh_stage1_tests COMMAND cartmesh_stage1_tests)
    add_executable(cartmesh_stage2_tests tests/stage2_test.cpp)
    target_link_libraries(cartmesh_stage2_tests PRIVATE cartmesh::cartmesh)
    target_compile_features(cartmesh_stage2_tests PRIVATE cxx_std_20)
    cartmesh_enable_warnings(cartmesh_stage2_tests)
    add_test(NAME cartmesh_stage2_core_tests COMMAND cartmesh_stage2_tests)
    add_executable(cartmesh_stage3_tests tests/stage3_test.cpp)
    target_link_libraries(cartmesh_stage3_tests PRIVATE cartmesh::cartmesh)
    target_compile_features(cartmesh_stage3_tests PRIVATE cxx_std_20)
    cartmesh_enable_warnings(cartmesh_stage3_tests)
    add_test(NAME cartmesh_stage3_geometry_tests COMMAND cartmesh_stage3_tests)
    add_executable(cartmesh_stage4_tests tests/stage4_test.cpp)
    target_link_libraries(cartmesh_stage4_tests PRIVATE cartmesh::cartmesh)
    target_compile_features(cartmesh_stage4_tests PRIVATE cxx_std_20)
    cartmesh_enable_warnings(cartmesh_stage4_tests)
    add_test(NAME cartmesh_stage4_industrial_geometry_tests COMMAND cartmesh_stage4_tests)
    add_executable(cartmesh_stage5_tests tests/stage5_test.cpp)
    target_link_libraries(cartmesh_stage5_tests PRIVATE cartmesh::cartmesh)
    target_compile_features(cartmesh_stage5_tests PRIVATE cxx_std_20)
    target_compile_definitions(cartmesh_stage5_tests PRIVATE
        CARTMESH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    cartmesh_enable_warnings(cartmesh_stage5_tests)
    add_test(NAME cartmesh_stage5_incremental_tests COMMAND cartmesh_stage5_tests)

    add_executable(cartmesh_stage6_tests tests/stage6_test.cpp)
    target_link_libraries(cartmesh_stage6_tests PRIVATE cartmesh::cartmesh)
    target_compile_features(cartmesh_stage6_tests PRIVATE cxx_std_20)
    cartmesh_enable_warnings(cartmesh_stage6_tests)
    add_test(NAME cartmesh_stage6_scalable_tests COMMAND cartmesh_stage6_tests)
    add_test(
        NAME cartmesh_cli_stl_smoke
        COMMAND cartmesh_cli --stl ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/closed_unit_cube_ascii.stl
                --resolution 12 --no-vtk
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage1_cube_smoke.json
    )
    add_test(
        NAME cartmesh_cli_rejects_open_stl
        COMMAND cartmesh_cli
                --stl ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/open_unit_cube_missing_triangle_ascii.stl
                --resolution 8 --no-vtk
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage1_open_rejected.json
    )
    set_tests_properties(cartmesh_cli_rejects_open_stl PROPERTIES WILL_FAIL TRUE)
    add_test(
        NAME cartmesh_cli_rejects_nonmanifold_vertex_stl
        COMMAND cartmesh_cli
                --stl ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nonmanifold_vertex_two_tetrahedra_ascii.stl
                --resolution 8 --no-vtk
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage1_nonmanifold_vertex_rejected.json
    )
    set_tests_properties(cartmesh_cli_rejects_nonmanifold_vertex_stl PROPERTIES WILL_FAIL TRUE)
    add_test(
        NAME cartmesh_cli_cube_smoke
        COMMAND cartmesh_cli --case cube --resolution 4 --no-vtk
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/cube_smoke.json
    )
    add_test(NAME cartmesh_cli_rejects_zero COMMAND cartmesh_cli --resolution 0 --no-vtk)
    set_tests_properties(cartmesh_cli_rejects_zero PROPERTIES WILL_FAIL TRUE)
    add_test(
        NAME cartmesh_octree_cli_smoke
        COMMAND cartmesh_octree_cli
                --stl ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/closed_unit_cube_ascii.stl
                --base-level 2 --max-level 4 --surface-level 4 --distance 0.1:3
                --no-vtk
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage2_cube_smoke.json
    )
    add_test(
        NAME cartmesh_cutcell_cli_convex_cube_smoke
        COMMAND cartmesh_cutcell_cli
                --stl ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/closed_unit_cube_ascii.stl
                --resolution 8 --no-vtk
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage3_cube_smoke.json
    )
    add_test(
        NAME cartmesh_cutcell_cli_nonconvex_l_prism_smoke
        COMMAND cartmesh_cutcell_cli
                --stl ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nonconvex_l_prism_ascii.stl
                --resolution 4 --padding-fraction 0.1 --no-vtk
                --geometry-output ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage3_l_prism_geometry.json
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage3_l_prism_smoke.json
    )
    add_test(
        NAME cartmesh_cutcell_cli_adaptive_cube_smoke
        COMMAND cartmesh_cutcell_cli
                --stl ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/closed_unit_cube_ascii.stl
                --adaptive --base-level 2 --max-level 4 --surface-level 4
                --distance 0.1:3 --no-vtk
                --geometry-output ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage4_adaptive_cube_geometry.json
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage4_adaptive_cube.json
    )
    add_test(
        NAME cartmesh_cutcell_cli_openfoam_cube_smoke
        COMMAND cartmesh_cutcell_cli
                --stl ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/closed_unit_cube_ascii.stl
                --resolution 4 --padding-fraction 0.1 --no-vtk
                --geometry-output ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage3_openfoam_cube_geometry.json
                --openfoam-case ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage3_openfoam_cube_case
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage3_openfoam_cube.json
    )
    add_test(
        NAME cartmesh_cutcell_cli_rejects_invalid_surface
        COMMAND cartmesh_cutcell_cli
                --stl ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/analytic/stage4/partial_overlap_ascii.stl
                --resolution 8 --no-vtk
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage4_overlap_rejected.json
    )
    set_tests_properties(cartmesh_cutcell_cli_rejects_invalid_surface
                         PROPERTIES WILL_FAIL TRUE)
    add_test(
        NAME cartmesh_octree_cli_rejects_open_stl
        COMMAND cartmesh_octree_cli
                --stl ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/open_unit_cube_missing_triangle_ascii.stl
                --base-level 2 --max-level 4 --no-vtk
                --report ${CMAKE_CURRENT_BINARY_DIR}/test-artifacts/stage2_open_rejected.json
    )
    set_tests_properties(cartmesh_octree_cli_rejects_open_stl PROPERTIES WILL_FAIL TRUE)
    if(CARTMESH_BUILD_BENCHMARKS)
        add_test(NAME cartmesh_benchmark_smoke COMMAND cartmesh_benchmark 8)
        add_test(NAME cartmesh_stage2_benchmark_smoke COMMAND cartmesh_stage2_benchmark 10000)
    endif()
endif()

install(TARGETS cartmesh cartmesh_cli cartmesh_octree_cli cartmesh_cutcell_cli
                cartmesh_incremental_cli)
install(DIRECTORY include/ DESTINATION include)
