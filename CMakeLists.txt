cmake_minimum_required(VERSION 3.20)

project(cartmesh2d LANGUAGES CXX)

add_library(cartmesh2d
    src/boundary_layer/BoundaryLayer2D.cpp
    src/geometry/Geometry2D.cpp
    src/geometry/ConstructionIdentity2D.cpp
    src/geometry/ConstructionRecovery2D.cpp
    src/geometry/IntersectionRegistry2D.cpp
    src/geometry/IntersectionConstruction2D.cpp
    src/geometry/BoundarySimplification2D.cpp
    src/grid/CartesianGrid2D.cpp
    src/hybrid/HybridMesh2D.cpp
    src/hybrid/TransitionCanonicalization2D.cpp
    src/spatial/BoundarySegmentIndex2D.cpp
    src/quadtree/Quadtree2D.cpp
    src/sizing/SizeField2D.cpp
    src/cutcell/CutCell2D.cpp
    src/topology/Topology2D.cpp
    src/topology/EdgeIncidence2D.cpp
    src/topology/PatchTransaction2D.cpp
    src/topology/SharedEdgePartition2D.cpp
    src/stabilization/SmallCell2D.cpp
    src/stabilization/Agglomeration2D.cpp
    src/quality/Quality2D.cpp
    src/quality/QualityContract2D.cpp
    src/quality/SolverQuality2D.cpp
    src/quality/PatchLocalQuality2D.cpp
    src/quality/SolverTopology2D.cpp
    src/io/MeshIO2D.cpp
    src/io/BoundaryMetadata2D.cpp
    src/io/OpenFoam2D.cpp
    src/io/Dxf2D.cpp
)
add_library(cartmesh2d::cartmesh2d ALIAS cartmesh2d)

target_include_directories(cartmesh2d
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_compile_features(cartmesh2d PUBLIC cxx_std_20)

if(COMMAND cartmesh_enable_warnings)
    cartmesh_enable_warnings(cartmesh2d)
elseif(MSVC)
    target_compile_options(cartmesh2d PRIVATE /W4 /permissive-)
else()
    target_compile_options(cartmesh2d PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

add_executable(cartmesh2d_cli apps/cartmesh2d_cli.cpp)
target_link_libraries(cartmesh2d_cli PRIVATE cartmesh2d::cartmesh2d)
target_compile_features(cartmesh2d_cli PRIVATE cxx_std_20)
if(COMMAND cartmesh_enable_warnings)
    cartmesh_enable_warnings(cartmesh2d_cli)
elseif(MSVC)
    target_compile_options(cartmesh2d_cli PRIVATE /W4 /permissive-)
else()
    target_compile_options(cartmesh2d_cli PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

add_executable(cartmesh2d_dxf_cli apps/cartmesh2d_dxf_cli.cpp)
target_link_libraries(cartmesh2d_dxf_cli PRIVATE cartmesh2d::cartmesh2d)
target_compile_features(cartmesh2d_dxf_cli PRIVATE cxx_std_20)
if(COMMAND cartmesh_enable_warnings)
    cartmesh_enable_warnings(cartmesh2d_dxf_cli)
elseif(MSVC)
    target_compile_options(cartmesh2d_dxf_cli PRIVATE /W4 /permissive-)
else()
    target_compile_options(cartmesh2d_dxf_cli PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

add_executable(cartmesh2d_boundary_layer_cli apps/cartmesh2d_boundary_layer_cli.cpp)
target_link_libraries(cartmesh2d_boundary_layer_cli PRIVATE cartmesh2d::cartmesh2d)
target_compile_features(cartmesh2d_boundary_layer_cli PRIVATE cxx_std_20)
if(COMMAND cartmesh_enable_warnings)
    cartmesh_enable_warnings(cartmesh2d_boundary_layer_cli)
elseif(MSVC)
    target_compile_options(cartmesh2d_boundary_layer_cli PRIVATE /W4 /permissive-)
else()
    target_compile_options(cartmesh2d_boundary_layer_cli PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

add_executable(cartmesh2d_hybrid_cli apps/cartmesh2d_hybrid_cli.cpp)
target_link_libraries(cartmesh2d_hybrid_cli PRIVATE cartmesh2d::cartmesh2d)
target_compile_features(cartmesh2d_hybrid_cli PRIVATE cxx_std_20)
if(COMMAND cartmesh_enable_warnings)
    cartmesh_enable_warnings(cartmesh2d_hybrid_cli)
elseif(MSVC)
    target_compile_options(cartmesh2d_hybrid_cli PRIVATE /W4 /permissive-)
else()
    target_compile_options(cartmesh2d_hybrid_cli PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

option(CARTMESH2D_BUILD_TESTS "Build cartmesh2d tests" ON)

if(CARTMESH2D_BUILD_TESTS)
    include(CTest)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_executable(cartmesh2d_geometry_tests tests/geometry_test.cpp)
    target_link_libraries(cartmesh2d_geometry_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_geometry_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_grid_tests tests/grid_test.cpp)
    target_link_libraries(cartmesh2d_grid_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_grid_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_quadtree_tests tests/quadtree_test.cpp)
    target_link_libraries(cartmesh2d_quadtree_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_quadtree_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_spatial_index_tests tests/spatial_index_test.cpp)
    target_link_libraries(cartmesh2d_spatial_index_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_spatial_index_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_size_field_tests tests/size_field_test.cpp)
    target_link_libraries(cartmesh2d_size_field_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_size_field_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_multiloop_tests tests/multiloop_test.cpp)
    target_link_libraries(cartmesh2d_multiloop_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_multiloop_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_solver_export_tests tests/solver_export_test.cpp)
    target_link_libraries(cartmesh2d_solver_export_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_solver_export_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_cutcell_tests tests/cutcell_test.cpp)
    target_link_libraries(cartmesh2d_cutcell_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_cutcell_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_cutcell_components_tests tests/cutcell_components_test.cpp)
    target_link_libraries(cartmesh2d_cutcell_components_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_cutcell_components_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_cutcell_quadtree_audit_tests tests/cutcell_quadtree_audit_test.cpp)
    target_link_libraries(cartmesh2d_cutcell_quadtree_audit_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_cutcell_quadtree_audit_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_topology_tests tests/topology_test.cpp)
    target_link_libraries(cartmesh2d_topology_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_topology_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_topology_quadtree_audit_tests tests/topology_quadtree_audit_test.cpp)
    target_link_libraries(cartmesh2d_topology_quadtree_audit_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_topology_quadtree_audit_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_small_cell_tests tests/small_cell_test.cpp)
    target_link_libraries(cartmesh2d_small_cell_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_small_cell_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_agglomeration_tests tests/agglomeration_test.cpp)
    target_link_libraries(cartmesh2d_agglomeration_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_agglomeration_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_quality_tests tests/quality_test.cpp)
    target_link_libraries(cartmesh2d_quality_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_quality_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_quality_contract_tests tests/quality_contract_test.cpp)
    target_link_libraries(cartmesh2d_quality_contract_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_quality_contract_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_io_tests tests/io_test.cpp)
    target_link_libraries(cartmesh2d_io_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_io_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_dxf_tests tests/dxf_test.cpp)
    target_link_libraries(cartmesh2d_dxf_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_dxf_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_boundary_layer_tests tests/boundary_layer_test.cpp)
    target_link_libraries(cartmesh2d_boundary_layer_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_boundary_layer_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_hybrid_mesh_tests tests/hybrid_mesh_test.cpp)
    target_link_libraries(cartmesh2d_hybrid_mesh_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_hybrid_mesh_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_h4_local_termination_tests
                   tests/h4_local_termination_test.cpp)
    target_link_libraries(cartmesh2d_h4_local_termination_tests
                          PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_h4_local_termination_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_weld_budget_tests tests/weld_budget_test.cpp)
    target_link_libraries(cartmesh2d_weld_budget_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_weld_budget_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_refinement_attribution_tests
                   tests/refinement_attribution_test.cpp)
    target_link_libraries(cartmesh2d_refinement_attribution_tests
                          PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_refinement_attribution_tests
                            PRIVATE cxx_std_20)

    foreach(target
            cartmesh2d_geometry_tests
            cartmesh2d_grid_tests
            cartmesh2d_quadtree_tests
            cartmesh2d_spatial_index_tests
            cartmesh2d_size_field_tests
            cartmesh2d_multiloop_tests
            cartmesh2d_solver_export_tests
            cartmesh2d_cutcell_tests
            cartmesh2d_cutcell_components_tests
            cartmesh2d_cutcell_quadtree_audit_tests
            cartmesh2d_topology_tests
            cartmesh2d_topology_quadtree_audit_tests
            cartmesh2d_small_cell_tests
            cartmesh2d_agglomeration_tests
            cartmesh2d_quality_tests
            cartmesh2d_quality_contract_tests
            cartmesh2d_io_tests
            cartmesh2d_dxf_tests
            cartmesh2d_boundary_layer_tests
            cartmesh2d_hybrid_mesh_tests
            cartmesh2d_h4_local_termination_tests
            cartmesh2d_refinement_attribution_tests
            cartmesh2d_weld_budget_tests)
        if(COMMAND cartmesh_enable_warnings)
            cartmesh_enable_warnings(${target})
        elseif(MSVC)
            target_compile_options(${target} PRIVATE /W4 /permissive-)
        else()
            target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
        endif()
    endforeach()

    add_test(NAME cartmesh2d_stage0_geometry_tests COMMAND cartmesh2d_geometry_tests)
    add_executable(cartmesh2d_intersection_registry_tests tests/intersection_registry_test.cpp)
    target_link_libraries(cartmesh2d_intersection_registry_tests PRIVATE cartmesh2d::cartmesh2d)
    add_test(NAME cartmesh2d_stageQ2_intersection_registry_tests COMMAND cartmesh2d_intersection_registry_tests)
    add_executable(cartmesh2d_shared_construction_tests tests/shared_construction_test.cpp)
    target_link_libraries(cartmesh2d_shared_construction_tests PRIVATE cartmesh2d::cartmesh2d)
    add_test(NAME cartmesh2d_stageQ2A_shared_construction_tests COMMAND cartmesh2d_shared_construction_tests)
    add_test(NAME cartmesh2d_stage1_grid_tests COMMAND cartmesh2d_grid_tests)
    add_test(NAME cartmesh2d_stage2_quadtree_tests COMMAND cartmesh2d_quadtree_tests)
    add_test(NAME cartmesh2d_stageP1_spatial_index_tests COMMAND cartmesh2d_spatial_index_tests)
    add_test(NAME cartmesh2d_stageWB_size_field_tests COMMAND cartmesh2d_size_field_tests)
    add_test(NAME cartmesh2d_stageD1_multiloop_tests COMMAND cartmesh2d_multiloop_tests)
    add_test(NAME cartmesh2d_stageS1_solver_export_tests COMMAND cartmesh2d_solver_export_tests)
    add_test(NAME cartmesh2d_stage3_cutcell_tests COMMAND cartmesh2d_cutcell_tests)
    add_test(NAME cartmesh2d_stage3_cutcell_components_tests COMMAND cartmesh2d_cutcell_components_tests)
    add_test(NAME cartmesh2d_stage3_quadtree_cutcell_audit_tests COMMAND cartmesh2d_cutcell_quadtree_audit_tests)
    add_test(NAME cartmesh2d_stage4_topology_tests COMMAND cartmesh2d_topology_tests)
    add_test(NAME cartmesh2d_stage4_adaptive_topology_audit_tests COMMAND cartmesh2d_topology_quadtree_audit_tests)
    add_test(NAME cartmesh2d_stage5a_small_cell_tests COMMAND cartmesh2d_small_cell_tests)
    add_test(NAME cartmesh2d_stage5b_agglomeration_tests COMMAND cartmesh2d_agglomeration_tests)
    add_test(NAME cartmesh2d_stage6a_quality_tests COMMAND cartmesh2d_quality_tests)
    add_test(NAME cartmesh2d_stageQ1_quality_contract_tests
             COMMAND cartmesh2d_quality_contract_tests)
    add_test(NAME cartmesh2d_stage6b_io_tests COMMAND cartmesh2d_io_tests)
    add_test(NAME cartmesh2d_stageDXF2_import_tests COMMAND cartmesh2d_dxf_tests)
    add_test(NAME cartmesh2d_stageH4_1_boundary_layer_tests
             COMMAND cartmesh2d_boundary_layer_tests)
    add_test(NAME cartmesh2d_stageH4_2_conformal_hybrid_tests
             COMMAND cartmesh2d_hybrid_mesh_tests)
    add_test(NAME cartmesh2d_stageH4_3_local_termination_tests
             COMMAND cartmesh2d_h4_local_termination_tests)
    add_test(NAME cartmesh2d_stageR2_weld_budget_tests
             COMMAND cartmesh2d_weld_budget_tests)
    add_test(NAME cartmesh2d_stageR2_refinement_attribution_tests
             COMMAND cartmesh2d_refinement_attribution_tests)
    add_test(NAME cartmesh2d_stageR2_refinement_ladder_gate_tests
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/refinement_ladder_gate_test.py)
    add_test(NAME cartmesh2d_product_cutcell_help
        COMMAND cartmesh2d_cli --help)
    add_test(NAME cartmesh2d_product_hybrid_help
        COMMAND cartmesh2d_hybrid_cli --help)
    add_test(NAME cartmesh2d_product_dxf_help
        COMMAND cartmesh2d_dxf_cli --help)

    set(CARTMESH2D_H4_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/boundary_layer)
    add_test(NAME cartmesh2d_stageH4_1_circle_e2e
        COMMAND cartmesh2d_boundary_layer_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy
                ${CARTMESH2D_H4_OUTPUT}/circle 4 first 0.02 1.2 exterior)
    add_test(NAME cartmesh2d_stageH4_1_circle_independent_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_boundary_layer2d.py
                ${CARTMESH2D_H4_OUTPUT}/circle.layer.vtk
                ${CARTMESH2D_H4_OUTPUT}/circle.layer.json
                --report ${CARTMESH2D_H4_OUTPUT}/circle.independent.json)
    set_tests_properties(cartmesh2d_stageH4_1_circle_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageH4_1_circle_e2e)

    add_test(NAME cartmesh2d_stageH4_1_rectangle_e2e
        COMMAND cartmesh2d_boundary_layer_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/rectangle.xy
                ${CARTMESH2D_H4_OUTPUT}/rectangle 3 first 0.04 1.25 exterior)
    add_test(NAME cartmesh2d_stageH4_1_rectangle_independent_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_boundary_layer2d.py
                ${CARTMESH2D_H4_OUTPUT}/rectangle.layer.vtk
                ${CARTMESH2D_H4_OUTPUT}/rectangle.layer.json
                --report ${CARTMESH2D_H4_OUTPUT}/rectangle.independent.json)
    set_tests_properties(cartmesh2d_stageH4_1_rectangle_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageH4_1_rectangle_e2e)

    add_test(NAME cartmesh2d_stageH4_1_concave_fails_closed
        COMMAND cartmesh2d_boundary_layer_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/concave.xy
                ${CARTMESH2D_H4_OUTPUT}/concave 3 first 0.02 1.1 exterior)
    set_tests_properties(cartmesh2d_stageH4_1_concave_fails_closed PROPERTIES
        WILL_FAIL TRUE)
    add_test(NAME cartmesh2d_stageH4_1_concave_failure_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_boundary_layer_failure.py
                ${CARTMESH2D_H4_OUTPUT}/concave.layer.json
                --reason concave_corner)
    set_tests_properties(cartmesh2d_stageH4_1_concave_failure_reader PROPERTIES
        DEPENDS cartmesh2d_stageH4_1_concave_fails_closed)

    add_test(NAME cartmesh2d_stageH4_1_circle_repeat_e2e
        COMMAND cartmesh2d_boundary_layer_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy
                ${CARTMESH2D_H4_OUTPUT}/circle_repeat 4 first 0.02 1.2 exterior)
    add_test(NAME cartmesh2d_stageH4_1_circle_vtk_deterministic
        COMMAND ${CMAKE_COMMAND} -E compare_files
                ${CARTMESH2D_H4_OUTPUT}/circle.layer.vtk
                ${CARTMESH2D_H4_OUTPUT}/circle_repeat.layer.vtk)
    set_tests_properties(cartmesh2d_stageH4_1_circle_vtk_deterministic PROPERTIES
        DEPENDS "cartmesh2d_stageH4_1_circle_e2e;cartmesh2d_stageH4_1_circle_repeat_e2e")
    add_test(NAME cartmesh2d_stageH4_1_circle_report_deterministic
        COMMAND ${CMAKE_COMMAND} -E compare_files
                ${CARTMESH2D_H4_OUTPUT}/circle.layer.json
                ${CARTMESH2D_H4_OUTPUT}/circle_repeat.layer.json)
    set_tests_properties(cartmesh2d_stageH4_1_circle_report_deterministic PROPERTIES
        DEPENDS "cartmesh2d_stageH4_1_circle_e2e;cartmesh2d_stageH4_1_circle_repeat_e2e")

    set(CARTMESH2D_H4_2_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hybrid_mesh)
    add_test(NAME cartmesh2d_stageH4_2_circle_e2e
        COMMAND cartmesh2d_hybrid_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy
                ${CARTMESH2D_H4_2_OUTPUT}/circle 6 3 6 4 0.02 1.2 1.0)
    add_test(NAME cartmesh2d_stageH4_2_circle_independent_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_hybrid_mesh2d.py
                ${CARTMESH2D_H4_2_OUTPUT}/circle.hybrid.vtk
                ${CARTMESH2D_H4_2_OUTPUT}/circle.hybrid.json
                --output ${CARTMESH2D_H4_2_OUTPUT}/circle.independent.json)
    set_tests_properties(cartmesh2d_stageH4_2_circle_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageH4_2_circle_e2e)

    add_test(NAME cartmesh2d_stageH4_2_superellipse_e2e
        COMMAND cartmesh2d_hybrid_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/superellipse_24.xy
                ${CARTMESH2D_H4_2_OUTPUT}/superellipse 6 3 6 3 0.015 1.15 1.0)
    add_test(NAME cartmesh2d_stageH4_2_superellipse_independent_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_hybrid_mesh2d.py
                ${CARTMESH2D_H4_2_OUTPUT}/superellipse.hybrid.vtk
                ${CARTMESH2D_H4_2_OUTPUT}/superellipse.hybrid.json
                --output ${CARTMESH2D_H4_2_OUTPUT}/superellipse.independent.json)
    set_tests_properties(cartmesh2d_stageH4_2_superellipse_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageH4_2_superellipse_e2e)

    add_test(NAME cartmesh2d_stageH4_2_circle_repeat_e2e
        COMMAND cartmesh2d_hybrid_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy
                ${CARTMESH2D_H4_2_OUTPUT}/circle_repeat 6 3 6 4 0.02 1.2 1.0)
    add_test(NAME cartmesh2d_stageH4_2_circle_topology_deterministic
        COMMAND ${CMAKE_COMMAND} -E compare_files
                ${CARTMESH2D_H4_2_OUTPUT}/circle.hybrid.cm2d
                ${CARTMESH2D_H4_2_OUTPUT}/circle_repeat.hybrid.cm2d)
    set_tests_properties(cartmesh2d_stageH4_2_circle_topology_deterministic PROPERTIES
        DEPENDS "cartmesh2d_stageH4_2_circle_e2e;cartmesh2d_stageH4_2_circle_repeat_e2e")
    add_test(NAME cartmesh2d_stageH4_2_circle_report_deterministic
        COMMAND ${CMAKE_COMMAND} -E compare_files
                ${CARTMESH2D_H4_2_OUTPUT}/circle.hybrid.json
                ${CARTMESH2D_H4_2_OUTPUT}/circle_repeat.hybrid.json)
    set_tests_properties(cartmesh2d_stageH4_2_circle_report_deterministic PROPERTIES
        DEPENDS "cartmesh2d_stageH4_2_circle_e2e;cartmesh2d_stageH4_2_circle_repeat_e2e")

    set(CARTMESH2D_DXF_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/dxf)
    add_test(NAME cartmesh2d_stageDXF1_convert_e2e
        COMMAND cartmesh2d_dxf_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/dxf/rounded_obstacle.dxf
                ${CARTMESH2D_DXF_OUTPUT}/rounded_obstacle.xy 0.002
                ${CARTMESH2D_DXF_OUTPUT}/rounded_obstacle.dxf.json)
    add_test(NAME cartmesh2d_stageDXF1_mesh_e2e
        COMMAND cartmesh2d_cli
                ${CARTMESH2D_DXF_OUTPUT}/airfoil_like.xy
                ${CARTMESH2D_DXF_OUTPUT}/airfoil_like 7 0.30 0.02 exterior
                ${CARTMESH2D_DXF_OUTPUT}/airfoil_like-case 0 0)
    set_tests_properties(cartmesh2d_stageDXF1_mesh_e2e PROPERTIES
        DEPENDS cartmesh2d_stageDXF1_convert_airfoil_e2e)
    add_test(NAME cartmesh2d_stageDXF1_independent_openfoam_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_DXF_OUTPUT}/airfoil_like-case
                --report ${CARTMESH2D_DXF_OUTPUT}/airfoil_like-case/independent_check.json)
    set_tests_properties(cartmesh2d_stageDXF1_independent_openfoam_reader PROPERTIES
        DEPENDS cartmesh2d_stageDXF1_mesh_e2e)
    add_test(NAME cartmesh2d_stageDXF1_convert_airfoil_e2e
        COMMAND cartmesh2d_dxf_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/dxf/airfoil_like.dxf
                ${CARTMESH2D_DXF_OUTPUT}/airfoil_like.xy 0.001
                ${CARTMESH2D_DXF_OUTPUT}/airfoil_like.dxf.json)

    add_test(NAME cartmesh2d_stageDXF2_convert_spline_mm_e2e
        COMMAND cartmesh2d_dxf_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/dxf/spline_circle_mm.dxf
                ${CARTMESH2D_DXF_OUTPUT}/spline_circle.xy 0.001
                ${CARTMESH2D_DXF_OUTPUT}/spline_circle.dxf.json)
    add_test(NAME cartmesh2d_stageDXF2_mesh_spline_mm_e2e
        COMMAND cartmesh2d_cli
                ${CARTMESH2D_DXF_OUTPUT}/spline_circle.xy
                ${CARTMESH2D_DXF_OUTPUT}/spline_circle 6 0.25 0.05 exterior
                ${CARTMESH2D_DXF_OUTPUT}/spline_circle-case 0 0)
    set_tests_properties(cartmesh2d_stageDXF2_mesh_spline_mm_e2e PROPERTIES
        DEPENDS cartmesh2d_stageDXF2_convert_spline_mm_e2e)
    add_test(NAME cartmesh2d_stageDXF2_independent_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_DXF_OUTPUT}/spline_circle-case
                --require-patch wall_spline:wall
                --report ${CARTMESH2D_DXF_OUTPUT}/spline_circle-case/independent_check.json)
    set_tests_properties(cartmesh2d_stageDXF2_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageDXF2_mesh_spline_mm_e2e)

    add_test(NAME cartmesh2d_stageDXF2_convert_symmetry_e2e
        COMMAND cartmesh2d_dxf_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/dxf/symmetry_circle_m.dxf
                ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle.xy 0.01
                ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle.dxf.json)
    add_test(NAME cartmesh2d_stageDXF2_mesh_symmetry_e2e
        COMMAND cartmesh2d_cli
                ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle.xy
                ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle 6 0.25 0.05 exterior
                ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle-case 0 0)
    set_tests_properties(cartmesh2d_stageDXF2_mesh_symmetry_e2e PROPERTIES
        DEPENDS cartmesh2d_stageDXF2_convert_symmetry_e2e)
    add_test(NAME cartmesh2d_stageDXF2_symmetry_patch_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle-case
                --require-patch symmetry_curve:symmetryPlane
                --require-field-boundary symmetry_curve:symmetryPlane)
    set_tests_properties(cartmesh2d_stageDXF2_symmetry_patch_reader PROPERTIES
        DEPENDS cartmesh2d_stageDXF2_mesh_symmetry_e2e)

    set(CARTMESH2D_ACCEPTANCE_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/acceptance)
    add_test(NAME cartmesh2d_stage6c_rectangle_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/rectangle.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/rectangle 5 0.25 0.05)
    add_test(NAME cartmesh2d_stage6c_circle_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/circle 5 0.25 0.05)
    add_test(NAME cartmesh2d_stage6c_concave_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/concave.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/concave 5 0.25 0.05)
    add_test(NAME cartmesh2d_stage6c_airfoil_like_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/airfoil_like 6 0.30 0.02)

    # R2/W1 regression: the shared construction path must remove the level-10
    # grid-corner spur, preserve wall patch identity through stabilization and
    # solver partitioning, and emit independently readable OpenFOAM topology.
    set(CARTMESH2D_W1_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/r2_w1)
    add_test(NAME cartmesh2d_stageR2_w1_circle_l10_e2e
        COMMAND cartmesh2d_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy
                ${CARTMESH2D_W1_OUTPUT}/circle_l10 10 0.25 0.10 exterior
                ${CARTMESH2D_W1_OUTPUT}/circle_l10-case 8)
    add_test(NAME cartmesh2d_stageR2_w1_circle_l10_independent_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_W1_OUTPUT}/circle_l10-case
                --report ${CARTMESH2D_W1_OUTPUT}/circle_l10-case/independent_check.json)
    set_tests_properties(cartmesh2d_stageR2_w1_circle_l10_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageR2_w1_circle_l10_e2e)

    set(CARTMESH2D_SIZING_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/sizing_field)
    add_test(NAME cartmesh2d_stageH1_sizing_field_e2e
        COMMAND cartmesh2d_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy
                ${CARTMESH2D_SIZING_OUTPUT}/airfoil 7 0.30 0.02 exterior - 3 0
                --distance-band 0.05 6
                --refine-box 0.80 -0.10 1.30 0.10 6)
    add_test(NAME cartmesh2d_stageH1_rejects_outside_box
        COMMAND cartmesh2d_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy
                ${CARTMESH2D_SIZING_OUTPUT}/invalid 6 0.30 0.02 exterior - 2 0
                --refine-box 10 10 11 11 5)
    set_tests_properties(cartmesh2d_stageH1_rejects_outside_box PROPERTIES WILL_FAIL TRUE)

    # The size field owns the domain and the tree depth, so it needs its own e2e
    # coverage: the positional max-level/padding-fraction above cannot exercise it.
    add_test(NAME cartmesh2d_stageWB_size_field_defaults_e2e
        COMMAND cartmesh2d_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy
                ${CARTMESH2D_SIZING_OUTPUT}/sf_circle 8 0.25 0.15 exterior
                ${CARTMESH2D_SIZING_OUTPUT}/sf_circle-case 0 0
                --size-field --cells-per-level 3)
    add_test(NAME cartmesh2d_stageWB_size_field_wake_e2e
        COMMAND cartmesh2d_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy
                ${CARTMESH2D_SIZING_OUTPUT}/sf_airfoil 8 0.25 0.15 exterior
                ${CARTMESH2D_SIZING_OUTPUT}/sf_airfoil-case 0 0
                --size-field --cells-per-level 3 --wake 4 6 0.6 4)
    # A wall request past the measured construction ceiling must fail closed rather
    # than silently emit a mesh whose gate outcome depends on grid alignment.
    add_test(NAME cartmesh2d_stageWB_rejects_unsafe_wall_level
        COMMAND cartmesh2d_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy
                ${CARTMESH2D_SIZING_OUTPUT}/sf_unsafe 8 0.25 0.15 exterior - 0 0
                --size-field --wall-cells-per-span 4096)
    set_tests_properties(cartmesh2d_stageWB_rejects_unsafe_wall_level
        PROPERTIES WILL_FAIL TRUE)
    # The two spellings of one field must not be combinable, or the reported sizing
    # would not describe the mesh that was built.
    add_test(NAME cartmesh2d_stageWB_rejects_mixed_sizing_spellings
        COMMAND cartmesh2d_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy
                ${CARTMESH2D_SIZING_OUTPUT}/sf_mixed 8 0.25 0.15 exterior - 0 0
                --size-field --distance-band 0.05 6)
    set_tests_properties(cartmesh2d_stageWB_rejects_mixed_sizing_spellings
        PROPERTIES WILL_FAIL TRUE)

    set(CARTMESH2D_ROBUSTNESS_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/robustness)
    add_test(NAME cartmesh2d_stageR_gear_star_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/gear_star.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/gear_star 8 0.25 0.01)
    add_test(NAME cartmesh2d_stageR_serpentine_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/serpentine_body.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/serpentine 7 0.20 0.01)
    add_test(NAME cartmesh2d_stageR_nozzle_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/nozzle_profile.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/nozzle 7 0.20 0.01)
    add_test(NAME cartmesh2d_stageR_naca2412_dense_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/naca2412_dense.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/naca2412_dense 8 0.30 0.005)
    add_test(NAME cartmesh2d_stageR_superellipse_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/superellipse_24.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/superellipse 7 0.25 0.01)
    add_test(NAME cartmesh2d_stageD1_two_obstacles_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/two_obstacles.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/two_obstacles 7 0.25 0.05 exterior)
    add_test(NAME cartmesh2d_stageD1_annulus_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/annulus.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/annulus 7 0.25 0.05 interior)
    add_test(NAME cartmesh2d_stageD1_nested_island_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/nested_island.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/nested_island 7 0.25 0.05 interior)
    set(CARTMESH2D_SOLVER_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/solver_ready)
    add_test(NAME cartmesh2d_stageS1_rectangle_openfoam_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/rectangle.xy
                ${CARTMESH2D_SOLVER_OUTPUT}/rectangle 5 0.25 0.05 exterior
                ${CARTMESH2D_SOLVER_OUTPUT}/rectangle-case)
    add_test(NAME cartmesh2d_stageS1_independent_openfoam_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_SOLVER_OUTPUT}/rectangle-case
                --report ${CARTMESH2D_SOLVER_OUTPUT}/rectangle-case/independent_check.json)
    set_tests_properties(cartmesh2d_stageS1_independent_openfoam_reader PROPERTIES
        DEPENDS cartmesh2d_stageS1_rectangle_openfoam_e2e)
    add_test(NAME cartmesh2d_stageV1_annulus_openfoam_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/annulus.xy
                ${CARTMESH2D_SOLVER_OUTPUT}/annulus 7 0.25 0.05 interior
                ${CARTMESH2D_SOLVER_OUTPUT}/annulus-case)
    add_test(NAME cartmesh2d_stageV1_annulus_independent_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_SOLVER_OUTPUT}/annulus-case)
    set_tests_properties(cartmesh2d_stageV1_annulus_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageV1_annulus_openfoam_e2e)
    add_test(NAME cartmesh2d_stageV1_nested_openfoam_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/nested_island.xy
                ${CARTMESH2D_SOLVER_OUTPUT}/nested 7 0.25 0.05 interior
                ${CARTMESH2D_SOLVER_OUTPUT}/nested-case)
    add_test(NAME cartmesh2d_stageV1_nested_independent_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_SOLVER_OUTPUT}/nested-case)
    set_tests_properties(cartmesh2d_stageV1_nested_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageV1_nested_openfoam_e2e)
    add_test(NAME cartmesh2d_stageV1b_circle_openfoam_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy
                ${CARTMESH2D_SOLVER_OUTPUT}/circle 6 0.25 0.05 exterior
                ${CARTMESH2D_SOLVER_OUTPUT}/circle-case)
    add_test(NAME cartmesh2d_stageV1b_circle_independent_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_SOLVER_OUTPUT}/circle-case)
    set_tests_properties(cartmesh2d_stageV1b_circle_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageV1b_circle_openfoam_e2e)
    add_test(NAME cartmesh2d_stageV1b_airfoil_openfoam_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy
                ${CARTMESH2D_SOLVER_OUTPUT}/airfoil 7 0.30 0.02 exterior
                ${CARTMESH2D_SOLVER_OUTPUT}/airfoil-case)
    add_test(NAME cartmesh2d_stageV1b_airfoil_independent_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_SOLVER_OUTPUT}/airfoil-case)
    set_tests_properties(cartmesh2d_stageV1b_airfoil_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageV1b_airfoil_openfoam_e2e)
    add_test(NAME cartmesh2d_stageV1b_two_obstacles_openfoam_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/two_obstacles.xy
                ${CARTMESH2D_SOLVER_OUTPUT}/two_obstacles 8 0.25 0.05 exterior
                ${CARTMESH2D_SOLVER_OUTPUT}/two_obstacles-case)
    add_test(NAME cartmesh2d_stageV1b_two_obstacles_independent_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CARTMESH2D_SOLVER_OUTPUT}/two_obstacles-case)
    set_tests_properties(cartmesh2d_stageV1b_two_obstacles_independent_reader PROPERTIES
        DEPENDS cartmesh2d_stageV1b_two_obstacles_openfoam_e2e)
    add_test(NAME cartmesh2d_stageD1_rejects_intersecting_loops
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/invalid/intersecting_loops.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/invalid_intersecting_loops 6 0.25 0.05 interior)
    set_tests_properties(cartmesh2d_stageD1_rejects_intersecting_loops PROPERTIES WILL_FAIL TRUE)
endif()
