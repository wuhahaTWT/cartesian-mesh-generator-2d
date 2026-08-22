cmake_minimum_required(VERSION 3.20)

project(cartmesh2d LANGUAGES CXX)

add_library(cartmesh2d
    src/geometry/Geometry2D.cpp
    src/grid/CartesianGrid2D.cpp
    src/quadtree/Quadtree2D.cpp
    src/cutcell/CutCell2D.cpp
    src/topology/Topology2D.cpp
    src/stabilization/SmallCell2D.cpp
    src/stabilization/Agglomeration2D.cpp
    src/quality/Quality2D.cpp
    src/io/MeshIO2D.cpp
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

option(CARTMESH2D_BUILD_TESTS "Build cartmesh2d tests" ON)

if(CARTMESH2D_BUILD_TESTS)
    include(CTest)
    add_executable(cartmesh2d_geometry_tests tests/geometry_test.cpp)
    target_link_libraries(cartmesh2d_geometry_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_geometry_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_grid_tests tests/grid_test.cpp)
    target_link_libraries(cartmesh2d_grid_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_grid_tests PRIVATE cxx_std_20)

    add_executable(cartmesh2d_quadtree_tests tests/quadtree_test.cpp)
    target_link_libraries(cartmesh2d_quadtree_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_quadtree_tests PRIVATE cxx_std_20)

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

    add_executable(cartmesh2d_io_tests tests/io_test.cpp)
    target_link_libraries(cartmesh2d_io_tests PRIVATE cartmesh2d::cartmesh2d)
    target_compile_features(cartmesh2d_io_tests PRIVATE cxx_std_20)

    foreach(target
            cartmesh2d_geometry_tests
            cartmesh2d_grid_tests
            cartmesh2d_quadtree_tests
            cartmesh2d_cutcell_tests
            cartmesh2d_cutcell_components_tests
            cartmesh2d_cutcell_quadtree_audit_tests
            cartmesh2d_topology_tests
            cartmesh2d_topology_quadtree_audit_tests
            cartmesh2d_small_cell_tests
            cartmesh2d_agglomeration_tests
            cartmesh2d_quality_tests
            cartmesh2d_io_tests)
        if(COMMAND cartmesh_enable_warnings)
            cartmesh_enable_warnings(${target})
        elseif(MSVC)
            target_compile_options(${target} PRIVATE /W4 /permissive-)
        else()
            target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
        endif()
    endforeach()

    add_test(NAME cartmesh2d_stage0_geometry_tests COMMAND cartmesh2d_geometry_tests)
    add_test(NAME cartmesh2d_stage1_grid_tests COMMAND cartmesh2d_grid_tests)
    add_test(NAME cartmesh2d_stage2_quadtree_tests COMMAND cartmesh2d_quadtree_tests)
    add_test(NAME cartmesh2d_stage3_cutcell_tests COMMAND cartmesh2d_cutcell_tests)
    add_test(NAME cartmesh2d_stage3_cutcell_components_tests COMMAND cartmesh2d_cutcell_components_tests)
    add_test(NAME cartmesh2d_stage3_quadtree_cutcell_audit_tests COMMAND cartmesh2d_cutcell_quadtree_audit_tests)
    add_test(NAME cartmesh2d_stage4_topology_tests COMMAND cartmesh2d_topology_tests)
    add_test(NAME cartmesh2d_stage4_adaptive_topology_audit_tests COMMAND cartmesh2d_topology_quadtree_audit_tests)
    add_test(NAME cartmesh2d_stage5a_small_cell_tests COMMAND cartmesh2d_small_cell_tests)
    add_test(NAME cartmesh2d_stage5b_agglomeration_tests COMMAND cartmesh2d_agglomeration_tests)
    add_test(NAME cartmesh2d_stage6a_quality_tests COMMAND cartmesh2d_quality_tests)
    add_test(NAME cartmesh2d_stage6b_io_tests COMMAND cartmesh2d_io_tests)

    set(CARTMESH2D_ACCEPTANCE_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/acceptance)
    add_test(NAME cartmesh2d_stage6c_rectangle_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/rectangle.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/rectangle 5 0.25 0.05)
    add_test(NAME cartmesh2d_stage6c_circle_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/circle 5 0.25 0.05)
    add_test(NAME cartmesh2d_stage6c_concave_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/concave.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/concave 5 0.25 0.05)
    add_test(NAME cartmesh2d_stage6c_airfoil_like_e2e
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/airfoil_like 6 0.30 0.02)

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
endif()
