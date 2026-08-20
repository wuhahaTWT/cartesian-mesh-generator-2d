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

    foreach(target
            cartmesh2d_geometry_tests
            cartmesh2d_grid_tests
            cartmesh2d_quadtree_tests
            cartmesh2d_cutcell_tests
            cartmesh2d_cutcell_quadtree_audit_tests
            cartmesh2d_topology_tests
            cartmesh2d_topology_quadtree_audit_tests
            cartmesh2d_small_cell_tests
            cartmesh2d_agglomeration_tests
            cartmesh2d_quality_tests)
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
    add_test(NAME cartmesh2d_stage3_quadtree_cutcell_audit_tests
             COMMAND cartmesh2d_cutcell_quadtree_audit_tests)
    add_test(NAME cartmesh2d_stage4_topology_tests COMMAND cartmesh2d_topology_tests)
    add_test(NAME cartmesh2d_stage4_adaptive_topology_audit_tests
             COMMAND cartmesh2d_topology_quadtree_audit_tests)
    add_test(NAME cartmesh2d_stage5a_small_cell_tests COMMAND cartmesh2d_small_cell_tests)
    add_test(NAME cartmesh2d_stage5b_agglomeration_tests COMMAND cartmesh2d_agglomeration_tests)
    add_test(NAME cartmesh2d_stage6a_quality_tests COMMAND cartmesh2d_quality_tests)
endif()
