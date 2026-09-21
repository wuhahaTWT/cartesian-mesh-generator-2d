cmake_minimum_required(VERSION 3.20)

project(cartmesh2d LANGUAGES CXX)

if(MSVC)
    # Ship self-contained command-line tools and make narrow source/argv text
    # consistently UTF-8 on supported Windows 10 and newer systems.
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    add_compile_options(/utf-8)
    add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN)
endif()

add_library(cartmesh2d
    src/boundary_layer/BoundaryLayer2D.cpp
    src/geometry/Geometry2D.cpp
    src/geometry/ConstructionIdentity2D.cpp
    src/geometry/ConstructionRecovery2D.cpp
    src/geometry/IntersectionRegistry2D.cpp
    src/geometry/IntersectionConstruction2D.cpp
    src/geometry/BoundarySimplification2D.cpp
    src/grid/CartesianGrid2D.cpp
    src/grid/RectilinearMesh2D.cpp
    src/hybrid/HybridMesh2D.cpp
    src/hybrid/TransitionCanonicalization2D.cpp
    src/spatial/BoundarySegmentIndex2D.cpp
    src/quadtree/Quadtree2D.cpp
    src/sizing/SizeField2D.cpp
    src/sizing/MeshResolution2D.cpp
    src/sizing/BoundaryLayerResolution2D.cpp
    src/cutcell/CutCell2D.cpp
    src/topology/Topology2D.cpp
    src/topology/EdgeIncidence2D.cpp
    src/topology/PatchTransaction2D.cpp
    src/topology/SharedEdgePartition2D.cpp
    src/stabilization/SmallCell2D.cpp
    src/stabilization/Agglomeration2D.cpp
    src/quality/Quality2D.cpp
    src/quality/SolverQuality2D.cpp
    src/quality/PatchLocalQuality2D.cpp
    src/quality/SolverTopology2D.cpp
    src/io/MeshIO2D.cpp
    src/io/BackgroundGridIO2D.cpp
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

if(MSVC)
    foreach(target
            cartmesh2d_cli
            cartmesh2d_dxf_cli
            cartmesh2d_boundary_layer_cli
            cartmesh2d_hybrid_cli)
        target_sources(${target} PRIVATE apps/cartmesh2d.exe.manifest)
    endforeach()
endif()

# Native two-dimensional finite-volume foundation, independently executable.
add_library(cartmesh2d_fv src/fv/FvMesh2D.cpp src/fv/Diffusion2D.cpp src/fv/Incompressible2D.cpp src/fv/Euler2D.cpp src/fv/FlowBoundaryIO2D.cpp src/fv/ScalarTransport2D.cpp src/fv/ThermalFlow2D.cpp src/fv/Sst2003m2D.cpp src/fv/WallDistance2D.cpp src/fv/SstRans2D.cpp)
target_link_libraries(cartmesh2d_fv PUBLIC cartmesh2d::cartmesh2d)
target_compile_features(cartmesh2d_fv PUBLIC cxx_std_20)
add_executable(cartmesh2d_fv_cli apps/cartmesh2d_fv_cli.cpp)
target_link_libraries(cartmesh2d_fv_cli PRIVATE cartmesh2d_fv)
if(MSVC)
    target_sources(cartmesh2d_fv_cli PRIVATE apps/cartmesh2d.exe.manifest)
endif()
if(COMMAND cartmesh_enable_warnings)
    cartmesh_enable_warnings(cartmesh2d_fv)
    cartmesh_enable_warnings(cartmesh2d_fv_cli)
elseif(MSVC)
    target_compile_options(cartmesh2d_fv PRIVATE /W4 /permissive-)
    target_compile_options(cartmesh2d_fv_cli PRIVATE /W4 /permissive-)
else()
    target_compile_options(cartmesh2d_fv PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    target_compile_options(cartmesh2d_fv_cli PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

add_executable(cartmesh2d_flow_cli apps/cartmesh2d_flow_cli.cpp)
add_executable(cartmesh2d_euler_cli apps/cartmesh2d_euler_cli.cpp)
target_link_libraries(cartmesh2d_euler_cli PRIVATE cartmesh2d_fv)
if(COMMAND cartmesh_enable_warnings)
    cartmesh_enable_warnings(cartmesh2d_euler_cli)
elseif(MSVC)
    target_compile_options(cartmesh2d_euler_cli PRIVATE /W4 /permissive-)
else()
    target_compile_options(cartmesh2d_euler_cli PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()
target_link_libraries(cartmesh2d_flow_cli PRIVATE cartmesh2d_fv)
if(MSVC)
    target_sources(cartmesh2d_flow_cli PRIVATE apps/cartmesh2d.exe.manifest)
endif()
if(COMMAND cartmesh_enable_warnings)
    cartmesh_enable_warnings(cartmesh2d_flow_cli)
elseif(MSVC)
    target_compile_options(cartmesh2d_flow_cli PRIVATE /W4 /permissive-)
else()
    target_compile_options(cartmesh2d_flow_cli PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

add_executable(cartmesh2d_transport_cli apps/cartmesh2d_transport_cli.cpp)
target_link_libraries(cartmesh2d_transport_cli PRIVATE cartmesh2d_fv)
if(MSVC)
    target_sources(cartmesh2d_transport_cli PRIVATE apps/cartmesh2d.exe.manifest)
endif()
if(COMMAND cartmesh_enable_warnings)
    cartmesh_enable_warnings(cartmesh2d_transport_cli)
elseif(MSVC)
    target_compile_options(cartmesh2d_transport_cli PRIVATE /W4 /permissive-)
else()
    target_compile_options(cartmesh2d_transport_cli PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

option(CARTMESH2D_BUILD_TESTS "Build cartmesh2d tests" ON)

if(CARTMESH2D_BUILD_TESTS)
    include(CTest)
    add_executable(cartmesh2d_rectilinear_tests tests/rectilinear_test.cpp)
    target_link_libraries(cartmesh2d_rectilinear_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_rectilinear COMMAND cartmesh2d_rectilinear_tests)
    add_executable(cartmesh2d_rectilinear_probe tests/rectilinear_probe.cpp)
    target_link_libraries(cartmesh2d_rectilinear_probe PRIVATE cartmesh2d_fv)
    add_executable(cartmesh2d_sst_rans_tests tests/sst_rans_test.cpp)
    target_link_libraries(cartmesh2d_sst_rans_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_sst_rans COMMAND cartmesh2d_sst_rans_tests)
    add_executable(cartmesh2d_sst_rans_probe tests/sst_rans_probe.cpp)
    target_link_libraries(cartmesh2d_sst_rans_probe PRIVATE cartmesh2d_fv)
    add_executable(cartmesh2d_sst_nonlinear_unit_tests tests/sst_nonlinear_unit_test.cpp)
    target_link_libraries(cartmesh2d_sst_nonlinear_unit_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_sst_nonlinear_unit COMMAND cartmesh2d_sst_nonlinear_unit_tests)
    add_executable(cartmesh2d_sst_nonlinear_probe tests/sst_nonlinear_probe.cpp)
    target_link_libraries(cartmesh2d_sst_nonlinear_probe PRIVATE cartmesh2d_fv)
    add_executable(cartmesh2d_wall_distance_tests tests/wall_distance_test.cpp)
    target_link_libraries(cartmesh2d_wall_distance_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_wall_distance COMMAND cartmesh2d_wall_distance_tests)
    add_executable(cartmesh2d_sst_spatial_tests tests/sst_spatial_test.cpp)
    target_link_libraries(cartmesh2d_sst_spatial_tests PRIVATE cartmesh2d_fv)
    add_executable(cartmesh2d_sst_transport_tests tests/sst_transport_test.cpp)
    target_link_libraries(cartmesh2d_sst_transport_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_sst_transport COMMAND cartmesh2d_sst_transport_tests)
    add_executable(cartmesh2d_thermal_flow_tests tests/thermal_flow_test.cpp)
    target_link_libraries(cartmesh2d_thermal_flow_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_thermal_flow COMMAND cartmesh2d_thermal_flow_tests)
    add_executable(cartmesh2d_scalar_transport_tests tests/scalar_transport_test.cpp)
    target_link_libraries(cartmesh2d_scalar_transport_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_scalar_transport COMMAND cartmesh2d_scalar_transport_tests)
    add_executable(cartmesh2d_scalar_workspace_tests tests/scalar_workspace_test.cpp)
    target_link_libraries(cartmesh2d_scalar_workspace_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_scalar_workspace COMMAND cartmesh2d_scalar_workspace_tests)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME cartmesh2d_background_grid
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/test_background_grid.py
                $<TARGET_FILE:cartmesh2d_cli>)
    set_tests_properties(cartmesh2d_background_grid PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_sst_flatplate_analysis COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/sst_flatplate_analysis_test.py)
    add_test(NAME cartmesh2d_rectilinear_refinement_verifier COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/rectilinear_refinement_verifier_test.py)
    add_test(NAME cartmesh2d_sst_rans_verifier COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/sst_rans_verifier_test.py
        --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
        --probe $<TARGET_FILE:cartmesh2d_sst_rans_probe>
        --rect-probe $<TARGET_FILE:cartmesh2d_rectilinear_probe>)
    add_test(NAME cartmesh2d_sst_nonlinear_verifier COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/sst_nonlinear_verifier_test.py
        --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
        --probe $<TARGET_FILE:cartmesh2d_sst_nonlinear_probe>)
    add_test(NAME cartmesh2d_sst_spatial_verifier COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/sst_spatial_verifier_test.py
        --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
        --probe $<TARGET_FILE:cartmesh2d_sst_spatial_tests>
        --wall-probe $<TARGET_FILE:cartmesh2d_wall_distance_tests>)
    add_test(NAME cartmesh2d_sst_decay_verifier COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/sst_decay_verifier_test.py
        --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
        --probe $<TARGET_FILE:cartmesh2d_sst_transport_tests>)
    add_test(NAME cartmesh2d_thermal_flow_cli COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/thermal_flow_cli_test.py
        --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
        --transport-cli $<TARGET_FILE:cartmesh2d_transport_cli>)
    add_test(NAME cartmesh2d_thermal_scale_verifier COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/thermal_scale_verifier_test.py
        --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
        --flow-cli $<TARGET_FILE:cartmesh2d_flow_cli>
        --transport-cli $<TARGET_FILE:cartmesh2d_transport_cli>)
    set_tests_properties(cartmesh2d_thermal_scale_verifier PROPERTIES TIMEOUT 60)
    add_test(NAME cartmesh2d_thermal_time_verifier COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/thermal_time_verifier_test.py)
    add_test(NAME cartmesh2d_scalar_transport_cli COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/scalar_transport_cli_test.py
        --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
        --flow-cli $<TARGET_FILE:cartmesh2d_flow_cli>
        --transport-cli $<TARGET_FILE:cartmesh2d_transport_cli>)
    add_test(NAME cartmesh2d_variable_viscosity_cli COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/variable_viscosity_cli_test.py
        --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
        --flow-cli $<TARGET_FILE:cartmesh2d_flow_cli>
        --transport-cli $<TARGET_FILE:cartmesh2d_transport_cli>)
    set_tests_properties(cartmesh2d_variable_viscosity_cli PROPERTIES TIMEOUT 60)
    add_executable(cartmesh2d_fv_tests tests/fv_test.cpp)
    add_executable(cartmesh2d_euler_tests tests/euler_test.cpp)
    target_link_libraries(cartmesh2d_euler_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_euler_core COMMAND cartmesh2d_euler_tests)
    set_tests_properties(cartmesh2d_euler_core PROPERTIES TIMEOUT 60)
    add_test(NAME cartmesh2d_euler_cli COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/euler_cli_test.py --cli $<TARGET_FILE:cartmesh2d_euler_cli>)
    set_tests_properties(cartmesh2d_euler_cli PROPERTIES TIMEOUT 120)
    target_link_libraries(cartmesh2d_fv_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_fv_foundation COMMAND cartmesh2d_fv_tests)
    add_executable(cartmesh2d_flow_linear_tests tests/flow_linear_test.cpp)
    target_link_libraries(cartmesh2d_flow_linear_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_flow_linear COMMAND cartmesh2d_flow_linear_tests)
    add_executable(cartmesh2d_flow_aggregation_tests tests/flow_aggregation_test.cpp)
    target_link_libraries(cartmesh2d_flow_aggregation_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_flow_aggregation COMMAND cartmesh2d_flow_aggregation_tests)
    add_executable(cartmesh2d_flow_face_tests tests/flow_face_test.cpp)
    target_link_libraries(cartmesh2d_flow_face_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_flow_face COMMAND cartmesh2d_flow_face_tests)
    add_executable(cartmesh2d_flow_boundary_tests tests/flow_boundary_test.cpp)
    target_link_libraries(cartmesh2d_flow_boundary_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_flow_boundary COMMAND cartmesh2d_flow_boundary_tests)
    set_tests_properties(cartmesh2d_flow_boundary PROPERTIES TIMEOUT 60)
    add_test(NAME cartmesh2d_transient_flow_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/transient_flow_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    add_test(NAME cartmesh2d_transient_verifier
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/transient_verifier_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    add_test(NAME cartmesh2d_transient_step_comparison
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/transient_step_comparison_test.py)
    set_tests_properties(cartmesh2d_transient_step_comparison PROPERTIES TIMEOUT 30)
    set_tests_properties(cartmesh2d_transient_flow_cli cartmesh2d_transient_verifier PROPERTIES TIMEOUT 90)
    add_executable(cartmesh2d_taylor_green_tests tests/taylor_green_test.cpp)
    target_link_libraries(cartmesh2d_taylor_green_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_taylor_green COMMAND cartmesh2d_taylor_green_tests)
    add_executable(cartmesh2d_flow_checkpoint_tests tests/flow_checkpoint_test.cpp)
    target_link_libraries(cartmesh2d_flow_checkpoint_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_flow_checkpoint COMMAND cartmesh2d_flow_checkpoint_tests)
    add_executable(cartmesh2d_flow_initialization_tests tests/flow_initialization_test.cpp)
    target_link_libraries(cartmesh2d_flow_initialization_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_flow_initialization COMMAND cartmesh2d_flow_initialization_tests)
    add_executable(cartmesh2d_manufactured_flow_tests tests/manufactured_flow_test.cpp)
    target_link_libraries(cartmesh2d_manufactured_flow_tests PRIVATE cartmesh2d_fv)
    add_test(NAME cartmesh2d_manufactured_flow COMMAND cartmesh2d_manufactured_flow_tests)
    add_test(NAME cartmesh2d_fv_cli_validation
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/fv_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_fv_cli>)
    add_test(NAME cartmesh2d_flow_cli_validation
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/flow_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>
            --mesh-cli $<TARGET_FILE:cartmesh2d_cli>)
    set_tests_properties(cartmesh2d_flow_cli_validation PROPERTIES TIMEOUT 120)
    add_test(NAME cartmesh2d_duct_flow_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/duct_flow_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>
            --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
            --transport-cli $<TARGET_FILE:cartmesh2d_transport_cli>)
    set_tests_properties(cartmesh2d_duct_flow_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_frame_limited_duct_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/duct_flow_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>
            --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
            --transport-cli $<TARGET_FILE:cartmesh2d_transport_cli>
            --convection face-limited-linear)
    set_tests_properties(cartmesh2d_frame_limited_duct_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_custom_flow_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/custom_flow_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>
            --mesh-cli $<TARGET_FILE:cartmesh2d_cli>)
    set_tests_properties(cartmesh2d_custom_flow_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_frame_limited_flow_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/custom_flow_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>
            --mesh-cli $<TARGET_FILE:cartmesh2d_cli>
            --convection face-limited-linear)
    set_tests_properties(cartmesh2d_frame_limited_flow_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_symmetry_flow_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/symmetry_flow_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    set_tests_properties(cartmesh2d_symmetry_flow_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_steady_acceleration_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/steady_acceleration_cli_test.py
        --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    set_tests_properties(cartmesh2d_steady_acceleration_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_adaptive_linear_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/adaptive_linear_cli_test.py
        --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    set_tests_properties(cartmesh2d_adaptive_linear_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_wall_split_batch_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/wall_split_batch_cli_test.py
        --cli $<TARGET_FILE:cartmesh2d_cli>)
    set_tests_properties(cartmesh2d_wall_split_batch_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_transient_linear_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/transient_linear_cli_test.py
        --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    add_test(NAME cartmesh2d_adaptive_time_linear_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/adaptive_flow_cli_test.py
        --cli $<TARGET_FILE:cartmesh2d_flow_cli> --linear-policy adaptive)
    set_tests_properties(cartmesh2d_transient_linear_cli cartmesh2d_adaptive_time_linear_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_steady_initial_guess_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/steady_initial_guess_cli_test.py
        --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    set_tests_properties(cartmesh2d_steady_initial_guess_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_pressure_opening_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/pressure_opening_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    set_tests_properties(cartmesh2d_pressure_opening_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_pressure_stencil_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/pressure_stencil_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>
            --mesh-cli $<TARGET_FILE:cartmesh2d_cli>)
    set_tests_properties(cartmesh2d_pressure_stencil_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_adaptive_flow_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/adaptive_flow_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    set_tests_properties(cartmesh2d_adaptive_flow_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_flow_initialization_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/flow_initialization_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    set_tests_properties(cartmesh2d_flow_initialization_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_flow_trajectory
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/flow_trajectory_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli> --mesh-cli $<TARGET_FILE:cartmesh2d_cli>)
    set_tests_properties(cartmesh2d_flow_trajectory PROPERTIES TIMEOUT 60)
    add_test(NAME cartmesh2d_rotating_annulus_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/rotating_annulus_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    set_tests_properties(cartmesh2d_rotating_annulus_cli PROPERTIES TIMEOUT 90)
    add_test(NAME cartmesh2d_flow_verifier
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/flow_verifier_test.py
            --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
    set_tests_properties(cartmesh2d_flow_verifier PROPERTIES TIMEOUT 30)
    add_test(NAME cartmesh2d_cavity_reference
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/cavity_reference_test.py)
    set_tests_properties(cartmesh2d_cavity_reference PROPERTIES TIMEOUT 30)


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
    add_executable(cartmesh2d_boundary_layer_resolution_tests tests/boundary_layer_resolution_test.cpp)
    target_link_libraries(cartmesh2d_boundary_layer_resolution_tests PRIVATE cartmesh2d::cartmesh2d)
    add_test(NAME cartmesh2d_boundary_layer_resolution_tests COMMAND cartmesh2d_boundary_layer_resolution_tests)
    add_test(NAME cartmesh2d_relative_resolution_cli
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/mesh_resolution_cli_test.py
            --cli $<TARGET_FILE:cartmesh2d_cli>
            --hybrid-cli $<TARGET_FILE:cartmesh2d_hybrid_cli>
            --log ${CMAKE_CURRENT_BINARY_DIR}/relative-resolution/cli.log
            --output-dir ${CMAKE_CURRENT_BINARY_DIR}/relative-resolution)
    set_tests_properties(cartmesh2d_relative_resolution_cli PROPERTIES TIMEOUT 150)
    add_test(NAME cartmesh2d_layer_reader_stopped_columns
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/layer_reader_test.py)
    set_tests_properties(cartmesh2d_layer_reader_stopped_columns PROPERTIES TIMEOUT 15)
    add_test(NAME cartmesh2d_directional_connectivity_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/directional_connectivity_test.py)
    set_tests_properties(cartmesh2d_directional_connectivity_reader PROPERTIES TIMEOUT 15)
    add_test(NAME cartmesh2d_face_plane_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/face_plane_audit_test.py)
    set_tests_properties(cartmesh2d_face_plane_reader PROPERTIES TIMEOUT 15)
    add_test(NAME cartmesh2d_nozzle_flow_verification
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/nozzle_flow_verification_test.py)
    set_tests_properties(cartmesh2d_nozzle_flow_verification PROPERTIES TIMEOUT 15)
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

    add_test(NAME cartmesh2d_interior_hybrid_nozzle_e2e
        COMMAND cartmesh2d_hybrid_cli
                ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/nozzle_profile.xy
                ${CMAKE_CURRENT_BINARY_DIR}/interior-hybrid/nozzle
                7 6 7 3 0.015 1.15 0.5
                ${CMAKE_CURRENT_BINARY_DIR}/interior-hybrid/nozzle-foam 0.06
                --fluid-region=interior --small-alpha=0.65)
    set_tests_properties(cartmesh2d_interior_hybrid_nozzle_e2e PROPERTIES
        PASS_REGULAR_EXPRESSION "hybrid_status=success" TIMEOUT 60)
    add_test(NAME cartmesh2d_interior_hybrid_nozzle_reader
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
                ${CMAKE_CURRENT_BINARY_DIR}/interior-hybrid/nozzle-foam
                --require-patch wall_0:wall)
    set_tests_properties(cartmesh2d_interior_hybrid_nozzle_reader PROPERTIES
        DEPENDS cartmesh2d_interior_hybrid_nozzle_e2e)

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
    set(CARTMESH2D_CURVE_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/curve_acceptance)
    file(MAKE_DIRECTORY ${CARTMESH2D_CURVE_OUTPUT})
    add_test(NAME cartmesh2d_curve_dense_cutcell
        COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/serpentine_body.xy
            ${CARTMESH2D_CURVE_OUTPUT}/pure 8 0.25 0.45 exterior
            ${CARTMESH2D_CURVE_OUTPUT}/pure-case 0 0 --size-field
            --far-field-spans 8 --wall-cells-per-span 48 --cells-per-level 3
            --far-level 0 --max-safe-wall-level 11)
    add_test(NAME cartmesh2d_curve_dense_cutcell_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_openfoam2d.py
            ${CARTMESH2D_CURVE_OUTPUT}/pure-case)
    set_tests_properties(cartmesh2d_curve_dense_cutcell_reader PROPERTIES DEPENDS cartmesh2d_curve_dense_cutcell)
    add_test(NAME cartmesh2d_curve_dense_hybrid
        COMMAND cartmesh2d_hybrid_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/serpentine_body.xy
            ${CARTMESH2D_CURVE_OUTPUT}/hybrid 7 3 7 4 0.039 1.2 3
            ${CARTMESH2D_CURVE_OUTPUT}/hybrid-case 0.06)
    set_tests_properties(cartmesh2d_curve_dense_hybrid PROPERTIES
        PASS_REGULAR_EXPRESSION "hybrid_status=success"
        FAIL_REGULAR_EXPRESSION "mesh_mode=pure_cutcell_fallback")
    add_test(NAME cartmesh2d_curve_dense_hybrid_reader
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/verification/check_hybrid_mesh2d.py
            ${CARTMESH2D_CURVE_OUTPUT}/hybrid.hybrid.vtk
            ${CARTMESH2D_CURVE_OUTPUT}/hybrid.hybrid.json)
    set_tests_properties(cartmesh2d_curve_dense_hybrid_reader PROPERTIES DEPENDS cartmesh2d_curve_dense_hybrid)
endif()
