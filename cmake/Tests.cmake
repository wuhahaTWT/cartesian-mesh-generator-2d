find_package(Python3 REQUIRED COMPONENTS Interpreter)

function(cartmesh2d_add_unit_test target test_name source link_lib)
    cmake_parse_arguments(ARG "WARNINGS" "TIMEOUT" "INCLUDE_DIRS" ${ARGN})
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE ${link_lib})
    target_compile_features(${target} PRIVATE cxx_std_20)
    if(ARG_INCLUDE_DIRS)
        target_include_directories(${target} PRIVATE ${ARG_INCLUDE_DIRS})
    endif()
    if(ARG_WARNINGS)
        cartmesh2d_enable_warnings(${target})
    endif()
    add_test(NAME ${test_name} COMMAND ${target})
    if(ARG_TIMEOUT)
        set_tests_properties(${test_name} PROPERTIES TIMEOUT ${ARG_TIMEOUT})
    endif()
endfunction()

# ---- Mesh core unit tests ----
cartmesh2d_add_unit_test(cartmesh2d_geometry_tests cartmesh2d_stage0_geometry_tests tests/geometry_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_intersection_registry_tests cartmesh2d_stageQ2_intersection_registry_tests tests/intersection_registry_test.cpp cartmesh2d::cartmesh2d)
cartmesh2d_add_unit_test(cartmesh2d_shared_construction_tests cartmesh2d_stageQ2A_shared_construction_tests tests/shared_construction_test.cpp cartmesh2d::cartmesh2d)
cartmesh2d_add_unit_test(cartmesh2d_grid_tests cartmesh2d_stage1_grid_tests tests/grid_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_quadtree_tests cartmesh2d_stage2_quadtree_tests tests/quadtree_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_spatial_index_tests cartmesh2d_stageP1_spatial_index_tests tests/spatial_index_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_size_field_tests cartmesh2d_stageWB_size_field_tests tests/size_field_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_boundary_layer_resolution_tests cartmesh2d_boundary_layer_resolution_tests tests/boundary_layer_resolution_test.cpp cartmesh2d::cartmesh2d)
cartmesh2d_add_unit_test(cartmesh2d_multiloop_tests cartmesh2d_stageD1_multiloop_tests tests/multiloop_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_solver_export_tests cartmesh2d_stageS1_solver_export_tests tests/solver_export_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_shape_agglomeration_tests cartmesh2d_shape_agglomeration_tests tests/shape_agglomeration_test.cpp cartmesh2d::cartmesh2d)
cartmesh2d_add_unit_test(cartmesh2d_cutcell_tests cartmesh2d_stage3_cutcell_tests tests/cutcell_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_cutcell_components_tests cartmesh2d_stage3_cutcell_components_tests tests/cutcell_components_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_cutcell_quadtree_audit_tests cartmesh2d_stage3_quadtree_cutcell_audit_tests tests/cutcell_quadtree_audit_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_topology_tests cartmesh2d_stage4_topology_tests tests/topology_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_topology_quadtree_audit_tests cartmesh2d_stage4_adaptive_topology_audit_tests tests/topology_quadtree_audit_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_small_cell_tests cartmesh2d_stage5a_small_cell_tests tests/small_cell_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_agglomeration_tests cartmesh2d_stage5b_agglomeration_tests tests/agglomeration_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_quality_tests cartmesh2d_stage6a_quality_tests tests/quality_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_io_tests cartmesh2d_stage6b_io_tests tests/io_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_fluent_msh_tests cartmesh2d_fluent_msh_export tests/fluent_msh_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_dxf_tests cartmesh2d_stageDXF2_import_tests tests/dxf_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_boundary_layer_tests cartmesh2d_stageH4_1_boundary_layer_tests tests/boundary_layer_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_hybrid_mesh_tests cartmesh2d_stageH4_2_conformal_hybrid_tests tests/hybrid_mesh_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_h4_local_termination_tests cartmesh2d_stageH4_3_local_termination_tests tests/h4_local_termination_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_weld_budget_tests cartmesh2d_stageR2_weld_budget_tests tests/weld_budget_test.cpp cartmesh2d::cartmesh2d WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_refinement_attribution_tests cartmesh2d_stageR2_refinement_attribution_tests tests/refinement_attribution_test.cpp cartmesh2d::cartmesh2d WARNINGS)

# ---- Built-in laminar flow & thermal transport unit tests ----
cartmesh2d_add_unit_test(cartmesh2d_rectilinear_tests cartmesh2d_rectilinear tests/rectilinear_test.cpp cartmesh2d_fv)
add_executable(cartmesh2d_rectilinear_probe tests/rectilinear_probe.cpp)
target_link_libraries(cartmesh2d_rectilinear_probe PRIVATE cartmesh2d_fv)

cartmesh2d_add_unit_test(cartmesh2d_fv_tests cartmesh2d_fv_foundation tests/fv_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_flow_adaptive_sizing_tests cartmesh2d_flow_adaptive_sizing tests/flow_adaptive_sizing_test.cpp cartmesh2d_fv WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_thermal_flow_tests cartmesh2d_thermal_flow tests/thermal_flow_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_scalar_transport_tests cartmesh2d_scalar_transport tests/scalar_transport_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_scalar_workspace_tests cartmesh2d_scalar_workspace tests/scalar_workspace_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_flow_efficiency_tests cartmesh2d_flow_efficiency tests/flow_efficiency_test.cpp cartmesh2d_fv WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_flow_simple_extensions_tests cartmesh2d_flow_simple_extensions tests/flow_simple_extensions_test.cpp cartmesh2d_fv WARNINGS TIMEOUT 60)
cartmesh2d_add_unit_test(cartmesh2d_flow_coupled_tests cartmesh2d_flow_coupled tests/flow_coupled_test.cpp cartmesh2d_fv WARNINGS TIMEOUT 60)
cartmesh2d_add_unit_test(cartmesh2d_flow_equation_tests cartmesh2d_flow_equation tests/flow_equation_test.cpp cartmesh2d_fv WARNINGS INCLUDE_DIRS src/fv)
cartmesh2d_add_unit_test(cartmesh2d_transport_precision_tests cartmesh2d_transport_precision tests/transport_precision_test.cpp cartmesh2d_fv WARNINGS TIMEOUT 120)
cartmesh2d_add_unit_test(cartmesh2d_wall_gradient_tests cartmesh2d_wall_gradient tests/wall_gradient_test.cpp cartmesh2d_fv WARNINGS)
cartmesh2d_add_unit_test(cartmesh2d_viscous_tests cartmesh2d_viscous_core tests/viscous_stress_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_flow_linear_tests cartmesh2d_flow_linear tests/flow_linear_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_flow_aggregation_tests cartmesh2d_flow_aggregation tests/flow_aggregation_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_flow_face_tests cartmesh2d_flow_face tests/flow_face_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_flow_cholesky_tests cartmesh2d_flow_cholesky tests/flow_cholesky_test.cpp cartmesh2d_fv TIMEOUT 30)
cartmesh2d_add_unit_test(cartmesh2d_flow_boundary_tests cartmesh2d_flow_boundary tests/flow_boundary_test.cpp cartmesh2d_fv TIMEOUT 60)
cartmesh2d_add_unit_test(cartmesh2d_taylor_green_tests cartmesh2d_taylor_green tests/taylor_green_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_flow_checkpoint_tests cartmesh2d_flow_checkpoint tests/flow_checkpoint_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_flow_initialization_tests cartmesh2d_flow_initialization tests/flow_initialization_test.cpp cartmesh2d_fv)
cartmesh2d_add_unit_test(cartmesh2d_manufactured_flow_tests cartmesh2d_manufactured_flow tests/manufactured_flow_test.cpp cartmesh2d_fv)

add_executable(cartmesh2d_flow_stability_tests tests/flow_stability_test.cpp)
target_link_libraries(cartmesh2d_flow_stability_tests PRIVATE cartmesh2d_fv)
set(laminar_fixture_dir "${CMAKE_CURRENT_BINARY_DIR}/laminar-fixtures")
file(MAKE_DIRECTORY "${laminar_fixture_dir}")
add_test(NAME cartmesh2d_laminar_fixtures COMMAND ${CMAKE_COMMAND} -E tar xzf
    "${CMAKE_CURRENT_SOURCE_DIR}/artifacts/current/laminar-foundation/failure-inputs.tar.gz")
set_tests_properties(cartmesh2d_laminar_fixtures PROPERTIES
    WORKING_DIRECTORY "${laminar_fixture_dir}" FIXTURES_SETUP laminar_inputs)
add_test(NAME cartmesh2d_flow_stability COMMAND cartmesh2d_flow_stability_tests "${laminar_fixture_dir}")
set_tests_properties(cartmesh2d_flow_stability PROPERTIES FIXTURES_REQUIRED laminar_inputs TIMEOUT 180)

# ---- CLI & end-to-end mesh/flow tests ----
add_test(NAME cartmesh2d_transient_flow_cli
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/transient_flow_cli_test.py
        --cli $<TARGET_FILE:cartmesh2d_flow_cli>)
set_tests_properties(cartmesh2d_transient_flow_cli PROPERTIES TIMEOUT 90)

add_test(NAME cartmesh2d_fv_cli_validation
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/fv_cli_test.py
        --cli $<TARGET_FILE:cartmesh2d_fv_cli>)
add_test(NAME cartmesh2d_flow_cli_validation
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/flow_cli_test.py
        --cli $<TARGET_FILE:cartmesh2d_flow_cli>
        --mesh-cli $<TARGET_FILE:cartmesh2d_cli>)
set_tests_properties(cartmesh2d_flow_cli_validation PROPERTIES TIMEOUT 120)

add_test(NAME cartmesh2d_product_cutcell_help COMMAND cartmesh2d_cli --help)
add_test(NAME cartmesh2d_product_hybrid_help COMMAND cartmesh2d_hybrid_cli --help)
add_test(NAME cartmesh2d_product_dxf_help COMMAND cartmesh2d_dxf_cli --help)

# Boundary layer e2e & determinism
set(CARTMESH2D_H4_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/boundary_layer)
add_test(NAME cartmesh2d_stageH4_1_circle_e2e
    COMMAND cartmesh2d_boundary_layer_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_H4_OUTPUT}/circle 4 first 0.02 1.2 exterior)
add_test(NAME cartmesh2d_stageH4_1_rectangle_e2e
    COMMAND cartmesh2d_boundary_layer_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/rectangle.xy ${CARTMESH2D_H4_OUTPUT}/rectangle 3 first 0.04 1.25 exterior)
add_test(NAME cartmesh2d_stageH4_1_concave_fails_closed
    COMMAND cartmesh2d_boundary_layer_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/concave.xy ${CARTMESH2D_H4_OUTPUT}/concave 3 first 0.02 1.1 exterior)
set_tests_properties(cartmesh2d_stageH4_1_concave_fails_closed PROPERTIES WILL_FAIL TRUE)
add_test(NAME cartmesh2d_stageH4_1_circle_repeat_e2e
    COMMAND cartmesh2d_boundary_layer_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_H4_OUTPUT}/circle_repeat 4 first 0.02 1.2 exterior)
add_test(NAME cartmesh2d_stageH4_1_circle_vtk_deterministic
    COMMAND ${CMAKE_COMMAND} -E compare_files ${CARTMESH2D_H4_OUTPUT}/circle.layer.vtk ${CARTMESH2D_H4_OUTPUT}/circle_repeat.layer.vtk)
set_tests_properties(cartmesh2d_stageH4_1_circle_vtk_deterministic PROPERTIES DEPENDS "cartmesh2d_stageH4_1_circle_e2e;cartmesh2d_stageH4_1_circle_repeat_e2e")
add_test(NAME cartmesh2d_stageH4_1_circle_report_deterministic
    COMMAND ${CMAKE_COMMAND} -E compare_files ${CARTMESH2D_H4_OUTPUT}/circle.layer.json ${CARTMESH2D_H4_OUTPUT}/circle_repeat.layer.json)
set_tests_properties(cartmesh2d_stageH4_1_circle_report_deterministic PROPERTIES DEPENDS "cartmesh2d_stageH4_1_circle_e2e;cartmesh2d_stageH4_1_circle_repeat_e2e")

# Hybrid mesh e2e & determinism
set(CARTMESH2D_H4_2_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hybrid_mesh)
add_test(NAME cartmesh2d_stageH4_2_circle_e2e
    COMMAND cartmesh2d_hybrid_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_H4_2_OUTPUT}/circle 6 3 6 4 0.02 1.2 1.0)
add_test(NAME cartmesh2d_stageH4_2_superellipse_e2e
    COMMAND cartmesh2d_hybrid_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/superellipse_24.xy ${CARTMESH2D_H4_2_OUTPUT}/superellipse 6 3 6 3 0.015 1.15 1.0)
add_test(NAME cartmesh2d_stageH4_2_circle_repeat_e2e
    COMMAND cartmesh2d_hybrid_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_H4_2_OUTPUT}/circle_repeat 6 3 6 4 0.02 1.2 1.0)
add_test(NAME cartmesh2d_stageH4_2_circle_topology_deterministic
    COMMAND ${CMAKE_COMMAND} -E compare_files ${CARTMESH2D_H4_2_OUTPUT}/circle.hybrid.cm2d ${CARTMESH2D_H4_2_OUTPUT}/circle_repeat.hybrid.cm2d)
set_tests_properties(cartmesh2d_stageH4_2_circle_topology_deterministic PROPERTIES DEPENDS "cartmesh2d_stageH4_2_circle_e2e;cartmesh2d_stageH4_2_circle_repeat_e2e")
add_test(NAME cartmesh2d_stageH4_2_circle_report_deterministic
    COMMAND ${CMAKE_COMMAND} -E compare_files ${CARTMESH2D_H4_2_OUTPUT}/circle.hybrid.json ${CARTMESH2D_H4_2_OUTPUT}/circle_repeat.hybrid.json)
set_tests_properties(cartmesh2d_stageH4_2_circle_report_deterministic PROPERTIES DEPENDS "cartmesh2d_stageH4_2_circle_e2e;cartmesh2d_stageH4_2_circle_repeat_e2e")

add_test(NAME cartmesh2d_interior_hybrid_nozzle_e2e
    COMMAND cartmesh2d_hybrid_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/nozzle_profile.xy
            ${CMAKE_CURRENT_BINARY_DIR}/interior-hybrid/nozzle 7 6 7 3 0.015 1.15 0.5
            ${CMAKE_CURRENT_BINARY_DIR}/interior-hybrid/nozzle-foam 0.06 --fluid-region=interior --small-alpha=0.65)
set_tests_properties(cartmesh2d_interior_hybrid_nozzle_e2e PROPERTIES PASS_REGULAR_EXPRESSION "hybrid_status=success" TIMEOUT 60)

# DXF e2e
set(CARTMESH2D_DXF_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/dxf)
add_test(NAME cartmesh2d_stageDXF1_convert_e2e
    COMMAND cartmesh2d_dxf_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/dxf/rounded_obstacle.dxf ${CARTMESH2D_DXF_OUTPUT}/rounded_obstacle.xy 0.002 ${CARTMESH2D_DXF_OUTPUT}/rounded_obstacle.dxf.json)
add_test(NAME cartmesh2d_stageDXF1_convert_airfoil_e2e
    COMMAND cartmesh2d_dxf_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/dxf/airfoil_like.dxf ${CARTMESH2D_DXF_OUTPUT}/airfoil_like.xy 0.001 ${CARTMESH2D_DXF_OUTPUT}/airfoil_like.dxf.json)
add_test(NAME cartmesh2d_stageDXF1_mesh_e2e
    COMMAND cartmesh2d_cli ${CARTMESH2D_DXF_OUTPUT}/airfoil_like.xy ${CARTMESH2D_DXF_OUTPUT}/airfoil_like 7 0.30 0.02 exterior ${CARTMESH2D_DXF_OUTPUT}/airfoil_like-case 0 0)
set_tests_properties(cartmesh2d_stageDXF1_mesh_e2e PROPERTIES DEPENDS cartmesh2d_stageDXF1_convert_airfoil_e2e)

add_test(NAME cartmesh2d_stageDXF2_convert_spline_mm_e2e
    COMMAND cartmesh2d_dxf_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/dxf/spline_circle_mm.dxf ${CARTMESH2D_DXF_OUTPUT}/spline_circle.xy 0.001 ${CARTMESH2D_DXF_OUTPUT}/spline_circle.dxf.json)
add_test(NAME cartmesh2d_stageDXF2_mesh_spline_mm_e2e
    COMMAND cartmesh2d_cli ${CARTMESH2D_DXF_OUTPUT}/spline_circle.xy ${CARTMESH2D_DXF_OUTPUT}/spline_circle 6 0.25 0.05 exterior ${CARTMESH2D_DXF_OUTPUT}/spline_circle-case 0 0)
set_tests_properties(cartmesh2d_stageDXF2_mesh_spline_mm_e2e PROPERTIES DEPENDS cartmesh2d_stageDXF2_convert_spline_mm_e2e)

add_test(NAME cartmesh2d_stageDXF2_convert_symmetry_e2e
    COMMAND cartmesh2d_dxf_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/dxf/symmetry_circle_m.dxf ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle.xy 0.01 ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle.dxf.json)
add_test(NAME cartmesh2d_stageDXF2_mesh_symmetry_e2e
    COMMAND cartmesh2d_cli ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle.xy ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle 6 0.25 0.05 exterior ${CARTMESH2D_DXF_OUTPUT}/symmetry_circle-case 0 0)
set_tests_properties(cartmesh2d_stageDXF2_mesh_symmetry_e2e PROPERTIES DEPENDS cartmesh2d_stageDXF2_convert_symmetry_e2e)

# Acceptance, sizing field, robustness & OpenFOAM export e2e
set(CARTMESH2D_ACCEPTANCE_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/acceptance)
add_test(NAME cartmesh2d_stage6c_rectangle_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/rectangle.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/rectangle 5 0.25 0.05)
add_test(NAME cartmesh2d_stage6c_circle_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/circle 5 0.25 0.05)
add_test(NAME cartmesh2d_stage6c_concave_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/concave.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/concave 5 0.25 0.05)
add_test(NAME cartmesh2d_stage6c_airfoil_like_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy ${CARTMESH2D_ACCEPTANCE_OUTPUT}/airfoil_like 6 0.30 0.02)

set(CARTMESH2D_W1_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/r2_w1)
add_test(NAME cartmesh2d_stageR2_w1_circle_l10_e2e
    COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_W1_OUTPUT}/circle_l10 10 0.25 0.10 exterior ${CARTMESH2D_W1_OUTPUT}/circle_l10-case 8)

set(CARTMESH2D_SIZING_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/sizing_field)
add_test(NAME cartmesh2d_stageH1_sizing_field_e2e
    COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy ${CARTMESH2D_SIZING_OUTPUT}/airfoil 7 0.30 0.02 exterior - 3 0 --distance-band 0.05 6 --refine-box 0.80 -0.10 1.30 0.10 6)
add_test(NAME cartmesh2d_stageH1_rejects_outside_box
    COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy ${CARTMESH2D_SIZING_OUTPUT}/invalid 6 0.30 0.02 exterior - 2 0 --refine-box 10 10 11 11 5)
set_tests_properties(cartmesh2d_stageH1_rejects_outside_box PROPERTIES WILL_FAIL TRUE)

add_test(NAME cartmesh2d_stageWB_size_field_defaults_e2e
    COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_SIZING_OUTPUT}/sf_circle 8 0.25 0.15 exterior ${CARTMESH2D_SIZING_OUTPUT}/sf_circle-case 0 0 --size-field --cells-per-level 3)
add_test(NAME cartmesh2d_stageWB_size_field_wake_e2e
    COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy ${CARTMESH2D_SIZING_OUTPUT}/sf_airfoil 8 0.25 0.15 exterior ${CARTMESH2D_SIZING_OUTPUT}/sf_airfoil-case 0 0 --size-field --cells-per-level 3 --wake 4 6 0.6 4)
add_test(NAME cartmesh2d_stageWB_rejects_unsafe_wall_level
    COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_SIZING_OUTPUT}/sf_unsafe 8 0.25 0.15 exterior - 0 0 --size-field --wall-cells-per-span 4096)
set_tests_properties(cartmesh2d_stageWB_rejects_unsafe_wall_level PROPERTIES WILL_FAIL TRUE)
add_test(NAME cartmesh2d_stageWB_rejects_mixed_sizing_spellings
    COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_SIZING_OUTPUT}/sf_mixed 8 0.25 0.15 exterior - 0 0 --size-field --distance-band 0.05 6)
set_tests_properties(cartmesh2d_stageWB_rejects_mixed_sizing_spellings PROPERTIES WILL_FAIL TRUE)

set(CARTMESH2D_ROBUSTNESS_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/robustness)
add_test(NAME cartmesh2d_stageR_gear_star_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/gear_star.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/gear_star 8 0.25 0.01)
add_test(NAME cartmesh2d_stageR_serpentine_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/serpentine_body.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/serpentine 7 0.20 0.01)
add_test(NAME cartmesh2d_stageR_nozzle_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/nozzle_profile.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/nozzle 7 0.20 0.01)
add_test(NAME cartmesh2d_stageR_naca2412_dense_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/naca2412_dense.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/naca2412_dense 8 0.30 0.005)
add_test(NAME cartmesh2d_stageR_superellipse_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/superellipse_24.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/superellipse 7 0.25 0.01)
add_test(NAME cartmesh2d_stageD1_two_obstacles_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/two_obstacles.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/two_obstacles 7 0.25 0.05 exterior)
add_test(NAME cartmesh2d_stageD1_annulus_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/annulus.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/annulus 7 0.25 0.05 interior)
add_test(NAME cartmesh2d_stageD1_nested_island_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/nested_island.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/nested_island 7 0.25 0.05 interior)
add_test(NAME cartmesh2d_stageD1_rejects_intersecting_loops COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/invalid/intersecting_loops.xy ${CARTMESH2D_ROBUSTNESS_OUTPUT}/invalid_intersecting_loops 6 0.25 0.05 interior)
set_tests_properties(cartmesh2d_stageD1_rejects_intersecting_loops PROPERTIES WILL_FAIL TRUE)

set(CARTMESH2D_SOLVER_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/solver_ready)
add_test(NAME cartmesh2d_stageS1_rectangle_openfoam_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/rectangle.xy ${CARTMESH2D_SOLVER_OUTPUT}/rectangle 5 0.25 0.05 exterior ${CARTMESH2D_SOLVER_OUTPUT}/rectangle-case)
add_test(NAME cartmesh2d_stageV1_annulus_openfoam_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/annulus.xy ${CARTMESH2D_SOLVER_OUTPUT}/annulus 7 0.25 0.05 interior ${CARTMESH2D_SOLVER_OUTPUT}/annulus-case)
add_test(NAME cartmesh2d_stageV1_nested_openfoam_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/nested_island.xy ${CARTMESH2D_SOLVER_OUTPUT}/nested 7 0.25 0.05 interior ${CARTMESH2D_SOLVER_OUTPUT}/nested-case)
add_test(NAME cartmesh2d_stageV1b_circle_openfoam_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/circle.xy ${CARTMESH2D_SOLVER_OUTPUT}/circle 6 0.25 0.05 exterior ${CARTMESH2D_SOLVER_OUTPUT}/circle-case)
add_test(NAME cartmesh2d_stageV1b_airfoil_openfoam_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/acceptance/airfoil_like.xy ${CARTMESH2D_SOLVER_OUTPUT}/airfoil 7 0.30 0.02 exterior ${CARTMESH2D_SOLVER_OUTPUT}/airfoil-case)
add_test(NAME cartmesh2d_stageV1b_two_obstacles_openfoam_e2e COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/two_obstacles.xy ${CARTMESH2D_SOLVER_OUTPUT}/two_obstacles 8 0.25 0.05 exterior ${CARTMESH2D_SOLVER_OUTPUT}/two_obstacles-case)

set(CARTMESH2D_CURVE_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/curve_acceptance)
file(MAKE_DIRECTORY ${CARTMESH2D_CURVE_OUTPUT})
add_test(NAME cartmesh2d_curve_dense_cutcell
    COMMAND cartmesh2d_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/serpentine_body.xy
            ${CARTMESH2D_CURVE_OUTPUT}/pure 8 0.25 0.45 exterior ${CARTMESH2D_CURVE_OUTPUT}/pure-case 0 0
            --size-field --far-field-spans 8 --wall-cells-per-span 48 --cells-per-level 3 --far-level 0 --max-safe-wall-level 11)
add_test(NAME cartmesh2d_curve_dense_hybrid
    COMMAND cartmesh2d_hybrid_cli ${CMAKE_CURRENT_SOURCE_DIR}/examples/complex/serpentine_body.xy
            ${CARTMESH2D_CURVE_OUTPUT}/hybrid 7 3 7 4 0.039 1.2 3 ${CARTMESH2D_CURVE_OUTPUT}/hybrid-case 0.06)
set_tests_properties(cartmesh2d_curve_dense_hybrid PROPERTIES
    PASS_REGULAR_EXPRESSION "hybrid_status=success"
    FAIL_REGULAR_EXPRESSION "mesh_mode=pure_cutcell_fallback")
