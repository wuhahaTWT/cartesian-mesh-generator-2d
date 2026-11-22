# Native 2D incompressible laminar flow & passive thermal transport
# (serves built-in quick preview and mesh-validation workflows).
add_library(cartmesh2d_fv
    src/fv/FvMesh2D.cpp
    src/fv/Diffusion2D.cpp
    src/fv/Incompressible2D.cpp
    src/fv/FlowCase2D.cpp
    src/fv/FlowMomentum2D.cpp
    src/fv/FlowFlux2D.cpp
    src/fv/FlowForces2D.cpp
    src/fv/FlowConvergence2D.cpp
    src/fv/FlowCoupled2D.cpp
    src/fv/FlowEquation2D.cpp
    src/fv/HeatConduction2D.cpp
    src/fv/ViscousStress2D.cpp
    src/fv/WallGradient2D.cpp
    src/fv/FlowBoundaryIO2D.cpp
    src/fv/FlowAdaptiveSizing2D.cpp
    src/fv/ScalarTransport2D.cpp
    src/fv/ThermalFlow2D.cpp
)
target_link_libraries(cartmesh2d_fv PUBLIC cartmesh2d::cartmesh2d)
if(APPLE)
    target_link_libraries(cartmesh2d_fv PUBLIC "-framework Accelerate")
endif()
target_compile_features(cartmesh2d_fv PUBLIC cxx_std_20)
cartmesh2d_enable_warnings(cartmesh2d_fv)

foreach(cli IN ITEMS
        cartmesh2d_fv_cli
        cartmesh2d_flow_cli
        cartmesh2d_transport_cli)
    add_executable(${cli} apps/${cli}.cpp)
    cartmesh2d_configure_cli(${cli} cartmesh2d_fv)
endforeach()
