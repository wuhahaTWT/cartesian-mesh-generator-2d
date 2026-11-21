# Core native 2D geometry, Cut-cell, conformal boundary layer, hybrid mesh,
# shape agglomeration/repartitioning, quality diagnostics, and mesh I/O.
add_library(cartmesh2d
    src/boundary_layer/BoundaryLayer2D.cpp
    src/boundary_layer/BoundaryLayerIO2D.cpp
    src/geometry/Geometry2D.cpp
    src/geometry/ConstructionIdentity2D.cpp
    src/geometry/ConstructionRecovery2D.cpp
    src/geometry/IntersectionRegistry2D.cpp
    src/geometry/IntersectionConstruction2D.cpp
    src/geometry/BoundarySimplification2D.cpp
    src/grid/CartesianGrid2D.cpp
    src/grid/RectilinearMesh2D.cpp
    src/hybrid/HybridMesh2D.cpp
    src/hybrid/HybridMeshFallback2D.cpp
    src/hybrid/HybridMeshIO2D.cpp
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
    src/quality/SolverTopologyCommon2D.cpp
    src/quality/SolverTopologyQualityRepair2D.cpp
    src/quality/SolverTopologyShortFaceRepair2D.cpp
    src/quality/SolverTopologyTerminationRepair2D.cpp
    src/quality/SolverShapeAgglomeration2D.cpp
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
cartmesh2d_enable_warnings(cartmesh2d)

foreach(cli IN ITEMS
        cartmesh2d_cli
        cartmesh2d_dxf_cli
        cartmesh2d_boundary_layer_cli
        cartmesh2d_hybrid_cli)
    add_executable(${cli} apps/${cli}.cpp)
    cartmesh2d_configure_cli(${cli} cartmesh2d::cartmesh2d)
endforeach()
