cmake_minimum_required(VERSION 3.20)

# This subproject intentionally owns a separate namespace and source tree.
# It must not link against cartmesh::cartmesh unless a future design review
# explicitly approves a stable, dimension-independent utility dependency.

add_library(cartmesh2d INTERFACE)
add_library(cartmesh2d::cartmesh2d ALIAS cartmesh2d)

target_include_directories(cartmesh2d INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(cartmesh2d INTERFACE cxx_std_20)

option(CARTMESH2D_BUILD_TESTS "Build cartmesh2d tests when implementations exist" OFF)

# Algorithm targets are intentionally not added in the architecture-baseline
# commit. 2D-0 should introduce the first compiled sources and tests only after
# the user explicitly starts Stage 2D-0.
