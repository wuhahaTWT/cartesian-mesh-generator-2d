cmake_minimum_required(VERSION 3.20)

project(cartmesh2d LANGUAGES CXX)

add_library(cartmesh2d
    src/geometry/Geometry2D.cpp
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
    if(COMMAND cartmesh_enable_warnings)
        cartmesh_enable_warnings(cartmesh2d_geometry_tests)
    elseif(MSVC)
        target_compile_options(cartmesh2d_geometry_tests PRIVATE /W4 /permissive-)
    else()
        target_compile_options(cartmesh2d_geometry_tests PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    endif()
    add_test(NAME cartmesh2d_stage0_geometry_tests COMMAND cartmesh2d_geometry_tests)
endif()
