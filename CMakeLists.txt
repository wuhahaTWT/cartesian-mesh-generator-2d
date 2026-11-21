cmake_minimum_required(VERSION 3.20)

project(cartmesh2d LANGUAGES CXX)

if(MSVC)
    # Ship self-contained command-line tools and make narrow source/argv text
    # consistently UTF-8 on supported Windows 10 and newer systems.
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    add_compile_options(/utf-8)
    add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN)
endif()

option(CARTMESH2D_BUILD_TESTS "Build cartmesh2d tests" ON)

function(cartmesh2d_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    endif()
endfunction()

function(cartmesh2d_configure_cli target link_lib)
    target_link_libraries(${target} PRIVATE ${link_lib})
    target_compile_features(${target} PRIVATE cxx_std_20)
    cartmesh2d_enable_warnings(${target})
    if(MSVC)
        target_sources(${target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/apps/cartmesh2d.exe.manifest")
    endif()
endfunction()

include(cmake/MeshTargets.cmake)
include(cmake/FvTargets.cmake)

if(CARTMESH2D_BUILD_TESTS)
    include(CTest)
    include(cmake/Tests.cmake)
endif()
