# Updated for manylinux2014 compatibility
cmake_minimum_required(VERSION 3.14.0)
project(solnp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use FetchContent to manage dependencies
include(FetchContent)

# Fetch dlib v20.0
FetchContent_Declare(
        dlib
        GIT_REPOSITORY https://github.com/davisking/dlib.git
        GIT_TAG v20.0
)

# Fetch Catch2 for testing
FetchContent_Declare(
        catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0
)

# Fetch pybind11 for Python bindings
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.13.6
)

FetchContent_MakeAvailable(dlib)
set_target_properties(dlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

IF(${BUILD_PYSOLNP})
    set(CMAKE_BUILD_TYPE Release)
    message("Building pysolnp")

    file(GLOB_RECURSE PYSOLNP_SOURCES "python_solnp/*.cpp")

    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/python_solnp)

    # Make pybind11 available when building Python bindings
    FetchContent_MakeAvailable(pybind11)

    set(PYBIND11_CPP_STANDARD -std=c++14) # Updated to C++14
    pybind11_add_module(pysolnp "${PYSOLNP_SOURCES}")
    target_link_libraries(pysolnp PRIVATE dlib::dlib)
ELSE()
    file(GLOB_RECURSE sources_test test/*.cpp)
    file(GLOB_RECURSE sources_test test/*.hpp)

    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    # Make Catch2 available for testing
    FetchContent_MakeAvailable(catch2)

    include_directories(src)

    add_library(solnp_lib STATIC src/solnp.hpp src/stdafx.h src/subnp.hpp src/utils.hpp)
    set_target_properties(solnp_lib PROPERTIES LINKER_LANGUAGE CXX)

    add_library(catch2_main STATIC test/test_main.cpp)
    add_executable(solnp_tests test/test_solnp.cpp)
    add_executable(utils_tests test/test_utils.cpp)

    # Link with Catch2::Catch2 instead of including the header directly
    target_link_libraries(
            catch2_main
            PUBLIC Catch2::Catch2WithMain
    )

    target_link_libraries(
            solnp_tests
            PUBLIC catch2_main
            PUBLIC solnp_lib
            PRIVATE dlib::dlib)
    target_link_libraries(
            utils_tests
            PUBLIC catch2_main
            PUBLIC solnp_lib
            PRIVATE dlib::dlib)

    IF(${RUN_CODECOV})
        # --- Coverage specific configuration --- start
        # required: Coverage or Debug
        # (alternatively via command line option: cmake -DCMAKE_BUILD_TYPE=Coverage ..)
        set(CMAKE_BUILD_TYPE Coverage)
        # add the project specific cmake modules directory to the cmake module path
        # (w.r.t. CI/CD builds this is not the best approach)
        set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
        # include the cmake code coverage module
        include(CodeCoverage)
        # test coverage build configuration for C++: debug build, no optimization, profiling enabled
        set(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
        # exclude /test and external libraries from coverage analysis
        set(LCOV_REMOVE_EXTRA "'test/*'" "'${catch2_SOURCE_DIR}/*'" "'${dlib_SOURCE_DIR}/*'")
        # cmake code coverage module target
        setup_target_for_coverage(solnp_tests_coverage solnp_tests coverage)
        setup_target_for_coverage(utils_tests_coverage utils_tests coverage)
        # --- Coverage specific configuration --- end
    ENDIF()

ENDIF()