# Note: This script builds pysolnp, for building the C++ code, see test/CMakeLists
# Note2: Version 2.8.12 is needed for compatability with anylinux docker
cmake_minimum_required(VERSION 2.8.12)
project(solnp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


IF(${BUILD_PYSOLNP})
    set(CMAKE_BUILD_TYPE Release)
    message("Building pysolnp")

    file(GLOB_RECURSE PYSOLNP_SOURCES "python_solnp/*.cpp")

    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/library/dlib)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/python_solnp)

    set(PYBIND11_CPP_STANDARD -std=c++1y)
    add_subdirectory(library/pybind11)
    pybind11_add_module(pysolnp "${PYSOLNP_SOURCES}")
    target_link_libraries(pysolnp PRIVATE dlib::dlib)
ELSE()
    IF(${RUN_CODECOV})
        # --- Coverage specific configuration --- start
        # required: Coverage or Debug
        # (alternativelly via command line option: cmake -DCMAKE_BUILD_TYPE=Coverage ..)
        set(CMAKE_BUILD_TYPE Coverage)
        # add the project specific cmake modules directory to the cmake module path
        # (w.r.t. CI/CD builds this is not the best approach)
        set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
        # include the cmake code coverage module
        include(CodeCoverage)
        # test coverage build configuration for C++: debug build, no optimization, profiling enabled
        set(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
        # exclude /test and /library directories from coverage analysis
        set(LCOV_REMOVE_EXTRA "'test/*'" "'library/*'")
        # cmake code coverage module target
        setup_target_for_coverage(solnp_tests_coverage solnp_tests coverage)
        setup_target_for_coverage(utils_tests_coverage utils_tests coverage)
        # --- Coverage specific configuration --- end
    ENDIF()
    file(GLOB_RECURSE sources_test test/*.cpp)
    file(GLOB_RECURSE sources_test test/*.hpp)

    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    add_subdirectory(library/dlib)
    include_directories(src)
    include_directories(library/Catch2/single_include)

    add_library(solnp_lib STATIC src/solnp.hpp src/stdafx.h src/subnp.hpp src/utils.hpp)
    set_target_properties(solnp_lib PROPERTIES LINKER_LANGUAGE CXX)

    add_library(catch2_main STATIC test/test_main.cpp)
    add_executable(solnp_tests test/test_solnp.cpp)
    add_executable(utils_tests test/test_utils.cpp)
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

ENDIF()
