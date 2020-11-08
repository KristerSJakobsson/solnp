# Note: This script builds pysolnp, for building the C++ code, see test/CMakeLists
# Note2: Version 2.8.12 is needed for compatability with anylinux docker
cmake_minimum_required(VERSION 2.8.12)
project(solnp)

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
    # --- Coverage specific configuration --- start
    # required: Coverage or Debug
    # (alternativelly via command line option: cmake -DCMAKE_BUILD_TYPE=Coverage ..)
    set(CMAKE_BUILD_TYPE Coverage)
    # add the project specific cmake modules directory to the cmake module path
    # (w.r.t. CI/CD builds this is not the best approach)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
    # include the cmake code coverage module
    include(CodeCoverage)
    # test coverage build configuration for C++: debug build, no optimization, profiling enabled
    set(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    # exclude /test and /ext directories from coverage analysis
    set(LCOV_REMOVE_EXTRA "'test/*'" "'ext/*'")
    # cmake code coverage module target
    setup_target_for_coverage(unit_tests_coverage unit_tests coverage)
    # --- Coverage specific configuration --- end

    file(GLOB_RECURSE sources_production src/*.cpp)
    file(GLOB_RECURSE sources_test test/*.cpp)
    file(GLOB_RECURSE sources_test test/*.hpp)

    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    add_subdirectory(library/dlib)
    include_directories(src)
    include_directories(library/Catch2/include)

    add_executable(unit_tests ${sources_test} ${sources_production})
    target_link_libraries(unit_tests PRIVATE dlib::dlib)

ENDIF()
