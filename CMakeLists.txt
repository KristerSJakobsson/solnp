# Note: CMake 3.12 is required for manylinux2010 compatibility
cmake_minimum_required(VERSION 3.12)
project(solnp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


IF(${BUILD_PYSOLNP})
    set(CMAKE_BUILD_TYPE Release)
    message("Building pysolnp")

    file(GLOB_RECURSE PYSOLNP_SOURCES "python_solnp/*.cpp")

    # Find Python3 and NumPy
    find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
    
    # Try to find dlib, fallback to subdirectory if not found
    find_package(dlib QUIET)
    if(NOT dlib_FOUND)
        message(STATUS "dlib not found, using subdirectory")
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/library/dlib)
    endif()
    
    # Try to find pybind11, fallback to subdirectory if not found
    find_package(pybind11 QUIET)
    if(NOT pybind11_FOUND)
        message(STATUS "pybind11 not found, using subdirectory")
        add_subdirectory(library/pybind11)
    endif()
    
    pybind11_add_module(pysolnp "${PYSOLNP_SOURCES}")
    target_include_directories(pysolnp PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/python_solnp
        ${Python3_NumPy_INCLUDE_DIRS}
    )
    target_link_libraries(pysolnp PRIVATE dlib::dlib Python3::NumPy)
ELSE()
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
        # exclude /test and /library directories from coverage analysis
        set(LCOV_REMOVE_EXTRA "'test/*'" "'library/*'")
        # cmake code coverage module target
        setup_target_for_coverage(solnp_tests_coverage solnp_tests coverage)
        setup_target_for_coverage(utils_tests_coverage utils_tests coverage)
        # --- Coverage specific configuration --- end
    ENDIF()
    file(GLOB_RECURSE sources_test test/*.cpp)
    file(GLOB_RECURSE sources_test test/*.hpp)

    # Try to find dlib, fallback to subdirectory if not found
    find_package(dlib QUIET)
    if(NOT dlib_FOUND)
        message(STATUS "dlib not found, using subdirectory")
        add_subdirectory(library/dlib)
    endif()

    add_library(solnp_lib STATIC src/solnp.hpp src/stdafx.h src/subnp.hpp src/utils.hpp)
    set_target_properties(solnp_lib PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(solnp_lib PUBLIC src)

    add_library(catch2_main STATIC test/test_main.cpp)
    target_include_directories(catch2_main PUBLIC library/Catch2/single_include)
    
    add_executable(solnp_tests test/test_solnp.cpp)
    target_include_directories(solnp_tests PRIVATE src library/Catch2/single_include)
    
    add_executable(utils_tests test/test_utils.cpp)
    target_include_directories(utils_tests PRIVATE src library/Catch2/single_include)
    
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
