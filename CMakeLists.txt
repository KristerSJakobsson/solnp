# Note: 2.8.12 is needed for compatability with anylinux docker
cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_BUILD_TYPE Release)
project(pysolnp)

file(GLOB_RECURSE PYSOLNP_SOURCES "python_solnp/*.cpp")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/library/dlib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/python_solnp)

set(PYBIND11_CPP_STANDARD -std=c++1y)
add_subdirectory(library/pybind11)
pybind11_add_module(pysolnp "${PYSOLNP_SOURCES}")
target_link_libraries(pysolnp PRIVATE dlib::dlib)

if (ENABLE_COVERAGE)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/library/cmake_codecov")
    list(APPEND LCOV_REMOTE_PATTERNS "${CMAKE_CURRENT_LIST_DIR}/library/*")
    list(APPEND LCOV_REMOTE_PATTERNS "${CMAKE_CURRENT_LIST_DIR}/python_solnp/*")

    find_package(codecov)
    add_coverage(pysolnp)
    list(APPEND LCOV_REMOVE_PATTERNS "/usr/")
    coverage_evaluate()
endif()

