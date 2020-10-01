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
    set(CMAKE_BUILD_TYPE Debug)
    message("Building solnp")

    add_subdirectory(library/dlib)
    find_package(dlib REQUIRED)
    message("Found dlib")

    add_subdirectory(library/Catch2)
    find_package(Catch2 REQUIRED)
    message("Found Catch2")

    add_subdirectory(src)
    find_package(solnp REQUIRED)
    message("Found solnp")

    add_subdirectory(test)

ENDIF()
