# Note: This script builds pysolnp, for building the C++ code, see test/CMakeLists
# Note2: Version 2.8.12 is needed for compatability with anylinux docker

# Note: This script builds pysolnp, for building the C++ code, see test/CMakeLists
# Note2: Version 2.8.12 is needed for compatability with anylinux docker
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



#cmake_minimum_required(VERSION 2.8.12)
#project(pysolnp)
#
#IF(${BUILD_PYSOLNP})
#    # This flag sets the build to build the python library Pysolnp
#    # This is run through setup.py, please don't run it manually unless you know what you are doing
#    set(CMAKE_BUILD_TYPE Release)
#
#    file(GLOB_RECURSE PYSOLNP_SOURCES "python_solnp/*.cpp")
#
#    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#
#    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/library/dlib)
#    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
#    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/python_solnp)
#
#    set(PYBIND11_CPP_STANDARD -std=c++1y)
#    add_subdirectory(library/pybind11)
#    pybind11_add_module(pysolnp "${PYSOLNP_SOURCES}")
#    target_link_libraries(pysolnp PRIVATE dlib::dlib)
#ELSE()
#    set(CMAKE_BUILD_TYPE Debug)
#    message("Hello world")
#ENDIF()
#
