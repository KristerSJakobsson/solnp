cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 11)

project(pysolnp)
include(FetchContent)

FetchContent_Declare(dlib
        GIT_REPOSITORY https://github.com/davisking/dlib.git
        GIT_TAG v19.18)
FetchContent_MakeAvailable(dlib)

FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.4.3)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(pysolnp src/pysolver.cpp)
target_link_libraries(pysolnp PRIVATE dlib::dlib)
target_link_libraries(pysolnp PRIVATE pybind11::pybind11)
add_subdirectory(test)