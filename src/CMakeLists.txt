cmake_minimum_required(VERSION 2.8.12)

set(CMAKE_CXX_STANDARD 11)

project(pysolnp)

pybind11_add_module(pysolnp pysolver.cpp)
target_link_libraries(pysolnp PRIVATE dlib::dlib)
target_link_libraries(pysolnp PRIVATE pybind11::pybind11)
