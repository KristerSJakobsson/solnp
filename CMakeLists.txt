cmake_minimum_required(VERSION 2.8.12) # Version 2.8.12 only used when building python package for compatability

project(solnp)

add_subdirectory(library/dlib dlib)
add_subdirectory(library/pybind11 pybind)
add_subdirectory(src)
