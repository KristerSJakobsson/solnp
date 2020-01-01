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

pybind11_add_module(pysolnp pysolnp.cpp)
