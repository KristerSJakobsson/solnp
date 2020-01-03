project(cppsolnp)

set(CPACK_PACKAGE_NAME "cppsolnp")

include(FetchContent)

FetchContent_Declare(dlib
        GIT_REPOSITORY https://github.com/davisking/dlib.git
        GIT_TAG v19.18)

FetchContent_MakeAvailable(dlib)

add_library(cppsolnp solver.hpp solnp.hpp stdafx.h subnp.hpp utils.hpp)

target_link_libraries(cppsolnp dlib::dlib)

add_library(cppsolnp::cppsolnp ALIAS cppsolnp)
