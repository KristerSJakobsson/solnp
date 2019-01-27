cmake_minimum_required(VERSION 3.10)
project(cpp_solnp)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(libraries/dlib-19.16/dlib dlib_build)
find_package(Catch2 REQUIRED)

add_executable(cpp_solnp test_main.cpp)

target_link_libraries(cpp_solnp dlib::dlib)
target_link_libraries(cpp_solnp Catch2::Catch2)

include(CTest)
include(Catch)
catch_discover_tests(cpp_solnp)