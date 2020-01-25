cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 11)

project(solnp_test)
include(FetchContent)

FetchContent_Declare(dlib
        GIT_REPOSITORY https://github.com/davisking/dlib.git
        GIT_TAG v19.18)

FetchContent_MakeAvailable(dlib)

FetchContent_Declare(
        catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2
        GIT_TAG v2.11.1)

FetchContent_GetProperties(catch2)
if (NOT catch2_POPULATED)
    FetchContent_Populate(catch2)
    add_subdirectory(${catch2_SOURCE_DIR} ${catch2_BINARY_DIR})
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)
    include(Catch)
endif()

add_executable(solnp_test test_main.cpp)

target_link_libraries(solnp_test dlib::dlib)
target_link_libraries(solnp_test Catch2::Catch2)

include(CTest)
include(Catch)
catch_discover_tests(solnp_test)
