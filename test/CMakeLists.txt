cmake_minimum_required(VERSION 2.8.12)

project(cpp_solnp_test)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

target_include_directories(Catch INTERFACE catch.hpp)

add_executable(cpp_solnp_test test_main.cpp)

target_link_libraries(cpp_solnp_test dlib::dlib)

include(CTest)
include(Catch)
catch_discover_tests(cpp_solnp_test)


if (ENABLE_COVERAGE)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../library/cmake_codecov")
    list(APPEND LCOV_REMOTE_PATTERNS "${CMAKE_CURRENT_LIST_DIR}/../library/*")
    list(APPEND LCOV_REMOTE_PATTERNS "${CMAKE_CURRENT_LIST_DIR}/../python_solnp/*")

    find_package(codecov)
    add_coverage(cpp_solnp_test)
    list(APPEND LCOV_REMOVE_PATTERNS "/usr/")
    coverage_evaluate()
endif()