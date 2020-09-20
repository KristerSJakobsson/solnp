project(cpp_solnp_test)
cmake_minimum_required(VERSION 3.12.4)

set(CMAKE_CXX_STANDARD 11)

project(solnp_test)
include(FetchContent)

FetchContent_Declare(dlib
        GIT_REPOSITORY https://github.com/davisking/dlib.git
        GIT_TAG v19.18)

FetchContent_MakeAvailable(dlib)
find_package(Catch2 REQUIRED)

add_executable(cpp_solnp_test test_main.cpp)
add_executable(solnp_test test_main.cpp)

target_link_libraries(cpp_solnp_test dlib::dlib)
target_link_libraries(cpp_solnp_test Catch2::Catch2)
target_link_libraries(solnp_test dlib::dlib)
target_link_libraries(solnp_test Catch2::Catch2)

set (USE_PRECOMPILED_DLIB true)

include(CTest)
include(Catch)
catch_discover_tests(cpp_solnp_test)
catch_discover_tests(solnp_test)


if (ENABLE_COVERAGE)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../library/cmake_codecov")
    list(APPEND LCOV_REMOTE_PATTERNS "${CMAKE_CURRENT_LIST_DIR}/../library/*")
    list(APPEND LCOV_REMOTE_PATTERNS "${CMAKE_CURRENT_LIST_DIR}/../python_solnp/*")

    find_package(codecov)
    add_coverage(pysolnp)
    list(APPEND LCOV_REMOVE_PATTERNS "/usr/")
    coverage_evaluate()
endif()