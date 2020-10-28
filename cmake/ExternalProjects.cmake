include(FindGit)
find_package(Git)

if (NOT Git_FOUND)
    message(FATAL_ERROR "Git not found!")
endif ()

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(ExternalProject)

message("Adding dlib")

ExternalProject_Add (
        dlib
        SOURCE_DIR        ${CMAKE_CURRENT_SOURCE_DIR}/library/dlib
)

message("Adding Catch")

#ExternalProject_Add(
#        catch
#        BUILD_IN_SOURCE 1
#        PREFIX ${CMAKE_BINARY_DIR}/catch
#        SOURCE_DIR        ${CMAKE_CURRENT_SOURCE_DIR}/library/Catch2
#)
#set (CATCH "Catch2")
#ExternalProject_Add(
#        catch
#        PREFIX ${CMAKE_BINARY_DIR}/catch
#        SOURCE_DIR ${source_dir}
#        GIT_REPOSITORY https://github.com/philsquared/Catch.git
#        TIMEOUT 10
#        UPDATE_COMMAND ${GIT_EXECUTABLE} pull
#        CONFIGURE_COMMAND ""
#        BUILD_COMMAND ""
#        INSTALL_COMMAND ""
#        LOG_DOWNLOAD ON
#)
#
## Expose required variable (CATCH_INCLUDE_DIR) to parent scope
#message("Source directory ${source_dir}")
#message("Binary directory 1 ${CMAKE_BINARY_DIR}/catch")
#ExternalProject_Get_Property(catch source_dir)
##set(CATCH_INCLUDE_DIR ${source_dir}/include CACHE INTERNAL "Path to include folder for Catch")
## CMAKE_MODULE_PATH
#LIST(APPEND CMAKE_MODULE_PATH source_dir)
#message(Source directory 2 ${source_dir})
#
#find_package(Catch2 MODULE)
#include(Catch)

find_package(Git REQUIRED)

ExternalProject_Add(
        catch
        PREFIX ${CMAKE_BINARY_DIR}/test/catch
        GIT_REPOSITORY https://github.com/philsquared/Catch.git
        TIMEOUT 10
        UPDATE_COMMAND ${GIT_EXECUTABLE} pull
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
)

# Expose required variable (CATCH_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(catch source_dir)
set(CATCH_INCLUDE_DIR ${source_dir}/include CACHE INTERNAL "Path to include folder for Catch")

# Add catch as an interface library
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})
