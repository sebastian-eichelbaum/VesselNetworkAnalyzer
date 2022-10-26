# ----------------------------------------------------------------------------------------
#
# Project: Analyzer
#
# Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
#
# You should have received a copy of the License along with this program.
#
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
#
# General CMake Setup
#
# ----------------------------------------------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.0)
project(analyzer)

# How to call the binary?
set(BinName ${PROJECT_NAME})

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please create a new directory and run CMake from there.")
endif()

# the build types
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Release
        CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo" FORCE
    )
endif(NOT CMAKE_BUILD_TYPE)

# guard against typos in build-type strings
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(NOT cmake_build_type_tolower STREQUAL "debug"
   AND NOT cmake_build_type_tolower STREQUAL "release"
   AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo"
   AND NOT cmake_build_type_tolower STREQUAL ""
)
    message(
        SEND_ERROR
            "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo  and \"\"."
    )
endif()

# ----------------------------------------------------------------------------------------------------------------------
#
# Compiler Setup
#
# ----------------------------------------------------------------------------------------------------------------------

# Create compile_commands.json - required by several tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Unfortunately, this variable is empty. We fill it here and hopefully this is fixed in the near future.
list(APPEND CMAKE_STANDARD_LIBRARIES "stdc++" "m")

# Force latest standard and activate nearly all warnings
set(CMAKE_CXX_FLAGS
    "-frtti -pedantic -std=c++14 -Wall -Wno-long-long -Wextra -pthread"
    CACHE STRING "" FORCE
)

# Darwin's ld isn't GNU and doesn't like the following
if(NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")
    # The following allows us to prevent cyclic dependencies even on linux
    set(CMAKE_SHARED_LINKER_FLAGS
        "-Wl,--no-undefined -Wl,--allow-shlib-undefined,--as-needed"
        CACHE STRING "" FORCE
    )
    set(CMAKE_EXE_LINKER_FLAGS
        "-Wl,--as-needed"
        CACHE STRING "" FORCE
    )
endif()

set(CMAKE_CXX_FLAGS_RELEASE
    "-O3"
    CACHE STRING "" FORCE
)

set(CMAKE_CXX_FLAGS_DEBUG
    "-g -DDEBUG -O0"
    CACHE STRING "" FORCE
)

# Supress compilation warnings from includes residing in system paths, see #230 for further details.
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX
    "-isystem"
    CACHE STRING "" FORCE
)

# Use src/ as include search path
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# ----------------------------------------------------------------------------------------------------------------------
#
# External Code Setup (/ext)
#
# ----------------------------------------------------------------------------------------------------------------------

# Add /ext directory to include search path
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ext/picojson)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ext/vigra/include)

# ----------------------------------------------------------------------------------------------------------------------
#
# Code Setup
#
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Collect everything to compile
# ----------------------------------------------------------------------------------------------------------------------

file(GLOB_RECURSE TARGET_CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE TARGET_H_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)

# ----------------------------------------------------------------------------------------------------------------------
# Build the binary
# ----------------------------------------------------------------------------------------------------------------------

# Setup the target
add_executable(${BinName} ${TARGET_CPP_FILES} ${TARGET_H_FILES})

# Linking
target_link_libraries(${BinName} ${CMAKE_STANDARD_LIBRARIES} ${ADDITIONAL_LINK_LIBRARIES})
