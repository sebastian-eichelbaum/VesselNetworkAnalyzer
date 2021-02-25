#----------------------------------------------------------------------------------------
#
# Project: Analyzer
#
# Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
#
# You should have received a copy of the License along with this program.
#
#----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------
#
# General CMake Setup
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------

CMAKE_MINIMUM_REQUIRED( VERSION 3.0 )
PROJECT( analyzer )

# How to call the binary?
SET( BinName ${PROJECT_NAME} )

# guard against in-source builds
IF( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
  MESSAGE( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there (you may need to remove CMakeCache.txt)." )
ENDIF()

# the build types
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE Release
         CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE
       )
ENDIF( NOT CMAKE_BUILD_TYPE )

# guard against typos in build-type strings
STRING( TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
IF( NOT cmake_build_type_tolower STREQUAL "debug" AND
    NOT cmake_build_type_tolower STREQUAL "release" AND
    NOT cmake_build_type_tolower STREQUAL "relwithdebinfo" AND
    NOT cmake_build_type_tolower STREQUAL "")
  MESSAGE( SEND_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo  and \"\" (case-insensitive).")
ENDIF()

# ---------------------------------------------------------------------------------------------------------------------------------------------------
#
# Compiler Setup
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# Unfortunately, this variable is empty. We fill it here and hopefully this is fixed in the near future.
LIST( APPEND CMAKE_STANDARD_LIBRARIES "stdc++" "m" )

# Force latest standard and activate nearly all warnings
SET( CMAKE_CXX_FLAGS "-frtti -pedantic -std=c++14 -Wall -Wno-long-long -Wextra -pthread" CACHE STRING "" FORCE )

# Darwin's ld isn't GNU and doesn't like the following
IF( NOT CMAKE_SYSTEM_NAME MATCHES "Darwin" )
    ## The following allows us to prevent cyclic dependencies even on linux
    SET( CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined -Wl,--allow-shlib-undefined,--as-needed" CACHE STRING "" FORCE )
    SET( CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed" CACHE STRING "" FORCE )
ENDIF()
SET( CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "" FORCE )
SET( CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG -O0" CACHE STRING "" FORCE )

# Supress compilation warnings from includes residing in system paths, see #230 for further details.
SET( CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem" CACHE STRING "" FORCE )

# Use src/ as include search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/src )

# ---------------------------------------------------------------------------------------------------------------------------------------------------
#
# External Code Setup (/ext)
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# Add /ext directory to include search path
INCLUDE_DIRECTORIES( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ext/picojson )
INCLUDE_DIRECTORIES( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ext/vigra/include )

# ---------------------------------------------------------------------------------------------------------------------------------------------------
#
# Code Setup
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Collect everything to compile
# ---------------------------------------------------------------------------------------------------------------------------------------------------

FILE( GLOB_RECURSE TARGET_CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp )
FILE( GLOB_RECURSE TARGET_H_FILES   ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h )

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Build the binary
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# Setup the target
ADD_EXECUTABLE( ${BinName} ${TARGET_CPP_FILES} ${TARGET_H_FILES} )

# Linking
TARGET_LINK_LIBRARIES( ${BinName} ${CMAKE_STANDARD_LIBRARIES}
                                  ${ADDITIONAL_LINK_LIBRARIES} )

