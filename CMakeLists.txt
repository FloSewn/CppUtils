#*********************************************************************
# CppUtils Library
#*********************************************************************
cmake_minimum_required(VERSION 3.10)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(CppUtils
        LANGUAGES CXX
        VERSION   1.0)

set( BIN ${CMAKE_SOURCE_DIR}/bin )

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  #set(CMAKE_BUILD_TYPE "Release")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# --- Compiler flags ---
if (CMAKE_BUILD_TYPE MATCHES "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
endif()

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

# Add config file
configure_file(aux/CppUtilsConfig.h.in ${CMAKE_BINARY_DIR}/CppUtilsConfig.h)
include_directories(${CMAKE_BINARY_DIR})

# Add subdirectories
add_subdirectory( src/cpputils )
add_subdirectory( src/test )
add_subdirectory( src/benchmark )

# Info
message(STATUS "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_COMPILER_ID is ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION is ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

#---------------------------------------------------------------------
# Testing
#---------------------------------------------------------------------
enable_testing()

# does the application run
add_test(NAME MathUtility COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "MathUtility")
add_test(NAME StringOps COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "StringOps")
add_test(NAME Vec2 COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "Vec2")
add_test(NAME Geometry COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "Geometry")
add_test(NAME QuadTree COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "QuadTree")
add_test(NAME Container COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "Container")
add_test(NAME ParaReader COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "ParaReader")
add_test(NAME VtkIO COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "VtkIO")
add_test(NAME Log COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "Log")
add_test(NAME Matrix COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils "Matrix")
