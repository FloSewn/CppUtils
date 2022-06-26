#*********************************************************************
# CppUtils Library
#*********************************************************************
cmake_minimum_required(VERSION 3.10)

project(CppUtils
        LANGUAGES CXX
        VERSION   1.0)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- Compiler flags ---
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable -std=c++17 -O3 -DNDEBUG -fconcepts")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable -std=c++17 -O0 -g -fconcepts")

set( BIN ${CMAKE_SOURCE_DIR}/bin )

# Add config file
configure_file(aux/CppUtilsConfig.h.in ${CMAKE_BINARY_DIR}/CppUtilsConfig.h)
include_directories(${CMAKE_BINARY_DIR})

# Add subdirectories
add_subdirectory( src/cpputils )
add_subdirectory( src/test )
add_subdirectory( src/benchmark )



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
