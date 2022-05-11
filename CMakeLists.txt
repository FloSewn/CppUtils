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
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable -std=c++17 -O3 -DNDEBUG -fconcepts")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable -std=c++17 -O0 -g -fconcepts")

set( BIN ${CMAKE_SOURCE_DIR}/bin )

add_subdirectory( src/cpputils )
add_subdirectory( src/test )



#---------------------------------------------------------------------
# Testing
#---------------------------------------------------------------------
enable_testing()

# does the application run
add_test(NAME run_tests COMMAND ${CMAKE_SOURCE_DIR}/bin/TestCppUtils )
