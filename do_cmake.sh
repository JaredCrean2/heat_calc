#!/bin/bash

# ASAN flags
#ASAN_FLAGS="-g -fsanitize=address -fno-omit-frame-pointer"
#CXXFLAGS="-O0 -Wall -g"
#BUILDTYPE="Debug"

CXXFLAGS="-O3 -Wall"
BUILDTYPE="Release"


cmake \
-D CMAKE_C_COMPILER=`which mpicc` \
-D CMAKE_EXPORT_COMPILE_COMMANDS=1 \
-D CMAKE_CXX_COMPILER=`which mpicxx` \
-D CMAKE_BUILD_TYPE=$BUILDTYPE \
-D CMAKE_CXX_FLAGS="$CXXFLAGS $ASAN_FLAGS" \
..


#cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS="-O3 -Wall -Werror" .. 
