#!/bin/bash

# ASAN flags
ASAN_FLAGS="-g -fsanitize=address -fno-omit-frame-pointer"
CXXFLAGS="-O0 -Wall -g -fdebug-default-version=4"
BUILDTYPE="Debug"

# Optimized build
#CXXFLAGS="-Ofast -march=native -mtune=native -ffast-math -Wall"
#BUILDTYPE="Release"

# VTune build
#CXXFLAGS="-Ofast -march=native -mtune=native -g -Rpass-analysis=loop-vectorize -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -ffast-math"
#BUILDTYPE="Release"

# set PETSC_DIR and PETSC_ARCH to guide which Petsc installation CMake finds

cmake \
-D CMAKE_C_COMPILER=`which clang` \
-D CMAKE_CXX_COMPILER=`which clang++` \
-D MPI_C_COMPILER=`which mpicc` \
-D MPI_CXX_COMPILER=`which mpicxx` \
-D CMAKE_EXPORT_COMPILE_COMMANDS=1 \
-D CMAKE_BUILD_TYPE=$BUILDTYPE \
-D CMAKE_CXX_FLAGS="$CXXFLAGS $ASAN_FLAGS" \
-D BLAS_ROOT="$HOME/build/OpenBLAS_install/" \
..


#cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS="-O3 -Wall -Werror" .. 
