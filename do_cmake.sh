#!/bin/bash

cmake \
-D CMAKE_C_COMPILER=`which mpicc` \
-D CMAKE_CXX_COMPILER=`which mpicxx` \
-D CMAKE_BUILD_TYPE=Debug \
-D CMAKE_CXX_FLAGS="-O0 -Wall -Werror -g" \
..


#cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS="-O3 -Wall -Werror" .. 
