cmake_minimum_required(VERSION 3.5.1)

# set project name
project(RetirementCalc)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(ALL_LIBS "")

add_subdirectory(src)
add_subdirectory(test)
