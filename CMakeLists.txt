cmake_minimum_required(VERSION 3.5.1)

# set project name
project(RetirementCalc)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(Boost_USE_MULTITHREADED, OFF)
find_package(Boost)
if (Boost_FOUND)
  message("Found Boost version ${Boost_VERSION}")
  message("Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "Unable to locate Boost")
endif()

# Allow the user to indicate where they installed SCOREC
# via "-DSCOREC_PREFIX=/home/somewhere" when calling `cmake`
set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")

# If SCOREC_PREFIX was specified, only link to that directory,
# i.e. don't link to another installation in /usr/lib by mistake
if (SCOREC_PREFIX)
  find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
else()
# IF SCOREC_PREFIX was not specified, look in typical system directories,
# and also in CMAKE_PREFIX_PATH (environment variable)
  find_package(
    SCOREC #package name, has to be SCOREC
    2.1.0  #version. can be omitted, and will match any installed version
    #greater than or equal to this one, as long as the major number
    #is the same
    REQUIRED #indicate that SCOREC is really needed to compile
    CONFIG   #skip the 'MODULE' search system, save some time and confusion
  )

endif()



set(ALL_LIBS "")

add_subdirectory(src)
add_subdirectory(test)
