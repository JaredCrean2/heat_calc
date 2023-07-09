cmake_minimum_required(VERSION 3.12.0)

# set project name
project(RetirementCalc)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

find_package(MPI REQUIRED)

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
message("Scorec libraries: ${SCOREC_LIBRARIES}")
endif()


# find BLAS
find_package(BLAS)
if (BLAS_FOUND)
  message("Found BLAS libraries: ${BLAS_LIBRARIES}")
  message("Linker flags: ${BLAS_LINKER_FLAGS}")
else()
  message(FATAL_ERROR, "Unable to locate BLAS")
endif()

if (BLA_STATIC)
  add_library(Blas STATIC IMPORTED GLOBAL)
else()
  add_library(Blas SHARED IMPORTED GLOBAL)
  target_link_options(Blas INTERFACE ${BLAS_LINKER_FLAGS})
endif()

get_filename_component(BLAS_INCLUDE_DIR ${BLAS_LIBRARIES} DIRECTORY)
set(BLAS_INCLUDE_DIR "${BLAS_INCLUDE_DIR}/../include")
message("BLAS_INCLUDE_DIR = ${BLAS_INCLUDE_DIR}")
target_include_directories(Blas INTERFACE ${BLAS_INCLUDE_DIR})
set_property(TARGET Blas PROPERTY
             IMPORTED_LOCATION ${BLAS_LIBRARIES})


# find LAPACK
find_package(LAPACK)
if (LAPACK_FOUND)
  message("Found LAPACK libraries: ${LAPACK_LIBRARIES}")
else()
  message(FATAL_ERROR, "Unable to locate LAPACK")
endif()

if (BLA_STATIC)
  add_library(Lapack STATIC IMPORTED GLOBAL)
else()
  add_library(Lapack SHARED IMPORTED GLOBAL)
  target_link_options(Lapack INTERFACE ${LAPACK_LINKER_FLAGS})
endif()

get_filename_component(LAPACK_INCLUDE_DIR ${LAPACK_LIBRARIES} DIRECTORY)
set(LAPACK_INCLUDE_DIR "${LAPACK_INCLUDE_DIR}/../include")
target_include_directories(Lapack INTERFACE ${LAPACK_INCLUDE_DIR})

set_property(TARGET Lapack PROPERTY
             IMPORTED_LOCATION ${LAPACK_LIBRARIES})

             

# find PETSc
find_package(PETSc)
if (PETSC_FOUND)
  add_library(Petsc SHARED IMPORTED GLOBAL)
  target_include_directories(Petsc INTERFACE ${PETSC_INCLUDES})
  set_property(TARGET Petsc PROPERTY IMPORTED_LOCATION ${PETSC_LIBRARIES})
else()
  message("Unable to locate PETSc")
endif()
#find_package(PkgConfig REQUIRED)
#message("about to search for petsc")
#pkg_check_modules(PETSc REQUIRED IMPORTED_TARGET)
##pkg_check_modules(PETSc REQUIRED)
#message("Found PETSc = ${PETSc_FOUND}")
#message("PETSc_LIBRARIES = ${PETSc_LIBRARIES}")
#message("PETSc_INCLUDE_DIRS = ${PETSc_INCLUDE_DIRS}")
#message("Found PETSc libraries: ${PETSc_LIBRARY_DIRS}")
#get_target_property(PETSC_INCLUDE_DIRS PkgConfig::PETSc INCLUDE_DIRECTORIES)
#message("include dirs = ${PETSC_INCLUDE_DIRS}")





set(ALL_LIBS "")

add_subdirectory(src)
add_subdirectory(test)
