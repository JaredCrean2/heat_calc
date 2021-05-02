
set(src_files EpwReader.cc utils/lagrange.cc)


set(Boost_USE_MULTITHREADED, OFF)
find_package(Boost)
if (Boost_FOUND)
  message("Found Boost version ${Boost_VERSION}")
  message("Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "Unable to locate Boost")
endif()


add_library(heat_lib "${src_files}")
target_include_directories(heat_lib PUBLIC
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
                          "${PROJECT_SOURCE_DIR}/include"
                          "${PROJECT_SOURCE_DIR}/include/utils"
                          "${Boost_INCLUDE_DIRS}"
                          )

add_executable(heat_calc main.cc)

target_include_directories(heat_calc PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/include"
                          )

target_link_libraries(heat_calc PUBLIC heat_lib)


add_subdirectory(mesh)
