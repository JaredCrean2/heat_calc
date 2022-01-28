message("Processing source directory")

set(src_files EpwReader.cc 
              utils/lagrange_basis.cc 
              utils/quadrature.cc 
              utils/legendre.cc
              utils/error_handling.cc
              discretization/volume_discretization.cc 
              discretization/surface_discretization.cc
              discretization/discretization.cc
              discretization/DirichletBC.cc
              discretization/dof_numbering.cc
              discretization/disc_vector.cc
              discretization/model.cc
              physics/PhysicsModel.cc
              physics/heat/HeatEquation.cc
              )

add_subdirectory(mesh)

add_library(heat_lib "${src_files}")
target_include_directories(heat_lib PUBLIC
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
                          "${PROJECT_SOURCE_DIR}/include"
                          "${PROJECT_SOURCE_DIR}/include/utils"
                          "${Boost_INCLUDE_DIRS}"
)

target_link_libraries(heat_lib PUBLIC mesh)

set(ALL_LIBS ${ALL_LIBS} heat_lib PARENT_SCOPE)
message("ALL_LIBS = ${ALL_LIBS}")

add_executable(heat_calc main.cc)

target_include_directories(heat_calc PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/include"
                          )

target_link_libraries(heat_calc PUBLIC heat_lib)


