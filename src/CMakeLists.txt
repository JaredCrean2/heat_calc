message("Processing source directory")

set(src_files EpwReader.cc 
              utils/lagrange_basis.cc 
              utils/quadrature.cc 
              utils/legendre.cc
              utils/error_handling.cc
              utils/initialization.cc
              discretization/volume_discretization.cc 
              discretization/surface_discretization.cc
              discretization/discretization.cc
              discretization/DirichletBC.cc
              discretization/dof_numbering.cc
              discretization/disc_vector.cc
              #discretization/model.cc
              physics/PhysicsModel.cc
              physics/post_processors.cc
              physics/post_processor_manager.cc
              physics/heat/HeatEquation.cc
              physics/heat/HeatEquationSolar.cc
              physics/heat/mass_matrix.cc
              physics/heat/volume_term.cc
              physics/heat/source_term.cc
              physics/heat/dirichlet_bc.cc
              physics/heat/neumann_bc.cc
              physics/heat/bc_defs.cc
              physics/heat/sky_radiation.cc
              physics/heat/tarp.cc
              physics/heat/dates.cc
              physics/heat/solar_position.cc
              physics/heat/basis_vals.cc
              physics/heat/AuxiliaryEquationsSolar.cc
              physics/heat/interior_temperature_update.cc
              linear_system/assembler.cc
              linear_system/bla_wrapper.cc
              linear_system/large_matrix.cc
              linear_system/large_matrix_dense.cc
              linear_system/large_matrix_petsc.cc
              linear_system/large_matrix_factory.cc
              linear_system/sparsity_pattern_mesh.cc
              time_solver/newton_result.cc
              time_solver/newton.cc
              time_solver/crank_nicolson.cc
              time_solver/crank_nicolson_function.cc
              time_solver/crank_nicolson_aux_equations.cc
              )

add_subdirectory(mesh)

add_library(heat_lib "${src_files}")
target_include_directories(heat_lib PUBLIC
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
                          "${PROJECT_SOURCE_DIR}/include"
                          "${PROJECT_SOURCE_DIR}/include/utils"
                          "${Boost_INCLUDE_DIRS}"
)

target_link_libraries(heat_lib PUBLIC mesh Blas Lapack Petsc MPI::MPI_CXX)

set(ALL_LIBS ${ALL_LIBS} heat_lib PARENT_SCOPE)
message("ALL_LIBS = ${ALL_LIBS}")

add_executable(heat_calc main.cc)

target_include_directories(heat_calc PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/include"
                          )
target_link_libraries(heat_calc PUBLIC heat_lib)


add_executable(simple_house simple_house.cc)

target_include_directories(simple_house PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/include"
                          )
target_link_libraries(simple_house PUBLIC heat_lib)



#add_executable(gauss_seidel gauss_seidel.cc)
#
#target_include_directories(gauss_seidel PUBLIC
#                          "${PROJECT_BINARY_DIR}"
#                          "${PROJECT_SOURCE_DIR}/include"
#                          )
#
#target_link_libraries(gauss_seidel PUBLIC heat_lib)
