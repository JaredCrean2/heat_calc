message("Processing mesh directory")

add_library(mesh mesh_generator.cc
                 mesh_generator_multi.cc
                 mesh_geometry_multi.cc
                 mesh_block.cc
                 mesh_generator_multi_block.cc
                 mesh_helper.cc
                 mesh.cc
                 apf_data.cc
                 dof_numbering.cc
                 reference_element.cc
                 reference_element_geometry.cc
                 reference_element_geometry_hex.cc
                 reference_element_geometry_interface.cc
                 reference_element_nodes.cc
                 reference_element_interface.cc
                 reference_element_apf.cc
                 gmiDataStructure.cc
                 apfShapeHex.cc
                 apfMDSField.cc
                 FieldDataManager.cc
                 dirichlet_update_map.cc
                 ghosting.cc)

target_include_directories(mesh PUBLIC
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
                          "${PROJECT_SOURCE_DIR}/include"
                          "${PROJECT_SOURCE_DIR}/include/utils"
                          "${Boost_INCLUDE_DIRS}"
                          )

target_link_libraries(mesh PUBLIC SCOREC::core SCOREC::mds SCOREC::pcu MPI::MPI_CXX)

#TODO: setup interface linking

set(ALL_LIBS ${ALL_LIBS} mesh PARENT_SCOPE)

message("ALL_LIBS = ${ALL_LIBS}")


add_executable(generate_mesh generate_mesh.cc)

target_include_directories(generate_mesh PUBLIC
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
                          "${PROJECT_SOURCE_DIR}/include"
                          "${PROJECT_SOURCE_DIR}/include/utils"
                          "${Boost_INCLUDE_DIRS}"
                          )
 
target_link_libraries(generate_mesh
                      PRIVATE mesh)
