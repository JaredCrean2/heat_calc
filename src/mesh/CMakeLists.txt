message("Processing mesh directory")

add_library(mesh mesh_generator.cc
                 mesh_helper.cc
                 mesh.cc
                 dof_numbering.cc)

target_include_directories(mesh PUBLIC
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
                          "${PROJECT_SOURCE_DIR}/include"
                          "${PROJECT_SOURCE_DIR}/include/utils"
                          "${Boost_INCLUDE_DIRS}"
                          )

target_link_libraries(mesh PUBLIC SCOREC::core)

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
