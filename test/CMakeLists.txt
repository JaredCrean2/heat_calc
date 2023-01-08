enable_testing()
find_package(GTest REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS})

file(GLOB src_files *.cc)

add_executable(unit_tests ${src_files})
target_include_directories(unit_tests PUBLIC
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}/include"
  "${Boost_INCLUDE_DIRS}"
)
target_link_libraries(unit_tests
                     ${GTEST_BOTH_LIBRARIES}
                     pthread 
                     mesh
                     heat_lib
                     SCOREC::core)

file(GLOB verification_src_files ./verification/*.cc)

add_executable(verification_tests ${verification_src_files} mesh_helper.cc gtest_main.cc test_helper.cc)
target_include_directories(verification_tests PUBLIC
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}/include"
  "${PROJECT_SOURCE_DIR}/test"
  "${Boost_INCLUDE_DIRS}"
)
target_link_libraries(verification_tests
                     ${GTEST_BOTH_LIBRARIES}
                     pthread 
                     mesh
                     heat_lib
                     SCOREC::core
                     MPI::MPI_CXX)

message("Linking to ALL_LIBS = ${ALL_LIBS}")
#target_link_libraries(unit_tests PUBLIC SCOREC::core)

file(COPY
    "${CMAKE_CURRENT_SOURCE_DIR}/data/USA_AZ_Phoenix-Sky.Harbor.Intl.AP.722780_TMY3.epw"
    DESTINATION
    "${CMAKE_CURRENT_BINARY_DIR}/data/")
