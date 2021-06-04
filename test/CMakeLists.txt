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
                     SCOREC::core)

message("Linking to ALL_LIBS = ${ALL_LIBS}")
#target_link_libraries(unit_tests PUBLIC SCOREC::core)
