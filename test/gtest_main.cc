#include "utils/initialization.h"
#include "gtest/gtest.h"
#include <iostream>
#include "mpi.h"

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  initialize(argc, argv);
  MPI_Barrier(MPI_COMM_WORLD);
  auto val = RUN_ALL_TESTS();
  MPI_Barrier(MPI_COMM_WORLD);
  finalize();
  return val;
}