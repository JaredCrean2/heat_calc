#include "utils/initialization.h"
#include "gtest/gtest.h"
#include <iostream>

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  initialize(argc, argv);
  auto val = RUN_ALL_TESTS();
  finalize();
  return val;
}