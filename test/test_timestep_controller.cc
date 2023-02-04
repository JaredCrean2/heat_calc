#include "gtest/gtest.h"
#include <ctime>
#include "time_solver/timestep_controller.h"

TEST(TimestepControllerConstant, All)
{
  Real delta_t = 0.1;
  timesolvers::TimestepControllerConstant controller(delta_t);

  EXPECT_EQ(controller.getNextTimestep(0), delta_t);
  EXPECT_EQ(controller.getNextTimestep(1), delta_t);
  EXPECT_EQ(controller.getNextTimestep(2), delta_t);

  EXPECT_NO_THROW(controller.recordLastIteration(1e-6));
}