#include "gtest/gtest.h"
#include "time_solver/timestep_controller.h"

TEST(TimestepControllerResidual, All)
{
  Real delta_t0 = 2;
  Real exponent = 2;
  timesolvers::TimestepControllerResidual controller(delta_t0, exponent);

  EXPECT_EQ(controller.getNextTimestep(0), delta_t0);

  controller.recordLastIteration(1);
  EXPECT_EQ(controller.getNextTimestep(0), delta_t0);

  controller.recordLastIteration(0.5);
  EXPECT_NEAR(controller.getNextTimestep(0), delta_t0*4, 1e-13);

  controller.recordLastIteration(0.25);
  EXPECT_NEAR(controller.getNextTimestep(0), delta_t0*16, 1e-13);  
}