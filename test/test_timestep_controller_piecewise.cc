#include "gtest/gtest.h"
#include "time_solver/timestep_controller_piecewise.h"

TEST(TimestepControllerPiecewise, TwoPoints)
{
  std::vector<timesolvers::TimestepControllerPiecewise::TimestepPoint> pts = {{0, 1}, {2, 11} };

  timesolvers::TimestepControllerPiecewise controller(pts, true);

  EXPECT_ANY_THROW(controller.getNextTimestep(-1));
  EXPECT_ANY_THROW(controller.getNextTimestep(3));

  EXPECT_NEAR(controller.getNextTimestep(0), 1, 1e-13);
  EXPECT_NEAR(controller.getNextTimestep(2), 11, 1e-13);
  EXPECT_NEAR(controller.getNextTimestep(1), 6, 1e-13);
}

TEST(TimestepControllerPiecewise, OnePointError)
{
  std::vector<timesolvers::TimestepControllerPiecewise::TimestepPoint> pts = {{0, 1}, {0, 11} };

  EXPECT_ANY_THROW(timesolvers::TimestepControllerPiecewise controller(pts, true));
}

TEST(TimestepControllerPiecewise, TwoPointsNonSorted)
{
  std::vector<timesolvers::TimestepControllerPiecewise::TimestepPoint> pts = {{2, 11}, {0, 1} };

  timesolvers::TimestepControllerPiecewise controller(pts, true);

  EXPECT_ANY_THROW(controller.getNextTimestep(-1));
  EXPECT_ANY_THROW(controller.getNextTimestep(3));

  EXPECT_NEAR(controller.getNextTimestep(0), 1, 1e-13);
  EXPECT_NEAR(controller.getNextTimestep(2), 11, 1e-13);
  EXPECT_NEAR(controller.getNextTimestep(1), 6, 1e-13);
}


// non-increasing points