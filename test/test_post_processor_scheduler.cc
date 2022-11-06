#include "gtest/gtest.h"
#include "physics/post_processor_scheduler.h"

TEST(PostProcessorScheduler, FixedInterval)
{
  physics::PostProcessorScheduleFixedInterval scheduler(3);

  double delta_t = 0.1;
  EXPECT_TRUE(scheduler.shouldOutput(0, 0));
  EXPECT_FALSE(scheduler.shouldOutput(1,   delta_t));
  EXPECT_FALSE(scheduler.shouldOutput(2, 2*delta_t));
  EXPECT_TRUE( scheduler.shouldOutput(3, 3*delta_t));

  EXPECT_FALSE(scheduler.shouldOutput(4, 4*delta_t));
  EXPECT_FALSE(scheduler.shouldOutput(5, 5*delta_t));
  EXPECT_TRUE( scheduler.shouldOutput(6, 6*delta_t));
}