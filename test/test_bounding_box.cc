#include "gtest/gtest.h"
#include "utils/bounding_box.h"

TEST(BoundingBox, Contains)
{
  utils::BoundingBox box({0, 0, 0}, {1, 1, 1});

  EXPECT_TRUE(box.contains({0.5, 0.5, 0.5}));
  EXPECT_FALSE(box.contains({-1, 0.5, 0.5}));
  EXPECT_FALSE(box.contains({0.5, -1, 0.5}));
  EXPECT_FALSE(box.contains({0.5, 0.5, -1}));
}