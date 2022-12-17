#include "gtest/gtest.h"
#include "utils/math.h"

TEST(Math, SmoothAbs)
{
  double delta = 0.1;
  EXPECT_EQ(smoothAbs(5, delta), std::abs(5));
  EXPECT_EQ(smoothAbs(-5, delta), std::abs(-5));
  EXPECT_NEAR(smoothAbs(0.05, delta), 0.05*0.05/delta, 1e-13);
  EXPECT_NEAR(smoothAbs(-0.05, delta), 0.05*0.05/delta, 1e-13);
}

TEST(Math, SmoothAbsDeriv)
{
  double delta = 0.1;
  EXPECT_EQ(smoothAbsDeriv(5, delta), 1);
  EXPECT_EQ(smoothAbsDeriv(-5, delta), -1);
  EXPECT_NEAR(smoothAbsDeriv(0.05, delta), 2*0.05/delta, 1e-13);
  EXPECT_NEAR(smoothAbsDeriv(-0.05, delta), -2*0.05/delta, 1e-13);
}