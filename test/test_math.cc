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

TEST(Math, ArrayOps)
{
  std::array<double, 3> a = {1, 2, 3}, b = {4, 5, 6};
  auto c = a + b;
  EXPECT_EQ(c[0], 5);
  EXPECT_EQ(c[1], 7);
  EXPECT_EQ(c[2], 9);
}