#include "gtest/gtest.h"
#include "utils/spline.h"

namespace {

void test_spline_rev(Real x_lower, Real x_upper, Real f_lower, Real fprime_lower, Real f_upper, Real fprime_upper)
{
  SingleCubicSpline spline;

  int nvals = 5;
  std::vector<Real> test_vals(nvals);
  for (int i=0; i < nvals; ++i)
    test_vals[i] = x_lower + i * (x_upper - x_lower)/(nvals-1);

  for (int i=0; i < nvals; ++i)
  {
    Real result_bar = 2;
    spline.setupSpline(x_lower, x_upper, f_lower, fprime_lower, f_upper, fprime_upper);
    Real f_val = spline.eval(test_vals[i]);

    Real f_lower_bar, fprime_lower_bar, f_upper_bar, fprime_upper_bar, x_bar;
    std::array<Real, 4> coeffs_bar;
    spline.evalRev(test_vals[i], result_bar, x_bar, coeffs_bar);
    spline.setupSplineRev(x_lower, x_upper, coeffs_bar, f_lower_bar, fprime_lower_bar, f_upper_bar, fprime_upper_bar);

    // compute FD
    double eps = 1e-7;

    spline.setupSpline(x_lower, x_upper, f_lower + eps, fprime_lower, f_upper, fprime_upper);
    Real f_lower_bar_fd = (spline.eval(test_vals[i]) - f_val)/eps;

    spline.setupSpline(x_lower, x_upper, f_lower, fprime_lower + eps, f_upper, fprime_upper);
    Real fprime_lower_bar_fd = (spline.eval(test_vals[i]) - f_val)/eps;

    spline.setupSpline(x_lower, x_upper, f_lower, fprime_lower, f_upper + eps, fprime_upper);
    Real f_upper_bar_fd = (spline.eval(test_vals[i]) - f_val)/eps;

    spline.setupSpline(x_lower, x_upper, f_lower, fprime_lower, f_upper, fprime_upper + eps);
    Real fprime_upper_bar_fd = (spline.eval(test_vals[i]) - f_val)/eps;

    spline.setupSpline(x_lower, x_upper, f_lower, fprime_lower, f_upper, fprime_upper);
    Real f_val_bar_fd = (spline.eval(test_vals[i] + eps) - f_val)/eps;

    EXPECT_NEAR(f_lower_bar,      result_bar*f_lower_bar_fd, 1e-6);
    EXPECT_NEAR(fprime_lower_bar, result_bar*fprime_lower_bar_fd, 1e-6);
    EXPECT_NEAR(f_upper_bar,      result_bar*f_upper_bar_fd, 1e-6);
    EXPECT_NEAR(fprime_upper_bar, result_bar*fprime_upper_bar_fd, 1e-6);
    EXPECT_NEAR(x_bar,            result_bar*f_val_bar_fd, 1e-6);
  }

}
}

TEST(SingleCubicSpline, Line)
{
  // line y=x
  SingleCubicSpline spline;

  spline.setupSpline(-1, 1, -1, 1, 1, 1);

  EXPECT_NEAR(spline.eval(-1), -1, 1e-13);
  EXPECT_NEAR(spline.eval(0),   0, 1e-13);
  EXPECT_NEAR(spline.eval(1),   1, 1e-13);
  EXPECT_NEAR(spline.eval(-2), -2, 1e-13);
  EXPECT_NEAR(spline.eval(2),   2, 1e-13);  

  Real dfdx;
  spline.evalDot(-1, 1, dfdx);
  EXPECT_NEAR(dfdx, 1, 1e-13);

  spline.evalDot(1, 1, dfdx);
  EXPECT_NEAR(dfdx, 1, 1e-13);

  test_spline_rev(-1, 1, -1, 1, 1, 1);
}

TEST(SingleCubicSpline, Curve)
{
  // line y=x
  SingleCubicSpline spline;

  spline.setupSpline(298, 302, -1, -0.1, 1, 0.1);

  EXPECT_NEAR(spline.eval(298), -1, 1e-13);
  EXPECT_NEAR(spline.eval(302),  1, 1e-13);

  Real dfdx;
  spline.evalDot(298, 1, dfdx);
  EXPECT_NEAR(dfdx, -0.1, 1e-13);

  spline.evalDot(302, 1, dfdx);
  EXPECT_NEAR(dfdx, 0.1, 1e-13);

  test_spline_rev(298, 302, -1, -0.1, 1, 0.1);
}