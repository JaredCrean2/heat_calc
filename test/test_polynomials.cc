#include <cmath>

#include "gtest/gtest.h"
#include "utils/legendre.h"
#include "utils/quadrature.h"

Real legendre0(const Real x)
{
  return 1;
}

Real legendre1(const Real x)
{
  return x;
}

Real legendre2(const Real x)
{
  return 0.5*(3*x*x -1);
}

Real legendre3(const Real x)
{
  return 0.5*(5*std::pow(x, 3) - 3*x);
}

Real legendre4(const Real x)
{
  return 0.125*(35*std::pow(x, 4) - 30*std::pow(x, 2) + 3);
}




TEST(Polynomials, Legendre)
{
  LegendrePoly legendre;


  std::vector<Real> xvals = {-1, -0.5, 0, 0.5, 1};
  int nvals = xvals.size();

  // p = 0
  for (int i=0; i < nvals; ++i)
    EXPECT_FLOAT_EQ(legendre.evalPoly(0, xvals[i]), legendre0(xvals[i]));

  // p = 1
  for (int i=0; i < nvals; ++i)
    EXPECT_FLOAT_EQ(legendre.evalPoly(1, xvals[i]), legendre1(xvals[i]));

  // p = 2
  for (int i=0; i < nvals; ++i)
    EXPECT_FLOAT_EQ(legendre.evalPoly(2, xvals[i]), legendre2(xvals[i]));

  // p = 3
  for (int i=0; i < nvals; ++i)
    EXPECT_FLOAT_EQ(legendre.evalPoly(3, xvals[i]), legendre3(xvals[i]));

  // p = 4
  for (int i=0; i < nvals; ++i)
    EXPECT_FLOAT_EQ(legendre.evalPoly(4, xvals[i]), legendre4(xvals[i]));
}

TEST(Polynomials, Quadrature)
{
  // test using orthogonality property of Legendre polnomials

  LegendrePoly legendre;
  std::vector<int> exactness = {1, 3, 5, 7};

  for (auto degree_max : exactness)
  {
    Quadrature quad = getGaussianQuadrature(degree_max);
    for (int degree=0; degree < degree_max/2; ++ degree)
    {
      // integrate
      Real val = 0;
      for (int i=0; i < quad.getNumPoints(); ++i)
      {
        Real integrand = std::pow(legendre.evalPoly(degree, quad.getPoint(i)), 2);
        val += quad.getWeight(i) * integrand;
      }

      EXPECT_FLOAT_EQ(val, 2.0/(2*degree + 1));
    }
  }

}
