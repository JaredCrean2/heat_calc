#include <cmath>

#include "gtest/gtest.h"
#include "utils/lagrange.h"
#include "utils/legendre.h"
#include "utils/quadrature.h"

Real legendre0(const Real x)
{
  return 1;
}

Real legendre0_dx(const Real x)
{
  return 0;
}

Real legendre1(const Real x)
{
  return x;
}

Real legendre1_dx(const Real x)
{
  return 1;
}

Real legendre2(const Real x)
{
  return 0.5*(3*x*x - 1);
}

Real legendre2_dx(const Real x)
{
  return 0.5*(6*x);
}

Real legendre3(const Real x)
{
  return 0.5*(5*std::pow(x, 3) - 3*x);
}

Real legendre3_dx(const Real x)
{
  return 0.5*(15*std::pow(x, 2) - 3);
}

Real legendre4(const Real x)
{
  return 0.125*(35*std::pow(x, 4) - 30*std::pow(x, 2) + 3);
}

Real legendre4_dx(const Real x)
{
  return 0.125*(4*35*std::pow(x, 3) - 2*30*x);
}



TEST(Polynomials, Legendre)
{
  LegendrePoly legendre;


  std::vector<Real> xvals = {-1, -0.5, 0, 0.5, 1};
  int nvals = xvals.size();

  // p = 0
  for (int i=0; i < nvals; ++i)
  {
    EXPECT_FLOAT_EQ(legendre.evalPoly(0, xvals[i]), legendre0(xvals[i]));
    EXPECT_FLOAT_EQ(legendre.evalPolyDeriv(0, xvals[i]), legendre0_dx(xvals[i]));
  }

  // p = 1
  for (int i=0; i < nvals; ++i)
  {
    EXPECT_FLOAT_EQ(legendre.evalPoly(1, xvals[i]), legendre1(xvals[i]));
    EXPECT_FLOAT_EQ(legendre.evalPolyDeriv(1, xvals[i]), legendre1_dx(xvals[i]));
  }

  // p = 2
  for (int i=0; i < nvals; ++i)
  {
    EXPECT_FLOAT_EQ(legendre.evalPoly(2, xvals[i]), legendre2(xvals[i]));
    EXPECT_FLOAT_EQ(legendre.evalPolyDeriv(2, xvals[i]), legendre2_dx(xvals[i]));
  }

  // p = 3
  for (int i=0; i < nvals; ++i)
  {
    EXPECT_FLOAT_EQ(legendre.evalPoly(3, xvals[i]), legendre3(xvals[i]));
    EXPECT_FLOAT_EQ(legendre.evalPolyDeriv(3, xvals[i]), legendre3_dx(xvals[i]));
  }

  // p = 4
  for (int i=0; i < nvals; ++i)
  {
    EXPECT_FLOAT_EQ(legendre.evalPoly(4, xvals[i]), legendre4(xvals[i]));
    EXPECT_FLOAT_EQ(legendre.evalPolyDeriv(4, xvals[i]), legendre4_dx(xvals[i]));
  }
}


const std::vector<Real> pts = {-1, -0.5, 0.5, 1};

Real lagrange0(const Real x)
{
  Real xi = pts[0];
  return (x - pts[1])*(x - pts[2])*(x - pts[3]) /( (xi - pts[1])*(xi - pts[2])*(xi - pts[3]) );
}

Real lagrange0_dx(const Real x)
{
  Real xi = pts[0];
  Real den = (xi - pts[1])*(xi - pts[2])*(xi - pts[3]);
  Real num = (x - pts[2])*(x - pts[3]) + (x - pts[1])*(x - pts[3]) + (x - pts[1])*(x - pts[2]);

  return num / den;
}

Real lagrange1(const Real x)
{
  Real xi = pts[1];
  return (x - pts[0])*(x - pts[2])*(x - pts[3]) /( (xi - pts[0])*(xi - pts[2])*(xi - pts[3]) );
}

Real lagrange1_dx(const Real x)
{
  Real xi = pts[1];
  Real den = (xi - pts[0])*(xi - pts[2])*(xi - pts[3]);
  Real num = (x - pts[2])*(x - pts[3]) + (x - pts[0])*(x - pts[3]) + (x - pts[0])*(x - pts[2]);

  return num / den;
}

Real lagrange2(const Real x)
{
  Real xi = pts[2];
  return (x - pts[0])*(x - pts[1])*(x - pts[3]) /( (xi - pts[0])*(xi - pts[1])*(xi - pts[3]) );
}

Real lagrange2_dx(const Real x)
{
  Real xi = pts[2];
  Real den = (xi - pts[0])*(xi - pts[1])*(xi - pts[3]);
  Real num = (x - pts[1])*(x - pts[3]) + (x - pts[0])*(x - pts[3]) + (x - pts[0])*(x - pts[1]);

  return num / den;
}

Real lagrange3(const Real x)
{
  Real xi = pts[3];
  return (x - pts[0])*(x - pts[1])*(x - pts[2]) /( (xi - pts[0])*(xi - pts[1])*(xi - pts[2]) );
}

Real lagrange3_dx(const Real x)
{
  Real xi = pts[3];
  Real den = (xi - pts[0])*(xi - pts[1])*(xi - pts[2]);
  Real num = (x - pts[1])*(x - pts[2]) + (x - pts[0])*(x - pts[2]) + (x - pts[0])*(x - pts[1]);

  return num / den;
}

TEST(Polynomials, Lagrange)
{
  std::vector<Real> pts_test = {-2, -1, -0.75, -0.5, 0.25, 0, 0.25, 0.5, 0.75, 1, 1};
  LagrangeBasis lagrange(pts);

  // j = 0
  for (auto pts_i : pts_test)
  {
    EXPECT_FLOAT_EQ(lagrange.evalPoly(0, pts_i), lagrange0(pts_i));
    EXPECT_FLOAT_EQ(lagrange.evalPolyDeriv(0, pts_i), lagrange0_dx(pts_i));
  }

  for (auto pts_i : pts_test)
  {
    EXPECT_FLOAT_EQ(lagrange.evalPoly(1, pts_i), lagrange1(pts_i));
    EXPECT_FLOAT_EQ(lagrange.evalPolyDeriv(1, pts_i), lagrange1_dx(pts_i));
  }

  for (auto pts_i : pts_test)
  {
    EXPECT_FLOAT_EQ(lagrange.evalPoly(2, pts_i), lagrange2(pts_i));
    EXPECT_FLOAT_EQ(lagrange.evalPolyDeriv(2, pts_i), lagrange2_dx(pts_i));
  }

  for (auto pts_i : pts_test)
  {
    EXPECT_FLOAT_EQ(lagrange.evalPoly(3, pts_i), lagrange3(pts_i));
    EXPECT_FLOAT_EQ(lagrange.evalPolyDeriv(3, pts_i), lagrange3_dx(pts_i));
  }
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
