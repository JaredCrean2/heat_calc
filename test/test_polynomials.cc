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

Real testPoly(const Real x, const Real y, const Real z)
{
  return x*x +  y*y + z*z + x*y + x*z + y*z;
}

Real testPoly_dx(const Real x, const Real y, const Real z)
{
  return 2*x + y + z;
}

Real testPoly_dy(const Real x, const Real y, const Real z)
{
  return 2*y + + x + z;
}

Real testPoly_dz(const Real x, const Real y, const Real z)
{
  return 2*z + x + y;
}

TEST(Polynomials, LagrangeTP)
{
  std::vector<Real> pts_in = {-1.0, 0.0, 1.0};
  std::vector<Real> pts_out = {-0.75, -0.25, 0.25, 0.75};
  LagrangeEvaluatorTP basis(pts_in, pts_out);

  {
    ArrayType<Real, 3> vals_in(boost::extents[3][3][3]);
    ArrayType<Real, 3> vals_out(boost::extents[4][4][4]);
    ArrayType<Real, 4> derivs_out(boost::extents[4][4][4][3]);

    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
        for (int k=0; k < 3; ++k)
          vals_in[i][j][k] = testPoly(pts_in[i], pts_in[j], pts_in[k]);

    // test interpolation functions
    basis.interpolateVals(vals_in, vals_out);
    basis.interpolateDerivs(vals_in, derivs_out);

    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
        for (int k=0; k < 4; ++k)
        {
          EXPECT_FLOAT_EQ(vals_out[i][j][k], testPoly(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out[i][j][k][0], testPoly_dx(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out[i][j][k][1], testPoly_dy(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out[i][j][k][2], testPoly_dz(pts_out[i], pts_out[j], pts_out[k]));
        }

    // test getting interpolant values
    for (int i_out=0; i_out < 4; ++i_out)
      for (int j_out=0; j_out < 4; ++j_out)
        for (int k_out=0; k_out < 4; ++k_out)
        {
          Real val = 0;
          Real val_x = 0;
          Real val_y = 0;
          Real val_z = 0;

          for (int i_in=0; i_in < 3; ++i_in)
            for (int j_in=0; j_in < 3; ++j_in)
              for (int k_in=0; k_in <3; ++k_in)
              {
                val += vals_in[i_in][j_in][k_in] *
                  basis.getInterpolantValue(i_in, j_in, k_in, i_out, j_out, k_out);
                val_x += vals_in[i_in][j_in][k_in] * 
                  basis.getInterpolantDx(i_in, j_in, k_in, i_out, j_out, k_out);
                val_y += vals_in[i_in][j_in][k_in] * 
                  basis.getInterpolantDy(i_in, j_in, k_in, i_out, j_out, k_out);
                val_z += vals_in[i_in][j_in][k_in] * 
                  basis.getInterpolantDz(i_in, j_in, k_in, i_out, j_out, k_out);
              }

          EXPECT_FLOAT_EQ(val, testPoly(pts_out[i_out], pts_out[j_out], pts_out[k_out]));
          EXPECT_FLOAT_EQ(val_x, testPoly_dx(pts_out[i_out], pts_out[j_out], pts_out[k_out]));
          EXPECT_FLOAT_EQ(val_y, testPoly_dy(pts_out[i_out], pts_out[j_out], pts_out[k_out]));
          EXPECT_FLOAT_EQ(val_z, testPoly_dz(pts_out[i_out], pts_out[j_out], pts_out[k_out]));
        }
     

  }

}


TEST(Polynomials, LagrangeTPIn)
{
  std::vector<Real> pts_in = {-1.0, 0.0, 1.0};
  ArrayType<Real, 2> pts_out(boost::extents[4][3]);
  pts_out[0][0] = -0.75; pts_out[0][1] = -0.25; pts_out[0][2] = -0.5;
  pts_out[1][0] = -0.5;  pts_out[1][1] = 0.0;   pts_out[1][2] = -0.25;
  pts_out[2][0] = 0.0;   pts_out[2][1] = 0.25;  pts_out[2][2] = 0.1;
  pts_out[3][0] = 0.5;   pts_out[3][1] = 0.75;  pts_out[3][2] = 0.5;
  LagrangeEvaluatorTPIn basis(pts_in, pts_out);

  EXPECT_EQ(basis.getNumPointsIn(), static_cast<unsigned int>(3));
  EXPECT_EQ(basis.getNumPointsOut(), static_cast<unsigned int>(4));

  {
    ArrayType<Real, 3> vals_in(boost::extents[3][3][3]);
    ArrayType<Real, 1> vals_out(boost::extents[4]);
    //ArrayType<Real, 4> derivs_out(boost::extents[4][3]);

    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
        for (int k=0; k < 3; ++k)
          vals_in[i][j][k] = testPoly(pts_in[i], pts_in[j], pts_in[k]);

    // test interpolation functions
    basis.interpolateVals(vals_in, vals_out);
    //basis.interpolateDerivs(vals_in, derivs_out);

    for (int i=0; i < 4; ++i)
      EXPECT_FLOAT_EQ(vals_out[i], testPoly(pts_out[i][0], pts_out[i][1],
                                            pts_out[i][2]));
  }
}
