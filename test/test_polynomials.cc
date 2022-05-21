#include <cmath>

#include "gtest/gtest.h"
#include "test_helper.h"
#include "utils/lagrange.h"
#include "utils/lagrange2d.h"
#include "utils/legendre.h"
#include "utils/quadrature.h"

#include "mesh/mesh.h"
#include "physics/heat/basis_vals.h"

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
  SERIAL_ONLY();

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
  SERIAL_ONLY();

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
  SERIAL_ONLY();

  // test using orthogonality property of Legendre polnomials

  LegendrePoly legendre;
  std::vector<int> exactness = {1, 3, 5, 7, 9};

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

    {
      // test setDomain
      quad.setDomain({1.0, 5.0});
      Real val = 0;
      for (int i=0; i < quad.getNumPoints(); ++i)
      {
        Real integrand = std::pow(quad.getPoint(i), degree_max);
        val += quad.getWeight(i) * integrand;
      }

      Real val_ex = (std::pow(5.0, degree_max + 1) -
                     std::pow(1.0, degree_max + 1)) / (degree_max + 1);
      EXPECT_FLOAT_EQ(val, val_ex);
    }
  }

}


TEST(Polynomials, LagrangeMemoizer)
{
  SERIAL_ONLY();

  std::vector<Real> pts_in1 = {-1, 0, 1};
  std::vector<Real> pts_out1 = {0.5, 0, 0.5};

  std::vector<Real> pts_in2 = {-0.75, 0, 0.75};
  std::vector<Real> pts_out2 = {0.55, 0, 0.55};

  const auto& vals1 = lagrange_memoizer.getValues(pts_in1, pts_out1);
  const auto& vals1a = lagrange_memoizer.getValues(pts_in1, pts_out1);

  EXPECT_EQ(&vals1, &vals1a);

  const auto& vals2 = lagrange_memoizer.getValues(pts_in2, pts_out2);
  const auto& vals2a = lagrange_memoizer.getValues(pts_in2, pts_out2);

  EXPECT_EQ(&vals2, &vals2a);
  EXPECT_NE(&vals1, &vals2);

  // test mixed
  const auto& vals3 = lagrange_memoizer.getValues(pts_in1, pts_out2);
  const auto& vals3a = lagrange_memoizer.getValues(pts_in1, pts_out2);

  EXPECT_EQ(&vals3, &vals3a);
  EXPECT_NE(&vals1, &vals3);
  EXPECT_NE(&vals2, &vals3);
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
  SERIAL_ONLY();

  std::vector<Real> pts_in = {-1.0, 0.0, 1.0};
  std::vector<Real> pts_out = {-0.75, -0.25, 0.25, 0.75};
  ArrayType<Real, 2> pts_nontp(boost::extents[64][3]);
  ArrayType<LocalIndex, 3> nodemap_in(boost::extents[3][3][3]);
  ArrayType<LocalIndex, 3> nodemap_out(boost::extents[4][4][4]);

  {
    ArrayType<Real, 3> vals_in_tp(boost::extents[3][3][3]);
    ArrayType<Real, 3> vals_out_tp(boost::extents[4][4][4]);
    ArrayType<Real, 4> derivs_out_tp(boost::extents[4][4][4][3]);

    ArrayType<Real, 1> vals_in_flat(boost::extents[nodemap_in.num_elements()]);
    ArrayType<Real, 1> vals_out_flat(boost::extents[nodemap_out.num_elements()]);
    ArrayType<Real, 2> derivs_out_flat(boost::extents[nodemap_out.num_elements()][3]);

    ArrayType<Real, 1> vals_out_nontp(boost::extents[64]);
    ArrayType<Real, 2> derivs_out_nontp(boost::extents[64][3]);


    int idx = 0;
    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
        for (int k=0; k < 3; ++k)
        {
          nodemap_in[i][j][k] = idx++;
          vals_in_tp[i][j][k] = testPoly(pts_in[i], pts_in[j], pts_in[k]);
          vals_in_flat[nodemap_in[i][j][k]] = vals_in_tp[i][j][k];
        }

    idx = 0;
    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
        for (int k=0; k < 4; ++k)
        {
          nodemap_out[i][j][k] = idx;
          pts_nontp[idx][0] = pts_out[i];
          pts_nontp[idx][1] = pts_out[j];
          pts_nontp[idx][2] = pts_out[k];
          ++idx;
        }


    LagrangeEvaluatorTPToTP tp_to_tp(pts_in, pts_out);
    LagrangeEvaluatorTPFlatToTP flat_to_tp(pts_in, pts_out, nodemap_in);
    LagrangeEvaluatorTPToTPFlat tp_to_flat(pts_in, pts_out, nodemap_out);
    LagrangeEvaluatorTPFlatToTPFlat flat_to_flat(pts_in, pts_out, nodemap_in, nodemap_out);
    LagrangeEvaluatorTPToNonTP tp_to_nontp(pts_in, pts_nontp);
    LagrangeEvaluatorTPFlatToNonTP flat_to_nontp(pts_in, pts_nontp, nodemap_in);

    // test sizes
    EXPECT_EQ(static_cast<Index>(tp_to_tp.getNumPointsIn()), 3);
    EXPECT_EQ(static_cast<Index>(tp_to_tp.getNumPointsOut()), 4);

    EXPECT_EQ(static_cast<Index>(flat_to_tp.getNumPointsIn()), 27);
    EXPECT_EQ(static_cast<Index>(flat_to_tp.getNumPointsOut()), 4);

    EXPECT_EQ(static_cast<Index>(tp_to_flat.getNumPointsIn()), 3);
    EXPECT_EQ(static_cast<Index>(tp_to_flat.getNumPointsOut()), 64);

    EXPECT_EQ(static_cast<Index>(flat_to_flat.getNumPointsIn()), 27);
    EXPECT_EQ(static_cast<Index>(flat_to_flat.getNumPointsOut()), 64);

    EXPECT_EQ(static_cast<Index>(tp_to_nontp.getNumPointsIn()), 3);
    EXPECT_EQ(static_cast<Index>(tp_to_nontp.getNumPointsOut()), 64);

    EXPECT_EQ(static_cast<Index>(flat_to_nontp.getNumPointsIn()), 27);
    EXPECT_EQ(static_cast<Index>(flat_to_nontp.getNumPointsOut()), 64);

    // test interpolating from tensor product
    tp_to_tp.interpolateVals(vals_in_tp, vals_out_tp);
    tp_to_tp.interpolateDerivs(vals_in_tp, derivs_out_tp);

    tp_to_flat.interpolateVals(vals_in_tp, vals_out_flat);
    tp_to_flat.interpolateDerivs(vals_in_tp, derivs_out_flat);

    tp_to_nontp.interpolateVals(vals_in_tp, vals_out_nontp);
    tp_to_nontp.interpolateDerivs(vals_in_tp, derivs_out_nontp);



    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
        for (int k=0; k < 4; ++k)
        {
          EXPECT_FLOAT_EQ(vals_out_tp[i][j][k], testPoly(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_tp[i][j][k][0], testPoly_dx(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_tp[i][j][k][1], testPoly_dy(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_tp[i][j][k][2], testPoly_dz(pts_out[i], pts_out[j], pts_out[k]));

          auto node = nodemap_out[i][j][k];
          EXPECT_FLOAT_EQ(vals_out_flat[node], testPoly(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_flat[node][0], testPoly_dx(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_flat[node][1], testPoly_dy(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_flat[node][2], testPoly_dz(pts_out[i], pts_out[j], pts_out[k]));

          EXPECT_FLOAT_EQ(vals_out_nontp[node], testPoly(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_nontp[node][0], testPoly_dx(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_nontp[node][1], testPoly_dy(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_nontp[node][2], testPoly_dz(pts_out[i], pts_out[j], pts_out[k]));
        }

    // test interpolating from flat
    flat_to_tp.interpolateVals(vals_in_flat, vals_out_tp);
    flat_to_tp.interpolateDerivs(vals_in_flat, derivs_out_tp);

    flat_to_flat.interpolateVals(vals_in_flat, vals_out_flat);
    flat_to_flat.interpolateDerivs(vals_in_flat, derivs_out_flat);

    flat_to_nontp.interpolateVals(vals_in_flat, vals_out_nontp);
    flat_to_nontp.interpolateDerivs(vals_in_flat, derivs_out_nontp);


    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
        for (int k=0; k < 4; ++k)
        {
          EXPECT_FLOAT_EQ(vals_out_tp[i][j][k], testPoly(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_tp[i][j][k][0], testPoly_dx(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_tp[i][j][k][1], testPoly_dy(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_tp[i][j][k][2], testPoly_dz(pts_out[i], pts_out[j], pts_out[k]));

          auto node = nodemap_out[i][j][k];
          EXPECT_FLOAT_EQ(vals_out_flat[node], testPoly(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_flat[node][0], testPoly_dx(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_flat[node][1], testPoly_dy(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_flat[node][2], testPoly_dz(pts_out[i], pts_out[j], pts_out[k]));

          EXPECT_FLOAT_EQ(vals_out_nontp[node], testPoly(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_nontp[node][0], testPoly_dx(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_nontp[node][1], testPoly_dy(pts_out[i], pts_out[j], pts_out[k]));
          EXPECT_FLOAT_EQ(derivs_out_nontp[node][2], testPoly_dz(pts_out[i], pts_out[j], pts_out[k]));
        }


/*
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
                val += vals_in_tp[i_in][j_in][k_in] *
                  tp_to_tp.getInterpolantValue(i_in, j_in, k_in, i_out, j_out, k_out);
                val_x += vals_in_tp[i_in][j_in][k_in] * 
                  tp_to_tp.getInterpolantDx(i_in, j_in, k_in, i_out, j_out, k_out);
                val_y += vals_in_tp[i_in][j_in][k_in] * 
                  tp_to_tp.getInterpolantDy(i_in, j_in, k_in, i_out, j_out, k_out);
                val_z += vals_in_tp[i_in][j_in][k_in] * 
                  tp_to_tp.getInterpolantDz(i_in, j_in, k_in, i_out, j_out, k_out);
              }

          EXPECT_FLOAT_EQ(val, testPoly(pts_out[i_out], pts_out[j_out], pts_out[k_out]));
          EXPECT_FLOAT_EQ(val_x, testPoly_dx(pts_out[i_out], pts_out[j_out], pts_out[k_out]));
          EXPECT_FLOAT_EQ(val_y, testPoly_dy(pts_out[i_out], pts_out[j_out], pts_out[k_out]));
          EXPECT_FLOAT_EQ(val_z, testPoly_dz(pts_out[i_out], pts_out[j_out], pts_out[k_out]));
        }
     
*/
  }

}

TEST(Polynomials, LagrangeTP2D)
{
  SERIAL_ONLY();

  std::vector<Real> pts_in = {-1.0, 0.0, 1.0};
  std::vector<Real> pts_out = {-0.75, -0.25, 0.25, 0.75};
  ArrayType<Real, 2> pts_nontp(boost::extents[16][3]);
  ArrayType<LocalIndex, 2> nodemap_in(boost::extents[3][3]);
  ArrayType<LocalIndex, 2> nodemap_out(boost::extents[4][4]);

  {
    ArrayType<Real, 2> vals_in_tp(boost::extents[3][3]);
    ArrayType<Real, 2> vals_out_tp(boost::extents[4][4]);
    ArrayType<Real, 3> derivs_out_tp(boost::extents[4][4][2]);

    ArrayType<Real, 1> vals_in_flat(boost::extents[nodemap_in.num_elements()]);
    ArrayType<Real, 1> vals_out_flat(boost::extents[nodemap_out.num_elements()]);
    ArrayType<Real, 2> derivs_out_flat(boost::extents[nodemap_out.num_elements()][2]);

    ArrayType<Real, 1> vals_out_nontp(boost::extents[16]);
    ArrayType<Real, 2> derivs_out_nontp(boost::extents[16][2]);


    int idx = 0;
    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
      {
        nodemap_in[i][j] = idx++;
        vals_in_tp[i][j] = testPoly(pts_in[i], pts_in[j], 1);
        vals_in_flat[nodemap_in[i][j]] = vals_in_tp[i][j];
      }

    idx = 0;
    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
      {
        nodemap_out[i][j] = idx;
        pts_nontp[idx][0] = pts_out[i];
        pts_nontp[idx][1] = pts_out[j];
        ++idx;
      }


    LagrangeEvaluator2DTPToTP tp_to_tp(pts_in, pts_out);
    LagrangeEvaluator2DTPFlatToTP flat_to_tp(pts_in, pts_out, nodemap_in);
    LagrangeEvaluator2DTPToTPFlat tp_to_flat(pts_in, pts_out, nodemap_out);
    LagrangeEvaluator2DTPFlatToTPFlat flat_to_flat(pts_in, pts_out, nodemap_in, nodemap_out);
    LagrangeEvaluator2DTPToNonTP tp_to_nontp(pts_in, pts_nontp);
    LagrangeEvaluator2DTPFlatToNonTP flat_to_nontp(pts_in, pts_nontp, nodemap_in);

    // test sizes
    EXPECT_EQ(static_cast<Index>(tp_to_tp.getNumPointsIn()), 3);
    EXPECT_EQ(static_cast<Index>(tp_to_tp.getNumPointsOut()), 4);

    EXPECT_EQ(static_cast<Index>(flat_to_tp.getNumPointsIn()), 9);
    EXPECT_EQ(static_cast<Index>(flat_to_tp.getNumPointsOut()), 4);

    EXPECT_EQ(static_cast<Index>(tp_to_flat.getNumPointsIn()), 3);
    EXPECT_EQ(static_cast<Index>(tp_to_flat.getNumPointsOut()), 16);

    EXPECT_EQ(static_cast<Index>(flat_to_flat.getNumPointsIn()), 9);
    EXPECT_EQ(static_cast<Index>(flat_to_flat.getNumPointsOut()), 16);

    EXPECT_EQ(static_cast<Index>(tp_to_nontp.getNumPointsIn()), 3);
    EXPECT_EQ(static_cast<Index>(tp_to_nontp.getNumPointsOut()), 16);

    EXPECT_EQ(static_cast<Index>(flat_to_nontp.getNumPointsIn()), 9);
    EXPECT_EQ(static_cast<Index>(flat_to_nontp.getNumPointsOut()), 16);

    // test interpolating from tensor product
    tp_to_tp.interpolateVals(vals_in_tp, vals_out_tp);
    tp_to_tp.interpolateDerivs(vals_in_tp, derivs_out_tp);

    tp_to_flat.interpolateVals(vals_in_tp, vals_out_flat);
    tp_to_flat.interpolateDerivs(vals_in_tp, derivs_out_flat);

    tp_to_nontp.interpolateVals(vals_in_tp, vals_out_nontp);
    tp_to_nontp.interpolateDerivs(vals_in_tp, derivs_out_nontp);



    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
      {
        EXPECT_FLOAT_EQ(vals_out_tp[i][j], testPoly(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_tp[i][j][0], testPoly_dx(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_tp[i][j][1], testPoly_dy(pts_out[i], pts_out[j], 1));

        auto node = nodemap_out[i][j];
        EXPECT_FLOAT_EQ(vals_out_flat[node], testPoly(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_flat[node][0], testPoly_dx(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_flat[node][1], testPoly_dy(pts_out[i], pts_out[j], 1));

        EXPECT_FLOAT_EQ(vals_out_nontp[node], testPoly(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_nontp[node][0], testPoly_dx(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_nontp[node][1], testPoly_dy(pts_out[i], pts_out[j], 1));
      }

    // test interpolating from flat
    flat_to_tp.interpolateVals(vals_in_flat, vals_out_tp);
    flat_to_tp.interpolateDerivs(vals_in_flat, derivs_out_tp);

    flat_to_flat.interpolateVals(vals_in_flat, vals_out_flat);
    flat_to_flat.interpolateDerivs(vals_in_flat, derivs_out_flat);

    flat_to_nontp.interpolateVals(vals_in_flat, vals_out_nontp);
    flat_to_nontp.interpolateDerivs(vals_in_flat, derivs_out_nontp);


    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
      {
        EXPECT_FLOAT_EQ(vals_out_tp[i][j], testPoly(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_tp[i][j][0], testPoly_dx(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_tp[i][j][1], testPoly_dy(pts_out[i], pts_out[j], 1));

        auto node = nodemap_out[i][j];
        EXPECT_FLOAT_EQ(vals_out_flat[node], testPoly(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_flat[node][0], testPoly_dx(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_flat[node][1], testPoly_dy(pts_out[i], pts_out[j], 1));

        EXPECT_FLOAT_EQ(vals_out_nontp[node], testPoly(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_nontp[node][0], testPoly_dx(pts_out[i], pts_out[j], 1));
        EXPECT_FLOAT_EQ(derivs_out_nontp[node][1], testPoly_dy(pts_out[i], pts_out[j], 1));
      }

  }

}

TEST(Polynomials, HeatBasisVals)
{
  SERIAL_ONLY();

  std::vector<Real> pts_in = {-1.0, 0.0, 1.0};
  std::vector<Real> pts_out = {-0.75, -0.25, 0.25, 0.75};
  Mesh::TensorProductMapper mapper_in(pts_in), mapper_out(pts_out);
  Heat::BasisVals basis_vals(mapper_in, mapper_out);

  auto& nodemap_in  = mapper_in.getNodemap();
  auto& nodemap_out = mapper_out.getNodemap();
  int npts_in       = nodemap_in.num_elements();
  int npts_out      = nodemap_out.num_elements();
  ArrayType<Real, 1> vals_in_flat(boost::extents[npts_in]);
  ArrayType<Real, 1> vals_out_flat(boost::extents[npts_out]);
  ArrayType<Real, 2> derivs_out_flat(boost::extents[npts_out][3]);

  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      for (int k=0; k < 3; ++k)
        vals_in_flat[nodemap_in[i][j][k]] = testPoly(pts_in[i], pts_in[j], pts_in[k]);

  for (int i=0; i < npts_out; ++i)
  {
    vals_out_flat[i] = 0;
    for (int d=0; d < 3; ++d)
      derivs_out_flat[i][d] = 0;
    for (int j=0; j < npts_in; ++j)
    {
      vals_out_flat[i] += basis_vals.getValue(j, i) * vals_in_flat[j];

      Real derivs[3];
      basis_vals.getDerivs(j, i, derivs);
      for (int d=0; d < 3; ++d)
        derivs_out_flat[i][d] += derivs[d] * vals_in_flat[j];
    }
  
  }

  for (int i=0; i < 4; ++i)
    for (int j=0; j < 4; ++j)
      for (int k=0; k < 4; ++k)
      {
        EXPECT_FLOAT_EQ(vals_out_flat[nodemap_out[i][j][k]], testPoly(pts_out[i], pts_out[j], pts_out[k]));
        EXPECT_FLOAT_EQ(derivs_out_flat[nodemap_out[i][j][k]][0], testPoly_dx(pts_out[i], pts_out[j], pts_out[k]));
        EXPECT_FLOAT_EQ(derivs_out_flat[nodemap_out[i][j][k]][1], testPoly_dy(pts_out[i], pts_out[j], pts_out[k]));
        EXPECT_FLOAT_EQ(derivs_out_flat[nodemap_out[i][j][k]][2], testPoly_dz(pts_out[i], pts_out[j], pts_out[k]));
      }
        
}
