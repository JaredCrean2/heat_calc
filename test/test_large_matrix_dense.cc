
#include "gtest/gtest.h"
#include "linear_system/large_matrix_dense.h"
#include "test_helper.h"


TEST(LargeMatrixDense, GeneralSolve)
{
  SERIAL_ONLY();

  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat(3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0, 1, 2};
  auto vals = make_mat(3, 3, {1, 2, 3,
                              4, 5, 6,
                              8, 8, 9});
  auto b = make_vec({3, 6, 9});
  auto x_ex = make_vec({0, 0, 1});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}

TEST(LargeMatrixDense, AssembleValuesAdditive)
{
  SERIAL_ONLY();

  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat(3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0, 1, 2};
  auto vals = make_mat(3, 3, {1, 2, 3,
                              4, 5, 3,
                              4, 4, 9});
  auto vals2 = make_mat(3, 3, {0, 0, 0,
                               0, 0, 3,
                               4, 4, 0});
  auto b = make_vec({3, 6, 9});
  auto x_ex = make_vec({0, 0, 1});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.assembleValues(dofs, vals2);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixDense, AssembleValuesIgnore)
{
  SERIAL_ONLY();

  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat(3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0, 1, 2, -1};
  auto vals = make_mat(4, 4, {1, 2, 3,       666,
                              4, 5, 6,       666,
                              8, 8, 9,       666,
                              666, 666, 666, 666});
  auto b = make_vec({3, 6, 9});
  auto x_ex = make_vec({0, 0, 1});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixDense, ZeroMatrix)
{
  SERIAL_ONLY();

  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat(3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0, 1, 2};
  auto vals = make_mat(3, 3, {1, 2, 3,
                              4, 5, 6,
                              8, 8, 9});
  auto b = make_vec({3, 6, 9});
  auto x_ex = make_vec({0, 0, 1});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);


  mat.zeroMatrix();
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      vals[i][j] *= 2;

  mat.assembleValues(dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], 0.5 * x_ex[i], 1e-13);
}


TEST(LargeMatrixDense, FactorInPlace)
{
  SERIAL_ONLY();

  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = true;

  linear_system::LargeMatrixDense mat(3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0, 1, 2};
  auto vals = make_mat(3, 3, {1, 2, 3,
                              4, 5, 6,
                              8, 8, 9});
  auto b = make_vec({3, 6, 9});
  auto x_ex = make_vec({0, 0, 1});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}

TEST(LargeMatrixDense, SPD)
{
  SERIAL_ONLY();

  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = true;
  opts.is_value_symmetric        = true;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat(3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0, 1, 2};
  auto vals = make_mat(3, 3, {10, 2, 3,
                              2, 12, 5,
                              3, 5, 20});
  auto b = make_vec({1, 2, 3});
  auto x_ex = make_vec({0.043027, 0.111276, 0.115727});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-5);
}


TEST(LargeMatrixDense, SPDInPlace)
{
  SERIAL_ONLY();
  
  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = true;
  opts.is_value_symmetric        = true;
  opts.factor_in_place           = true;

  linear_system::LargeMatrixDense mat(3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);
  
  std::vector<DofInt> dofs{0, 1, 2};
  auto vals = make_mat(3, 3, {10, 2, 3,
                              2, 12, 5,
                              3, 5, 20});
  auto b = make_vec({1, 2, 3});
  auto x_ex = make_vec({0.043027, 0.111276, 0.115727});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-5);
}