#include "gtest/gtest.h"
#include "linear_system/large_matrix_petsc.h"
#include "linear_system/sparsity_pattern.h"
#include "utils/initialization.h"
#include "test_helper.h"

namespace {

class SparsityPatternTest : public linear_system::SparsityPattern
{
  public:
    explicit SparsityPatternTest(int size) :
      m_local(size, size),
      m_offproc(size, size),
      m_owned_to_local_dofs(size)
    {
      for (int i=0; i < size; ++i)
        m_owned_to_local_dofs[i] = i;
    }

    const std::vector<PetscInt>& getDiagonalCounts() override { return m_local; }

    const std::vector<PetscInt>& getDiagonalCountsSym() override { return m_local; }

    const std::vector<PetscInt>& getOffProcCounts() override { return m_offproc; }
        
    const std::vector<PetscInt>& getOffProcCountsSym() override { return m_offproc; }

    const std::vector<PetscInt>& getGhostGlobalIndices() override { return m_ghost_dofs; }

    const std::vector<PetscInt>& getGhostLocalIndices() override { return m_ghost_dofs; }

    const std::vector<PetscInt>& getOwnedToLocalInfo() override { return m_owned_to_local_dofs; };

  private:
    std::vector<PetscInt> m_local;
    std::vector<PetscInt> m_offproc;
    std::vector<PetscInt> m_ghost_dofs;
    std::vector<PetscInt> m_owned_to_local_dofs;
};

linear_system::LargeMatrixOptsPetsc get_options()
{
  linear_system::LargeMatrixOptsPetsc opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;
  opts.petsc_opts["ksp_atol"] = "1e-15";
  opts.petsc_opts["ksp_rtol"] = "1e-50";
  opts.petsc_opts["ksp_monitor"] = "";

  return opts;
}

}


TEST(LargeMatrixPetsc, GeneralSolve)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat(3, 3, opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
  auto vals = make_mat(3, 3, {1, 2, 3,
                              4, 5, 6,
                              8, 8, 9});
  auto b = make_vec({3, 6, 9});
  auto x_ex = make_vec({0, 0, 1});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  std::cout << "finished setting values" << std::endl;
  mat.finishMatrixAssembly();
  std::cout << "finished matrix assembly" << std::endl;
  mat.factor();
  std::cout << "finished factoring matrix" << std::endl;
  mat.solve(b, x);
  std::cout << "finished solve" << std::endl;

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, AssembleValuesAdditive)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat(3, 3, opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
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
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, AssembleValuesIgnore)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat(3, 3, opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank, -1};
  auto vals = make_mat(4, 4, {1, 2, 3,       666,
                              4, 5, 6,       666,
                              8, 8, 9,       666,
                              666, 666, 666, 666});
  auto b = make_vec({3, 6, 9});
  auto x_ex = make_vec({0, 0, 1});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, ZeroMatrix)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat(3, 3, opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
  auto vals = make_mat(3, 3, {1, 2, 3,
                              4, 5, 6,
                              8, 8, 9});
  auto b = make_vec({3, 6, 9});
  auto x_ex = make_vec({0, 0, 1});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);


  mat.zeroMatrix();
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      vals[i][j] *= 2;

  mat.assembleValues(dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], 0.5 * x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, FactorInPlace)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  opts.factor_in_place  = true;
  opts.petsc_opts["ksp_type"] = "preonly";
  opts.petsc_opts.erase("ksp_monitor");
  opts.petsc_opts.erase("ksp_atol");
  opts.petsc_opts.erase("ksp_rtol");
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat(3, 3, opts, sparsity_pattern);


  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
  auto vals = make_mat(3, 3, {1, 2, 3,
                              4, 5, 6,
                              8, 8, 9});
  auto b = make_vec({3, 6, 9});
  auto x_ex = make_vec({0, 0, 1});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, SPD)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  
  auto opts = get_options();
  opts.is_structurally_symmetric = true;
  opts.is_value_symmetric        = true;
  opts.petsc_opts["pc_type"] = "bjacobi";
  opts.petsc_opts["pc_sub_type"] = "icc";
  opts.petsc_opts["mat_type"] = "aij";  //TODO: need to update assembly procedure to filter
                                        //      out below-diagonal entries
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat(3, 3, opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::vector<DofInt> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
  auto vals = make_mat(3, 3, {10, 2, 3,
                              2, 12, 5,
                              3, 5, 20});
  auto b = make_vec({1, 2, 3});
  auto x_ex = make_vec({0.043027, 0.111276, 0.115727});
  ArrayType<Real, 1> x(boost::extents[3]);

  mat.assembleValues(dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-5);
}