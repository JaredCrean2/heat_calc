#include "gtest/gtest.h"
#include "linear_system/augmented_assembler.h"
#include "linear_system/sparsity_pattern.h"
#include "linear_system/sparsity_pattern_augmented.h"
#include "linear_system/large_matrix_petsc.h"
#include "linear_system/sparsity_pattern_mesh.h"
#include "utils/mpi_utils.h"
#include "mesh_helper.h"

namespace {

class SparsityPatternTest : public linear_system::SparsityPattern
{
  public:
    SparsityPatternTest(int num_owned_dofs) :
      m_num_owned_dofs(num_owned_dofs)
    {
      setup();
    }

    PetscInt getNumOwnedDofs() const override { return m_num_owned_dofs; }

    // returns vector giving number of local dofs connected to each dof
    const std::vector<PetscInt>& getDiagonalCounts() override { return m_onproc_dofs; }

    // similar to the above, but returns the count for only the matrix diagonal + upper triangle
    const std::vector<PetscInt>& getDiagonalCountsSym() override 
    { 
      throw std::runtime_error("not supported");

      return m_onproc_dofs;
    }

    // returns vector giving number of remote dofs connected to each dof
    const std::vector<PetscInt>& getOffProcCounts() override
    {
      return m_remote_dofs;
    }
        
    // similar to the above, but returns the count for only the upper triangle
    const std::vector<PetscInt>& getOffProcCountsSym() override
    {
      throw std::runtime_error("not supported");
    }

    const std::vector<PetscInt>& getGhostGlobalIndices() override
    {
      return m_ghost_global_dofs;
    }

    const std::vector<PetscInt>& getGhostLocalIndices() override
    {
      return m_ghost_onproc_dofs;
    }

    const std::vector<PetscInt>& getOwnedToLocalInfo() override
    {
      return m_owned_dof_to_local;
    }

  private:
    void setup()
    {
      for (int i=0; i < m_num_owned_dofs; ++i)
      {
        m_onproc_dofs.push_back(1);
        m_remote_dofs.push_back(2);
        m_owned_dof_to_local.push_back(i);
      }

      for (int i=0; i < 2; ++i)
      {
        m_ghost_global_dofs.push_back(666);
        m_ghost_onproc_dofs.push_back(m_num_owned_dofs + i);
      }
    }

    int m_num_owned_dofs;
    std::vector<PetscInt> m_onproc_dofs;
    std::vector<PetscInt> m_remote_dofs;
    std::vector<PetscInt> m_ghost_global_dofs;
    std::vector<PetscInt> m_ghost_onproc_dofs;
    std::vector<PetscInt> m_owned_dof_to_local;       
};


class AugmentedAssemblerTester : public StandardDiscSetup,
                                 public ::testing::Test
{
  public:
    AugmentedAssemblerTester()
    {
      setup(3, 1, getStandardMeshSpec());
    }
};

}


TEST_F(AugmentedAssemblerTester, Values)
{
  int num_augmented_dofs = 2;
  bool am_i_last_rank = commRank(MPI_COMM_WORLD) == (commSize(MPI_COMM_WORLD) - 1);
  auto base_pattern = std::make_shared<linear_system::SparsityPatternMesh>(mesh);
  auto augmented_pattern = std::make_shared<linear_system::SparsityPatternAugmented>(base_pattern, num_augmented_dofs, MPI_COMM_WORLD);
  linear_system::LargeMatrixOptsPetsc opts;
  auto mat = std::make_shared<linear_system::LargeMatrixPetsc>(opts, augmented_pattern);
  auto assembler = std::make_shared<linear_system::AugmentedAssembler>(disc, mat, num_augmented_dofs);

  std::vector<DofInt> local_dofs_to_global;
  mesh->getLocalToGlobalDofs(local_dofs_to_global);

  PetscInt num_owned_dofs = mesh->getNumOwnedDofs() + (am_i_last_rank ? num_augmented_dofs : 0);
  PetscInt num_dofs_total;
  MPI_Allreduce(&num_owned_dofs, &num_dofs_total, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);  

  for (int augmented_row=0; augmented_row < num_augmented_dofs; ++augmented_row)
  {
    std::cout << "number of owned dofs = " << mesh->getNumOwnedDofs() << std::endl;
    std::cout << "size of localToGlobal dofs = " << local_dofs_to_global.size() << std::endl;
    std::vector<DofInt> dofs(mesh->getNumOwnedDofs());
    std::vector<Real> vals(mesh->getNumOwnedDofs());
    for (size_t i=0; i < local_dofs_to_global.size(); ++i)
    {
      dofs[i] = i;
      vals[i] = local_dofs_to_global[i] + augmented_row;
    }

    std::cout << "assembling augmented row of size " << dofs.size() << std::endl;

    assembler->assembleValuesRow(augmented_row, dofs, vals);
  }

  {
    std::vector<DofInt> augmented_rows = {0, 1}, augmented_cols = {0, 1};
    ArrayType<Real, 2> vals(boost::extents[num_augmented_dofs][num_augmented_dofs]);
    vals[0][0] = num_dofs_total - num_augmented_dofs;     vals[0][1] = num_dofs_total - num_augmented_dofs + 1;
    vals[1][0] = num_dofs_total - num_augmented_dofs + 1; vals[1][1] = num_dofs_total - num_augmented_dofs + 2;

    assembler->assembleAugmentedValuesDiag(augmented_rows, augmented_cols, vals);
  }

  assembler->startAssembly();
  assembler->finishAssembly();
  mat->finishMatrixAssembly();

  mat->printToStdout();


  int num_dofs_local = mesh->getNumDofs();
  if (am_i_last_rank)
    num_dofs_local += num_augmented_dofs;

  ArrayType<Real, 1> x(boost::extents[num_dofs_local]), b(boost::extents[num_dofs_local]);

  for (PetscInt i=0; i < num_dofs_total; ++i)
  {
    std::cout << "\ni = " << i << std::endl;
    std::fill(x.begin(), x.end(), 0);
    std::fill(b.begin(), b.end(), 0);

    auto it = std::find(local_dofs_to_global.begin(), local_dofs_to_global.begin() + mesh->getNumOwnedDofs(), i);
    if (it != local_dofs_to_global.end())
    {
      std::cout << "x[" << *it << "] = 1" << std::endl;
      x[*it] = 1;
    }

    mat->matVec(x, b);

    for (size_t i=0; i < num_dofs_local; ++i)
      if (b[i] != 0)
        std::cout << "b[" << i << "] = " << b[i] << std::endl;

    for (int augmented_dof=0; augmented_dof < num_augmented_dofs; ++augmented_dof)
      if (am_i_last_rank)
      {
        std::cout << "checking b idx " << num_dofs_local - num_augmented_dofs + augmented_dof << std::endl;
        EXPECT_NEAR(b[num_dofs_local - num_augmented_dofs + augmented_dof], i + augmented_dof, 1e-12);
      }
  }
}