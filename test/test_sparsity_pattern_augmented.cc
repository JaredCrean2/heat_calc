#include "gtest/gtest.h"
#include "linear_system/sparsity_pattern.h"
#include "linear_system/sparsity_pattern_augmented.h"
#include "utils/mpi_utils.h"

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

}


TEST(SparsityPatternAugmented, Values)
{
  int mesh_dofs = 3, augmented_dofs = 2, mesh_dofs_total = commSize(MPI_COMM_WORLD) * mesh_dofs;
  bool am_i_last_rank = commRank(MPI_COMM_WORLD) == (commSize(MPI_COMM_WORLD) - 1);
  auto base_pattern = std::make_shared<SparsityPatternTest>(mesh_dofs);
  auto augmented_pattern = std::make_shared<linear_system::SparsityPatternAugmented>(base_pattern, augmented_dofs, MPI_COMM_WORLD);


  if (am_i_last_rank)
  {
    EXPECT_EQ(augmented_pattern->getNumOwnedDofs(), mesh_dofs + augmented_dofs);
    std::vector<PetscInt> expected_diagonal_counts = {1 + augmented_dofs, 
                                                      1 + augmented_dofs,
                                                      1 + augmented_dofs,
                                                      mesh_dofs + augmented_dofs,
                                                      mesh_dofs + augmented_dofs};
    EXPECT_EQ(augmented_pattern->getDiagonalCounts(),  expected_diagonal_counts);

    int num_dofs_on_other_procs = mesh_dofs * (commSize(MPI_COMM_WORLD) - 1);
    std::vector<PetscInt> expected_off_proc_counts = {2, 2, 2, num_dofs_on_other_procs, num_dofs_on_other_procs};
    EXPECT_EQ(augmented_pattern->getOffProcCounts(), expected_off_proc_counts);

    std::vector<PetscInt> ghost_dofs_to_local = base_pattern->getGhostLocalIndices();
    for (auto& val : ghost_dofs_to_local)
    { 
      val += augmented_dofs;
    }
    EXPECT_EQ(augmented_pattern->getGhostLocalIndices(), ghost_dofs_to_local);

    std::vector<PetscInt> ghost_dofs_to_global = {666, 666};
    EXPECT_EQ(augmented_pattern->getGhostGlobalIndices(), ghost_dofs_to_global);

    std::vector<PetscInt> owned_to_local_dofs = {0, 1, 2, mesh_dofs + 2, mesh_dofs + 2 + 1};
    EXPECT_EQ(augmented_pattern->getOwnedToLocalInfo(), owned_to_local_dofs);

  } else
  {
    EXPECT_EQ(augmented_pattern->getNumOwnedDofs(), mesh_dofs);
    EXPECT_EQ(augmented_pattern->getDiagonalCounts(),  std::vector<PetscInt>(3, 1));
    EXPECT_EQ(augmented_pattern->getOffProcCounts(), std::vector<PetscInt>(3, 2 + augmented_dofs));

    std::vector<PetscInt> expected_ghost_local_dofs = {mesh_dofs, mesh_dofs + 1, mesh_dofs + 2, mesh_dofs + 3};
    EXPECT_EQ(augmented_pattern->getGhostLocalIndices(), expected_ghost_local_dofs);

    std::vector<PetscInt> expected_ghost_global_dofs = {666, 666, mesh_dofs_total, mesh_dofs_total + 1};
    EXPECT_EQ(augmented_pattern->getGhostGlobalIndices(), expected_ghost_global_dofs);

    std::vector<PetscInt> expected_owned_to_local_dofs = {0, 1, 2};   
    EXPECT_EQ(augmented_pattern->getOwnedToLocalInfo(), expected_owned_to_local_dofs);
  }
}
