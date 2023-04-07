#ifndef SPARSITY_PATTERN_AUGMENTED_H
#define SPARSITY_PATTERN_AUGMENTED_H

#include "linear_system/sparsity_pattern.h"

namespace linear_system {


// adds n dense rows and columns to the bottom/right of the matrix
// The local vector dof layout is [local dofs from finite element mesh, augmented dofs]
class SparsityPatternAugmented : public SparsityPattern
{
  public:
    SparsityPatternAugmented(std::shared_ptr<SparsityPattern> base_pattern, int num_augmented_rows, MPI_Comm comm);

    PetscInt getNumOwnedDofs() const override;

    // returns vector giving number of local dofs connected to each dof
    const std::vector<PetscInt>& getDiagonalCounts() override;

    // similar to the above, but returns the count for only the matrix diagonal + upper triangle
    const std::vector<PetscInt>& getDiagonalCountsSym() override;

    // returns vector giving number of remote dofs connected to each dof
    const std::vector<PetscInt>& getOffProcCounts() override;
        
    // similar to the above, but returns the count for only the upper triangle
    const std::vector<PetscInt>& getOffProcCountsSym() override;

    const std::vector<PetscInt>& getGhostGlobalIndices() override;

    const std::vector<PetscInt>& getGhostLocalIndices() override;

    const std::vector<PetscInt>& getOwnedToLocalInfo() override;

  private:
    std::shared_ptr<SparsityPattern> m_base_pattern;
    int m_num_augmented_rows;
    MPI_Comm m_comm;
    bool m_am_i_last_rank;

    std::vector<PetscInt> m_onproc_dofs;
    std::vector<PetscInt> m_remote_dofs;
    std::vector<PetscInt> m_owned_dof_to_local;    
};

}  // namespace

#endif