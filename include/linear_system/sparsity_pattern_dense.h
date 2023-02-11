#ifndef SPARSITY_PATTERN_DENSE_H
#define SPARSITY_PATTERN_DENSE_H

#include "linear_system/sparsity_pattern.h"
#include "utils/error_handling.h"

namespace linear_system {

// SparsityPattern object to be used for a dense, single process matrix
class SparsityPatternDense : public SparsityPattern
{
  public:
    explicit SparsityPatternDense(int ndof) :
      m_ndof(ndof)
    {}

    PetscInt getNumOwnedDofs() const override { return m_ndof; }

    // returns vector giving number of local dofs connected to each dof
    const std::vector<PetscInt>& getDiagonalCounts() override
    {
      assertAlways(false, "not supported");
      return m_empty_vector;
    }

    // similar to the above, but returns the count for only the matrix diagonal + upper triangle
    const std::vector<PetscInt>& getDiagonalCountsSym() override
    {
      assertAlways(false, "not supported");
      return m_empty_vector;
    }

    // returns vector giving number of remote dofs connected to each dof
    const std::vector<PetscInt>& getOffProcCounts() override
    {
      assertAlways(false, "not supported");
      return m_empty_vector;
    }

    // similar to the above, but returns the count for only the upper triangle
    const std::vector<PetscInt>& getOffProcCountsSym() override
    {
      assertAlways(false, "not supported");
      return m_empty_vector;
    }

    const std::vector<PetscInt>& getGhostGlobalIndices() override
    {
      assertAlways(false, "not supported");
      return m_empty_vector;
    }

    const std::vector<PetscInt>& getGhostLocalIndices() override
    {
      assertAlways(false, "not supported");
      return m_empty_vector;
    }

    const std::vector<PetscInt>& getOwnedToLocalInfo() override
    {
      assertAlways(false, "not supported");
      return m_empty_vector;
    }    

  private:
    int m_ndof;
    std::vector<PetscInt> m_empty_vector;
};

} // namespace

#endif