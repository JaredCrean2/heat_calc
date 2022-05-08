#ifndef SPARSITY_PATTERN_H
#define SPARSITY_PATTERN_H

#include <vector>
#include "petscksp.h"

namespace linear_system {

class SparsityPattern
{
  public:
    virtual ~SparsityPattern() = default;

    // returns vector giving number of local dofs connected to each dof
    virtual const std::vector<PetscInt>& getDiagonalCounts() = 0;

    // similar to the above, but returns the count for only the matrix diagonal + upper triangle
    virtual const std::vector<PetscInt>& getDiagonalCountsSym() = 0;

    // returns vector giving number of remote dofs connected to each dof
    virtual const std::vector<PetscInt>& getOffProcCounts() = 0;
        
    // similar to the above, but returns the count for only the upper triangle
    virtual const std::vector<PetscInt>& getOffProcCountsSym() = 0;
};

} // namespace

#endif