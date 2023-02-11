#ifndef LINEAR_SYSTEM_LARGE_MATRIX_FACTORY_H
#define LINEAR_SYSTEM_LARGE_MATRIX_FACTORY_H

#include "large_matrix.h"
#include "linear_system/sparsity_pattern.h"


namespace linear_system
{

enum class LargeMatrixType
{
  Dense,
  Petsc,
  Unknown
};


std::shared_ptr<LargeMatrix> largeMatrixFactory(LargeMatrixType type,
                                               std::shared_ptr<LargeMatrixOpts> opts, 
                                               std::shared_ptr<SparsityPattern> sparsity);



}  // namespace

#endif