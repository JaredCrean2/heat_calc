#include "linear_system/large_matrix_factory.h"
#include "linear_system/large_matrix_dense.h"
#include "linear_system/large_matrix_petsc.h"

namespace linear_system {

std::shared_ptr<LargeMatrix> largeMatrixFactory(LargeMatrixType type, int mlocal, int nlocal, 
                                               std::shared_ptr<LargeMatrixOpts> opts, 
                                               std::shared_ptr<SparsityPattern> sparsity)
{
  switch(type)
  {
    case LargeMatrixType::Dense:
    {
      return std::make_shared<LargeMatrixDense>(mlocal, nlocal, *opts);
    }

    case LargeMatrixType::Petsc:
    {
      auto& opts_petsc = dynamic_cast<LargeMatrixOptsPetsc&>(*opts);
      return std::make_shared<LargeMatrixPetsc>(mlocal, nlocal, opts_petsc, sparsity);
    }

    default:
      throw std::runtime_error("unrecognized type");
  }
}

}