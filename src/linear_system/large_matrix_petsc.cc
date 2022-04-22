#include "linear_system/large_matrix_petsc.h"

namespace linear_system {

  
void setPetscOptions(const LargeMatrixOptsPetsc& opts)
{
  for (auto& p : opts.petsc_opts)
  {
    std::string full_key = std::string("-") + opts.opts_prefix + p.first;
    PetscOptionsSetValue(NULL, full_key.c_str(), p.second.c_str());
  }
}

}