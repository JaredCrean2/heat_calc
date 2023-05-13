#include "physics/AuxiliaryEquations.h"

void splitVector(const ArrayType<Real, 1>& combined_vec, DiscVectorPtr sol_vec,
                 AuxiliaryEquationsStoragePtr sol_aux)
{
  std::cout << "splitting vector" << std::endl;
  if (!sol_vec->isVectorCurrent())
    sol_vec->set(0);

  auto& sol_vec_vec = sol_vec->getVector();
  int dof=0;
  for (int i=0; i < sol_vec_vec.shape()[0]; ++i)
    sol_vec_vec[i] = combined_vec[dof++];

  std::cout << "assigned to dofs 0 - " << dof << std::endl;
  for (int block=1; block < sol_aux->getNumBlocks(); ++block)
  {
    auto& vec = sol_aux->getVector(block);
    for (int i=0; i < vec.shape()[0]; ++i)
      vec[i] = combined_vec[dof++];
  }

  sol_vec->markVectorModified();
}

void combineVector(DiscVectorPtr res_vec, AuxiliaryEquationsStoragePtr res_aux,
                   ArrayType<Real, 1>& res_combined)
{
  if (!res_vec->isVectorCurrent())
    res_vec->syncArrayToVector();

  auto& res_vec_vec = res_vec->getVector();
  int dof = 0;
  for (int i=0; i < res_vec_vec.shape()[0]; ++i)
    res_combined[dof++] = res_vec_vec[i];

  for (int block=1; block < res_aux->getNumBlocks(); ++block)
  {
    const auto& vec = res_aux->getVector(block);
    for (int i=0; i < vec.shape()[0]; ++i)
      res_combined[dof++] = vec[i];
  }
}