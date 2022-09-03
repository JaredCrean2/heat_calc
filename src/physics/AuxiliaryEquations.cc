#include "physics/AuxiliaryEquations.h"

void AuxiliaryEquations::setAuxiliaryBlockSolution(int block, const ArrayType<Real, 1>& vals)
{
  auto& solution = m_solutions->getVector(block);
  for (int i=0; i < vals.shape()[0]; ++i)
    solution[i] = vals[i];
}

ArrayType<Real, 1>& AuxiliaryEquations::getAuxiliaryBlockSolution(int block)
{
  return m_solutions->getVector(block);
}