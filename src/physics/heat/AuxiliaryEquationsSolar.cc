#include "physics/heat/AuxiliaryEquationsSolar.h"

namespace Heat {

void AuxiliaryEquationsSolar::setAuxiliaryBlockSolution(int block, const ArrayType<Real, 1>& vals)
{
  auto& solution = m_solutions->getVector(block+1);
  for (int i=0; i < vals.shape()[0]; ++i)
    solution[i] = vals[i];
}

ArrayType<Real, 1>& AuxiliaryEquationsSolar::getAuxiliaryBlockSolution(int block)
{
  return m_solutions->getVector(block+1);
}

}