#include "physics/heat/AuxiliaryEquationsSolar.h"

namespace Heat {

void AuxiliaryEquationsSolar::computeAuxiliaryRhs(int block, DiscVectorPtr u_vec, Real t, ArrayType<Real, 1>& rhs)
{
  Real interior_temp = getAuxiliaryBlockSolution(block)[0];
  rhs[0] = m_air_temp->computeNetFlux(u_vec, interior_temp, t);
}

void AuxiliaryEquationsSolar::computeAuxiliaryMassMatrix(int block, Real t, linear_system::SimpleAssemblerPtr assembler)
{
  std::vector<DofInt> dofs = {0};
  ArrayType<Real, 2> mat(boost::extents[1][1]);
  mat[0][0] = m_air_temp->getThermalMass();
  assembler->assembleEntry(dofs, mat);
}

void AuxiliaryEquationsSolar::multiplyAuxiliaryMassMatrix(int block, Real t,
                                const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
{
  Real val = m_air_temp->getThermalMass();
  b[0] = x[0] * val;
}

void AuxiliaryEquationsSolar::computeAuxiliaryJacobian(int block, DiscVectorPtr u_vec, Real t, 
                                                       linear_system::SimpleAssemblerPtr mat)
{
  Real interior_temp = getAuxiliaryBlockSolution(0)[0];
  Real val = m_air_temp->computeNetFluxJacobian(u_vec, interior_temp, t);

  std::cout << "val = " << val << std::endl;
  std::vector<DofInt> dofs = {0};
  ArrayType<Real, 2> vals(boost::extents[1][1]);
  vals[0][0] = val;
  mat->assembleEntry(dofs, vals);
}

void AuxiliaryEquationsSolar::computeAuxiliaryJacobianVectorProduct(int iblock, int jblock, 
      DiscVectorPtr u_vec, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
{
  Real t_interior = getAuxiliaryBlockSolution(0)[0];
  auto u_bar = makeDiscVector(m_heat_eqn.getDiscretization());  //TODO: cache this
  u_bar->set(0);
  m_air_temp->computeNetFlux_rev(u_vec, t_interior, t, u_bar, 1);

  if (!u_bar->isVectorCurrent())
    u_bar->syncArrayToVector();

  auto& u_bar_vec = u_bar->getVector();
  b[0] = 0;
  for (size_t i=0; i < u_bar_vec.shape()[0]; ++i)
    b[0] += u_bar_vec[i] * x[i];
}


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