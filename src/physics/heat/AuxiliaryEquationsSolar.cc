#include "physics/heat/AuxiliaryEquationsSolar.h"

namespace Heat {

void AuxiliaryEquationsSolar::computeAuxiliaryRhs(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, ArrayType<Real, 1>& rhs)
{
  Real interior_temp = u_aux_vec->getVector(1)[0];
  rhs[0] = m_air_temp->computeNetFlux(u_vec, interior_temp, t);
}

void AuxiliaryEquationsSolar::computeAuxiliaryMassMatrix(int block, Real t, linear_system::SimpleAssemblerPtr assembler)
{
  std::vector<DofInt> dofs = {0};
  ArrayType<Real, 2> mat(boost::extents[1][1]);
  mat[0][0] = m_air_temp->getThermalMass();
  std::cout << "AuxiliaryEquationsSolar Mass matrix value = " << mat[0][0] << std::endl;
  assembler->assembleEntry(dofs, mat);
}

void AuxiliaryEquationsSolar::multiplyAuxiliaryMassMatrix(int block, Real t,
                                const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
{
  Real val = m_air_temp->getThermalMass();
  b[0] = x[0] * val;
}

void AuxiliaryEquationsSolar::computeAuxiliaryJacobian(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, 
                                                       linear_system::SimpleAssemblerPtr mat)
{
  Real interior_temp = u_aux_vec->getVector(1)[0];
  Real val = m_air_temp->computeNetFluxJacobian(u_vec, interior_temp, t);

  std::vector<DofInt> dofs = {0};
  ArrayType<Real, 2> vals(boost::extents[1][1]);
  vals[0][0] = val;
  mat->assembleEntry(dofs, vals);
}

void AuxiliaryEquationsSolar::computeFiniteElementJacobianVectorProduct(int jblock, DiscVectorPtr u_vec,
                                                AuxiliaryEquationsStoragePtr u_aux_vec, Real t, 
                                                const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
{
  Real interior_temp = u_aux_vec->getVector(1)[0];
  m_heat_eqn.computedRdTinterior_airProduct(u_vec, interior_temp, t, x[0], b);
}

void AuxiliaryEquationsSolar::computeAuxiliaryJacobianVectorProduct(int iblock, int jblock, 
      DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
{
  Real t_interior = u_aux_vec->getVector(1)[0];
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

}