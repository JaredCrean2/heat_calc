#include "time_solver/crank_nicolson_aux_equations.h"
#include "discretization/disc_vector.h"

namespace timesolvers {

AuxiliaryEquationsStoragePtr CrankNicolsonAuxiliaryEquations::createStorage()
{
  return makeAuxiliaryEquationsStorage(m_aux_eqns);
}


Real CrankNicolsonAuxiliaryEquations::computeRhs(int block, const ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, bool compute_norm, ArrayType<Real, 1>& rhs)
{
  int num_vars = getBlockSize(block);

  // compute M * (u_np1 - u_n)/delta_t
  ArrayType<Real, 1> delta_u(boost::extents[num_vars]);
  auto& u_np1 = u_aux_vec->getVector(block);
  auto& u_n   = m_aux_un->getVector(block);
  Real delta_t = m_tnp1 - m_tn;
  for (int i=0; i < num_vars; ++i)
    delta_u[i] = (u_np1[i] - u_n[i])/delta_t;
  m_aux_eqns->multiplyMassMatrix(block, m_tnp1, delta_u, rhs);
  //for (int i=0; i < rhs.shape()[0]; ++i)
  //  rhs[i] = 0;

  // compute 1/2(f(u_np1, t_np1) + f(u_n, t_n))
  auto u_disc_vec = makeDiscVector(m_physics_model->getDiscretization());
  copyToVector(u_vec, u_disc_vec);
  ArrayType<Real, 1> rhs_tmp(boost::extents[num_vars]), rhs_tmp2(boost::extents[num_vars]);
  m_aux_eqns->computeRhs(block, m_un, m_aux_un, m_tn, rhs_tmp);
  m_aux_eqns->computeRhs(block, u_disc_vec, u_aux_vec, m_tnp1, rhs_tmp2);

  Real norm = 0;
  for (int i=0; i < num_vars; ++i)
  {
    rhs[i] -= 0.5*(rhs_tmp[i] + rhs_tmp2[i]);
    norm += rhs[i]*rhs[i];
  }

  return std::sqrt(norm);
}

void CrankNicolsonAuxiliaryEquations::computeJacobian(int block, const ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, linear_system::LargeMatrixPtr mat)
{
  auto assembler = std::make_shared<linear_system::SimpleAssembler>(mat);
  Real delta_t = m_tnp1 - m_tn;
  assembler->setAlpha(1.0/delta_t);
  m_aux_eqns->computeMassMatrix(block, m_tnp1, assembler);

  auto u_disc_vec = makeDiscVector(m_physics_model->getDiscretization());
  copyToVector(u_vec, u_disc_vec);
  assembler->setAlpha(-0.5);
  m_aux_eqns->computeJacobian(block, u_disc_vec, u_aux_vec, m_tnp1, assembler);
}

void CrankNicolsonAuxiliaryEquations::multiplyOffDiagonal(int iblock, int jblock, const ArrayType<Real, 1>& u_vec,
                                                          AuxiliaryEquationsStoragePtr u_aux_vec, 
                                                          const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
{

  auto u_disc_vec = makeDiscVector(m_physics_model->getDiscretization());
  copyToVector(u_vec, u_disc_vec);

  m_aux_eqns->multiplyOffDiagonal(iblock, jblock, u_disc_vec, u_aux_vec, m_tnp1, x, b);
  int num_vars = getBlockSize(iblock);
  for (int i=0; i < num_vars; ++i)
    b[i] *= -0.5;
}

void CrankNicolsonAuxiliaryEquations::setTnp1(DiscVectorPtr u_n, AuxiliaryEquationsStoragePtr u_aux_vec, Real t_np1)
{
  m_tn   = m_tnp1;
  m_tnp1 = t_np1;

  if (m_use_aux_eqns)
  {
    *m_un  = *u_n;

    for (int iblock=1; iblock < m_aux_eqns->getNumBlocks(); ++iblock)
    {
      auto& aux_u_np1 = u_aux_vec->getVector(iblock); //m_aux_unp1.getVector(iblock);
      auto& aux_u_n = m_aux_un->getVector(iblock);
      for (int i=0; i < m_aux_eqns->getBlockSize(iblock); ++i)
        aux_u_n[i] = aux_u_np1[i];
    }
  }
}

}  // namespace