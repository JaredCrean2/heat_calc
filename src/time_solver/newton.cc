#include "time_solver/newton.h"
#include "discretization/disc_vector.h"
#include "linear_system/large_matrix.h"
#include "linear_system/large_matrix_factory.h"

namespace timesolvers {

NewtonSolver::NewtonSolver(NewtonFunctionPtr func, linear_system::LargeMatrixPtr jac) :
  m_func(func),
  m_jac(jac),
  m_aux_jacs(func->getAuxiliaryEquations()->getJacobians()),
  m_f(func->createVector()),
  m_delta_u(func->createVector()),
  m_aux_u(func->getAuxiliaryEquations()->createStorage()),
  m_aux_delta_u(func->getAuxiliaryEquations()->createStorage()),
  m_aux_rhs(func->getAuxiliaryEquations()->createStorage())
{
  //TODO: need to set initial values for auxiliary variables
}

NewtonResult NewtonSolver::solve(DiscVectorPtr u, Real abs_tol, Real rel_tol, int itermax)
{
  setupForSolve(u, abs_tol, rel_tol, itermax);

  Real norm = m_func->computeFunc(u, true, m_f);
  Real norm0 = norm;

  if (norm < m_abs_tol)
    return NewtonResult(norm0, norm, 0, abs_tol, rel_tol, itermax);

  for (int i=0; i < itermax; ++i)
  {

    solveStep(u);

    norm = m_func->computeFunc(u, true, m_f);
    
    if (norm < m_abs_tol || norm/norm0 < rel_tol)
      return NewtonResult(norm0, norm, i, abs_tol, rel_tol, itermax);
  }

  return NewtonResult(norm0, norm, itermax + 1, abs_tol, rel_tol, itermax);
}


void NewtonSolver::setupForSolve(DiscVectorPtr u, Real abs_tol, Real rel_tol, int itermax)
{
  assertAlways(u->getNumDofs() == m_f->getNumDofs(), "vectors u and f must be the same size");
  m_abs_tol = abs_tol;
  m_rel_tol = rel_tol;
  m_itermax = itermax;
  m_func->resetForNewSolve();
}


void NewtonSolver::solveStep(DiscVectorPtr u)
{
  if (!m_f->isVectorCurrent())
    m_f->syncArrayToVector();
  m_delta_u->set(0);

  computeJacobians(u);
  gaussSeidelStep(u);

  auto& u_vec = u->getVector();
  auto& delta_u_vec = m_delta_u->getVector();
  for (DofInt i=0; i < u->getNumDofs(); ++i)
    u_vec[i] -= delta_u_vec[i];

  // update auxiliary quantities
  auto aux_eqns = m_func->getAuxiliaryEquations();
  for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
  {
    int num_vars = aux_eqns->getBlockSize(block);
    auto& u_block = m_aux_u->getVector(block);
    auto& delta_u_block = m_aux_delta_u->getVector(block);
    for (int i=0; i < num_vars; ++i)
      u_block[i] -= delta_u_block[i];

    aux_eqns->setBlockSolution(block, u_block);
  }

  u->markVectorModified();

  //m_func->updateDependentQuantities(u);
}

void NewtonSolver::computeJacobians(DiscVectorPtr u)
{
  m_jac->zeroMatrix();
  m_func->computeJacobian(u, m_jac);

  //TODO: split startMatrixAssembly/finishMatrixAssembly
  m_jac->finishMatrixAssembly();
  m_jac->factor();

  for (int iblock=1; iblock < m_func->getAuxiliaryEquations()->getNumBlocks(); ++iblock)
  {
    auto jac = m_aux_jacs->getMatrix(iblock);
    m_func->getAuxiliaryEquations()->computeJacobian(iblock, u, jac);
  }
}

void NewtonSolver::gaussSeidelStep(DiscVectorPtr u)
{
  //TODO: zero out vectors?
  auto aux_eqns = m_func->getAuxiliaryEquations();
  // compute first row
  {
    ArrayType<Real, 1> rhs(boost::extents[aux_eqns->getBlockSize(0)]);
    auto& rhs0 = m_f->getVector();
    for (int i=0; i < aux_eqns->getBlockSize(0); ++i)
      rhs[i] = rhs0[i];
      
    ArrayType<Real, 1> rhs_tmp(boost::extents[aux_eqns->getBlockSize(0)]);
    for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
    {
      aux_eqns->multiplyOffDiagonal(0, block, u, m_aux_delta_u->getVector(block), rhs_tmp);
      for (size_t j=0; j < rhs.shape()[0]; ++j)
        rhs[j] -= rhs_tmp[j];
    }

    m_jac->solve(rhs, m_delta_u->getVector());
  }

  //TODO: there is a more efficient way to do this.  when going from one iblock to the
  //      next, only 1 matrix-vector product needs to be updated

  // do all other rows
  for (int iblock=1; iblock < aux_eqns->getNumBlocks(); ++iblock)
  {
    auto& rhs = m_aux_rhs->getVector(iblock);
    aux_eqns->computeRhs(iblock, u, rhs);
    ArrayType<Real, 1> rhs_tmp(boost::extents[aux_eqns->getBlockSize(iblock)]);
    for (int jblock=0; jblock < aux_eqns->getNumBlocks(); ++jblock)
    {
      if (iblock == jblock)
        continue;

      auto& delta_u_j = jblock == 0 ? m_delta_u->getVector() : m_aux_delta_u->getVector(jblock);
      aux_eqns->multiplyOffDiagonal(iblock, jblock, u, delta_u_j, rhs_tmp);
      for (int i=0; i < aux_eqns->getBlockSize(iblock); ++i)
        rhs[i] -= rhs_tmp[i];
    }

    auto jac = m_aux_jacs->getMatrix(iblock);
    jac->solve(rhs, m_aux_delta_u->getVector(iblock));
  }
}


}  // namespace