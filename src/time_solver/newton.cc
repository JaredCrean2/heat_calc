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
  m_aux_delta_u(func->getAuxiliaryEquations()->createStorage()),
  m_aux_rhs(func->getAuxiliaryEquations()->createStorage())
{
}

NewtonResult NewtonSolver::solve(DiscVectorPtr u, NewtonOpts opts /*Real abs_tol, Real rel_tol, int itermax*/)
{
  setupForSolve(u, opts);

  //Real norm = m_func->computeFunc(u, true, m_f);
  Real norm = computeRhsAndNorm(u);
  Real norm0 = norm;

  if (norm < opts.nonlinear_abs_tol)
    return NewtonResult(norm0, norm, 0, opts);

  for (int i=0; i < opts.nonlinear_itermax; ++i)
  {
    std::cout << "\nnewton iteration " << i << std::endl;
    solveStep(u);

    std::cout << "after solveStep" << std::endl;

    Real norm = computeRhsAndNorm(u);
    
    if (norm < opts.nonlinear_abs_tol || norm/norm0 < opts.nonlinear_rel_tol)
      return NewtonResult(norm0, norm, i, opts);
  }

  return NewtonResult(norm0, norm, opts.nonlinear_itermax + 1, opts);
}


void NewtonSolver::setupForSolve(DiscVectorPtr u, NewtonOpts opts)
{
  assertAlways(u->getNumDofs() == m_f->getNumDofs(), "vectors u and f must be the same size");
  //m_abs_tol = abs_tol;
  //m_rel_tol = rel_tol;
  //m_itermax = itermax;
  m_opts = opts;
  m_func->resetForNewSolve();
}

Real NewtonSolver::computeRhsAndNorm(DiscVectorPtr u)
{
  Real norm = m_func->computeFunc(u, true, m_f);
  Real norm_squared = norm*norm;

  std::cout << "after first block, norm = " << std::sqrt(norm) << std::endl;


  auto aux_eqns = m_func->getAuxiliaryEquations();
  for (int iblock=1; iblock < aux_eqns->getNumBlocks(); ++iblock)
  {
    auto& rhs = m_aux_rhs->getVector(iblock);
    norm = aux_eqns->computeRhs(iblock, u, true, rhs);
    norm_squared += norm*norm;
    std::cout << "after " << iblock << " block, norm = " << std::sqrt(norm) << std::endl;

  }

  return std::sqrt(norm_squared);
}


void NewtonSolver::solveStep(DiscVectorPtr u)
{
  if (!m_f->isVectorCurrent())
    m_f->syncArrayToVector();
  m_delta_u->set(0);


  computeJacobians(u);

  for (int i=0; i < m_opts.linear_itermax; ++i)
  {
    std::cout << "Gauss Seidel iteration " << i << std::endl;
    Real delta_u_relative_norm = gaussSeidelStep(u);
    std::cout << "relative norm = " << delta_u_relative_norm << std::endl;

    if (i > 0 && delta_u_relative_norm < m_opts.linear_delta_u_tol)
      break;
  }

  auto& u_vec = u->getVector();
  auto& delta_u_vec = m_delta_u->getVector();
  for (DofInt i=0; i < u->getNumDofs(); ++i)
  {
    u_vec[i] -= delta_u_vec[i];
    delta_u_vec[i] = 0;
  }

  // update auxiliary quantities
  auto aux_eqns = m_func->getAuxiliaryEquations();
  for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
  {
    int num_vars = aux_eqns->getBlockSize(block);
    //auto& u_block = m_aux_u->getVector(block);
    auto& u_block = aux_eqns->getBlockSolution(block);
    auto& delta_u_block = m_aux_delta_u->getVector(block);
    for (int i=0; i < num_vars; ++i)
    {
      u_block[i] -= delta_u_block[i];
      delta_u_block[i] = 0;
    }
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
    jac->zeroMatrix();
    m_func->getAuxiliaryEquations()->computeJacobian(iblock, u, jac);
    jac->factor();
  }
}

Real NewtonSolver::gaussSeidelStep(DiscVectorPtr u)
{
  //TODO: zero out vectors?
  auto aux_eqns = m_func->getAuxiliaryEquations();
  ArrayType<Real, 1> delta_u_tmp(boost::extents[u->getNumDofs()]);
  auto delta_u_aux_tmp = aux_eqns->createStorage();
  Real delta_u_relative_norm = 0;

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

    m_jac->solve(rhs, delta_u_tmp);

    auto& delta_u_vec = m_delta_u->getVector();
    for (int i=0; i < rhs.shape()[0]; ++i)
    {
      Real delta_delta_u = std::abs(delta_u_tmp[i] - delta_u_vec[i]);
      if (std::abs(delta_u_vec[i]) > 1e-13)
        delta_u_relative_norm = std::max(delta_u_relative_norm, std::abs(delta_delta_u/delta_u_vec[i]));
      delta_u_vec[i] = delta_u_tmp[i];
    }
  }

  //TODO: there is a more efficient way to do this.  when going from one iblock to the
  //      next, only 1 matrix-vector product needs to be updated
  //TODO: also, if doing a linear solve, only need to compute rhs once
  //TODO: no need to allocate a temporary vector each iteration
  
  // do all other rows
  for (int iblock=1; iblock < aux_eqns->getNumBlocks(); ++iblock)
  {
    auto& rhs0 = m_aux_rhs->getVector(iblock);
    ArrayType<Real, 1> rhs_tmp(boost::extents[aux_eqns->getBlockSize(iblock)]),
                       rhs(boost::extents[aux_eqns->getBlockSize(iblock)]);
    for (int i=0; i < aux_eqns->getBlockSize(iblock); ++i)
      rhs[i] = rhs0[i];
      
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
    auto& delta_u_tmp_vec = delta_u_aux_tmp->getVector(iblock);
    jac->solve(rhs, delta_u_tmp_vec);

    std::cout << "rhs = " << rhs[0] << ", delta_u = " << delta_u_tmp_vec[0] << std::endl;

    auto& delta_u_vec = m_aux_delta_u->getVector(iblock);
    for (int i=0; i < aux_eqns->getBlockSize(iblock); ++i)
    {
      Real delta_delta_u = std::abs(delta_u_tmp_vec[i] - delta_u_vec[i]);
      if (std::abs(delta_u_vec[i]) > 1e-13)
        delta_u_relative_norm = std::max(delta_u_relative_norm, std::abs(delta_delta_u/delta_u_vec[i]));
      delta_u_vec[i] = delta_u_tmp_vec[i];
    }
  }

  return delta_u_relative_norm;
}


}  // namespace