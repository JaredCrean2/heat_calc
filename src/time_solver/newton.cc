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

    norm = computeRhsAndNorm(u);
    std::cout << "newton iteration " << i << " norm = " << norm << ", ratio = " << norm/norm0 << std::endl;
    
    if (norm < opts.nonlinear_abs_tol || norm/norm0 < opts.nonlinear_rel_tol)
      return NewtonResult(norm0, norm, i, opts);
  }

  //std::cout << "norm0 = " << norm0 << ", norm = " << norm << ", ratio = " << norm/norm0 << std::endl;
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
  //std::cout << "computing rhs and norm" << std::endl;
  Real norm = m_func->computeFunc(u, true, m_f);
  Real norm_squared = norm*norm;

  //std::cout << "block 0 norm = " << norm << std::endl;
  //std::cout << "after block 0, norm = " << std::sqrt(norm_squared) << std::endl;


  auto aux_eqns = m_func->getAuxiliaryEquations();
  for (int iblock=1; iblock < aux_eqns->getNumBlocks(); ++iblock)
  {
    auto& rhs = m_aux_rhs->getVector(iblock);
    norm = aux_eqns->computeRhs(iblock, u, true, rhs);
    norm_squared += norm*norm;
    //std::cout << "block " << iblock << " norm = " << norm << std::endl;
    //std::cout << "after " << iblock << " block, norm = " << std::sqrt(norm_squared) << std::endl;
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
    std::cout << "\nGauss Seidel iteration " << i << std::endl;
    Real delta_u_relative_norm = gaussSeidelStep(u);

    //computeLinearResidual(u);
    std::cout << "relative norm = " << delta_u_relative_norm << std::endl;

    if (i > 0 && delta_u_relative_norm < m_opts.linear_delta_u_tol)
      break;
  }

  updateNonlinearSolution(u);
}

void NewtonSolver::updateNonlinearSolution(DiscVectorPtr u)
{
  Real damp_factor = 1;
  auto& u_vec = u->getVector();
  auto& delta_u_vec = m_delta_u->getVector();
  for (DofInt i=0; i < u->getNumDofs(); ++i)
  {
    u_vec[i] -= damp_factor * delta_u_vec[i];
    delta_u_vec[i] = 0;
  }
  
  u->markVectorModified();

  // update auxiliary quantities
  auto aux_eqns = m_func->getAuxiliaryEquations();
  for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
  {
    int num_vars = aux_eqns->getBlockSize(block);
    auto& u_block = aux_eqns->getBlockSolution(block);
    auto& delta_u_block = m_aux_delta_u->getVector(block);
    for (int i=0; i < num_vars; ++i)
    {
      u_block[i] -= damp_factor * delta_u_block[i];
      delta_u_block[i] = 0;
    }
  }
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
  auto delta_u_aux_tmp = aux_eqns->createStorage();

  // compute first row
  Real delta_u_relative_norm_first = gaussSeidelStepFirstRow(u);

  Real delta_u_relative_norm_second = gaussSeidelStepOtherRows(u);

  return std::max(delta_u_relative_norm_first, delta_u_relative_norm_second);

  //TODO: no need to allocate a temporary vector each iteration

}

Real NewtonSolver::gaussSeidelStepFirstRow(DiscVectorPtr u)
{
  auto aux_eqns = m_func->getAuxiliaryEquations();
  ArrayType<Real, 1> delta_u_tmp(boost::extents[u->getNumDofs()]);
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
  return updateLinearSolution(delta_u_tmp, delta_u_vec);
}


Real NewtonSolver::gaussSeidelStepOtherRows(DiscVectorPtr u)
{
  auto aux_eqns = m_func->getAuxiliaryEquations();
  auto delta_u_aux_tmp = aux_eqns->createStorage();

  Real delta_u_relative_norm = 0;
  for (int iblock=1; iblock < aux_eqns->getNumBlocks(); ++iblock)
  {
    ArrayType<Real, 1> rhs(boost::extents[aux_eqns->getBlockSize(iblock)]);
    gaussSeidelComputeRhs(iblock, u, rhs);

    auto jac = m_aux_jacs->getMatrix(iblock);
    auto& delta_u_tmp_vec = delta_u_aux_tmp->getVector(iblock);
    jac->solve(rhs, delta_u_tmp_vec);

    auto& delta_u_vec = m_aux_delta_u->getVector(iblock);
    Real delta_u_relative_norm_i = updateLinearSolution(delta_u_tmp_vec, delta_u_vec);
    delta_u_relative_norm = std::max(delta_u_relative_norm, delta_u_relative_norm_i);
  }

  return delta_u_relative_norm;
}

void NewtonSolver::gaussSeidelComputeRhs(int iblock, DiscVectorPtr u, ArrayType<Real, 1>& rhs)
{
  assertAlways(iblock > 0, "Cannot do first row");

  auto aux_eqns = m_func->getAuxiliaryEquations();
  auto& rhs0 = m_aux_rhs->getVector(iblock);
  ArrayType<Real, 1> rhs_tmp(boost::extents[aux_eqns->getBlockSize(iblock)]);
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
}

Real NewtonSolver::updateLinearSolution(const ArrayType<Real, 1>& delta_u_tmp, ArrayType<Real, 1>& delta_u)
{
  Real delta_u_relative_norm = 0;
  for (int i=0; i < delta_u.shape()[0]; ++i)
  {
    Real delta_delta_u = std::abs(delta_u_tmp[i] - delta_u[i]);
    if (std::abs(delta_u[i]) > 1e-13)
      delta_u_relative_norm = std::max(delta_u_relative_norm, std::abs(delta_delta_u/delta_u[i]));

    delta_u[i] = delta_u_tmp[i];
  }

  return delta_u_relative_norm;
}

//TODO: DEBUGGING
void NewtonSolver::computeLinearResidual(DiscVectorPtr u_vec)
{
  auto aux_eqns = m_func->getAuxiliaryEquations();
  ArrayType<Real, 1> residual(boost::extents[u_vec->getNumDofs()]), 
                     residual_tmp(boost::extents[u_vec->getNumDofs()]);
  auto aux_residual = aux_eqns->createStorage();

  // first row
  {
    m_jac->matVec(m_delta_u->getVector(), residual);
    for (int jblock=1; jblock < aux_eqns->getNumBlocks(); ++jblock)
    {
      aux_eqns->multiplyOffDiagonal(0, jblock, u_vec, m_aux_delta_u->getVector(jblock), residual_tmp);
      for (int i=0; i < u_vec->getNumDofs(); ++i)
        residual[i] += residual_tmp[i];
    }

    auto residual_tmp2 = m_func->createVector();
    m_func->computeFunc(u_vec, false, residual_tmp2);
    for (int i=0; i < u_vec->getNumDofs(); ++i)
      residual[i] -= residual_tmp2->getVector()[i];
  }

  // other rows
  {
    for (int iblock=1; iblock < aux_eqns->getNumBlocks(); ++iblock)
    {
      //std::cout << "computing linear residual for block " << iblock << std::endl;
      auto& aux_residual_i = aux_residual->getVector(iblock);
      ArrayType<Real, 1> aux_residual_tmp(boost::extents[aux_eqns->getBlockSize(iblock)]);
      auto jac = m_aux_jacs->getMatrix(iblock);
      jac->zeroMatrix();
      aux_eqns->computeJacobian(iblock, u_vec, jac);
      jac->matVec(m_aux_delta_u->getVector(iblock), aux_residual_i);
      //std::cout << "diagonal contribution = " << aux_residual_i[0] << std::endl;

      for (int jblock=0; jblock < aux_eqns->getNumBlocks(); ++jblock)
      {
        if (iblock == jblock)
          continue;

        auto& delta_u = jblock == 0 ? m_delta_u->getVector() : m_aux_delta_u->getVector(jblock);
        aux_eqns->multiplyOffDiagonal(iblock, jblock, u_vec, delta_u, aux_residual_tmp);
        //std::cout << "block " << jblock << " contribution = " << aux_residual_tmp[0] << std::endl;
        for (int i=0; i < aux_eqns->getBlockSize(iblock); ++i)
          aux_residual_i[i] += aux_residual_tmp[i];
      }

      aux_eqns->computeRhs(iblock, u_vec, false, aux_residual_tmp);
      //std::cout << "rhs = " << aux_residual_tmp[0] << std::endl;
      for (int i=0; i < aux_eqns->getBlockSize(iblock); ++i)
          aux_residual_i[i] -= aux_residual_tmp[i];
      //std::cout << "final linear residual = " << aux_residual_i[0] << std::endl;     
    }
  }

  std::vector<Real> residuals(aux_eqns->getNumBlocks(), 0);
  //std::cout << "block 0 linear residual: " << std::endl;
  for (int i=0; i < u_vec->getNumDofs(); ++i)
  {
    //std::cout << "dof " << i << " residual = " << residual[i] << std::endl;
    residuals[0] += residual[i]*residual[i];
  }
  residuals[0] = std::sqrt(residuals[0]);

  for (int iblock=1; iblock < aux_eqns->getNumBlocks(); ++iblock)
  {
    //std::cout << "block " << iblock << " linear residual: " << std::endl;
    auto& aux_residual_i = aux_residual->getVector(iblock);
    for (int i=0; i < aux_residual_i.shape()[0]; ++i)
    {
      //std::cout << "dof " << i << " residual = " << aux_residual_i[i] << std::endl;
      residuals[iblock] += aux_residual_i[i]*aux_residual_i[i];
    }

    residuals[iblock] = std::sqrt(residuals[iblock]);
  }

  //for (int i=0; i < aux_eqns->getNumBlocks(); ++i)
  //  std::cout << "block " << i << " linear residual norm = " << residuals[i] << std::endl;

  computeJacobians(u_vec); // basically to factor the matrix again
}


}  // namespace