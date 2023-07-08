#include "time_solver/newton.h"
#include "discretization/disc_vector.h"
#include "linear_system/large_matrix.h"
#include "linear_system/large_matrix_factory.h"
#include <iomanip>

namespace timesolvers {

NewtonSolver::NewtonSolver(NewtonFunctionPtr func, linear_system::LargeMatrixPtr jac, MPI_Comm comm) :
  m_func(func),
  m_jac(jac),
  m_aux_jacs(func->getAuxiliaryEquations()->getJacobians()),
  m_f(boost::extents[jac->getSparsityPattern()->getNumLocalDofs()] /*boost::extents[jac->getNLocal()]*/),
  m_delta_u(boost::extents[jac->getSparsityPattern()->getNumLocalDofs()] /*boost::extents[jac->getNLocal()]*/),
  m_aux_delta_u(func->getAuxiliaryEquations()->createStorage()),
  m_aux_rhs(func->getAuxiliaryEquations()->createStorage()),
  //m_comm(comm),
  m_am_i_root(commRank(comm) == 0)
{
}

NewtonResult NewtonSolver::solve(ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, NewtonOpts opts /*Real abs_tol, Real rel_tol, int itermax*/)
{
  setupForSolve(u, u_aux_vec, opts);

  //Real norm = m_func->computeFunc(u, true, m_f);
  Real norm = computeRhsAndNorm(u, u_aux_vec);
  Real norm0 = norm;

  if (m_am_i_root)
    std::cout << "at start of Newton iteration, norm0 = " << norm0 << ", norm = " << norm << ", ratio = " << norm/norm0 << std::endl;
  //checkJacobianFiniteDifference(u, u_aux_vec);

  if (norm < opts.nonlinear_abs_tol)
    return NewtonResult(norm0, norm, 0, opts);

  for (int i=0; i < opts.nonlinear_itermax; ++i)
  {
    if (m_am_i_root)
      std::cout << "\nnewton iteration " << i << std::endl;
    solveStep(u, u_aux_vec);

    norm = computeRhsAndNorm(u, u_aux_vec);
    if (m_am_i_root)
      std::cout << "newton iteration " << i << " norm = " << norm << ", ratio = " << norm/norm0 << std::endl;
    
    if (norm < opts.nonlinear_abs_tol || norm/norm0 < opts.nonlinear_rel_tol)
      return NewtonResult(norm0, norm, i, opts);
  }

  if (m_am_i_root)
    std::cout << "norm0 = " << norm0 << ", norm = " << norm << ", ratio = " << norm/norm0 << std::endl;
  return NewtonResult(norm0, norm, opts.nonlinear_itermax + 1, opts);
}


void NewtonSolver::setupForSolve(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, NewtonOpts opts)
{
  assertAlways(u.shape()[0] == m_f.shape()[0], "vectors u and f must be the same size");
  //m_abs_tol = abs_tol;
  //m_rel_tol = rel_tol;
  //m_itermax = itermax;
  m_opts = opts;
  m_func->resetForNewSolve();
}

Real NewtonSolver::computeRhsAndNorm(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec)
{
  //std::cout << "computing rhs and norm" << std::endl;
  Real norm = m_func->computeFunc(u, u_aux_vec, true, m_f);
  Real norm_squared = norm*norm;

  //std::cout << "block 0 norm = " << norm << std::endl;
  //std::cout << "after block 0, norm = " << std::sqrt(norm_squared) << std::endl;


  auto aux_eqns = m_func->getAuxiliaryEquations();
  for (int iblock=1; iblock < aux_eqns->getNumBlocks(); ++iblock)
  {
    auto& rhs = m_aux_rhs->getVector(iblock);
    norm = aux_eqns->computeRhs(iblock, u, u_aux_vec, true, rhs);
    norm_squared += norm*norm;
    //std::cout << "block " << iblock << " norm = " << norm << std::endl;
    //std::cout << "after " << iblock << " block, norm = " << std::sqrt(norm_squared) << std::endl;
  }

  return std::sqrt(norm_squared);
}


void NewtonSolver::solveStep(ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec)
{
  fill(m_delta_u, 0);


  computeJacobians(u, u_aux_vec);

  // if there are no auxiliary blocks, no need to do more than one Gauss-Seidel step
  int linear_itermax = u_aux_vec->getNumBlocks() > 1 ? m_opts.linear_itermax : 1;
  for (int i=0; i < linear_itermax; ++i)
  {
    if (m_am_i_root)
      std::cout << "\nGauss Seidel iteration " << i << std::endl;
    Real delta_u_relative_norm = gaussSeidelStep(u, u_aux_vec);

    //computeLinearResidual(u);
    if (m_am_i_root)
      std::cout << "relative norm = " << delta_u_relative_norm << std::endl;

    if (i > 0 && delta_u_relative_norm < m_opts.linear_delta_u_tol)
      break;
  }

  updateNonlinearSolution(u, u_aux_vec);
}

void NewtonSolver::updateNonlinearSolution(ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec)
{
  Real damp_factor = 1;
  for (DofInt i=0; i < u_vec.shape()[0]; ++i)
  {
    u_vec[i] -= damp_factor * m_delta_u[i];
    m_delta_u[i] = 0;
  }
  
  // update auxiliary quantities
  auto aux_eqns = m_func->getAuxiliaryEquations();
  for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
  {
    int num_vars = aux_eqns->getBlockSize(block);
    auto& u_block = u_aux_vec->getVector(block);
    auto& delta_u_block = m_aux_delta_u->getVector(block);
    for (int i=0; i < num_vars; ++i)
    {
      u_block[i] -= damp_factor * delta_u_block[i];
      delta_u_block[i] = 0;
    }
  }
}

void NewtonSolver::computeJacobians(ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec)
{
  //checkJacobianFiniteDifference(u, u_aux_vec);
  m_jac->zeroMatrix();
  m_func->computeJacobian(u, u_aux_vec, m_jac);

  //TODO: split startMatrixAssembly/finishMatrixAssembly
  m_jac->finishMatrixAssembly();
  m_jac->factor();

  for (int iblock=1; iblock < m_func->getAuxiliaryEquations()->getNumBlocks(); ++iblock)
  {
    auto jac = m_aux_jacs->getMatrix(iblock);
    jac->zeroMatrix();
    m_func->getAuxiliaryEquations()->computeJacobian(iblock, u, u_aux_vec, jac);
    jac->factor();
  }
}
namespace {
Real computeError(Real val1, Real val2, Real eps, Real tol)
{
  Real error = std::abs((val1 - val2)/(std::max(std::abs(val1), std::abs(val2))));
  if (std::abs(val1) < 5*eps && std::abs(val2) < 5*eps)
    error = 0;

  return error;
}
}

void NewtonSolver::checkJacobianFiniteDifference(ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec)
{
  std::cout << std::setprecision(16);
  auto aux_eqns = m_func->getAuxiliaryEquations();

  ArrayType<Real, 1> x_u(boost::extents[u_vec.shape()[0]]);
  ArrayType<Real, 1> b_u(boost::extents[u_vec.shape()[0]]);
  ArrayType<Real, 1> b1_u(boost::extents[u_vec.shape()[0]]);
  ArrayType<Real, 1> b2_u(boost::extents[u_vec.shape()[0]]);
  ArrayType<Real, 1> b3_u(boost::extents[u_vec.shape()[0]]);
  ArrayType<Real, 1> x_T(boost::extents[1]);
  ArrayType<Real, 1> b1_T(boost::extents[1]);
  ArrayType<Real, 1> b2_T(boost::extents[1]);
  ArrayType<Real, 1> b3_T(boost::extents[1]);

  m_jac->zeroMatrix();
  m_func->computeJacobian(u_vec, u_aux_vec, m_jac);
  m_jac->finishMatrixAssembly();

  //auto jac = m_aux_jacs->getMatrix(1);
  //jac->zeroMatrix();
  //aux_eqns->computeJacobian(1, u_vec, u_aux_vec, jac);


  m_func->computeFunc(u_vec, u_aux_vec, false, b1_u);
  if (m_am_i_root)
    std::cout << "computing auxiliary equations rhs initially" << std::endl;
  //aux_eqns->computeRhs(1, u_vec, u_aux_vec, false, b1_T);
  const int nvectors = 10;
  const Real eps = 1e-7;
  const Real tol = 1e-2;

  // finite difference u
  if (m_am_i_root)
    std::cout << "\ndoing finite difference check of u" << std::endl;
  for (int v=0; v < nvectors; ++v)
  {
    if (m_am_i_root)
      std::cout << "\nvector " << v << std::endl;
    for (int i=0; i < x_u.shape()[0]; ++i)
    {
      x_u[i] = (i % nvectors) == v ? 1 : 0;
      u_vec[i] += x_u[i]*eps;

      //if (i == x_u.shape()[0] - 1)
      //{
      //  x_u[i] = 1;
      //  u_vec[i] += x_u[i]*eps;
      //}

    }

    m_func->computeFunc(u_vec, u_aux_vec, false, b2_u);
    //aux_eqns->computeRhs(1, u_vec, u_aux_vec, false, b2_T);


    m_jac->matVec(x_u, b3_u);
    //aux_eqns->multiplyOffDiagonal(1, 0, u_vec, u_aux_vec, x_u, b3_T);

    for (int i=0; i < x_u.shape()[0]; ++i)
    {
      Real val_fd = (b2_u[i] - b1_u[i])/eps;
      //Real error = std::abs((val_fd - b3_u[i])/(std::max(std::abs(val_fd), std::abs(b3_u[i]))));
      Real error = computeError(val_fd, b3_u[i], eps, tol);
      if (m_am_i_root)
        std::cout << "dof " << i << ", val_fd = " << val_fd << ", val_matvec = " << b3_u[i] << ", diff = " << error << std::endl;
      if (error > tol)
        throw std::runtime_error("finite difference test failed");
    }

    {
      Real val_fd = (b2_T[0] - b1_T[0])/eps;
      //Real error = std::abs((val_fd - b3_T[0])/(std::max(std::abs(val_fd), std::abs(b3_T[0]))));
      Real error = computeError(val_fd, b3_T[0], eps, tol);

      if (m_am_i_root)
        std::cout << "for aux equations, val_fd = " << val_fd << ", val_matvec = " << b3_T[0] << ", diff = " << error << std::endl;
      if (error > tol)
        throw std::runtime_error("finite difference test failed");
    }

    for (int i=0; i < x_u.shape()[0]; ++i)
      u_vec[i] -= x_u[i]*eps;
  }
/*
  // finite difference T
  std::cout << "\ndoing finite difference check of T" << std::endl;
  { 
    x_T[0] = 1;

    u_aux_vec->getVector(1)[0] += x_T[0]*eps;

    m_func->computeFunc(u_vec, u_aux_vec, false, b2_u);
    std::cout << "computing auxiliary equation rhs at perturbed state" << std::endl;
    aux_eqns->computeRhs(1, u_vec, u_aux_vec, false, b2_T);
    
    u_aux_vec->getVector(1)[0] -= x_T[0]*eps;

    aux_eqns->multiplyOffDiagonal(0, 1, u_vec, u_aux_vec, x_T, b3_u);
    jac->matVec(x_T, b3_T);

    for (int i=0; i < x_u.shape()[0]; ++i)
    {
      Real val_fd = (b2_u[i] - b1_u[i])/eps;
      //Real error = std::abs((val_fd - b3_u[i])/(std::max(std::abs(val_fd), std::abs(b3_u[i]))));
      Real error = computeError(val_fd, b3_u[i], eps, tol);

      std::cout << "dof " << i << ", val_fd = " << val_fd << ", val_matvec = " << b3_u[i] << ", diff = " << error << std::endl;
      if (error > tol)
        throw std::runtime_error("finite difference test failed");
    }    

    {
      Real val_fd = (b2_T[0] - b1_T[0])/eps;
      //Real error = std::abs((val_fd - b3_T[0])/(std::max(std::abs(val_fd), std::abs(b3_T[0]))));
      Real error = computeError(val_fd, b3_T[0], eps, tol);

      std::cout << "for aux equations, val_fd = " << val_fd << ", val_matvec = " << b3_T[0] << ", diff = " << error << std::endl;
      if (error > tol)
        throw std::runtime_error("finite difference test failed");
    }

  }
*/
}

Real NewtonSolver::gaussSeidelStep(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec)
{
  //TODO: zero out vectors?
  auto aux_eqns = m_func->getAuxiliaryEquations();
  auto delta_u_aux_tmp = aux_eqns->createStorage();

  // compute first row
  Real delta_u_relative_norm_first = gaussSeidelStepFirstRow(u, u_aux_vec);

  Real delta_u_relative_norm_second = gaussSeidelStepOtherRows(u, u_aux_vec);

  if (m_am_i_root)
  {
    std::cout << "delta_u_relative_norm_first = " << delta_u_relative_norm_first << std::endl;
    std::cout << "delta_u_relative_norm_second = " << delta_u_relative_norm_second << std::endl;
  }
  //computeLinearResidual(u, u_aux_vec);
  return std::max(delta_u_relative_norm_first, delta_u_relative_norm_second);

  //TODO: no need to allocate a temporary vector each iteration

}

Real NewtonSolver::gaussSeidelStepFirstRow(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec)
{
  auto aux_eqns = m_func->getAuxiliaryEquations();
  ArrayType<Real, 1> delta_u_tmp(boost::extents[u.shape()[0]]);
  ArrayType<Real, 1> rhs(boost::extents[u.shape()[0]]);
  
  for (int i=0; i < u.shape()[0]; ++i)
    rhs[i] = m_f[i];
    
  ArrayType<Real, 1> rhs_tmp(boost::extents[aux_eqns->getBlockSize(0)]);
  for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
  {
    aux_eqns->multiplyOffDiagonal(0, block, u, u_aux_vec, m_aux_delta_u->getVector(block), rhs_tmp);
    for (size_t j=0; j < rhs.shape()[0]; ++j)
    {
      rhs[j] -= rhs_tmp[j];
    }
  }

  m_jac->solve(rhs, delta_u_tmp);

  return updateLinearSolution(delta_u_tmp, m_delta_u);
}


Real NewtonSolver::gaussSeidelStepOtherRows(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec)
{
  auto aux_eqns = m_func->getAuxiliaryEquations();
  auto delta_u_aux_tmp = aux_eqns->createStorage();

  Real delta_u_relative_norm = 0;
  for (int iblock=1; iblock < aux_eqns->getNumBlocks(); ++iblock)
  {
    ArrayType<Real, 1> rhs(boost::extents[aux_eqns->getBlockSize(iblock)]);
    gaussSeidelComputeRhs(iblock, u, u_aux_vec, rhs);

    auto jac = m_aux_jacs->getMatrix(iblock);
    auto& delta_u_tmp_vec = delta_u_aux_tmp->getVector(iblock);
    jac->solve(rhs, delta_u_tmp_vec);

    auto& delta_u_vec = m_aux_delta_u->getVector(iblock);
    Real delta_u_relative_norm_i = updateLinearSolution(delta_u_tmp_vec, delta_u_vec);
    delta_u_relative_norm = std::max(delta_u_relative_norm, delta_u_relative_norm_i);
  }

  return delta_u_relative_norm;
}

void NewtonSolver::gaussSeidelComputeRhs(int iblock, const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, ArrayType<Real, 1>& rhs)
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

    auto& delta_u_j = jblock == 0 ? m_delta_u : m_aux_delta_u->getVector(jblock);
    aux_eqns->multiplyOffDiagonal(iblock, jblock, u, u_aux_vec, delta_u_j, rhs_tmp);
    for (int i=0; i < aux_eqns->getBlockSize(iblock); ++i)
      rhs[i] -= rhs_tmp[i];
  }
}

Real NewtonSolver::updateLinearSolution(const ArrayType<Real, 1>& delta_u_tmp, ArrayType<Real, 1>& delta_u)
{
  Real delta_u_relative_norm = 0;
  Real max_delta_delta_u = 0;
  for (int i=0; i < delta_u.shape()[0]; ++i)
  {
    Real delta_delta_u = std::abs(delta_u_tmp[i] - delta_u[i]);
    max_delta_delta_u = std::max(max_delta_delta_u, std::abs(delta_delta_u));
    if (std::abs(delta_u[i]) > 1e-13)
      delta_u_relative_norm = std::max(delta_u_relative_norm, std::abs(delta_delta_u/delta_u[i]));

    delta_u[i] = delta_u_tmp[i];
  }

  if (m_am_i_root)
    std::cout << "max_delta_delta_u = " << max_delta_delta_u << std::endl;
  return delta_u_relative_norm;
}

//TODO: DEBUGGING
void NewtonSolver::computeLinearResidual(ArrayType<Real, 1>& u_vec, AuxiliaryEquationsStoragePtr u_aux_vec)
{
  auto aux_eqns = m_func->getAuxiliaryEquations();
  ArrayType<Real, 1> residual(boost::extents[u_vec.shape()[0]]), 
                     residual_tmp(boost::extents[u_vec.shape()[0]]),
                     residual_tmp2(boost::extents[u_vec.shape()[0]]);
  auto aux_residual = aux_eqns->createStorage();

  // first row
  {
    m_jac->matVec(m_delta_u, residual);
    for (int jblock=1; jblock < aux_eqns->getNumBlocks(); ++jblock)
    {
      aux_eqns->multiplyOffDiagonal(0, jblock, u_vec, u_aux_vec, m_aux_delta_u->getVector(jblock), residual_tmp);
      for (int i=0; i < u_vec.shape()[0]; ++i)
        residual[i] += residual_tmp[i];
    }

    m_func->computeFunc(u_vec, u_aux_vec, false, residual_tmp2);
    for (int i=0; i < u_vec.shape()[0]; ++i)
      residual[i] -= residual_tmp2[i];
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
      aux_eqns->computeJacobian(iblock, u_vec, u_aux_vec, jac);
      jac->matVec(m_aux_delta_u->getVector(iblock), aux_residual_i);
      //std::cout << "diagonal contribution = " << aux_residual_i[0] << std::endl;

      for (int jblock=0; jblock < aux_eqns->getNumBlocks(); ++jblock)
      {
        if (iblock == jblock)
          continue;

        auto& delta_u = jblock == 0 ? m_delta_u : m_aux_delta_u->getVector(jblock);
        aux_eqns->multiplyOffDiagonal(iblock, jblock, u_vec, u_aux_vec, delta_u, aux_residual_tmp);
        //std::cout << "block " << jblock << " contribution = " << aux_residual_tmp[0] << std::endl;
        for (int i=0; i < aux_eqns->getBlockSize(iblock); ++i)
          aux_residual_i[i] += aux_residual_tmp[i];
      }

      aux_eqns->computeRhs(iblock, u_vec, u_aux_vec, false, aux_residual_tmp);
      //std::cout << "rhs = " << aux_residual_tmp[0] << std::endl;
      for (int i=0; i < aux_eqns->getBlockSize(iblock); ++i)
          aux_residual_i[i] -= aux_residual_tmp[i];
      //std::cout << "final linear residual = " << aux_residual_i[0] << std::endl;     
    }
  }

  std::vector<Real> residuals(aux_eqns->getNumBlocks(), 0);
  //std::cout << "block 0 linear residual: " << std::endl;
  for (int i=0; i < u_vec.shape()[0]; ++i)
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

  if (m_am_i_root)
    for (int i=0; i < aux_eqns->getNumBlocks(); ++i)
      std::cout << "block " << i << " linear residual norm = " << residuals[i] << std::endl;

  computeJacobians(u_vec, u_aux_vec); // basically to factor the matrix again
}


}  // namespace