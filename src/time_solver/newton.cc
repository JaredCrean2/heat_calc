#include "time_solver/newton.h"
#include "discretization/disc_vector.h"
#include "linear_system/large_matrix.h"

namespace timesolvers {

NewtonSolver::NewtonSolver(NewtonFunctionPtr func, linear_system::LargeMatrixPtr jac) :
  m_func(func),
  m_jac(jac),
  m_f(func->createVector()),
  m_delta_u(func->createVector())
{}

NewtonResult NewtonSolver::solve(DiscVectorPtr u, Real abs_tol, Real rel_tol, int itermax)
{
  setupForSolve(u, abs_tol, rel_tol, itermax);

  Real norm = m_func->computeFunc(u, true, m_f);
  Real norm0 = norm;

  if (norm < m_abs_tol)
    return NewtonResult(norm0, norm, 0, abs_tol, rel_tol, itermax);

  for (int i=0; i < itermax; ++i)
  {
    std::cout << "\nNewton iteration " << i << std::endl;
    m_jac->zeroMatrix();
    m_func->computeJacobian(u, m_jac);

    solveStep(u);

    norm = m_func->computeFunc(u, true, m_f);

    std::cout << "at end of Newton iteration " << i << ", absolute norm = " << norm << ", relative norm = " << norm/norm0 << std::endl;
    
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
  m_jac->finishMatrixAssembly();
  m_jac->factor();

  if (!m_f->isVectorCurrent())
    m_f->syncArrayToVector();
  m_delta_u->set(0);

  m_jac->solve(m_f->getVector(), m_delta_u->getVector());

  auto& u_vec = u->getVector();
  auto& delta_u_vec = m_delta_u->getVector();
  for (DofInt i=0; i < u->getNumDofs(); ++i)
  {
    //auto u_orig = u_vec[i];
    u_vec[i] -= delta_u_vec[i];
    //std::cout << "dof " << i << ", original u value = " << u_orig << ", new u value = " << u_vec[i] << ", delta_u = " << delta_u_vec[i] << std::endl;
  }

  u->markVectorModified();
}


}  // namespace