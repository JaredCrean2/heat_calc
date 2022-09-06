#include "time_solver/crank_nicolson.h"
#include "linear_system/sparsity_pattern_mesh.h"

namespace timesolvers {


AuxiliaryEquationsStoragePtr CrankNicolsonAuxiliaryEquations::createStorage()
{
  return std::make_shared<AuxiliaryEquationStorage>(*m_aux_eqns);
}

void CrankNicolsonAuxiliaryEquations::multiplyOffDiagonal(int iblock, int jblock, DiscVectorPtr u_vec, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
{
  std::cout << "\nEntered multiplyOffDiagonal for blocks " << iblock << ", " << jblock << std::endl;
  m_aux_eqns->multiplyOffDiagonal(iblock, jblock, u_vec, m_tnp1, x, b);
  int num_vars = getBlockSize(iblock);
  for (int i=0; i < num_vars; ++i)
    b[i] *= -0.5;
}


CrankNicolsonFunction::CrankNicolsonFunction(std::shared_ptr<PhysicsModel> physics_model, linear_system::LargeMatrixPtr mat, Real t0) :
  m_physics_model(physics_model),
  m_aux_eqns(std::make_shared<CrankNicolsonAuxiliaryEquations>(physics_model, t0)),
  m_assembler(std::make_shared<linear_system::Assembler>(physics_model->getDiscretization(), mat)),
  m_tn(t0),
  m_un(makeDiscVector(physics_model->getDiscretization())),
  m_fn(makeDiscVector(physics_model->getDiscretization())),
  m_tnp1(t0),
  m_delta_u(makeDiscVector(physics_model->getDiscretization())),
  m_Mdelta_u(makeDiscVector(physics_model->getDiscretization()))
{
  m_physics_model->getDiscretization()->getMesh()->getOwnedLocalDofInfo(m_owned_dof_to_local);
  m_aux_eqns->setCNFunc(this);
}

void CrankNicolsonFunction::resetForNewSolve()
{
  m_fn->set(0);
  m_physics_model->computeRhs(m_un, m_tn, m_fn);
  if (!m_fn->isVectorCurrent())
    m_fn->syncArrayToVector();
}

Real CrankNicolsonFunction::computeFunc(const DiscVectorPtr u_np1, bool compute_norm, DiscVectorPtr f_np1)
{
  //TODO: add flag for when u_np1 == un, avoid computing M * (u_np1 - u_n) on first iteration
  assertAlways(m_tnp1 - m_tn > 1e-12, "delta_t must be > 1e-12");

  f_np1->set(0);
  m_physics_model->computeRhs(u_np1, m_tnp1, f_np1);
  if (!f_np1->isVectorCurrent())
    f_np1->syncArrayToVector();

  if (!u_np1->isVectorCurrent())
    u_np1->syncArrayToVector();

  if (!m_un->isVectorCurrent())
    m_un->syncArrayToVector();

  // compute M * (u_np1 - u_n)
  m_delta_u->set(0);
  auto& u_np1_vec   = u_np1->getVector();
  auto& u_n_vec     = m_un->getVector();
  auto& delta_u_vec = m_delta_u->getVector();
  for (int i=0; i < m_delta_u->getNumDofs(); ++i)
    delta_u_vec[i] = u_np1_vec[i] - u_n_vec[i];
  m_delta_u->markVectorModified();

  m_physics_model->applyMassMatrix(m_delta_u, m_Mdelta_u);
  if (!m_Mdelta_u->isVectorCurrent())
    m_Mdelta_u->syncArrayToVector();

  // compute f_np1 = M * (u_np1 - u_n) / delta_t - 0.5 * f(u_n, t_n) - 0.5 * f(u_np1, t_np1)
  auto& Mdelta_u_vec = m_Mdelta_u->getVector();
  auto& f_np1_vec    = f_np1->getVector();
  auto& f_n_vec      = m_fn->getVector();
  Real delta_t_inv    = 1.0/(m_tnp1 - m_tn);
  for (int i=0; i < f_np1->getNumDofs(); ++i)
    f_np1_vec[i] = delta_t_inv * Mdelta_u_vec[i] - 0.5*f_np1_vec[i] - 0.5*f_n_vec[i];
  f_np1->markVectorModified();

  Real norm = 0, norm_global = 0;
  if (compute_norm)
  {
    //TODO: not sure this is the right norm to use
    for (auto dof : m_owned_dof_to_local)
      norm += f_np1_vec[dof] * f_np1_vec[dof];
    
    MPI_Allreduce(&norm, &norm_global, 1, REAL_MPI_DATATYPE, MPI_SUM, MPI_COMM_WORLD);
    norm = std::sqrt(norm_global);
  }

  return norm;
}

// compute jac = df/du, overwriting jac
void CrankNicolsonFunction::computeJacobian(const DiscVectorPtr u, linear_system::LargeMatrixPtr jac)
{
  m_assembler->setAlpha(-0.5);
  m_physics_model->computeJacobian(u, m_tnp1, m_assembler);

  m_assembler->setAlpha(1.0/(m_tnp1 - m_tn));
  m_physics_model->computeMassMatrix(m_assembler);

  m_assembler->setAlpha(1);
}

/*
void CrankNicolsonFunction::updateDependentQuantities(DiscVectorPtr u)
{
  m_physics_model->updateDependentQuantities(u, m_tnp1);
}
*/
/*
void CrankNicolsonFunction::completeTimestep(DiscVectorPtr u)
{
  m_physics_model->completeTimestep(u, m_tnp1);
}
*/

// create an empty vector
DiscVectorPtr CrankNicolsonFunction::createVector()
{ 
  return makeDiscVector(m_physics_model->getDiscretization()); 
}

void CrankNicolsonFunction::setTnp1(DiscVectorPtr u_n, Real t_np1)
{
  m_tn = m_tnp1;
  m_tnp1 = t_np1;
  *m_un = *u_n;

  m_aux_eqns->setTnp1(u_n, t_np1);
}


void checkTimeStepperOpts(const TimeStepperOpts& opts, bool check_implicit)
{
  assertAlways(opts.t_end > opts.t_start, "t_end must be > t_start");
  assertAlways(opts.delta_t > 0, "delta_t must be > 0");
  if (check_implicit)
  {
    assertAlways(opts.mat_type != linear_system::LargeMatrixType::Unknown, "Matrix type cannot be unknown");
    assertAlways(!!opts.matrix_opts, "matrix_opts must be defined");
    assertAlways(opts.nonlinear_abs_tol > 0 || opts.nonlinear_rel_tol > 0, "Either nonlinear_abs_tol or nonlinear_rel_tol must be > 0");
    assertAlways(opts.nonlinear_itermax > 0, "nonlinear_itermax must be > 0");
  }
}


CrankNicolson::CrankNicolson(std::shared_ptr<PhysicsModel> physics_model, DiscVectorPtr u, TimeStepperOpts opts) :
  m_physics_model(physics_model),
  m_u(u),
  m_opts(opts)
{
  checkTimeStepperOpts(opts);
  auto mesh     = physics_model->getDiscretization()->getMesh();
  auto sparsity = std::make_shared<linear_system::SparsityPatternMesh>(mesh);
  m_matrix      = largeMatrixFactory(opts.mat_type, mesh->getNumOwnedDofs(), mesh->getNumOwnedDofs(), opts.matrix_opts, sparsity);
  m_func        = std::make_shared<CrankNicolsonFunction>(physics_model, m_matrix, opts.t_start);
  m_newton      = std::make_shared<NewtonSolver>(m_func, m_matrix);
}


void CrankNicolson::solve()
{
  int nsteps = numWholeSteps();
  Real t = m_opts.t_start;
  for (int i=0; i < nsteps; ++i)
  {
    advanceTimestep(t + m_opts.delta_t, m_opts.delta_t);
    t += m_opts.delta_t;
  }

  Real delta_t_final = finalStepSize();
  // Need to be careful with this: the Newton problem has a term
  // (u_np1 - u_n)/delta_t, so we don't want delta_t too small
  if (delta_t_final > 1e-10)
    advanceTimestep(t + delta_t_final, delta_t_final);
}


void CrankNicolson::advanceTimestep(Real t_new, Real delta_t)
{
  std::cout << "CN advancing to time " << t_new << std::endl;
  m_func->setTnp1(m_u, t_new);

  NewtonOpts opts;
  opts.nonlinear_abs_tol = m_opts.nonlinear_abs_tol;
  opts.nonlinear_rel_tol = m_opts.nonlinear_rel_tol;
  opts.nonlinear_itermax = m_opts.nonlinear_itermax;

  NewtonResult result = m_newton->solve(m_u, opts);
  //m_func->completeTimestep(m_u);

  if (!result.isConverged())
  {
    std::stringstream ss;
    ss << result;
    throw std::runtime_error(ss.str());
  }
}


int CrankNicolson::numWholeSteps()
{
  return std::floor((m_opts.t_end - m_opts.t_start)/m_opts.delta_t);
}


double CrankNicolson::finalStepSize()
{
  int num_whole_steps = numWholeSteps();
  Real t_whole = m_opts.delta_t * num_whole_steps;
  Real t_range = m_opts.t_end - m_opts.t_start;
  return std::max(t_range - t_whole, 0.0);
}


}  // namespace