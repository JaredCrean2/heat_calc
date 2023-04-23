#include "time_solver/crank_nicolson.h"
#include "linear_system/sparsity_pattern.h"
#include "linear_system/sparsity_pattern_mesh.h"
#include "linear_system/sparsity_pattern_augmented.h"

namespace timesolvers {

void checkTimeStepperOpts(const TimeStepperOpts& opts, bool check_implicit)
{
  assertAlways(opts.t_end > opts.t_start, "t_end must be > t_start");
  if (check_implicit)
  {
    assertAlways(opts.mat_type != linear_system::LargeMatrixType::Unknown, "Matrix type cannot be unknown");
    assertAlways(!!opts.matrix_opts, "matrix_opts must be defined");
    assertAlways(opts.nonlinear_abs_tol > 0 || opts.nonlinear_rel_tol > 0, "Either nonlinear_abs_tol or nonlinear_rel_tol must be > 0");
    assertAlways(opts.nonlinear_itermax > 0, "nonlinear_itermax must be > 0");
  }
}


CrankNicolson::CrankNicolson(std::shared_ptr<PhysicsModel> physics_model, DiscVectorPtr u,
                             AuxiliaryEquationsStoragePtr u_aux_vec, 
                             TimeStepperOpts opts) :
  m_physics_model(physics_model),
  m_aux_eqns(physics_model->getAuxEquations()),
  m_u(u),
  m_u_aux(u_aux_vec),
  m_opts(opts)
{
  checkTimeStepperOpts(opts);
  auto mesh     = physics_model->getDiscretization()->getMesh();
  std::shared_ptr<linear_system::SparsityPattern> sparsity = std::make_shared<linear_system::SparsityPatternMesh>(mesh);
  if (opts.solve_auxiliary_equations_combined_system)
  {
    int num_augmented = 0;
    for (int i=1; i < physics_model->getAuxEquations()->getNumBlocks(); ++i)
      num_augmented += physics_model->getAuxEquations()->getBlockSize(i);
    
    sparsity = std::make_shared<linear_system::SparsityPatternAugmented>(sparsity, num_augmented, MPI_COMM_WORLD);
  }
  m_matrix      = largeMatrixFactory(opts.mat_type, opts.matrix_opts, sparsity);
  m_func        = std::make_shared<CrankNicolsonFunction>(physics_model, m_matrix, opts.t_start, opts.solve_auxiliary_equations_combined_system);
  m_newton      = std::make_shared<NewtonSolver>(m_func, m_matrix);
}


void CrankNicolson::solve()
{
  //int nsteps = numWholeSteps();
  Real t = m_opts.t_start;

  //TODO: reevaluate this: now that we pass in u_aux, everything should be initialized
  // this is unfortunate: it would be useful to have the initial state in the log
  // file, but some of the BCs haven't been initialized yet
  //m_physics_model->runPostProcessors(0, m_u, m_aux_aux, t);

  int iter = 0;
  while (t < m_opts.t_end)
  {
    Real delta_t = m_opts.timestep_controller->getNextTimestep(t);
    delta_t = std::min(delta_t, m_opts.t_end - t);
    std::cout << "delta_t = " << delta_t << std::endl;
    assertAlways(delta_t > 0, "delta_t must be > 0");

    advanceTimestep(t + delta_t, delta_t);
    t += delta_t;

    m_physics_model->runPostProcessors(iter, m_u, m_u_aux, t);
    if (iter % 10 == 0)
      m_physics_model->getDiscretization()->getMesh()->writeVtkFiles(std::string("mesh") + std::to_string(iter));

    m_opts.timestep_controller->recordLastIteration(m_func->getLastPhysicsRhsNorm());

    iter++;

  }
/*
  Real delta_t_final = finalStepSize();
  // Need to be careful with this: the Newton problem has a term
  // (u_np1 - u_n)/delta_t, so we don't want delta_t too small
  if (delta_t_final > 1e-10)
  {
    advanceTimestep(t + delta_t_final, delta_t_final);
    t += delta_t_final;
    m_physics_model->runPostProcessors(nsteps, m_u, m_u_aux, t);
  }
*/
}


void CrankNicolson::advanceTimestep(Real t_new, Real delta_t)
{
  std::cout << "CN advancing to time " << t_new << std::endl;
  m_func->setTnp1(m_u->getVector(), m_u_aux, t_new);

  NewtonOpts opts;
  opts.nonlinear_abs_tol = m_opts.nonlinear_abs_tol;
  opts.nonlinear_rel_tol = m_opts.nonlinear_rel_tol;
  opts.nonlinear_itermax = m_opts.nonlinear_itermax;

  NewtonResult result = m_newton->solve(m_u->getVector(), m_u_aux, opts);
  m_u->markVectorModified();
  //m_func->completeTimestep(m_u);

  if (!result.isConverged())
  {
    std::stringstream ss;
    ss << result;
    throw std::runtime_error(ss.str());
  }
}

/*
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
*/

}  // namespace