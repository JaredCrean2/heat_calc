#include "time_solver/crank_nicolson_function.h"
#include "mesh/mesh.h"

namespace timesolvers {

CrankNicolsonFunction::CrankNicolsonFunction(std::shared_ptr<PhysicsModel> physics_model, linear_system::LargeMatrixPtr mat, Real t0) :
  m_physics_model(physics_model),
  m_aux_eqns(std::make_shared<CrankNicolsonAuxiliaryEquations>(physics_model, t0)),
  m_assembler(std::make_shared<linear_system::Assembler>(physics_model->getDiscretization(), mat)),
  m_tn(t0),
  m_un(makeDiscVector(physics_model->getDiscretization())),
  m_u_aux_n(makeAuxiliaryEquationsStorage(physics_model->getAuxEquations())),
  m_fn(makeDiscVector(physics_model->getDiscretization())),
  m_tnp1(t0),
  m_delta_u(makeDiscVector(physics_model->getDiscretization())),
  m_Mdelta_u(makeDiscVector(physics_model->getDiscretization()))
{
  m_physics_model->getDiscretization()->getMesh()->getOwnedLocalDofInfo(m_owned_dof_to_local);
}


void CrankNicolsonFunction::resetForNewSolve()
{
  m_fn->set(0);
  m_physics_model->computeRhs(m_un, m_u_aux_n, m_tn, m_fn);
  if (!m_fn->isVectorCurrent())
    m_fn->syncArrayToVector();
}


Real CrankNicolsonFunction::computeFunc(const DiscVectorPtr u_np1, AuxiliaryEquationsStoragePtr u_aux_np1, bool compute_norm, DiscVectorPtr f_np1)
{
  //TODO: add flag for when u_np1 == un, avoid computing M * (u_np1 - u_n) on first iteration
  assertAlways(m_tnp1 - m_tn > 1e-12, "delta_t must be > 1e-12");

  //std::cout << "evaluating CN function" << std::endl;

  f_np1->set(0);
  m_physics_model->computeRhs(u_np1, u_aux_np1, m_tnp1, f_np1);
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
  {
    //std::cout << "i = " << i << ", delta_t_inv = " << delta_t_inv << ", Mdelta_u = "
    //         << Mdelta_u_vec[i] << ", f_np1 = " << f_np1_vec[i] << ", f_n_vec = " << f_n_vec[i] << std::endl;
    f_np1_vec[i] = delta_t_inv * Mdelta_u_vec[i] - 0.5*f_np1_vec[i] - 0.5*f_n_vec[i];
  }
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
void CrankNicolsonFunction::computeJacobian(const DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux_vec, linear_system::LargeMatrixPtr jac)
{
  m_assembler->setAlpha(-0.5);
  m_physics_model->computeJacobian(u, u_aux_vec, m_tnp1, m_assembler);

  m_assembler->setAlpha(1.0/(m_tnp1 - m_tn));
  m_physics_model->computeMassMatrix(m_assembler);

  m_assembler->setAlpha(1);
}


// create an empty vector
DiscVectorPtr CrankNicolsonFunction::createVector()
{ 
  return makeDiscVector(m_physics_model->getDiscretization()); 
}

void CrankNicolsonFunction::setTnp1(DiscVectorPtr u_n, AuxiliaryEquationsStoragePtr u_aux_n, Real t_np1)
{
  m_tn = m_tnp1;
  m_tnp1 = t_np1;
  *m_un = *u_n;
  *m_u_aux_n = *u_aux_n;

  m_aux_eqns->setTnp1(u_n, u_aux_n, t_np1);
}

}