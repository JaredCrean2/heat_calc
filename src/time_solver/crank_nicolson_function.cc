#include "time_solver/crank_nicolson_function.h"
#include "discretization/disc_vector.h"
#include "linear_system/assembler.h"
#include "linear_system/augmented_assembler.h"
#include "linear_system/large_matrix.h"
#include "linear_system/large_matrix_dense.h"
#include "mesh/mesh.h"
#include "physics/AuxiliaryEquations.h"

namespace timesolvers {

CrankNicolsonFunction::CrankNicolsonFunction(std::shared_ptr<PhysicsModel> physics_model, linear_system::LargeMatrixPtr mat,
                                             Real t0, const TimeStepperOpts& opts) :
  m_physics_model(physics_model),
  m_aux_eqns(std::make_shared<CrankNicolsonAuxiliaryEquations>(physics_model, t0, !opts.solve_auxiliary_equations_combined_system)),
  m_assembler(std::make_shared<linear_system::Assembler>(physics_model->getDiscretization(), mat)),
  m_augmented_assembler(std::make_shared<linear_system::AugmentedAssembler>(physics_model->getDiscretization(), mat, 
                                                        getNumAuxiliaryVariables(physics_model->getAuxEquations()))),
  m_tn(t0),
  m_aux_eqns_combined_system(opts.solve_auxiliary_equations_combined_system),
  m_precompute_linear_jacobian(opts.precompute_linear_jacobian),
  m_un(makeDiscVector(physics_model->getDiscretization())),
  m_u_aux_n(makeAuxiliaryEquationsStorage(physics_model->getAuxEquations())),
  m_fn(makeDiscVector(physics_model->getDiscretization())),
  m_fn_aux(makeAuxiliaryEquationsStorage(physics_model->getAuxEquations())),
  m_tnp1(t0),
  m_delta_u(makeDiscVector(physics_model->getDiscretization())),
  m_delta_u_aux(makeAuxiliaryEquationsStorage(physics_model->getAuxEquations())),
  m_Mdelta_u(makeDiscVector(physics_model->getDiscretization())),
  m_Mdelta_u_aux(makeAuxiliaryEquationsStorage(physics_model->getAuxEquations()))

{
  m_physics_model->getDiscretization()->getMesh()->getOwnedLocalDofInfo(m_owned_dof_to_local);
}


void CrankNicolsonFunction::resetForNewSolve()
{
  m_fn->set(0);
  m_physics_model->computeRhs(m_un, m_u_aux_n, m_tn, m_fn);
  auto aux_eqns = m_physics_model->getAuxEquations();

  if (m_aux_eqns_combined_system)
    for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
      aux_eqns->computeRhs(block, m_un, m_u_aux_n, m_tn, m_fn_aux->getVector(block));

  if (!m_fn->isVectorCurrent())
    m_fn->syncArrayToVector();
}

void CrankNicolsonFunction::computePhysicsRhs(DiscVectorPtr u_disc_vec, AuxiliaryEquationsStoragePtr u_aux,
                                              DiscVectorPtr f_disc_vec, AuxiliaryEquationsStoragePtr f_aux)
{
  m_physics_model->computeRhs(u_disc_vec, u_aux, m_tnp1, f_disc_vec);
  if (!f_disc_vec->isVectorCurrent())
    f_disc_vec->syncArrayToVector();

  if (m_aux_eqns_combined_system)
  {
    auto aux_eqns = m_physics_model->getAuxEquations();
    for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
      aux_eqns->computeRhs(block, u_disc_vec, u_aux, m_tnp1, f_aux->getVector(block));
  }
}

Real CrankNicolsonFunction::computeNorm(DiscVectorPtr f_disc_vec, AuxiliaryEquationsStoragePtr f_aux)
{
  bool am_i_last_rank = commRank(MPI_COMM_WORLD) == (commSize(MPI_COMM_WORLD) - 1);
  Real physics_rhs_norm = 0;

  if (!f_disc_vec->isVectorCurrent())
    f_disc_vec->syncArrayToVector();

  auto& f_vec = f_disc_vec->getVector();
  for (auto dof : m_owned_dof_to_local)
    physics_rhs_norm += f_vec[dof] * f_vec[dof];

  if (m_aux_eqns_combined_system && am_i_last_rank)
  {
    auto aux_eqns = m_physics_model->getAuxEquations();
    for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
      for (auto& val : f_aux->getVector(block))
        physics_rhs_norm += val*val;
  }

  return physics_rhs_norm;
}

void CrankNicolsonFunction::computeTimeTerm(DiscVectorPtr u_disc_vec, AuxiliaryEquationsStoragePtr u_aux,
                                            DiscVectorPtr Mdelta_u, AuxiliaryEquationsStoragePtr Mdelta_u_aux)
{
  auto aux_eqns = m_physics_model->getAuxEquations();

  // compute M * (u_np1 - u_n)
  auto& u_n_vec     = m_un->getVector();
  auto& u_np1_vec   = u_disc_vec->getVector();
  auto& delta_u_vec = m_delta_u->getVector();
  for (int i=0; i < delta_u_vec.shape()[0]; ++i)
    delta_u_vec[i] = u_np1_vec[i] - u_n_vec[i];
  m_delta_u->markVectorModified();

  if (m_aux_eqns_combined_system)
    for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
    {
      auto& u_n_aux_vec     = m_u_aux_n->getVector(block);
      auto& u_np1_aux_vec   = u_aux->getVector(block);
      auto& delta_u_aux_vec = m_delta_u_aux->getVector(block);
      for (int i=0; i < aux_eqns->getBlockSize(block); ++i)
        delta_u_aux_vec[i] = u_np1_aux_vec[i] - u_n_aux_vec[i];
    }

  m_physics_model->applyMassMatrix(m_delta_u, Mdelta_u);
  if (!Mdelta_u->isVectorCurrent())
    Mdelta_u->syncArrayToVector();

  if (m_aux_eqns_combined_system)
    for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
      aux_eqns->multiplyMassMatrix(block, m_tnp1, m_delta_u_aux->getVector(block), Mdelta_u_aux->getVector(block));

}

void CrankNicolsonFunction::combineTerms(DiscVectorPtr Mdelta_u, AuxiliaryEquationsStoragePtr Mdelta_u_aux,
                                         DiscVectorPtr f_disc_vec, AuxiliaryEquationsStoragePtr f_aux)
{
  auto aux_eqns = m_physics_model->getAuxEquations();

  // compute f_np1 = M * (u_np1 - u_n) / delta_t - 0.5 * f(u_n, t_n) - 0.5 * f(u_np1, t_np1)
  auto& Mdelta_u_vec  = Mdelta_u->getVector();
  auto& f_np1_vec     = f_disc_vec->getVector();
  auto& f_n_vec       = m_fn->getVector();
  Real delta_t_inv    = 1.0/(m_tnp1 - m_tn);
  for (int i=0; i < f_np1_vec.shape()[0]; ++i)
  {
    //std::cout << "i = " << i << ", delta_t_inv = " << delta_t_inv << ", Mdelta_u = "
    //         << Mdelta_u_vec[i] << ", f_np1 = " << f_np1_vec[i] << ", f_n_vec = " << f_n_vec[i] << std::endl;
    f_np1_vec[i] = delta_t_inv * Mdelta_u_vec[i] - 0.5*f_np1_vec[i] - 0.5*f_n_vec[i];
  }

  if (m_aux_eqns_combined_system)
    for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
    {
      auto& Mdelta_u_vec = Mdelta_u_aux->getVector(block);
      auto& f_np1_vec    = f_aux->getVector(block);
      auto& f_n_vec      = m_fn_aux->getVector(block);
      for (int i=0; i < aux_eqns->getBlockSize(block); ++i)
        f_np1_vec[i] = delta_t_inv * Mdelta_u_vec[i] - 0.5 * f_np1_vec[i] - 0.5*f_n_vec[i];
    }  
}

Real CrankNicolsonFunction::computeFunc(const ArrayType<Real, 1>& u_np1, AuxiliaryEquationsStoragePtr u_aux_np1,
                                        bool compute_norm, ArrayType<Real, 1>& f_np1)
{
  //TODO: add flag for when u_np1 == un, avoid computing M * (u_np1 - u_n) on first iteration
  assertAlways(m_tnp1 - m_tn > 1e-12, "delta_t must be > 1e-12");

  //std::cout << "evaluating CN function" << std::endl;
  //TODO: cache these
  //bool am_i_last_rank = commRank(MPI_COMM_WORLD) == (commSize(MPI_COMM_WORLD) - 1);
  auto aux_eqns = m_aux_eqns_combined_system ? m_physics_model->getAuxEquations() : nullptr;
  auto u_disc_vec = makeDiscVector(m_physics_model->getDiscretization());
  auto f_disc_vec = makeDiscVector(m_physics_model->getDiscretization());
  AuxiliaryEquationsStoragePtr u_aux, f_aux;
  if (m_aux_eqns_combined_system)
  {
    u_aux = makeAuxiliaryEquationsStorage(aux_eqns);
    f_aux = makeAuxiliaryEquationsStorage(aux_eqns);
    splitVector(u_np1, u_disc_vec, u_aux);
  } else
  {
    copyToVector(u_np1, u_disc_vec);
    u_aux = u_aux_np1;
  }

  computePhysicsRhs(u_disc_vec, u_aux, f_disc_vec, f_aux);

  if (!m_un->isVectorCurrent())
    m_un->syncArrayToVector();

  Real physics_rhs_norm = 0;
  if (compute_norm)
  {
    physics_rhs_norm = computeNorm(f_disc_vec, f_aux);
  }

  computeTimeTerm(u_disc_vec, u_aux, m_Mdelta_u, m_Mdelta_u_aux);

  combineTerms(m_Mdelta_u, m_Mdelta_u_aux, f_disc_vec, f_aux);

  if (m_aux_eqns_combined_system)
    combineVector(f_disc_vec, f_aux, f_np1);
  else
    copyFromVector(f_disc_vec, f_np1);

  Real norm = 0;
  if (compute_norm)
  {
    norm = computeNorm(f_disc_vec, f_aux);

    
    std::array<Real, 2> norms_local = {norm, physics_rhs_norm}, norms_global;
    MPI_Allreduce(norms_local.data(), norms_global.data(), 2, REAL_MPI_DATATYPE, MPI_SUM, MPI_COMM_WORLD);
    norm = std::sqrt(norms_global[0]);
    m_last_physics_rhs_norm = std::sqrt(norms_global[1]);
  }

  return norm;
}

// compute jac = df/du, overwriting jac
void CrankNicolsonFunction::computeJacobian(const ArrayType<Real, 1>& u, AuxiliaryEquationsStoragePtr u_aux_vec, linear_system::LargeMatrixPtr jac)
{
  auto aux_eqns = m_physics_model->getAuxEquations();
  auto u_disc_vec = makeDiscVector(m_physics_model->getDiscretization());  //TODO
  AuxiliaryEquationsStoragePtr u_aux      = makeAuxiliaryEquationsStorage(aux_eqns);

  if (m_precompute_linear_jacobian && !m_is_first_jacobian && !m_mass_matrix)
    precomputeLinearJacobian(jac);

  bool use_precomputed_jacobian = m_precompute_linear_jacobian && !m_is_first_jacobian;

  if (m_aux_eqns_combined_system)
  {
    u_aux = makeAuxiliaryEquationsStorage(aux_eqns);    
    splitVector(u, u_disc_vec, u_aux);
  } else
  {
    copyToVector(u, u_disc_vec);
    u_aux = u_aux_vec;
  }

  if (use_precomputed_jacobian)
  {
    jac->axpy(1.0/(m_tnp1 - m_tn), m_mass_matrix);
    jac->axpy(-0.5, m_linear_stiffness_matrix);
  }

  if (use_precomputed_jacobian)
  {
    m_assembler->setAlpha(-0.5);
    m_physics_model->computeJacobian(u_disc_vec, u_aux, m_tnp1, m_assembler, JacobianTerms::Nonlinear);
  } else
  {
    m_assembler->setAlpha(-0.5);
    m_physics_model->computeJacobian(u_disc_vec, u_aux, m_tnp1, m_assembler);
  }

  if (m_aux_eqns_combined_system)
  {
    m_augmented_assembler->setAlpha(-0.5);
    aux_eqns->computeJacobian(u_disc_vec, u_aux, m_tnp1, m_augmented_assembler);

    m_augmented_assembler->startAssembly();
    m_augmented_assembler->finishAssembly();
  }


  if (!use_precomputed_jacobian)
  {
    m_assembler->setAlpha(1.0/(m_tnp1 - m_tn));
    m_physics_model->computeMassMatrix(m_assembler);
  }

  if (m_aux_eqns_combined_system)
  {
    m_augmented_assembler->setAlpha(1.0/(m_tnp1 - m_tn));

    int start_dof = 0;
    for (int block=1; block < aux_eqns->getNumBlocks(); ++block)
    {
      int block_size = aux_eqns->getBlockSize(block);
      auto jac = std::make_shared<linear_system::LargeMatrixDense>(block_size, block_size, linear_system::LargeMatrixOpts());
      auto assembler = std::make_shared<linear_system::SimpleAssembler>(jac);
      aux_eqns->computeMassMatrix(block, m_tnp1, assembler);

      ArrayType<Real, 2> vals(boost::extents[block_size][block_size]);
      for (int i=0; i < block_size; ++i)
        for (int j=0; j < block_size; ++j)
          vals[i][j] = (*jac)(i, j);


      std::vector<DofInt> augmented_rows(block_size);
      for (int i=0; i < block_size; ++i)
        augmented_rows[i] = i + start_dof;

      if (m_augmented_assembler->amILastRank())
        m_augmented_assembler->assembleAugmentedValuesDiag(augmented_rows, augmented_rows, vals);
      start_dof += block_size;
    }
  }


  m_assembler->setAlpha(1);
  m_augmented_assembler->setAlpha(1);
  m_is_first_jacobian = false;
}

void CrankNicolsonFunction::setTnp1(const ArrayType<Real, 1>& u_n, AuxiliaryEquationsStoragePtr u_aux_n, Real t_np1)
{
  m_tn = m_tnp1;
  m_tnp1 = t_np1;

  if (m_aux_eqns_combined_system)
  {
    splitVector(u_n, m_un, m_u_aux_n);
  } else
  {
    copyToVector(u_n, m_un);
    *m_u_aux_n = *u_aux_n;
  }

  m_aux_eqns->setTnp1(m_un, u_aux_n, t_np1);
} 

void CrankNicolsonFunction::precomputeLinearJacobian(linear_system::LargeMatrixPtr jac)
{
  m_mass_matrix = jac->clone();
  m_linear_stiffness_matrix = jac->clone();

  auto assembler_mass = std::make_shared<linear_system::Assembler>(m_physics_model->getDiscretization(), m_mass_matrix);
  m_physics_model->computeMassMatrix(assembler_mass);
  
  auto assembler_jac = std::make_shared<linear_system::Assembler>(m_physics_model->getDiscretization(), m_linear_stiffness_matrix);
  // values in m_un and m_u_aux_n are irrelevent because we are only computing the linear part of the Jacobian
  m_physics_model->computeJacobian(m_un, m_u_aux_n, 0.0, assembler_jac, JacobianTerms::Linear);

  m_mass_matrix->finishMatrixAssembly();
  m_linear_stiffness_matrix->finishMatrixAssembly();
}

}