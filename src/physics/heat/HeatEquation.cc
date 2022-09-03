#include "physics/heat/HeatEquation.h"
#include "discretization/NeumannBC.h"
#include "discretization/disc_vector.h"
#include "discretization/dof_numbering.h"
#include "discretization/volume_discretization.h"
#include "physics/heat/basis_vals.h"
#include "physics/heat/helper_funcs.h"
#include "physics/heat/interior_temperature_update.h"
#include "mesh/mesh.h"

#include "linear_system/assembler.h"
#include <petscmat.h>

#include "physics/heat/mass_matrix.h"
#include "physics/heat/volume_term.h"
#include "physics/heat/source_term.h"
#include "physics/heat/dirichlet_bc.h"
#include "physics/heat/neumann_bc.h"

namespace Heat {

//-----------------------------------------------------------------------------
// HeatEquation

void HeatEquation::computeRhs(DiscVectorPtr u, const Real t, DiscVectorPtr rhs)
{
  rhs->set(0);
  if (!(u->isArrayCurrent()))
  {
    //std::cout << "syncing u vector to array" << std::endl;
    u->syncVectorToArray();
  }

  //std::cout << "u array = " << std::endl;
  //printArray(u);

  //std::cout << "before setting dirichlet values, u = " << std::endl;
  //printArray(u);

  // apply Dirichlet values to array
  applyDirichletValues(*this, t, u);
  //std::cout << "after setting dirichlet values, u = " << std::endl;
  //printArray(u);

  // compute volume terms
  computeVolumeTerm(*this, u, rhs);

  //std::cout << "\nafter volume term " << std::endl;
  //printArray(rhs);

  // compute Neumann BC terms
  computeNeumannBC(*this, t, u, rhs);

  // compute unsteady Dirichlet BC term
  //computeUnsteadyDirichletBC(*this, t, rhs);

  //std::cout << "\nafter Neumann term " << std::endl;
  //printArray(rhs);

  // compute source term
  computeSourceTerm(*this, t, rhs);

  //std::cout << "\nafter source term " << std::endl;
  //printArray(rhs);
  //printVector(rhs);
}


void HeatEquation::computeJacobian(DiscVectorPtr u, const Real t, linear_system::AssemblerPtr assembler)
{
  if (!u->isArrayCurrent())
    u->syncArrayToVector();

  applyDirichletValues(*this, t, u);

  computeVolumeJacobian(*this, u, assembler);

  // typical Neumann and source terms don't contribute to the Jacobian
  computeNeumannBCJacobian(*this, u, t, assembler);
}

void HeatEquation::applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out)
{
  if (!vec_in->isArrayCurrent())
    vec_in->syncVectorToArray();

  ::Heat::applyMassMatrix(*this, vec_in, vec_out);
}

void HeatEquation::computeMassMatrix(linear_system::AssemblerPtr assembler)
{
  ::Heat::computeMassMatrix(*this, assembler);
}


void HeatEquation::checkInitialization()
{
  PhysicsModel::checkInitialization();
  
  if (m_params.size() != getDiscretization()->getNumVolDiscs())
    throw std::runtime_error("Incorrect number of Heat::VolumeGroupParams.  Should be equal to number of volume groups");
}


}  // namespace 
