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

//-----------------------------------------------------------------------------
// HeatEquationSolar

void HeatEquationSolar::initialize(DiscVectorPtr sol_vec, Real t_start) 
{ 
  HeatEquation::initialize();
  
  std::vector<NeumannBCPtr> interior_bcs;
  for (size_t i=0; i < getNeumannBCs().size(); ++i)
    if (!m_is_neumann_bc_exterior[i])
      interior_bcs.push_back(getNeumannBCs()[i]);

  m_air_temp->initialize(this, interior_bcs, sol_vec, t_start); 
}


void HeatEquationSolar::addNeumannBC(NeumannBCPtr bc, bool is_exterior)
{
  addNeumannBC(bc); 
  m_is_neumann_bc_exterior.push_back(is_exterior);
}


void HeatEquationSolar::setTimeParameters(Real t)
{
  DirectionCosines solar_dir = m_solar_position.computePosition(t);
  EnvironmentData env_data   = m_environment->getEnvironmentData(t);
  const auto& neumann_bcs    = getNeumannBCs();
  for (size_t i=0; i < neumann_bcs.size(); ++i)
  {
    auto bc_air_wind_sky = std::dynamic_pointer_cast<AirWindSkyNeumannBC>(neumann_bcs[i]);
    if (bc_air_wind_sky)
    {
      if (m_is_neumann_bc_exterior[i])
      {
        bc_air_wind_sky->setAirTemperature(env_data.air_temp);
        bc_air_wind_sky->setAirSpeed(env_data.air_speed);
        bc_air_wind_sky->setAirDirection(env_data.air_direction);
        bc_air_wind_sky->setIRHorizontalRadiation(env_data.ir_horizontal_radiation);
        bc_air_wind_sky->setDirectNormalRadiation(env_data.direct_normal_radiation);
        bc_air_wind_sky->setDiffuseRadiation(env_data.diffuse_radiation);
        bc_air_wind_sky->setSolarDirection(solar_dir);
      } else
      {
        bc_air_wind_sky->setAirTemperature(m_air_temp->getTemperature());
        bc_air_wind_sky->setAirSpeed(0);
        bc_air_wind_sky->setAirDirection(std::array<Real, 3>{1, 0, 0});
        bc_air_wind_sky->setIRHorizontalRadiation(0);
        bc_air_wind_sky->setDirectNormalRadiation(0);
        bc_air_wind_sky->setDiffuseRadiation(0);
        bc_air_wind_sky->setSolarDirection(solar_dir);         
      }
    }
  }

  m_env_data = env_data;
}

void HeatEquationSolar::updateDependentQuantities(DiscVectorPtr u, Real t)
{
  m_air_temp->updateTemperature(u, t);
}

void HeatEquationSolar::completeTimestep(DiscVectorPtr u, Real t)
{
  m_air_temp->startNewTimestep(u, t);
}

}  // namespace 
