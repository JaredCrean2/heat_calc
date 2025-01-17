#ifndef PHYSICS_HEAT_HEAT_EQN_SOLAR_H
#define PHYSICS_HEAT_HEAT_EQN_SOLAR_H

#include "physics/PhysicsModel.h"
#include "physics/heat/HeatEquation.h"
#include "physics/heat/source_terms_def.h"


namespace Heat {

class InteriorAirTemperatureUpdator;

class AuxiliaryEquationsSolar;
using AuxiliaryEquationsSolarPtr = std::shared_ptr<AuxiliaryEquationsSolar>;

class HeatEquationSolar : public HeatEquation
{
  public:

    HeatEquationSolar(DiscPtr disc, std::shared_ptr<SolarPositionCalculator> solar_position,
                      std::shared_ptr<EnvironmentInterface> environment_interface,
                      std::shared_ptr<InteriorAirTemperatureUpdator> air_temp_updator,
                      MPI_Comm comm=MPI_COMM_WORLD);

    using HeatEquation::initialize; 
    using HeatEquation::addNeumannBC;

    void initialize() override;

    // is_exterior: true if this BC uses the environment air temperature and is exposed to
    //              thermal radiation
    void addNeumannBC(NeumannBCPtr bc, bool is_exterior);

    void setTimeParameters(Real t, Real interior_air_temp);

    const EnvironmentData& getEnvironmentData() const { return m_env_data; }

    void computeRhs(DiscVectorPtr u,  AuxiliaryEquationsStoragePtr u_aux, const Real t, DiscVectorPtr rhs) override
    {
      setTimeParameters(t, u_aux->getVector(1)[0]);
      HeatEquation::computeRhs(u, u_aux, t, rhs);
    }
        
    void computeJacobian(DiscVectorPtr u,  AuxiliaryEquationsStoragePtr u_aux, const Real t, linear_system::AssemblerPtr assembler,
                         JacobianTerms terms=JacobianTerms::All) override
    {
      setTimeParameters(t, u_aux->getVector(1)[0]);
      HeatEquation::computeJacobian(u, u_aux, t, assembler, terms);
    }

    // compute Jacobian-vector product dR/dT, where T is the air temperature
    void computedRdTinterior_airProduct(DiscVectorPtr u, Real interior_temp, Real t, Real x, ArrayType<Real, 1>& b);

    AuxiliaryEquationsPtr getAuxEquations() override;

    AuxiliaryEquationsSolarPtr getAuxEquationsSolar();

    bool isNeumannBCExterior(int i) const { return m_is_neumann_bc_exterior[i]; }

  private:
    std::vector<bool> m_is_neumann_bc_exterior;
    std::shared_ptr<SolarPositionCalculator> m_solar_position;
    std::shared_ptr<EnvironmentInterface> m_environment;
    std::shared_ptr<InteriorAirTemperatureUpdator> m_air_temp;
    EnvironmentData m_env_data;
    AuxiliaryEquationsSolarPtr m_aux_equations;
};

void computeNeumannBC_dotTair(const HeatEquationSolar& physics, const Real t, DiscVectorPtr u,
                              Real t_interior, Real t_interior_dot, 
                              DiscVectorPtr rhs, DiscVectorPtr rhs_dot);

void computeNeumannBC_dotTair(NeumannBCPtr bc, DiscVectorPtr u,
                              Real t_interior, Real t_interior_dot,
                              const Real t, DiscVectorPtr rhs, DiscVectorPtr rhs_dot);

void computeSourceTerm_dotTair(const HeatEquation& physics, Real t, Real t_interior_dot,
                               DiscVectorPtr rhs, DiscVectorPtr rhs_dot);

void computeSourceTerm_dotTair(const VolDiscPtr vol_disc, SourceTermAirWindSkyPtr src, Real t, Real t_interior_dot,
                               ArrayType<Real, 2>& rhs_arr, ArrayType<Real, 2>& rhs_arr_dot);                              

}

#endif