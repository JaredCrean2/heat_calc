#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

#include "lagrange_basis.h"
#include "physics/PhysicsModel.h"
#include "physics/heat/air_leakage.h"
#include "physics/heat/bc_defs.h"
#include "physics/heat/solar_position.h"
#include "quadrature.h"
#include <memory>

namespace Heat {

class VolumeGroupParams
{
  public:
    VolumeGroupParams (Real kappa, Real rho, Real Cp) :
      kappa(kappa),
      rho(rho),
      Cp(Cp)
    {}

    Real kappa; // thermal conductivity, W/(m K)
    Real rho;   // density (kg/m^2)
    Real Cp;    // specific heat capacityh (J/(kg K))
};

// solves heat equation rho Cp dT/dt = d/dx (kappa dT/dx) + S
// and S is a source term.  This module implements the
// two terms on the right hand side
// In weak form:
//  \int w rho Cp dT/dt = \int w d/dx (kappa dT/dx) + \int w S dOmega
//  \int w rho Cp dT/dt dOmega = - \int dw/dx kappa dT/dx dOmega + \int w kappa dT/dx n_i dGamma + \int w S dOmega = 0
//  \int w rho Cp dT/dt dOmega = - \int dw/dx kappa dT/dx dOmega + \int w h dGamma_h + \int w S dOmega = 0
// where h is the prescibes heat flux, Gamma is the boundary, and Gamma_h is the portion of the boundary
// where the heat flux is prescribed
// Note that w = 0 on Gamma_d, the portion of the boundary where the solution is specified, and
// Gamma_d \union Gamma_h = Gamma
class HeatEquation : public PhysicsModel
{
  public:
    HeatEquation(DiscPtr disc)
    : PhysicsModel(disc)
    {}

    void computeRhs(DiscVectorPtr u, const Real t, DiscVectorPtr rhs) override;
        
    void computeJacobian(DiscVectorPtr u, const Real t, linear_system::AssemblerPtr assembler) override;

    void applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out) override;

    void computeMassMatrix(linear_system::AssemblerPtr assembler) override;

    void addVolumeGroupParams(const VolumeGroupParams& params) { m_params.push_back(params); }

    const VolumeGroupParams& getVolumeGroupParams(int idx) const { return m_params.at(idx); }

  private:
    void checkInitialization() override;  // check everything was initialized


    std::vector<VolumeGroupParams> m_params;
};

class SolarPositionCalculator
{
  public:
    SolarPositionCalculator(Real julian_date_start, int time_zone, Real longitude, Real latitude) :
      m_time_zone(time_zone),
      m_longitude(longitude),
      m_latitude(latitude)
    {}

    DirectionCosines computePosition(int t)
    {
      Real julian_date = m_julian_date_start + t;
      int julian_date_whole = static_cast<int>(julian_date);
      Real time_hours = julian_date_whole - julian_date;
      return solar::computeDirectionCosines(julian_date_whole, time_hours, m_time_zone, m_longitude, m_latitude);
    }

  private:
    Real m_julian_date_start;
    int m_time_zone;
    Real m_longitude;
    Real m_latitude;
};

struct EnvironmentData
{
  Real air_temp;  // Kelvin
  Real air_speed;  // m/s
  std::array<Real, 3> air_direction;  // nondimensional unit vector
  Real ir_horizontal_radiation;  // W/m^2
  Real direct_normal_radiation;  // W/m^2
  Real diffuse_radiation;        // W/m^2
};

class EnvironmentInterface
{
  public:
    virtual EnvironmentData getEnvironmentData(Real t) = 0;
};


class InteriorAirTemperatureUpdator;

class HeatEquationSolar : public HeatEquation
{
  public:

    HeatEquationSolar(DiscPtr disc, SolarPositionCalculator solar_position, std::shared_ptr<EnvironmentInterface> environment_interface,
                      std::shared_ptr<InteriorAirTemperatureUpdator> air_temp_updator)
    : HeatEquation(disc),
      m_solar_position(solar_position),
      m_environment(environment_interface),
      m_air_temp(air_temp_updator)
    {}

    using HeatEquation::addNeumannBC;

    // is_exterior: true if this BC uses the environment air temperature and is exposed to
    //              thermal radiation
    void addNeumannBC(NeumannBCPtr bc, bool is_exterior);

    void setTimeParameters(Real t);

    const EnvironmentData& getEnvironmentData() const { return m_env_data; }

    void computeRhs(DiscVectorPtr u, const Real t, DiscVectorPtr rhs) override
    {
      setTimeParameters(t);
      HeatEquation::computeRhs(u, t, rhs);
    }
        
    void computeJacobian(DiscVectorPtr u, const Real t, linear_system::AssemblerPtr assembler) override
    {
      setTimeParameters(t);
      HeatEquation::computeJacobian(u, t, assembler);
    }

  private:
    std::vector<bool> m_is_neumann_bc_exterior;
    SolarPositionCalculator m_solar_position;
    std::shared_ptr<EnvironmentInterface> m_environment;
    std::shared_ptr<InteriorAirTemperatureUpdator> m_air_temp;
    EnvironmentData m_env_data;
};


} // namespace

#endif
