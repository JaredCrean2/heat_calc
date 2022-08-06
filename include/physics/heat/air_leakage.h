#ifndef PHYSICS_HEAT_AIR_LEAKAGE_H
#define PHYSICS_HEAT_AIR_LEAKAGE_H

#include "ProjectDefs.h"

namespace Heat {

class AirLeakageModel
{
  public:

    virtual ~AirLeakageModel() {}

    // compute the power of the air leakage (W).  Positive values mean
    // energy is removed from the interior air
    virtual Real computeAirLeakagePower(Real t_interior, Real t_exterior) = 0;
};

// implements a simple model: given ach50, compute ach natural as natural pressure / 50;
class AirLeakageModelPressure : public AirLeakageModel
{
  public:
    AirLeakageModelPressure(Real ach50, Real expected_pressure, Real volume, Real cp, Real rho) :
      m_ach50(ach50),
      m_expected_pressure(expected_pressure),
      m_volume(volume),
      m_cp(cp),
      m_rho(rho)
    {}

    Real computeAirLeakagePower(Real t_interior, Real t_exterior) override
    {
      Real ach_natural         = m_ach50 * m_expected_pressure / 50;
      Real air_volume_per_hour = ach_natural * m_volume;
      return m_rho * m_cp * air_volume_per_hour * (t_interior - t_exterior) / 3600;
    }

  private:
    Real m_ach50;
    Real m_expected_pressure;
    Real m_volume;
    Real m_cp;
    Real m_rho;
};


class HRVModel : public AirLeakageModel
{
  public:
    // flow rate in m^3/s
    HRVModel(Real flow_rate, Real efficiency, Real cp, Real rho) :
      m_flow_rate(flow_rate),
      m_efficiency(efficiency),
      m_cp(cp),
      m_rho(rho)
    {}

    Real computeAirLeakagePower(Real t_interior, Real t_exterior) override
    {
      return m_efficiency * m_rho * m_cp * m_flow_rate * (t_exterior - t_interior);
    }

  private:
    Real m_flow_rate;
    Real m_efficiency; 
    Real m_cp;
    Real m_rho;
};


}  // namespace
#endif