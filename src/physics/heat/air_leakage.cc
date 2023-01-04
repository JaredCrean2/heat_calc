#include "physics/heat/air_leakage.h"

namespace Heat {

Real AirLeakageModelPressure::computeAirLeakagePower(Real t_interior)
{
  Real ach_natural         = m_ach50 * m_expected_pressure / 50;
  Real air_volume_per_hour = ach_natural * m_volume;
  return m_rho * m_cp * air_volume_per_hour * (t_interior - m_t_exterior) / 3600;
}

Real AirLeakageModelPressure::computeAirLeakagePowerDot(Real t_interior, Real& flux_dot)
{
  Real ach_natural         = m_ach50 * m_expected_pressure / 50;
  Real air_volume_per_hour = ach_natural * m_volume;
  flux_dot =  m_rho * m_cp * air_volume_per_hour/ 3600;      
  return m_rho * m_cp * air_volume_per_hour * (t_interior - m_t_exterior) / 3600;
}

Real HRVModel::computeAirLeakagePower(Real t_interior)
{
  return (1 - m_efficiency) * m_rho * m_cp * m_flow_rate * (t_interior - m_t_exterior);
}

Real HRVModel::computeAirLeakagePowerDot(Real t_interior, Real& flux_dot)
{
  flux_dot = (1 - m_efficiency) * m_rho * m_cp * m_flow_rate;
  return (1 - m_efficiency) * m_rho * m_cp * m_flow_rate * (t_interior - m_t_exterior); 
}

}