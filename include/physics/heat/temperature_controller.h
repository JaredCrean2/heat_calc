#ifndef HEAT_CALC_PHYSICS_HEAT_TEMPERATURE_CONTROLLER_H
#define HEAT_CALC_PHYSICS_HEAT_TEMPERATURE_CONTROLLER_H

#include "ProjectDefs.h"
#include "error_handling.h"

namespace Heat {

class TemperatureController
{
  public:
    virtual ~TemperatureController() = default;
    
    virtual Real get_value(Real temp) = 0;

    virtual Real get_value_dot(Real temp, Real temp_dot) = 0;

};

class TemperatureControllerConstant : public TemperatureController
{
  public:
    Real get_value(Real temp) override { return 1; }

    Real get_value_dot(Real temp, Real temp_dot) override { return 0; }

};

class TemperatureControllerHeatQuadratic : public TemperatureController
{
  public:
    TemperatureControllerHeatQuadratic(Real t_lower, Real t_upper) :
      m_t_mid((t_upper + t_lower)/2),
      m_delta_t(t_upper - t_lower)
    {
      assertAlways(t_upper >= t_lower, "t_upper must be >= t_lower");
    }

    Real get_value(Real temp) override
    {
      if (temp >= m_t_mid)
        return 0;
      else if (temp < m_t_mid - m_delta_t)
        return 1;
      else
        return std::pow(temp - m_t_mid, 2)/(m_delta_t * m_delta_t);
    }

    Real get_value_dot(Real temp, Real temp_dot) override
    {
      if (temp >= m_t_mid)
        return 0;
      else if (temp < m_t_mid - m_delta_t)
        return 0;
      else
        return 2 * (temp - m_t_mid) * temp_dot / (m_delta_t * m_delta_t);
    }


  private:
    Real m_t_mid;
    Real m_delta_t;

};
}

#endif