#ifndef PHYSICS_HEAT_ENVIRONMENT_INTERFACE_H
#define PHYSICS_HEAT_ENVIRONMENT_INTERFACE_H

#include "ProjectDefs.h"

namespace Heat {

struct EnvironmentData
{
  Real air_temp;  // Kelvin
  Real air_speed;  // m/s
  std::array<Real, 3> air_direction;  // nondimensional unit vector
  Real ir_horizontal_radiation;  // W/m^2
  Real direct_normal_radiation;  // W/m^2
  Real diffuse_radiation;        // W/m^2
};

EnvironmentData operator+(const EnvironmentData& a, const EnvironmentData& b);

EnvironmentData operator-(const EnvironmentData& a, const EnvironmentData& b);

EnvironmentData operator*(const EnvironmentData& a, Real b);

EnvironmentData operator*(Real a, const EnvironmentData& b);

EnvironmentData operator/(const EnvironmentData& a, Real b);


class EnvironmentInterface
{
  public:
    virtual ~EnvironmentInterface() {}
    
    // t is in seconds since the starting time
    virtual EnvironmentData getEnvironmentData(Real t) = 0;
};


class EnvironmentInterfaceConstant : public EnvironmentInterface
{
  public:
    explicit EnvironmentInterfaceConstant(const EnvironmentData& data) :
      m_data(data)
    {}
    
    EnvironmentData getEnvironmentData(Real t) override { return m_data; }

  private:
    EnvironmentData m_data;
};

}

#endif