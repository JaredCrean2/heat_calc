#include "physics/heat/environment_interface.h"

#include "utils/math.h"


namespace Heat {

EnvironmentData operator+(const EnvironmentData& a, const EnvironmentData& b)
{
  EnvironmentData c;
  c.air_temp                = a.air_temp + b.air_temp;
  c.air_speed               = a.air_speed + b.air_speed;
  c.air_direction           = a.air_direction + b.air_direction;
  c.ir_horizontal_radiation = a.ir_horizontal_radiation + b.ir_horizontal_radiation;
  c.direct_normal_radiation = a.direct_normal_radiation + b.direct_normal_radiation;
  c.diffuse_radiation       = a.diffuse_radiation + b.diffuse_radiation;

  return c;
}

EnvironmentData operator-(const EnvironmentData& a, const EnvironmentData& b)
{
  EnvironmentData c;
  c.air_temp                = a.air_temp - b.air_temp;
  c.air_speed               = a.air_speed - b.air_speed;
  c.air_direction           = a.air_direction - b.air_direction;
  c.ir_horizontal_radiation = a.ir_horizontal_radiation - b.ir_horizontal_radiation;
  c.direct_normal_radiation = a.direct_normal_radiation - b.direct_normal_radiation;
  c.diffuse_radiation       = a.diffuse_radiation - b.diffuse_radiation;

  return c;
}

EnvironmentData operator*(const EnvironmentData& a, Real b)
{
  EnvironmentData c;
  c.air_temp                = a.air_temp*b;
  c.air_speed               = a.air_speed*b;
  c.air_direction           = a.air_direction*b;
  c.ir_horizontal_radiation = a.ir_horizontal_radiation*b;
  c.direct_normal_radiation = a.direct_normal_radiation*b;
  c.diffuse_radiation       = a.diffuse_radiation*b;

  return c;
}

EnvironmentData operator*(Real a, const EnvironmentData& b)
{
  return b * a;
}

EnvironmentData operator/(const EnvironmentData& a, Real b)
{
  return a * (1/b);
}

}