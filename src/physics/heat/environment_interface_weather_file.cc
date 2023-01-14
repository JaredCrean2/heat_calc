#include "physics/heat/environment_interface_weather_file.h"

namespace Heat {

// t is in seconds since the starting time
EnvironmentData EnvironmentInterfaceWeatherFile::getEnvironmentData(Real t)
{
  return interpolateData(t);
}

void EnvironmentInterfaceWeatherFile::computeSpacing()
{
  assertAlways(m_data.size() > 1, "Must have at least 2 data points");

  Real t0 = getJulianDate(m_data[0]);
  Real t1 = getJulianDate(m_data[1]);
  m_data_spacing = t1 - t0;

  checkUniformSpacing();
}

void EnvironmentInterfaceWeatherFile::checkUniformSpacing()
{
  const double tol = 1e-8;
  for (size_t i=2; i < m_data.size(); ++i)
  {
    Real t0 = getJulianDate(m_data[i-1]);
    Real t1 = getJulianDate(m_data[i]);
    Real spacing = t1 - t0;
    if (std::abs(spacing - m_data_spacing) > tol)
      throw std::runtime_error("data points are not uniformly spaced");
  }
}

int EnvironmentInterfaceWeatherFile::getIdx(Real t)
{
  Real t_days = getTDays(t);
  return std::floor(t_days/m_data_spacing);
}

EnvironmentData EnvironmentInterfaceWeatherFile::interpolateData(Real t)
{
  int idx = getIdx(t);

  if (idx == m_data.size()-1)
    return convertToEnvData(m_data[idx]);

  Real t_days = getTDays(t);
  Real xi = (t_days - idx*m_data_spacing)/m_data_spacing;

  return (1 - xi)*convertToEnvData(m_data[idx]) + xi*convertToEnvData(m_data[idx+1]);

}

Real EnvironmentInterfaceWeatherFile::getTDays(Real t)
{
  return t/(24*60*60);
}
}