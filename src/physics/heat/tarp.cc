#include "physics/heat/tarp.h"

namespace Heat {

Real TarpModel::computeHeatTransferCoeff(Real wall_temp, const std::array<Real, 3>& pt, const std::array<Real, 3>& unit_normal)
{
  Real delta_t = m_air_temp - wall_temp;
  Real air_speed = computeLocalWindSpeed(pt);
  bool is_windward = dot(unit_normal, m_air_direction) < 0;

  Real Wf = is_windward ? 1.0 : 0.5;

  Real h_f = 2.537 * Wf * m_roughness_fac * std::sqrt(m_surface_perimeter * air_speed/m_surface_area);
  Real cos_eps = computeCosTiltAngle(unit_normal);
  Real delta_t_one_third = std::pow(std::abs(delta_t), 1.0/3.0);

  Real h_n = 0.0;
  if ((delta_t < 0 && cos_eps > 0) || (delta_t > 0 && cos_eps < 0)) {
    h_n = 9.482 * delta_t_one_third / (7.283 - std::abs(cos_eps));
  } else if ( (delta_t > 0 && cos_eps > 0) || (delta_t < 0 && cos_eps < 0)) {
    h_n = 1.810 * delta_t_one_third/ (1.382 + std::abs(cos_eps));
  } else {
    h_n = 1.31 * delta_t_one_third;
  }

  return h_f + h_n;
}

Real TarpModel::computeHeatTransferCoeffDeriv(Real wall_temp, const std::array<Real, 3>& pt, const std::array<Real, 3>& unit_normal, Real& h_n_dot)
{
  Real delta_t = m_air_temp - wall_temp;
  Real delta_t_dot = -1;

  Real air_speed = computeLocalWindSpeed(pt);
  bool is_windward = dot(unit_normal, m_air_direction) < 0;

  Real Wf = is_windward ? 1.0 : 0.5;
  Real h_f = 2.537 * Wf * m_roughness_fac * std::sqrt(m_surface_perimeter * air_speed/m_surface_area);
  Real cos_eps = computeCosTiltAngle(unit_normal);

  Real abs_delta_t = std::abs(delta_t);  
  Real abs_delta_t_dot = delta_t > 0 ? delta_t_dot : -delta_t_dot;

  Real delta_t_one_third = std::pow(abs_delta_t, 1.0/3.0);
  Real delta_t_one_third_dot = (1.0/3.0) * std::pow(abs_delta_t, -2.0/3.0) * abs_delta_t_dot;

  Real h_n = 0.0;
  h_n_dot = 0.0;
  if ((delta_t < 0 && cos_eps > 0) || (delta_t > 0 && cos_eps < 0))
  {
    h_n     = 9.482 * delta_t_one_third     / (7.283 - std::abs(cos_eps));
    h_n_dot = 9.482 * delta_t_one_third_dot / (7.283 - std::abs(cos_eps));
  } else if ( (delta_t > 0 && cos_eps > 0) || (delta_t < 0 && cos_eps < 0))
  {
    h_n     = 1.810 * delta_t_one_third     / (1.382 + std::abs(cos_eps));
    h_n_dot = 1.810 * delta_t_one_third_dot / (1.382 + std::abs(cos_eps));
  } else
  {
    h_n     = 1.31 * delta_t_one_third;
    h_n_dot = 1.31 * delta_t_one_third_dot;
  }

  return h_f + h_n;
}


Real TarpModel::computeCosTiltAngle(const std::array<Real, 3>& unit_normal)
{
  return dot(unit_normal, m_vertical_vector);
}

Real TarpModel::computeLocalWindSpeed(const std::array<Real, 3>& pt)
{
  auto r = pt - m_point_at_zero_altitude;
  auto local_height = dot(r, m_vertical_vector);
  return m_local_wind_velocity_calc.computeLocalVelocity(m_air_speed, local_height);
}

}