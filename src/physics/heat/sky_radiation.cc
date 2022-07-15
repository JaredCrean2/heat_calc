#include "physics/heat/sky_radiation.h"

namespace Heat {

Real SkyRadiationModel::computeFlux(Real t_surf, const std::array<Real, 3>& unit_normal)
{
  Real cos_phi = computeCosTiltAngle(unit_normal);
  Real beta     = std::sqrt(0.5*(1 + cos_phi));

  Real f_ground = 0.5 * (1 - cos_phi);
  Real f_sky    = beta * 0.5 * (1 + cos_phi);
  Real f_air    = (1 - beta) * 0.5 * (1 + cos_phi);

  Real t_surf4   = std::pow(t_surf, 4);
  Real t_ground4 = m_t_air4;

  Real flux = m_emittance * m_sigma * ( f_ground * (t_ground4 - t_surf4) +
                                        f_sky    * (m_t_sky4  - t_surf4) +
                                        f_air    * (m_t_air4  - t_surf4) );
  return flux;
}

Real SkyRadiationModel::computeFluxDeriv(Real t_surf, const std::array<Real, 3>& unit_normal, Real& flux_dot)
{
  Real cos_phi = computeCosTiltAngle(unit_normal);
  Real beta     = std::sqrt(0.5*(1 + cos_phi));

  Real f_ground = 0.5 * (1 - cos_phi);
  Real f_sky    = beta * 0.5 * (1 + cos_phi);
  Real f_air    = (1 - beta) * 0.5 * (1 + cos_phi);

  Real t_surf4   = std::pow(t_surf, 4);
  Real t_surf4_dot = 4 * std::pow(t_surf, 3);
  Real t_ground4 = m_t_air4;

  Real flux = m_emittance * m_sigma * ( f_ground * (t_ground4 - t_surf4) +
                                        f_sky    * (m_t_sky4  - t_surf4) +
                                        f_air    * (m_t_air4  - t_surf4) );

  flux_dot   = m_emittance * m_sigma * ( f_ground * (-t_surf4_dot) +
                                          f_sky    * (-t_surf4_dot) +
                                          f_air    * (-t_surf4_dot) );
  return flux;
}

}