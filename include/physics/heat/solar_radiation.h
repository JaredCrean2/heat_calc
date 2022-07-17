#ifndef PHYSICS_HEAT_SOLAR_RADIATION_H
#define PHYSICS_HEAT_SOLAR_RADIATION_H

#include "ProjectDefs.h"
#include "solar_position.h"

namespace Heat {

class SolarRadiationModel
{
  public:

    SolarRadiationModel(Real absorptivity) :
      m_absorptivity(absorptivity)
    {}

    // sets the solar flux including both Direct Normal Radiation
    // and Diffuse Horzontal Radiation
    void setDirectNormalRadiation(Real flux) { m_direct_flux = flux; }

    void setDiffuseRadiation(Real flux) { m_diffuse_flux = flux; }

    void setSolarDirection(const AzimuthZenith& az)
    {
      Real cos_phi = az.cos_zenith;
      Real sin_phi = std::sqrt(1 - cos_phi*cos_phi);

      Real cos_psi = az.cos_azimuth;
      Real sin_psi = std::sqrt(1 - cos_psi*cos_psi);

      m_solar_direction_unit = std::array<Real, 3>{sin_phi * sin_psi, sin_phi * cos_psi, cos_phi};
    }

    Real computeFlux(const std::array<Real, 3>& unit_normal)
    {
      return m_absorptivity *(m_diffuse_flux * m_direct_flux * dot(unit_normal, m_solar_direction_unit));
    }

  private:
    Real m_absorptivity;
    Real m_direct_flux;
    Real m_diffuse_flux;
    std::array<Real, 3> m_solar_direction_unit;
};

}  // namespace

#endif