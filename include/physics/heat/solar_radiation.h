#ifndef PHYSICS_HEAT_SOLAR_RADIATION_H
#define PHYSICS_HEAT_SOLAR_RADIATION_H

#include "ProjectDefs.h"
#include "solar_position.h"

#include <iostream>

namespace Heat {

class SolarRadiationModel
{
  public:

    SolarRadiationModel(Real absorptivity) :
      m_absorptivity(absorptivity)
    {}

    void setDirectNormalRadiation(Real flux) { m_direct_flux = flux; }

    void setDiffuseRadiation(Real flux) { m_diffuse_flux = flux; }

    void setSolarDirection(const DirectionCosines& cosines)
    {
      m_solar_direction_unit = std::array<Real, 3>{cosines.cs1, cosines.cs2, cosines.cs3};
    }

    Real computeFlux(const std::array<Real, 3>& unit_normal)
    {
      //std::cout << "solar flux direction factor = " << std::max(dot(unit_normal, m_solar_direction_unit), 0.0) << ", flux = " <<  m_absorptivity * m_direct_flux * std::max(dot(unit_normal, m_solar_direction_unit), 0.0) << std::endl;
      return m_absorptivity * (m_diffuse_flux + m_direct_flux * std::max(dot(unit_normal, m_solar_direction_unit), 0.0));
    }

  private:
    Real m_absorptivity;
    Real m_direct_flux;
    Real m_diffuse_flux;
    std::array<Real, 3> m_solar_direction_unit;
};

}  // namespace

#endif