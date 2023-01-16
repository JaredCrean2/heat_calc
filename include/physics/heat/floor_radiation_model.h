#ifndef PHYSICS_HEAT_FLOOR_RADIATION_MODEL_H
#define PHYSICS_HEAT_FLOOR_RADIATION_MODEL_H

#include "ProjectDefs.h"
#include "solar_position.h"

#include <iostream>

namespace Heat
{

// models heating of the floor via solar radiation
// from windows
// This model takes the flux through the windows and distributes it evenly
// over the floor
class FloorRadiationModel
{
  public:
    FloorRadiationModel(Real window_area, std::array<Real, 3> window_normal, Real shgc,
                        Real floor_area, Real floor_absorptivity) :
      m_window_area(window_area),
      m_window_normal(window_normal),
      m_solar_heat_gain_coefficient(shgc),
      m_floor_area(floor_area),
      m_floor_absorptivity(floor_absorptivity)
    {
      m_window_normal = m_window_normal / std::sqrt(dot(window_normal, window_normal));
    }

    void setDirectNormalRadiation(Real flux) { m_direct_flux = flux; }

    void setDiffuseRadiation(Real flux)
    {
      m_diffuse_flux = flux;
    }

    void setSolarDirection(const DirectionCosines& cosines)
    {
      m_solar_direction_unit = std::array<Real, 3>{cosines.cs1, cosines.cs2, cosines.cs3};
    }

    Real computeFlux()
    {
      Real direction_val = std::max(dot(m_window_normal, m_solar_direction_unit), 0.0);
      Real net_flux = m_solar_heat_gain_coefficient * m_window_area * (direction_val * m_direct_flux + m_diffuse_flux);
      return m_floor_absorptivity * net_flux/m_floor_area;
    }

  private:
    Real m_window_area;
    std::array<Real, 3> m_window_normal;
    Real m_solar_heat_gain_coefficient;
    Real m_floor_area;
    Real m_floor_absorptivity;

    Real m_direct_flux = -1;
    Real m_diffuse_flux = -1;
    std::array<Real, 3> m_solar_direction_unit = {0, 0, 0};
};

}

#endif