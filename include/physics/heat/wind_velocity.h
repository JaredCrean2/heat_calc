#ifndef HEAT_WIND_VELOCITY_H
#define HEAT_WIND_VELOCITY_H

#include "ProjectDefs.h"
#include "utils/error_handling.h"

namespace Heat {

// Compute the local wind velocity from the meterological wind
// velocity measurement
// See Chapter 16 of the ASHRAE Handbook of Fundamentals, or
// the EnergyPlus Engineering Reference version 9.5, Section
// 3.5.4.2
class AshraeWindVelocity
{
  public:
    
    AshraeWindVelocity(int met_terrain_index, Real meterological_altitude, int local_terrain_index) :
      m_met_alpha(getAlpha(met_terrain_index)),
      m_met_delta(getBoundaryLayerThickness(met_terrain_index)),
      m_local_alpha(getAlpha(local_terrain_index)),
      m_local_delta(getBoundaryLayerThickness(local_terrain_index))
    {
      m_met_term = std::pow(m_met_delta / meterological_altitude, m_met_alpha);
    }

    Real computeLocalVelocity(Real met_velocity, Real local_altitude) const
    {
      return met_velocity * m_met_term * std::pow(local_altitude / m_local_delta, m_local_alpha);
    }


  private:
    // Gets the exponent for a given terrain type, see table 3.7 in EnergyPlus
    // Engineering reference.
    //   0: Flat, open country
    //   1: Rough, wooded country
    //   2: Towns and cities
    //   4: Ocean
    //   5: Urban, industrial, forest
    Real getAlpha(int terrain_index)
    {
      assertAlways(terrain_index >= 0 && terrain_index <= 4, "Terrain index must be in the range [0, 4]");
      static constexpr std::array<Real, 5> alpha_values{0.14, 0.22, 0.33, 0.10, 0.22};

      return alpha_values[terrain_index];
    }

    // Gets the boundary layer thickness for a given terrain type.
    // The returned value is in meters.
    Real getBoundaryLayerThickness(int terrain_index)
    {
      assertAlways(terrain_index >= 0 && terrain_index <= 4, "Terrain index must be in the range [0, 4]");
      static constexpr std::array<Real, 5> delta_values{270, 370, 460, 210, 370};

      return delta_values[terrain_index];
    }

    Real m_met_alpha;
    Real m_met_delta;
    Real m_local_alpha;
    Real m_local_delta;

    Real m_met_term; // (met_delta / z_met)^alpha_met


};
}

#endif