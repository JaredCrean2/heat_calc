#ifndef PHYSICS_HEAT_WINDOW_RADIATION_MODEL_H
#define PHYSICS_HEAT_WINDOW_RADIATION_MODEL_H

#include "ProjectDefs.h"
#include "error_handling.h"
#include "physics/heat/solar_position.h"

#include <iostream>  // TODO: DEBUGGING

namespace Heat
{

// takes the flux through the window and distributes it evenly over the floor
class WindowRadiationModel
{
  public:
    WindowRadiationModel(Real window_area, Real floor_area, Real solar_heat_gain_coefficient,
                         std::array<Real, 3> window_outward_normal,
                         Real shading_cos_zenith=-1) :
      m_window_area(window_area),
      m_floor_area(floor_area),
      m_solar_heat_gain_coefficient(solar_heat_gain_coefficient),
      m_shading_cos_zenith(shading_cos_zenith),
      m_window_outward_normal(window_outward_normal / std::sqrt(dot(window_outward_normal, window_outward_normal)))
    {
      assertAlways(solar_heat_gain_coefficient >= 0 && solar_heat_gain_coefficient <= 1, 
                   "solar heat gain coefficient must be in range [0, 1]");
    }

    void setIRHorizontalRadiation(Real flux) { m_ir_horizontal_radiation = flux; }

    void setDirectNormalRadiation(Real flux) { m_direct_normal_radiation = flux;}

    void setDiffuseRadiation(Real flux) { m_diffuse_radiation = flux; }

    void setSolarDirection(const DirectionCosines& cosines)
    { 
      m_solar_direction = {cosines.cs1, cosines.cs2, cosines.cs3};
      std::cout << "magnitude = " << std::sqrt(dot(m_solar_direction, m_solar_direction)); 
      assertAlways(std::sqrt(dot(m_solar_direction, m_solar_direction)) - 1 < 1e-12, "DirectionCosines should form a unit vector");
    }

    Real computeFlux()
    {
      //TODO: this function takes no arguments, cache the result?
      Real cos_theta = dot(m_window_outward_normal, m_solar_direction);
      //std::cout << "cos_theta = " << cos_theta, ", expected_val = " << 
      // this is a simplistic shading model: cos(angle between window normal and sun)
      // is negative, assume the window is in the shadow of the building.
      // For overhang shadowing, if sun zenith is above given angle, assume window is in shadow.
      // (this implicitly assumes the overhang is infinitely wide)
      bool in_building_shadow = cos_theta < 0;
      bool in_overhang_shadow = (m_shading_cos_zenith != -1 && m_solar_direction[2] > m_shading_cos_zenith);
      Real flux_direct = in_building_shadow || in_overhang_shadow ? 0 : m_direct_normal_radiation * cos_theta;
      Real flux_indirect = m_ir_horizontal_radiation + m_diffuse_radiation;
      Real effective_flux = flux_direct + flux_indirect;

      return m_solar_heat_gain_coefficient * m_window_area * effective_flux / m_floor_area;
    }


  private:
    Real m_window_area;
    Real m_floor_area;
    Real m_solar_heat_gain_coefficient;
    Real m_ir_horizontal_radiation;
    Real m_direct_normal_radiation;
    Real m_diffuse_radiation;
    Real m_shading_cos_zenith;  // at angles above this, the direct radiation 
                                // is blocked by the shading
    std::array<Real, 3> m_window_outward_normal;
    std::array<Real, 3> m_solar_direction;
};

}  // namespace

#endif