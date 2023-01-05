#ifndef HEAT_TARP_MODEL_H
#define HEAT_TARP_MODEL_H

#include "ProjectDefs.h"
#include "error_handling.h"
#include "utils/math.h"
#include "wind_velocity.h"

#include <iostream>

// TARP (Thermal Analysis Research Program) model for computing heat
// transfer coefficients
// See Section 3.5.5.2 of EnergyPlus Engineering reference Version 9.5

namespace Heat {

class TarpModel
{
  public:
    TarpModel(Real surface_area, Real perimeter, int roughness_index, const std::array<Real, 3>& vertical_vector,
              const std::array<Real, 3>& point_at_zero_altitude,
              int met_terrain_index, Real meterological_altitude, int local_terrain_index ) :
      m_surface_area(surface_area),
      m_surface_perimeter(perimeter),
      m_roughness_fac(getRoughnessFactor(roughness_index)),
      m_vertical_vector(vertical_vector / std::sqrt(dot(vertical_vector, vertical_vector))),
      m_point_at_zero_altitude(point_at_zero_altitude),
      m_local_wind_velocity_calc(met_terrain_index, meterological_altitude, local_terrain_index),
      m_alpha(std::pow(m_quadratic_intercept, -5.0/3.0))
    {}

    void setAirTemperature(Real temp) { m_air_temp = temp; }

    Real getAirTemperature() const { return m_air_temp; }

    void setAirSpeed(Real velocity) { m_air_speed = velocity; }

    Real getAirSpeed() const { return m_air_speed; }

    void setAirDirection(std::array<Real, 3> direction)
    {
      m_air_direction = direction / std::sqrt(dot(direction, direction));
    }

    // Get the roughness factor for the given roughness index, which identifies
    // how rough the surface of the material is.
    // From Table 3.7 in the EnergyPlus Engineering reference, example materials
    // for each roughness index are:
    //    0: Stucco
    //    1: Brick
    //    2: Concrete
    //    3: Clear pine
    //    4: Smooth Plaster
    //    5: Glass
    Real getRoughnessFactor(int roughness_index)
    {
      assertAlways(roughness_index >= 0 && roughness_index <= 5, "roughness index must be in range [0, 5]");
      static constexpr std::array<Real, 6> roughness_factors{2.17, 1.67, 1.52, 1.13, 1.11, 1.00};
      return roughness_factors[roughness_index];
    }

    Real computeHeatTransferCoeff(Real wall_temp, const std::array<Real, 3>& pt, const std::array<Real, 3>& unit_normal);

    // returns the heat transfer coefficient, and h_n_dot is overwritten with its derivative
    // with respect to the surface temperature
    Real computeHeatTransferCoeffdTwall(Real wall_temp, const std::array<Real, 3>& pt, const std::array<Real, 3>& unit_normal, Real& h_n_dot);

    Real computeHeatTransferCoeffdTair(Real wall_temp, const std::array<Real, 3>& pt, const std::array<Real, 3>& unit_normal, Real& h_n_dot);

    Real computeCosTiltAngle(const std::array<Real, 3>& unit_normal);

    Real computeLocalWindSpeed(const std::array<Real, 3>& pt);

  private:

    Real computeHeatTransferCoeffDeriv_impl(Real delta_t, Real delta_t_dot, const std::array<Real, 3>& pt, 
                                            const std::array<Real, 3>& unit_normal, Real& h_n_dot);

    Real m_air_temp     = std::numeric_limits<Real>::min();
    Real m_air_speed = std::numeric_limits<Real>::min();
    std::array<Real, 3>  m_air_direction{0, 0, 0};
    Real m_surface_area;
    Real m_surface_perimeter;
    Real m_roughness_fac;
    std::array<Real, 3> m_vertical_vector{0, 0, 1}; // vector that defines an upward facing surface, used in the
                                                    // tilt angle
    std::array<Real, 3> m_point_at_zero_altitude;
    AshraeWindVelocity m_local_wind_velocity_calc;
    const Real m_quadratic_intercept = 0.01;  // the derivative of delta_t^(1/3) does not exist
                                              // at delta_t = 0.  Use this parameter to fit a quadratic
                                              // at delta_t = 0 and delta_t = m_quadratic_intercept
    const Real m_alpha;

};

}

#endif