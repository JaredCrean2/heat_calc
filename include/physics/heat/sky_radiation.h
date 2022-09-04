#ifndef PHYSICS_HEAT_SKY_RADIATION_H
#define PHYSICS_HEAT_SKY_RADIATION_H

#include "ProjectDefs.h"
#include "utils/math.h"
#include <iostream>

namespace Heat {

// Accounts for radiation from environment (ground, sky), but not solar
// radiation
class SkyRadiationModel
{
  public:
    SkyRadiationModel(Real emittance, std::array<Real, 3> vertical_vector) :
      m_emittance(emittance),
      m_vertical_vector(vertical_vector / std::sqrt(dot(vertical_vector, vertical_vector)))
    {}

    // set the Horizontal Infrared Radiation intensity (W/m^2)
    void setIRHorizontalRadiation(Real flux) { m_t_sky4 = flux/m_sigma; }

    // set the air temperature (which the model uses as the temperature of the ground)
    void setAirTemperature(Real t_air) 
    {
      m_t_air = t_air;
      m_t_air4 = std::pow(t_air, 4);
    }

    Real computeFlux(Real t_surf, const std::array<Real, 3>& unit_normal);

    Real computeFluxdTwall(Real t_surf, const std::array<Real, 3>& unit_normal, Real& flux_dot);

    Real computeFluxdTair(Real t_surf, const std::array<Real, 3>& unit_normal, Real& flux_dot);

  private:

    Real computeCosTiltAngle(const std::array<Real, 3>& unit_normal)
    {
      return dot(unit_normal, m_vertical_vector);
    }

    Real m_sigma = 5.6697E-8;  // Stefan-Boltzmann constant
    Real m_emittance;
    Real m_t_air  = -1;
    Real m_t_air4 = -1;
    Real m_t_sky4 = -1;
    std::array<Real, 3> m_vertical_vector;

};

}  // namespace

#endif