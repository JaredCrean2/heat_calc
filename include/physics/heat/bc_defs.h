#ifndef HEAT_BC_DEFS_H
#define HEAT_BC_DEFS_H

#include "discretization/NeumannBC.h"
#include "discretization/surface_discretization.h"
#include "physics/heat/solar_position.h"
#include "tarp.h"
#include "sky_radiation.h"
#include "solar_radiation.h"

namespace Heat {

class NewtonCooling : public NeumannBC
{
  public:
    NewtonCooling(SurfDiscPtr surf, Real heat_transfer_coeff) :
      NeumannBC(surf, true),
      m_heat_transfer_coeff(heat_transfer_coeff)
    {}

    void setExternalTemperature(Real temp) { m_temp = temp; }

    Real getExternalTemperature() const { return m_temp; }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals);

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals_deriv);

  private:
    Real m_heat_transfer_coeff; // Units W/(m^2 K)
    Real m_temp = std::numeric_limits<Real>::min();
};


// Thermal Analysis Research Program (TARP) BC
class TarpBC : public NeumannBC
{
  public:
    TarpBC(SurfDiscPtr surf, Real surface_area, Real perimeter, int roughness_index, const std::array<Real, 3>& vertical_vector,
              const std::array<Real, 3>& point_at_zero_altitude,
              int met_terrain_index, Real meterological_altitude, int local_terrain_index ) : 
      NeumannBC(surf, true),
      m_tarp(surface_area, perimeter, roughness_index, vertical_vector, 
             point_at_zero_altitude, met_terrain_index, meterological_altitude, 
             local_terrain_index),
      m_quad_coords(boost::extents[m_surf->getNumQuadPtsPerFace()][3])
    {}

    void setAirTemperature(Real temp) { m_tarp.setAirTemperature(temp); }
    
    void setAirSpeed(Real velocity) { m_tarp.setAirSpeed(velocity); }
    
    void setAirDirection(std::array<Real, 3> direction) { m_tarp.setAirDirection(direction); }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override;

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) override;

  private:
    TarpModel m_tarp;
    ArrayType<Real, 2> m_quad_coords;
};

// Accounts for radiation from environment (ground, sky), but not solar
// radiation
class SkyRadiationBC : public NeumannBC
{
  public:
    SkyRadiationBC(SurfDiscPtr surf, Real emittance, std::array<Real, 3> vertical_vector) :
      NeumannBC(surf, true),
      m_model(emittance, vertical_vector)
    {}

    // set the Horizontal Infrared Radiation intensity (W/m^2)
    void setIRHorizontalRadiation(Real flux) { m_model.setIRHorizontalRadiation(flux); }

    // set the air temperature (which the model uses as the temperature of the ground)
    void setAirTemperature(Real t_air) { m_model.setAirTemperature(t_air); }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        //TODO: move this to base class
        std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
        auto unit_normal = normal / std::sqrt(dot(normal, normal));

        Real flux = m_model.computeFlux(sol_vals[i], unit_normal);
        for (int d=0; d < 3; ++d)
          flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
      }
    }

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        //TODO: move this to base class
        std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
        auto unit_normal = normal / std::sqrt(dot(normal, normal));

        Real flux_dot = 0;
        m_model.computeFluxDeriv(sol_vals[i], unit_normal, flux_dot);
        for (int d=0; d < 3; ++d)
          flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux_dot;
      }
    }

  private:
    SkyRadiationModel m_model;
};


class SolarRadiationBC : public NeumannBC
{
  public:
    SolarRadiationBC(SurfDiscPtr surf, Real absorbtivity) :
      NeumannBC(surf, false),
      m_model(absorbtivity)
    {}

    // sets the solar flux including both Direct Normal Radiation
    // and Diffuse Horzontal Radiation
    void setDirectNormalRadiation(Real flux) { m_model.setDirectNormalRadiation(flux); }

    void setDiffuseRadiation(Real flux) { m_model.setDiffuseRadiation(flux); }

    void setSolarDirection(const DirectionCosines& cosines) { m_model.setSolarDirection(cosines); }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        //TODO: move this to base class
        std::array<Real, 3> normal{m_surf->normals[face][i][0], m_surf->normals[face][i][1], m_surf->normals[face][i][2]};
        auto unit_normal = normal / std::sqrt(dot(normal, normal));

        Real flux = m_model.computeFlux(unit_normal);
        for (int d=0; d < 3; ++d)
          flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
      }
    }

  private:
    SolarRadiationModel m_model;
};

}  // namespace

#endif