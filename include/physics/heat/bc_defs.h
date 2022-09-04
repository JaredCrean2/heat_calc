#ifndef HEAT_BC_DEFS_H
#define HEAT_BC_DEFS_H

#include "discretization/NeumannBC.h"
#include "discretization/surface_discretization.h"
#include "error_handling.h"
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

// base class for BCs that depend on air, wind, and solar parameters
class AirWindSkyNeumannBC : public NeumannBC
{
  public:
    AirWindSkyNeumannBC(SurfDiscPtr surf, bool is_nonlinear) :
      NeumannBC(surf, is_nonlinear)
    {}

    virtual void setAirTemperature(Real temp) {}
    
    virtual void setAirSpeed(Real velocity) {}
    
    virtual void setAirDirection(std::array<Real, 3> direction) {}

    virtual void setIRHorizontalRadiation(Real flux) {}

    virtual void setDirectNormalRadiation(Real flux) {}

    virtual void setDiffuseRadiation(Real flux) {}

    virtual void setSolarDirection(const DirectionCosines& cosines) {}

    // compute derivative of flux_vals wrt air temperature
    virtual void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) = 0;

    virtual void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) = 0;
};


// Thermal Analysis Research Program (TARP) BC
class TarpBC : public AirWindSkyNeumannBC
{
  public:
    TarpBC(SurfDiscPtr surf, Real surface_area, Real perimeter, int roughness_index, const std::array<Real, 3>& vertical_vector,
              const std::array<Real, 3>& point_at_zero_altitude,
              int met_terrain_index, Real meterological_altitude, int local_terrain_index ) : 
      AirWindSkyNeumannBC(surf, true),
      m_tarp(surface_area, perimeter, roughness_index, vertical_vector, 
             point_at_zero_altitude, met_terrain_index, meterological_altitude, 
             local_terrain_index),
      m_quad_coords(boost::extents[m_surf->getNumQuadPtsPerFace()][3])
    {}

    // this constructor is for a Tarp BC on the interior of the structure (so meterological data is unneeded)
    TarpBC(SurfDiscPtr surf, Real surface_area, Real perimeter, int roughness_index, const std::array<Real, 3>& vertical_vector) : 
      AirWindSkyNeumannBC(surf, true),
      m_tarp(surface_area, perimeter, roughness_index, vertical_vector, 
             {0, 0, 0}, 0, 1, 0),
      m_quad_coords(boost::extents[m_surf->getNumQuadPtsPerFace()][3])
    {}

    void setAirTemperature(Real temp) override { m_tarp.setAirTemperature(temp); }
    
    void setAirSpeed(Real velocity) override { m_tarp.setAirSpeed(velocity); }
    
    void setAirDirection(std::array<Real, 3> direction) override { m_tarp.setAirDirection(direction); }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override;

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) override;

    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override;

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override;


  private:
    TarpModel m_tarp;
    ArrayType<Real, 2> m_quad_coords;
};

// Accounts for radiation from environment (ground, sky), but not solar
// radiation
class SkyRadiationBC : public AirWindSkyNeumannBC
{
  public:
    SkyRadiationBC(SurfDiscPtr surf, Real emittance, std::array<Real, 3> vertical_vector) :
      AirWindSkyNeumannBC(surf, true),
      m_model(emittance, vertical_vector)
    {}

    // set the Horizontal Infrared Radiation intensity (W/m^2)
    void setIRHorizontalRadiation(Real flux) override { m_model.setIRHorizontalRadiation(flux); }

    // set the air temperature (which the model uses as the temperature of the ground)
    void setAirTemperature(Real t_air) override { m_model.setAirTemperature(t_air); }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        auto unit_normal = getUnitNormal(face, i);
        Real flux = m_model.computeFlux(sol_vals[i], unit_normal);
        for (int d=0; d < 3; ++d)
          flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
      }
    }

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        auto unit_normal = getUnitNormal(face, i);
        Real flux_dot = 0;
        m_model.computeFluxdTwall(sol_vals[i], unit_normal, flux_dot);
        for (int d=0; d < 3; ++d)
          flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux_dot;
      }
    }

    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        auto unit_normal = getUnitNormal(face, i);
        Real flux_dot = 0;
        Real flux = m_model.computeFluxdTair(sol_vals[i], unit_normal, flux_dot);
        for (int d=0; d < 3; ++d)
        {
          flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
          flux_vals_deriv[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux_dot;
        }
      }
    }

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        auto unit_normal = getUnitNormal(face, i);
        Real flux_dot;
        m_model.computeFluxdTwall(sol_vals[i], unit_normal, flux_dot);
        Real flux_bar = 0;
        for (int d=0; d < 3; ++d)
        {
          //flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
          flux_bar += unit_normal[d] * flux_vals_bar[d * m_surf->getNumQuadPtsPerFace() + i];
        }

        sol_vals_bar[i] = flux_dot * flux_bar;
      }
    }

  private:
    SkyRadiationModel m_model;
};


class SolarRadiationBC : public AirWindSkyNeumannBC
{
  public:
    SolarRadiationBC(SurfDiscPtr surf, Real absorbtivity) :
      AirWindSkyNeumannBC(surf, false),
      m_model(absorbtivity)
    {}

    void setDirectNormalRadiation(Real flux) override { m_model.setDirectNormalRadiation(flux); }

    void setDiffuseRadiation(Real flux) override { m_model.setDiffuseRadiation(flux); }

    void setSolarDirection(const DirectionCosines& cosines) override { m_model.setSolarDirection(cosines); }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
      {
        auto unit_normal = getUnitNormal(face, i);
        Real flux = m_model.computeFlux(unit_normal);
        for (int d=0; d < 3; ++d)
          flux_vals[d * m_surf->getNumQuadPtsPerFace() + i] = unit_normal[d] * flux;
      }
    }

    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override
    {
      for (int i=0; i < 3*m_surf->getNumQuadPtsPerFace(); ++i)
        flux_vals_deriv[i] = 0;
    }

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
        sol_vals_bar[i] = 0;
    }

  private:
    SolarRadiationModel m_model;
};




}  // namespace

#endif