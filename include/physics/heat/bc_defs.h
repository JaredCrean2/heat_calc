#ifndef HEAT_BC_DEFS_H
#define HEAT_BC_DEFS_H

#include "discretization/NeumannBC.h"
#include "discretization/surface_discretization.h"
#include "error_handling.h"
#include "physics/heat/solar_position.h"
#include "tarp.h"
#include "sky_radiation.h"
#include "solar_radiation.h"
#include "floor_radiation_model.h"

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
    AirWindSkyNeumannBC(SurfDiscPtr surf, bool is_nonlinear, const std::string& name) :
      NeumannBC(surf, is_nonlinear),
      m_name(name)
    {}

    std::string getName() const { return m_name; }

    virtual void setAirTemperature(Real temp) {}

    virtual Real getAirTemperature() const = 0;
    
    virtual void setAirSpeed(Real velocity) {}
    
    virtual void setAirDirection(std::array<Real, 3> direction) {}

    virtual void setIRHorizontalRadiation(Real flux) {}

    virtual void setDirectNormalRadiation(Real flux) {}

    virtual void setDiffuseRadiation(Real flux) {}

    virtual void setSolarDirection(const DirectionCosines& cosines) {}

    // compute derivative of flux_vals wrt air temperature
    virtual void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) = 0;

    virtual void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) = 0;

  private:
    std::string m_name;
};

class NewtonCoolingFromAir : public AirWindSkyNeumannBC
{
  public:
    NewtonCoolingFromAir(SurfDiscPtr surf, Real heat_transfer_coeff) :
      AirWindSkyNeumannBC(surf, true, "newton_cooling"),
      m_heat_transfer_coeff(heat_transfer_coeff)
    {}

    virtual void setAirTemperature(Real temp) override { m_temp = temp;}

    virtual Real getAirTemperature() const override { return m_temp; }  

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override;

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals_deriv) override;

    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override;

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override;

  private:
    Real m_heat_transfer_coeff; // Units W/(m^2 K)
    Real m_temp = std::numeric_limits<Real>::min();
};


// Thermal Analysis Research Program (TARP) BC
class TarpBC : public AirWindSkyNeumannBC
{
  public:
    TarpBC(SurfDiscPtr surf, Real surface_area, Real perimeter, int roughness_index, const std::array<Real, 3>& vertical_vector,
              const std::array<Real, 3>& point_at_zero_altitude,
              int met_terrain_index, Real meterological_altitude, int local_terrain_index ) : 
      AirWindSkyNeumannBC(surf, true, "cond"),
      m_tarp(surface_area, perimeter, roughness_index, vertical_vector, 
             point_at_zero_altitude, met_terrain_index, meterological_altitude, 
             local_terrain_index),
      m_quad_coords(boost::extents[m_surf->getNumQuadPtsPerFace()][3])
    {}

    // this constructor is for a Tarp BC on the interior of the structure (so meterological data is unneeded)
    TarpBC(SurfDiscPtr surf, Real surface_area, Real perimeter, int roughness_index, const std::array<Real, 3>& vertical_vector) : 
      AirWindSkyNeumannBC(surf, true, "cond"),
      m_tarp(surface_area, perimeter, roughness_index, vertical_vector, 
             {0, 0, -10000}, 0, 1, 0),
      m_quad_coords(boost::extents[m_surf->getNumQuadPtsPerFace()][3])
    {}

    void setAirTemperature(Real temp) override { m_tarp.setAirTemperature(temp); }

    Real getAirTemperature() const override { return m_tarp.getAirTemperature(); }
    
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
      AirWindSkyNeumannBC(surf, true, "sky_rad"),
      m_model(emittance, vertical_vector)
    {}

    // set the Horizontal Infrared Radiation intensity (W/m^2)
    void setIRHorizontalRadiation(Real flux) override { m_model.setIRHorizontalRadiation(flux); }

    // set the air temperature (which the model uses as the temperature of the ground)
    void setAirTemperature(Real t_air) override { m_model.setAirTemperature(t_air); }

    Real getAirTemperature() const override { return m_model.getAirTemperature(); }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override;

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) override;

    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override;

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override;

  private:
    SkyRadiationModel m_model;
};


class SolarRadiationBC : public AirWindSkyNeumannBC
{
  public:
    SolarRadiationBC(SurfDiscPtr surf, Real absorbtivity) :
      AirWindSkyNeumannBC(surf, false, "solar_rad"),
      m_model(absorbtivity)
    {}

    void setAirTemperature(Real t_air) override { m_t_air = t_air; }

    Real getAirTemperature() const override { return m_t_air; }

    void setDirectNormalRadiation(Real flux) override { m_model.setDirectNormalRadiation(flux); }

    void setDiffuseRadiation(Real flux) override { m_model.setDiffuseRadiation(flux); }

    void setSolarDirection(const DirectionCosines& cosines) override { m_model.setSolarDirection(cosines); }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override;

    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override;

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override;

  private:
    SolarRadiationModel m_model;
    Real m_t_air;
};



class SimpleConvectionBC : public AirWindSkyNeumannBC
{
  public:
    SimpleConvectionBC(SurfDiscPtr surf, Real h) : 
      AirWindSkyNeumannBC(surf, true, "simple_convection"),
      m_quad_coords(boost::extents[m_surf->getNumQuadPtsPerFace()][3]),
      m_heat_transfer_coeff(h)
    {}

    void setAirTemperature(Real temp) override { m_air_temp = temp; }

    Real getAirTemperature() const override { return m_air_temp; }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override;

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) override;

    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override;

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override;

  private:
    ArrayType<Real, 2> m_quad_coords;
    Real m_heat_transfer_coeff;
    Real m_air_temp;
};


class FloorRadiationBC : public AirWindSkyNeumannBC
{
  public:
    FloorRadiationBC(SurfDiscPtr surf, Real window_area, std::array<Real, 3> window_normal,
                     Real shgc, Real floor_area, Real floor_absorbtivity, const std::string& name) :
      AirWindSkyNeumannBC(surf, false, name),
      m_model(window_area, window_normal, shgc, floor_area, floor_absorbtivity)
    {}

    void setAirTemperature(Real t_air) override { m_t_air = t_air; }

    Real getAirTemperature() const override { return m_t_air; }

    void setDirectNormalRadiation(Real flux) override { m_model.setDirectNormalRadiation(flux); }

    void setDiffuseRadiation(Real flux) override { m_model.setDiffuseRadiation(flux); }

    void setSolarDirection(const DirectionCosines& cosines) override { m_model.setSolarDirection(cosines); }

    void getValue(const Index face, const Real t, const Real* sol_vals,  Real* flux_vals) override;

    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override;

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override;

  private:
    FloorRadiationModel m_model;
    Real m_t_air;
};

class AirWindSkyZeroBC : public AirWindSkyNeumannBC
{
  public:
    AirWindSkyZeroBC(SurfDiscPtr surf) :
      AirWindSkyNeumannBC(surf, false, "zero_bc")
    {}

    void setAirTemperature(Real temp) override{ m_air_temp = temp;}

    Real getAirTemperature() const override { return m_air_temp; }
    
    void getValue(const Index face, const Real t, const Real* sol_vals, Real* flux_vals) override
    {
      for (int i=0; i < 3*m_surf->getNumQuadPtsPerFace(); ++i)
        flux_vals[i] = 0;
    }


    // compute derivative of flux_vals wrt air temperature
    void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override
    {
      for (int i=0; i < 3*m_surf->getNumQuadPtsPerFace(); ++i)
      {
        flux_vals[i] = 0;
        flux_vals_deriv[i] = 0;
      }
    }

    void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override
    {
      for (int i=0; i < m_surf->getNumQuadPtsPerFace(); ++i)
        sol_vals_bar[i] = 0;
    }

  private:
    Real m_air_temp;
};


class CombinedAirWindSkyNeumannBC : public AirWindSkyNeumannBC
{
  public:
    CombinedAirWindSkyNeumannBC(std::vector<std::shared_ptr<AirWindSkyNeumannBC>> bcs);

    int getNumBCs() const { return m_bcs.size(); }

    std::shared_ptr<AirWindSkyNeumannBC> getBC(int i) { return m_bcs.at(i); }

    void setAirTemperature(Real temp) override;

    Real getAirTemperature() const override;
    
    void setAirSpeed(Real velocity) override;
    
    void setAirDirection(std::array<Real, 3> direction) override;

    void setIRHorizontalRadiation(Real flux) override;

    void setDirectNormalRadiation(Real flux) override;

    void setDiffuseRadiation(Real flux) override;  

    void setSolarDirection(const DirectionCosines& cosines) override;

    void getValue(const Index face, const Real t, const Real* sol_vals, Real* flux_vals) override;

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) override;

    // compute derivative of flux_vals wrt air temperature
    virtual void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override;

    virtual void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override;

  private:
    void updateAndZero(Real* arr, Real* arr_tmp, int npts);

    std::vector<std::shared_ptr<AirWindSkyNeumannBC>> m_bcs;
};



}  // namespace

#endif