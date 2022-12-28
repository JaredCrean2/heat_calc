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
      AirWindSkyNeumannBC(surf, true),
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
      AirWindSkyNeumannBC(surf, false),
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
      AirWindSkyNeumannBC(surf, true),
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

namespace
{
  inline bool anyBcsNonlinear(const std::vector<std::shared_ptr<AirWindSkyNeumannBC>>& bcs)
  {
    bool is_nonlinear = false;
    for (auto& bc : bcs)
      is_nonlinear = is_nonlinear || bc->isNonlinear();

    return is_nonlinear;
  }
}

class CombinedAirWindSkyNeumannBC : public AirWindSkyNeumannBC
{
  public:
    CombinedAirWindSkyNeumannBC(std::vector<std::shared_ptr<AirWindSkyNeumannBC>> bcs) :
      AirWindSkyNeumannBC(bcs[0]->getSurfDisc(), anyBcsNonlinear(bcs)),
      m_bcs(bcs)
    {
      auto surf0 = bcs[0]->getSurfDisc();
      for (size_t i=1; i < bcs.size(); ++i)
        if (surf0 != bcs[i]->getSurfDisc())
          throw std::runtime_error("BCs have different surfaces");
    }

    void setAirTemperature(Real temp) override
    {
      for (auto& bc : m_bcs)
        bc->setAirTemperature(temp);
    }

    Real getAirTemperature() const override { return m_bcs[0]->getAirTemperature(); }
    
    void setAirSpeed(Real velocity) override
    {
      for (auto& bc : m_bcs)
        bc->setAirSpeed(velocity);
    }
    
    void setAirDirection(std::array<Real, 3> direction) override
    {
      for (auto& bc : m_bcs)
        bc->setAirDirection(direction);
    }

    void setIRHorizontalRadiation(Real flux) override
    {
      for (auto& bc : m_bcs)
        bc->setIRHorizontalRadiation(flux);
    }

    void setDirectNormalRadiation(Real flux) override
    {
      for (auto& bc : m_bcs)
        bc->setDirectNormalRadiation(flux);
    }

    void setDiffuseRadiation(Real flux) override
    {
      for (auto& bc : m_bcs)
        bc->setDiffuseRadiation(flux);
    }    

    void setSolarDirection(const DirectionCosines& cosines) override
    {
      for (auto& bc : m_bcs)
        bc->setSolarDirection(cosines);
    }

    void getValue(const Index face, const Real t, const Real* sol_vals, Real* flux_vals) override
    {
      int npts = getSurfDisc()->getNumQuadPtsPerFace()*3;
      std::vector<Real> flux_vals_tmp(npts, 0);

      for (int i=0; i < npts; ++i)
        flux_vals[i] = 0;

      //int i=0;
      for (auto& bc : m_bcs)
      {
        //std::cout << "doing inner bc " << i << std::endl;
        bc->getValue(face, t, sol_vals, flux_vals_tmp.data());
        updateAndZero(flux_vals, flux_vals_tmp.data(), npts);
        //++i;
      }
    }

    void getValueDeriv(const Index face, const Real t, const Real* sol_vals, Real* flux_vals_deriv) override
    {
      int npts = getSurfDisc()->getNumQuadPtsPerFace()*3;
      std::vector<Real> flux_vals_deriv_tmp(npts, 0);

      for (int i=0; i < npts; ++i)
        flux_vals_deriv[i] = 0;

      for (auto& bc : m_bcs)
      {
        bc->getValueDeriv(face, t, sol_vals, flux_vals_deriv_tmp.data());
        updateAndZero(flux_vals_deriv, flux_vals_deriv_tmp.data(), npts);
      }
    }

    // compute derivative of flux_vals wrt air temperature
    virtual void getValuedTair(const Index face, const Real t, const Real* sol_vals, Real* flux_vals, Real* flux_vals_deriv) override
    {
      int npts = getSurfDisc()->getNumQuadPtsPerFace()*3;
      std::vector<Real> flux_vals_tmp(npts, 0), flux_vals_deriv_tmp(npts, 0);

      for (int i=0; i < npts; ++i)
      {
        flux_vals[i] = 0;
        flux_vals_deriv[i] = 0;
      }

      for (auto& bc : m_bcs)
      {
        bc->getValuedTair(face, t, sol_vals, flux_vals_tmp.data(), flux_vals_deriv_tmp.data());

        updateAndZero(flux_vals, flux_vals_tmp.data(), npts);
        updateAndZero(flux_vals_deriv, flux_vals_deriv_tmp.data(), npts);
      }


    }

    virtual void getValue_rev(const Index face, const Real t, const Real* sol_vals, Real* sol_vals_bar, const Real* flux_vals_bar) override
    {
      int npts = getSurfDisc()->getNumQuadPtsPerFace()*3;
      std::vector<Real> sol_vals_bar_tmp(npts, 0);

      for (int i=0; i < npts; ++i)
        sol_vals_bar[i] = 0;

      for (auto& bc : m_bcs)
      {
        bc->getValue_rev(face, t, sol_vals, sol_vals_bar_tmp.data(), flux_vals_bar);
        updateAndZero(sol_vals_bar, sol_vals_bar_tmp.data(), npts);
      }
    }

  private:
    void updateAndZero(Real* arr, Real* arr_tmp, int npts)
    {
      for (int i=0; i < npts; ++i)
      {
        if (std::isnan(arr_tmp[i]))
          throw std::runtime_error("found nan");

        arr[i] += arr_tmp[i];
        arr_tmp[i] = 0;
      }
    }

    std::vector<std::shared_ptr<AirWindSkyNeumannBC>> m_bcs;
};





}  // namespace

#endif