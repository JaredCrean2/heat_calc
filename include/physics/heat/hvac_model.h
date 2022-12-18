#ifndef PHYSICS_HEAT_HVAC_MODEL
#define PHYSICS_HEAT_HVAC_MODEL

#include "ProjectDefs.h"

#include "error_handling.h"
#include "utils/spline.h"

namespace Heat {

class HVACModel
{
  public:
    virtual ~HVACModel() {}
    
    // computes the HVAC flux.  Negative values denote cooling (decreasing interior air 
    // temperature), positive values denote heating (increasing interior air temperature)
    virtual Real enforceTemperatureLimit(Real interior_temp, Real load_flux) = 0;

    // compute derivative of HVAC flux wrt interior_temp
    virtual Real enforceTemperatureLimitDotTair(Real interior_temp, Real load_flux, Real load_flux_dot) = 0;

    // reverse mode of enforceTemperatureLimit
    virtual void enforceTemperatureLimit_rev(Real interior_temp, Real& interior_temp_bar, Real load_flux, Real& load_flux_bar, Real hvac_flux_bar) = 0;
};


// simple on/off switch for HVAC
class HVACModelSwitch : public HVACModel
{
  public:
    HVACModelSwitch(Real min_temp, Real max_temp, Real rho_cp, Real air_volume, Real hvac_restore_time) :
      m_min_temp(min_temp),
      m_max_temp(max_temp),
      m_rho_cp(rho_cp),
      m_air_volume(air_volume),
      m_hvac_restore_time(hvac_restore_time)
    {
      assertAlways(max_temp >= min_temp, "max temperature must be >= min temperature");
    }

    Real enforceTemperatureLimit(Real interior_temp, Real load_flux) override
    {
      if (interior_temp > m_max_temp)
        return enforceTemperatureLimit(m_max_temp, interior_temp, load_flux);
      else if (interior_temp < m_min_temp)
        return enforceTemperatureLimit(m_min_temp, interior_temp, load_flux);
      else
        return 0;
    }

    Real enforceTemperatureLimitDotTair(Real interior_temp, Real load_flux, Real load_flux_dot) override
    {
      if (interior_temp > m_max_temp)
        return enforceTemperatureLimitDotTair(m_max_temp, interior_temp, load_flux, load_flux_dot);
      else if (interior_temp < m_min_temp)
        return enforceTemperatureLimitDotTair(m_min_temp, interior_temp, load_flux, load_flux_dot);
      else
        return 0;
    }

    void enforceTemperatureLimit_rev(Real interior_temp, Real& interior_temp_bar, 
                                     Real load_flux, Real& load_flux_bar, Real hvac_flux_bar) override
    {
      if (interior_temp > m_max_temp)
        enforceTemperatureLimit_rev(m_max_temp, interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);
      else if (interior_temp < m_min_temp)
        enforceTemperatureLimit_rev(m_min_temp, interior_temp, interior_temp_bar, load_flux, hvac_flux_bar);
    }
  
    void enforceTemperatureLimitDotTair_rev(Real interior_temp,
                                            Real load_flux,
                                            Real load_flux_dot, Real& load_flux_dot_bar,
                                            Real hvac_flux_dot_bar)
    {
      if (interior_temp > m_max_temp)
        return enforceTemperatureLimitDotTair_rev(m_max_temp, interior_temp, load_flux, load_flux_dot, load_flux_dot_bar, hvac_flux_dot_bar);
      else if (interior_temp < m_min_temp)
        return enforceTemperatureLimitDotTair_rev(m_min_temp, interior_temp, load_flux, load_flux_dot, load_flux_dot_bar, hvac_flux_dot_bar);
    }

  private:

    Real enforceTemperatureLimit(Real temp_limit, Real interior_temp, Real load_flux)
    {
      return m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
    }  

    Real enforceTemperatureLimitDotTair(Real temp_limit, Real interior_temp,
                                        Real load_flux, Real load_flux_dot)
    {
      //m_hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
      Real hvac_flux_dot = m_rho_cp * m_air_volume / m_hvac_restore_time - load_flux_dot;

      return hvac_flux_dot;
    }

    void enforceTemperatureLimit_rev(Real temp_limit, Real interior_temp,
                                     Real& interior_temp_bar, Real load_flux, 
                                     Real& load_flux_bar, Real hvac_flux_bar)
    {
      //Real hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
      interior_temp_bar += -m_rho_cp * m_air_volume * hvac_flux_bar / m_hvac_restore_time;
      load_flux_bar     -= hvac_flux_bar;
    }

    void enforceTemperatureLimitDotTair_rev(Real temp_limit, Real interior_temp,
                                            Real load_flux,
                                            Real load_flux_dot, Real& load_flux_dot_bar,
                                            Real hvac_flux_dot_bar)
    {
      //m_hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
      //Real hvac_flux_dot = m_rho_cp * m_air_volume / m_hvac_restore_time - load_flux_dot;

      load_flux_dot_bar -= hvac_flux_dot_bar;
    }

    Real m_min_temp;
    Real m_max_temp;
    Real m_rho_cp;
    Real m_air_volume;
    Real m_hvac_restore_time;  // (approximate) time for the HVAC system to restore the
                               // interior air temperature to within bands
};


// use a spline to compute HVAC flux when within the temperature range.
// This solves a problem with Newton's method convergence, but leads to
// the interior air temperature remaining constant
class HVACModelSpline : public HVACModel
{
  public:
    HVACModelSpline(Real min_temp, Real max_temp, Real rho_cp, Real air_volume, Real hvac_restore_time, Real expansion_factor=0.1) :
      m_switch_model(min_temp, max_temp, rho_cp, air_volume, hvac_restore_time),
      m_min_temp(min_temp),
      m_max_temp(max_temp),
      m_expansion_factor(expansion_factor)
    {}

    Real enforceTemperatureLimit(Real interior_temp, Real load_flux) override
    {
      std::array<Real, 6> params = getSplineParams(interior_temp, load_flux);

      if (interior_temp > params[1] || interior_temp < params[0])
        return m_switch_model.enforceTemperatureLimit(interior_temp, load_flux);
      else
      {
        SingleCubicSpline spline;
        spline.setupSpline(params[0], params[1], params[2], params[3], params[4], params[5]);
        return spline.eval(interior_temp);
      }
    }

    Real enforceTemperatureLimitDotTair(Real interior_temp, Real load_flux, Real load_flux_dot) override
    {
      std::array<Real, 6> params = getSplineParams(interior_temp, load_flux);

      if (interior_temp > params[1] || interior_temp < params[0])
        return m_switch_model.enforceTemperatureLimitDotTair(interior_temp, load_flux, load_flux_dot);
      else
      {
        SingleCubicSpline spline;
        spline.setupSpline(params[0], params[1], params[2], params[3], params[4], params[5]);
        Real hvac_flux_dot;
        spline.evalDot(interior_temp, 1, hvac_flux_dot);
        return hvac_flux_dot;
      }
    }

    void enforceTemperatureLimit_rev(Real interior_temp, Real& interior_temp_bar,
                                     Real load_flux, Real& load_flux_bar,
                                     Real hvac_flux_bar) override
    {
      std::array<Real, 6> params = getSplineParams(interior_temp, load_flux);

      if (interior_temp > params[1])
        m_switch_model.enforceTemperatureLimit_rev(interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);
      else
      {
        SingleCubicSpline spline;
        spline.setupSpline(params[0], params[1], params[2], params[3], params[4], params[5]);

        std::array<Real, 4> coeffs_bar;
        spline.evalRev(interior_temp, hvac_flux_bar, interior_temp_bar, coeffs_bar);

        std::array<Real, 6> params_bar = {0};
        spline.setupSplineRev(params[0], params[1], coeffs_bar, params_bar[2], params_bar[3], params_bar[4], params_bar[5]);

        getSplineParams_rev(interior_temp, interior_temp_bar, load_flux, load_flux_bar, params_bar);
      }
    }

  private:

    // setup spline to smooth the activation of the HVAC system, which causes
    // convergence problems for Newtons method
    std::array<Real, 6> getSplineParams(Real interior_temp, Real load_flux)
    {
      Real delta_t = (m_max_temp - m_min_temp);
      Real t_upper = m_max_temp + m_expansion_factor * delta_t;
      Real t_lower = m_min_temp - m_expansion_factor * delta_t;
      Real load_flux_dot = 0; // We are computing the spline so the HVAC flux is a continuous function of the air temperature,
                              // so don't consider the derivative wrt load flux

      Real flux_upper     = m_switch_model.enforceTemperatureLimit(t_upper, load_flux);
      Real flux_upper_dot = m_switch_model.enforceTemperatureLimitDotTair(t_upper, load_flux, load_flux_dot);

      Real flux_lower     = m_switch_model.enforceTemperatureLimit(t_lower, load_flux);
      Real flux_lower_dot = m_switch_model.enforceTemperatureLimitDotTair(t_lower, load_flux, load_flux_dot);

      return {t_lower, t_upper, flux_lower, flux_lower_dot, flux_upper, flux_upper_dot};
    }

    void getSplineParams_rev(Real interior_temp, Real& interior_temp_bar,
                             Real load_flux, Real& load_flux_bar,
                             const std::array<Real, 6> params_bar)
    {
      Real delta_t = (m_max_temp - m_min_temp);
      Real t_upper = m_max_temp + m_expansion_factor * delta_t;
      Real t_lower = m_min_temp - m_expansion_factor * delta_t;
      Real load_flux_dot = 0; // We are computing the spline so the HVAC flux is a continuous function of the air temperature,
                              // so don't consider the derivative wrt load flux

      //Real flux_upper     = enforceTemperatureLimitStraightLine(m_max_temp, t_upper, load_flux);
      //Real flux_upper_dot = enforceTemperatureLimitStraightLineDotTair(m_max_temp, t_upper, load_flux, load_flux_dot);

      //Real flux_lower     = enforceTemperatureLimitStraightLine(m_min_temp, t_lower, load_flux);
      //Real flux_lower_dot = enforceTemperatureLimitStraightLineDotTair(m_min_temp, t_lower, load_flux, load_flux_dot);

      Real load_flux_dot_bar = 0;
      m_switch_model.enforceTemperatureLimitDotTair_rev(t_lower, load_flux, load_flux_dot, load_flux_dot_bar, params_bar[3]);
      m_switch_model.enforceTemperatureLimit_rev(t_lower, interior_temp_bar, load_flux, load_flux_bar, params_bar[2]);

      m_switch_model.enforceTemperatureLimitDotTair_rev(t_upper, load_flux, load_flux_dot, load_flux_dot_bar, params_bar[5]);
      m_switch_model.enforceTemperatureLimit_rev(t_upper, interior_temp_bar, load_flux, load_flux_bar, params_bar[4]);
    }

    HVACModelSwitch m_switch_model;
    Real m_min_temp;
    Real m_max_temp;
    Real m_expansion_factor;  // used to determine the points used to construct the spline
};

}  // namespace

#endif