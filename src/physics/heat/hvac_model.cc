#include "physics/heat/hvac_model.h"

namespace Heat {

//-----------------------------------------------------------------------------
// HVACModelSwitch

HVACModelSwitch::HVACModelSwitch(Real min_temp, Real max_temp, Real rho_cp, Real air_volume, Real hvac_restore_time) :
      m_min_temp(min_temp),
      m_max_temp(max_temp),
      m_rho_cp(rho_cp),
      m_air_volume(air_volume),
      m_hvac_restore_time(hvac_restore_time)
{
  assertAlways(max_temp >= min_temp, "max temperature must be >= min temperature");
}

Real HVACModelSwitch::enforceTemperatureLimit(Real interior_temp, Real load_flux)
{
  if (interior_temp > m_max_temp)
    return enforceTemperatureLimit(m_max_temp, interior_temp, load_flux);
  else if (interior_temp < m_min_temp)
    return enforceTemperatureLimit(m_min_temp, interior_temp, load_flux);
  else
    return 0;
}

Real HVACModelSwitch::enforceTemperatureLimit_dot(Real interior_temp, Real interior_temp_dot, Real load_flux, Real load_flux_dot)
{
  if (interior_temp > m_max_temp)
    return enforceTemperatureLimit_dot(m_max_temp, interior_temp, interior_temp_dot, load_flux, load_flux_dot);
  else if (interior_temp < m_min_temp)
    return enforceTemperatureLimit_dot(m_min_temp, interior_temp, interior_temp_dot, load_flux, load_flux_dot);
  else
    return 0;
}

void HVACModelSwitch::enforceTemperatureLimit_rev(Real interior_temp, Real& interior_temp_bar, 
                                                  Real load_flux, Real& load_flux_bar, Real hvac_flux_bar)
{
  if (interior_temp > m_max_temp)
    enforceTemperatureLimit_rev(m_max_temp, interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);
  else if (interior_temp < m_min_temp)
    enforceTemperatureLimit_rev(m_min_temp, interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);
}

void HVACModelSwitch::enforceTemperatureLimit_dot_rev(Real interior_temp,
                                                      Real interior_temp_dot, Real& interior_temp_dot_bar,
                                                      Real load_flux,
                                                      Real load_flux_dot, Real& load_flux_dot_bar,
                                                      Real hvac_flux_dot_bar)
{
  if (interior_temp > m_max_temp)
    return enforceTemperatureLimit_dot_rev(m_max_temp, interior_temp, interior_temp_dot, interior_temp_dot_bar, 
                                          load_flux, load_flux_dot, load_flux_dot_bar, hvac_flux_dot_bar);
  else if (interior_temp < m_min_temp)
    return enforceTemperatureLimit_dot_rev(m_min_temp, interior_temp, interior_temp_dot, interior_temp_dot_bar,
                                            load_flux, load_flux_dot, load_flux_dot_bar, hvac_flux_dot_bar);
}

Real HVACModelSwitch::enforceTemperatureLimit(Real temp_limit, Real interior_temp, Real load_flux)
{
  return m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
}  

Real HVACModelSwitch::enforceTemperatureLimit_dot(Real temp_limit, Real interior_temp, Real interior_temp_dot,
                                                  Real load_flux, Real load_flux_dot)
{
  //m_hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
  Real hvac_flux_dot = -m_rho_cp * m_air_volume * interior_temp_dot / m_hvac_restore_time - load_flux_dot;

  return hvac_flux_dot;
}

void HVACModelSwitch::enforceTemperatureLimit_rev(Real temp_limit, Real interior_temp,
                                                  Real& interior_temp_bar, Real load_flux, 
                                                  Real& load_flux_bar, Real hvac_flux_bar)
{
  //Real hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
  std::cout << "original load_flux_bar value = " << load_flux_bar << std::endl;
  std::cout << "hvac_flux_bar = " << hvac_flux_bar << std::endl;
  
  interior_temp_bar += -m_rho_cp * m_air_volume * hvac_flux_bar / m_hvac_restore_time;
  load_flux_bar     -= hvac_flux_bar;
  std::cout << "final load_flux_bar value = " << load_flux_bar << std::endl;

}

void HVACModelSwitch::enforceTemperatureLimit_dot_rev(Real temp_limit, Real interior_temp,
                                                      Real interior_temp_dot, Real& interior_temp_dot_bar,
                                                      Real load_flux,
                                                      Real load_flux_dot, Real& load_flux_dot_bar,
                                                      Real hvac_flux_dot_bar)
{
  //m_hvac_flux = m_rho_cp * m_air_volume * (temp_limit - interior_temp)/m_hvac_restore_time - load_flux;
  //Real hvac_flux_dot = -m_rho_cp * m_air_volume * interior_temp_dot / m_hvac_restore_time - load_flux_dot;

  interior_temp_dot_bar -= -m_rho_cp * m_air_volume * hvac_flux_dot_bar / m_hvac_restore_time;
  load_flux_dot_bar     -= hvac_flux_dot_bar;
}


//-----------------------------------------------------------------------------
// HVACModelTempOnly

HVACModelTempOnly::HVACModelTempOnly(Real temp_lower, Real temp_upper, Real rho_cp, Real air_volume,
                                     Real hvac_restore_time, int poly_degree) :
  //m_hvac_restore_time(hvac_restore_time),
  m_temp_lower(temp_lower),
  m_temp_upper(temp_upper),
  m_poly_degree(poly_degree),
  m_A1(0),
  m_A2(0)
{
  assertAlways(poly_degree >= 1, "HVACModelTempOnly polynomial degree must be >= 1");
  assertAlways(std::abs(temp_upper - temp_lower) > 1e-3, "temperature range must be > 1e-3");
  assertAlways(temp_upper >= temp_lower, "upper temperature must be greater than lower temperature");
  Real temp_range_term = std::pow((temp_upper - temp_lower)/2, poly_degree-1);
  m_A1 = rho_cp*air_volume/(hvac_restore_time * temp_range_term);
  m_A2 = (poly_degree % 2 == 0) ? m_A1 : -m_A1;
}

// computes the HVAC flux.  Negative values denote cooling (decreasing interior air 
// temperature), positive values denote heating (increasing interior air temperature)
Real HVACModelTempOnly::enforceTemperatureLimit(Real interior_temp, Real load_flux)
{
  Real t_avg   = (m_temp_upper + m_temp_lower)/2;
  Real delta_t = interior_temp - t_avg;
  if (interior_temp >= t_avg)
    return -m_A1 * std::pow(delta_t, m_poly_degree);
  else  // interior_temp < t_avg
    return m_A2 * std::pow(delta_t, m_poly_degree);
}

// compute derivative of HVAC flux wrt interior_temp
Real HVACModelTempOnly::enforceTemperatureLimit_dot(Real interior_temp, Real interior_temp_dot, Real load_flux, Real load_flux_dot)
{
  Real t_avg       = (m_temp_upper + m_temp_lower)/2;
  Real delta_t     = interior_temp - t_avg;
  Real delta_t_dot = interior_temp_dot;
  if (interior_temp >= t_avg)
    return -m_poly_degree * m_A1 * std::pow(delta_t, m_poly_degree - 1) * delta_t_dot;
  else  // interior_temp < t_avg
    return m_poly_degree * m_A2 * std::pow(delta_t, m_poly_degree - 1) * delta_t_dot;      
}

// reverse mode of enforceTemperatureLimit
void HVACModelTempOnly::enforceTemperatureLimit_rev(Real interior_temp, Real& interior_temp_bar, Real load_flux, 
                                                    Real& load_flux_bar, Real hvac_flux_bar)
{
  Real t_avg   = (m_temp_upper + m_temp_lower)/2;
  Real delta_t = interior_temp - t_avg;

  Real delta_t_bar;
  if (interior_temp >= t_avg)
    delta_t_bar = -m_poly_degree * m_A1 * std::pow(delta_t, m_poly_degree - 1) * hvac_flux_bar;
  else  // interior_temp < t_avg
    delta_t_bar = m_poly_degree * m_A2 * std::pow(delta_t, m_poly_degree - 1) * hvac_flux_bar;

  interior_temp_bar = delta_t_bar;
}

//-----------------------------------------------------------------------------
// HVACModelTempOnlyCoolingOnly

HVACModelTempOnlyCoolingOnly::HVACModelTempOnlyCoolingOnly(Real temp_lower, Real temp_upper, Real rho_cp, Real air_volume,
                                     Real hvac_restore_time, int poly_degree) :
  m_temp_only_model(temp_lower, temp_upper, rho_cp, air_volume, hvac_restore_time, poly_degree)
{
}

Real HVACModelTempOnlyCoolingOnly::enforceTemperatureLimit(Real interior_temp, Real load_flux)
{
  Real flux = m_temp_only_model.enforceTemperatureLimit(interior_temp, load_flux);
  if (flux > 0)
    return 0;
  else
    return flux;
}

// compute derivative of HVAC flux wrt interior_temp
Real HVACModelTempOnlyCoolingOnly::enforceTemperatureLimit_dot(Real interior_temp, Real interior_temp_dot, Real load_flux, Real load_flux_dot)
{
  Real flux = m_temp_only_model.enforceTemperatureLimit(interior_temp, load_flux);

  if (flux > 0)
  {
    return 0;
  } else
  {
    return m_temp_only_model.enforceTemperatureLimit_dot(interior_temp, interior_temp_dot, load_flux, load_flux_dot);
  }
}

// reverse mode of enforceTemperatureLimit
void HVACModelTempOnlyCoolingOnly::enforceTemperatureLimit_rev(Real interior_temp, Real& interior_temp_bar, Real load_flux, 
                                                    Real& load_flux_bar, Real hvac_flux_bar)
{
  Real flux = m_temp_only_model.enforceTemperatureLimit(interior_temp, load_flux);
  if (flux > 0)
  {
    interior_temp_bar = 0;
  } else
  {
    m_temp_only_model.enforceTemperatureLimit_rev(interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);
  }
}

}