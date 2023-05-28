#include "physics/heat/steady_state_temp_calculator.h"

namespace Heat {

SteadyStateTempCaluclator::SteadyStateTempCaluclator(std::shared_ptr<EnvironmentInterface> env_interface,
                          std::shared_ptr<SolarPositionCalculator> solar_position,
                          const TarpModel& tarp_model,
                          const SkyRadiationModel& sky_model,
                          const SolarRadiationModel& solar_model) :
  m_env_interface(env_interface),
  m_solar_position(solar_position),
  m_tarp_model(tarp_model),
  m_sky_model(sky_model),
  m_solar_model(solar_model)
{}

double SteadyStateTempCaluclator::calculate(double t_start, double t_end, double delta_t)
{
  int nsteps_whole = std::floor((t_end - t_start)/delta_t);
  int rem = (nsteps_whole - 3) % 2;
  nsteps_whole -= rem;
  delta_t = (t_end - t_start)/(nsteps_whole - 1);
  std::vector<double> temp_vals, t_vals;

  double ground_temp = 0, time = 0;
  for (int i=0; i < nsteps_whole; ++i)
  {
    time = delta_t * i;
    double air_temp = setEnvironmentData(time);

    if (i == 0)
      ground_temp = air_temp;

    ground_temp = computeSteadyStateTemp(ground_temp, air_temp, 1e-11);

    temp_vals.push_back(ground_temp);
    t_vals.push_back(time);
  }

  return integrate(t_vals, temp_vals);
}
                        

double SteadyStateTempCaluclator::computeSteadyStateTemp(double ground_temp_guess, double air_temp, double tol)
{
  double ground_temp = ground_temp_guess;
  double rhs = computeRhs(ground_temp, air_temp);

  while (std::abs(rhs) > tol)
  {
    double deriv = computeDeriv(ground_temp, air_temp);
    double delta_temp = rhs/deriv;
    ground_temp = ground_temp - delta_temp;

    rhs = computeRhs(ground_temp, air_temp);
  }

  return ground_temp;
}

double SteadyStateTempCaluclator::computeRhs(double ground_temp, double air_temp)
{
  std::array<Real, 3> unit_normal = {0, 0, 1};
  double h = m_tarp_model.computeHeatTransferCoeff(ground_temp, {0, 0, 0}, unit_normal);
  double convective_flux = h * (air_temp - ground_temp); //TODO: check sign
  double sky_flux = m_sky_model.computeFlux(ground_temp, unit_normal);
  double solar_flux = m_solar_model.computeFlux(unit_normal);

  return convective_flux + sky_flux + solar_flux;
}

double SteadyStateTempCaluclator::computeDeriv(double ground_temp, double air_temp)
{
  std::array<Real, 3> unit_normal = {0, 0, 1};
  Real h_dot, sky_flux_dot;
  double h = m_tarp_model.computeHeatTransferCoeffdTwall(ground_temp, {0, 0, 0}, unit_normal, h_dot);

  double convective_flux_dot = -h + h_dot * (air_temp - ground_temp);
  m_sky_model.computeFluxdTwall(ground_temp, unit_normal, sky_flux_dot);
  
  return convective_flux_dot + sky_flux_dot;
}

double SteadyStateTempCaluclator::setEnvironmentData(double time)
{
  EnvironmentData edata = m_env_interface->getEnvironmentData(time);

  m_tarp_model.setAirSpeed(edata.air_speed);
  m_tarp_model.setAirDirection(edata.air_direction);

  m_sky_model.setAirTemperature(edata.air_temp);
  m_sky_model.setIRHorizontalRadiation(edata.ir_horizontal_radiation);

  m_solar_model.setDiffuseRadiation(edata.diffuse_radiation);
  m_solar_model.setDirectNormalRadiation(edata.direct_normal_radiation);
  m_solar_model.setSolarDirection(m_solar_position->computePosition(time));

  return edata.air_temp;
}

double SteadyStateTempCaluclator::integrate(const std::vector<double>& t_vals, const std::vector<double>& temp_vals)
{
  assertAlways((t_vals.size() - 3) % 2 == 0, "number of values must be divisible into groups of 3");

  double val = 0;
  for (size_t i=0; i < t_vals.size()-1; i+=2)
  {        
    double delta_t = t_vals[i+2] - t_vals[i];
    double val_i = delta_t * (temp_vals[i] + 4*temp_vals[i+1] + temp_vals[i+2])/6;
    val += val_i;
  }

  return val/(t_vals.back() - t_vals.front());
}

}  // namespace