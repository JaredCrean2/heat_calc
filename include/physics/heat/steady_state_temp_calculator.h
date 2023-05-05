#ifndef PHYSICS_HEAT_STEADY_STATE_TEMP_CALCULATOR_H
#define PHYSICS_HEAT_STEADY_STATE_TEMP_CALCULATOR_H

#include "physics/heat/sky_radiation.h"
#include "physics/heat/solar_position_calculator.h"
#include "physics/heat/solar_radiation.h"
#include "physics/heat/tarp.h"
#include "physics/heat/environment_interface.h"

namespace Heat {

class SteadyStateTempCaluclator
{
  public:
    SteadyStateTempCaluclator(std::shared_ptr<EnvironmentInterface> env_interface,
                              std::shared_ptr<SolarPositionCalculator> solar_position,
                              const TarpModel& tarp_model,
                              const SkyRadiationModel& sky_model,
                              const SolarRadiationModel& solar_model);

    double calculate(double t_start, double t_end, double delta_t);
                            

  private:
    double computeSteadyStateTemp(double ground_temp_guess, double air_temp, double tol);
    
    double computeRhs(double ground_temp, double air_temp);

    double computeDeriv(double ground_temp, double air_temp);

    double setEnvironmentData(double time);

    double integrate(const std::vector<double>& t_vals, const std::vector<double>& temp_vals);

    std::shared_ptr<EnvironmentInterface> m_env_interface;
    std::shared_ptr<SolarPositionCalculator> m_solar_position;
    TarpModel m_tarp_model;
    SkyRadiationModel m_sky_model;
    SolarRadiationModel m_solar_model;

};

}  // namespace

#endif