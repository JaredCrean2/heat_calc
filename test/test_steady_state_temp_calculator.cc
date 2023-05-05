#include "gtest/gtest.h"
#include "physics/heat/environment_interface.h"
#include "physics/heat/sky_radiation.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/solar_position_calculator.h"
#include "physics/heat/solar_radiation.h"
#include "physics/heat/steady_state_temp_calculator.h"
#include "physics/heat/tarp.h"

TEST(SteadyStateTempCalc, AirOnly)
{
  Heat::EnvironmentData edata;
  edata.air_temp = 300;
  edata.air_speed = 5;
  edata.air_direction = {1, 0, 0};
  edata.ir_horizontal_radiation = 5.6697E-8 * std::pow(edata.air_temp, 4);
  edata.direct_normal_radiation = 0;
  edata.diffuse_radiation = 0;

  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto solar_position_calc = std::make_shared<Heat::SolarPositionCalculatorConstant>(Heat::DirectionCosines{1, 0, 0});
  Heat::TarpModel tarp_model(1, 4, 0, {0, 0, 1}, {0, 0, -1}, 1, 100, 0);
  Heat::SkyRadiationModel sky_model(1, {0, 0, 1});
  Heat::SolarRadiationModel solar_model(1);

  Heat::SteadyStateTempCaluclator calc(env_interface, solar_position_calc, tarp_model, sky_model, solar_model);
  double val = calc.calculate(0, 2, 0.1);
  EXPECT_EQ(val, edata.air_temp);
}

TEST(SteadyStateTempCalc, AirPlusSunOnly)
{
  Heat::EnvironmentData edata;
  edata.air_temp = 300;
  edata.air_speed = 5;
  edata.air_direction = {1, 0, 0};
  edata.ir_horizontal_radiation = 5.6697E-8 * std::pow(edata.air_temp, 4);
  edata.direct_normal_radiation = 300;
  edata.diffuse_radiation = 0;

  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto solar_position_calc = std::make_shared<Heat::SolarPositionCalculatorConstant>(Heat::DirectionCosines{0.5, 0.5, 1.0/std::sqrt(8)});
  Heat::TarpModel tarp_model(1, 4, 0, {0, 0, 1}, {0, 0, -1}, 1, 100, 0);
  Heat::SkyRadiationModel sky_model(1, {0, 0, 1});
  Heat::SolarRadiationModel solar_model(1);

  Heat::SteadyStateTempCaluclator calc(env_interface, solar_position_calc, tarp_model, sky_model, solar_model);
  double val = calc.calculate(0, 2, 0.1);
  EXPECT_GT(val, edata.air_temp);
}