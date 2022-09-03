#include "gtest/gtest.h"
#include "discretization/disc_vector.h"
#include "mesh_helper.h"
#include "physics/heat/HeatEquation.h"
#include "physics/heat/HeatEquationSolar.h"
#include "physics/heat/air_leakage.h"
#include "physics/heat/bc_defs.h"
#include "physics/heat/interior_temperature_update.h"
#include "physics/heat/environment_interface.h"
#include "physics/heat/window_conduction_model.h"

namespace {

class TemperatureUpdatorTester : public StandardDiscSetup, public ::testing::Test
{
  protected:
    TemperatureUpdatorTester() :
      solar_position_calc(0, 0, 0, 0)
    {
      auto meshspec = Mesh::getMeshSpec(0, 2, 0, 2, 0, 2, 3, 3, 3);
      std::vector<bool> is_surf_dirichlet{false, false, false, false, false, false};
      setup(3, 1, meshspec, is_surf_dirichlet);
    }

  protected:
    Heat::SolarPositionCalculator solar_position_calc;
};
}

TEST_F(TemperatureUpdatorTester, ConstantInteriorLoad)
{
  Heat::EnvironmentData edata{298, 0, {1, 0, 0}, 0, 0, 0};
  Real hvac_restore_time = 5.0;
  Real min_temp = 0;
  Real max_temp = 10000;
  Real rho = 2;
  Real cp = 3;
  Real air_volume = 4;
  Real ach50 = 0;
  Real expected_pressure = 5;
  Real interior_load = 6;
  Real r_val = 7;
  Real window_area = 0;
  Real initial_air_temp = 298;
  Real t_start = 0;
  
  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
  auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
  auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
  auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

  auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
    air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp, hvac_restore_time);
  Heat::HeatEquationSolar heat_eqn(disc, solar_position_calc, env_interface, air_updator);

  auto sol_vec = makeDiscVector(disc);
  sol_vec->set(0);
  heat_eqn.initialize(sol_vec, t_start);

  EXPECT_NEAR(air_updator->computeNetFlux(sol_vec, initial_air_temp, 0.0), interior_load, 1e-13);
}

TEST_F(TemperatureUpdatorTester, ConstantInteriorLoad_UpperLimit)
{
  Heat::EnvironmentData edata{320, 0, {1, 0, 0}, 0, 0, 0};
  Real hvac_restore_time = 5.0;
  Real min_temp = 0;
  Real max_temp = 320;
  Real rho = 2;
  Real cp = 3;
  Real air_volume = 4;
  Real ach50 = 0;
  Real expected_pressure = 5;
  Real interior_load = 6;
  Real r_val = 7;
  Real window_area = 0;
  Real initial_air_temp = max_temp + 20;
  Real t_start = 0;
  
  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
  auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
  auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
  auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

  auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
    air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp, hvac_restore_time);
  Heat::HeatEquationSolar heat_eqn(disc, solar_position_calc, env_interface, air_updator);

  auto sol_vec = makeDiscVector(disc);
  sol_vec->set(0);
  heat_eqn.initialize(sol_vec, t_start);

  Real expected_flux = rho * cp * air_volume * (max_temp - initial_air_temp)/hvac_restore_time;
  EXPECT_NEAR(air_updator->computeNetFlux(sol_vec, initial_air_temp, 0.0), expected_flux ,1e-13);
}

TEST_F(TemperatureUpdatorTester, ConstantInteriorLoad_LowerLimit)
{
  Heat::EnvironmentData edata{298, 0, {1, 0, 0}, 0, 0, 0};
  Real hvac_restore_time = 5;
  Real min_temp = 298;
  Real max_temp = 320;
  Real rho = 2;
  Real cp = 3;
  Real air_volume = 4;
  Real ach50 = 0;
  Real expected_pressure = 5;
  Real interior_load = -6;
  Real r_val = 7;
  Real window_area = 0;
  Real initial_air_temp = min_temp - 20;
  Real t_start = 0;
  
  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
  auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
  auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
  auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

  auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
    air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp, hvac_restore_time);
  Heat::HeatEquationSolar heat_eqn(disc, solar_position_calc, env_interface, air_updator);

  auto sol_vec = makeDiscVector(disc);
  sol_vec->set(0);
  heat_eqn.initialize(sol_vec, t_start);

  Real expected_flux = rho * cp * air_volume * (min_temp - initial_air_temp)/hvac_restore_time;
  EXPECT_NEAR(air_updator->computeNetFlux(sol_vec, initial_air_temp, 0.0), expected_flux ,1e-13);
}


TEST_F(TemperatureUpdatorTester, AirLeakage)
{
  Heat::EnvironmentData edata{320, 0, {1, 0, 0}, 0, 0, 0};
  Real hvac_restore_time = 5;
  Real min_temp = 0;
  Real max_temp = 10000;
  Real rho = 2;
  Real cp = 3;
  Real air_volume = 4;
  Real ach50_air = 6;
  Real ach50_vent = 0;
  Real expected_pressure = 5;
  Real interior_load = 0;
  Real r_val = 7;
  Real window_area = 0;
  Real initial_air_temp = 298;
  Real t_start = 0;
  
  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_air, expected_pressure, air_volume, cp, rho);
  auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_vent, expected_pressure, air_volume, cp, rho);
  auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
  auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

  auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
    air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp, hvac_restore_time);
  Heat::HeatEquationSolar heat_eqn(disc, solar_position_calc, env_interface, air_updator);

  auto sol_vec = makeDiscVector(disc);
  sol_vec->set(0);
  heat_eqn.initialize(sol_vec, t_start);

  Real expected_flux =  ach50_air * (5.0/50) * (edata.air_temp - initial_air_temp) * rho * cp * air_volume / 3600;
  EXPECT_NEAR(air_updator->computeNetFlux(sol_vec, initial_air_temp, 0.0), expected_flux , 1e-13);
}


TEST_F(TemperatureUpdatorTester, Ventilation)
{
  Heat::EnvironmentData edata{320, 0, {1, 0, 0}, 0, 0, 0};
  Real hvac_restore_time = 5;
  Real min_temp = 0;
  Real max_temp = 10000;
  Real rho = 2;
  Real cp = 3;
  Real air_volume = 4;
  Real ach50_air = 0;
  Real ach50_vent = 6;
  Real expected_pressure = 5;
  Real interior_load = 0;
  Real r_val = 7;
  Real window_area = 0;
  Real initial_air_temp = 298;
  Real t_start = 0;
  
  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_air, expected_pressure, air_volume, cp, rho);
  auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_vent, expected_pressure, air_volume, cp, rho);
  auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
  auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

  auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
    air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp, hvac_restore_time);
  Heat::HeatEquationSolar heat_eqn(disc, solar_position_calc, env_interface, air_updator);

  auto sol_vec = makeDiscVector(disc);
  sol_vec->set(0);
  heat_eqn.initialize(sol_vec, t_start);

  Real expected_flux =  ach50_vent * (5.0/50) * (edata.air_temp - initial_air_temp) * rho * cp * air_volume  / 3600;
  EXPECT_NEAR(air_updator->computeNetFlux(sol_vec, initial_air_temp, 0.0), expected_flux , 1e-13);
}


TEST_F(TemperatureUpdatorTester, WindowConduction)
{
  Heat::EnvironmentData edata{320, 0, {1, 0, 0}, 0, 0, 0};
  Real hvac_restore_time = 5;
  Real min_temp = 0;
  Real max_temp = 10000;
  Real rho = 2;
  Real cp = 3;
  Real air_volume = 4;
  Real ach50_air = 0;
  Real ach50_vent = 0;
  Real expected_pressure = 5;
  Real interior_load = 0;
  Real r_val = 7;
  Real window_area = 6;
  Real initial_air_temp = 298;
  Real t_start = 0;
  
  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_air, expected_pressure, air_volume, cp, rho);
  auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_vent, expected_pressure, air_volume, cp, rho);
  auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
  auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

  auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
    air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp, hvac_restore_time);
  Heat::HeatEquationSolar heat_eqn(disc, solar_position_calc, env_interface, air_updator);

  auto sol_vec = makeDiscVector(disc);
  sol_vec->set(0);
  heat_eqn.initialize(sol_vec, t_start);

  Real expected_flux =  (window_area * (edata.air_temp - initial_air_temp) / r_val);
  EXPECT_NEAR(air_updator->computeNetFlux(sol_vec, initial_air_temp, 0.0), expected_flux , 1e-13);
}


TEST_F(TemperatureUpdatorTester, WallConduction)
{
  Heat::EnvironmentData edata{298, 0, {1, 0, 0}, 0, 0, 0};
  Real hvac_restore_time = 5;
  Real min_temp = -10000;
  Real max_temp = 10000;
  Real rho = 2;
  Real cp = 3;
  Real air_volume = 400;
  Real ach50 = 0;
  Real expected_pressure = 5;
  Real interior_load = 0;
  Real r_val = 7;
  Real window_area = 0;
  Real initial_air_temp = 298;
  Real wall_temp = 320;
  Real t_start = 0;
  
  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
  auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
  auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
  auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

  auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
    air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp, hvac_restore_time);
  Heat::HeatEquationSolar heat_eqn(disc, solar_position_calc, env_interface, air_updator);

  Real surface_area = 2;
  Real perimeter = 3;
  int roughness_index = 0;
  std::array<Real, 3> vertical_vector = {0, 0, 1};
  std::array<Real, 3> point_at_zero_altitude = {0, 0, 0};
  int met_terrain_index = 0;
  Real meterological_altitude = 1;
  int local_terrain_index = 0;
  auto wall_bc = std::make_shared<Heat::TarpBC>(disc->getSurfDisc(0), surface_area, perimeter, roughness_index, vertical_vector,
                                                point_at_zero_altitude, met_terrain_index, meterological_altitude, local_terrain_index);
  heat_eqn.addNeumannBC(wall_bc, false);
  heat_eqn.getAuxEquations()->getBlockSolution(1)[0] = initial_air_temp;

  Heat::TarpModel model(surface_area, perimeter, roughness_index, vertical_vector,
                        point_at_zero_altitude,
                        met_terrain_index, meterological_altitude, local_terrain_index );
  model.setAirDirection({0, 1, 0});
  model.setAirSpeed(0);
  model.setAirTemperature(initial_air_temp);
  Real wall_area = 4;
  Real heat_transfer_coeff = model.computeHeatTransferCoeff(wall_temp, {0, 0, 0}, {0, 1, 0});

  auto sol_vec = makeDiscVector(disc);
  sol_vec->set(wall_temp);
  heat_eqn.initialize(sol_vec, t_start);

  Real expected_flux = (wall_temp - initial_air_temp) * wall_area * heat_transfer_coeff;
  EXPECT_NEAR(air_updator->computeNetFlux(sol_vec, initial_air_temp, 0.0), expected_flux , 1e-12);
}



//TODO: test multiple iterations of same time step