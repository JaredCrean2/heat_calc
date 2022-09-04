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

class HeatEquationSolarTester : public StandardDiscSetup, public ::testing::Test
{
  protected:
    HeatEquationSolarTester()
    {
      auto meshspec = Mesh::getMeshSpec(0, 2, 0, 2, 0, 2, 3, 3, 3);
      std::vector<bool> is_surf_dirichlet{false, false, false, false, false, false};
      setup(3, 1, meshspec, is_surf_dirichlet);
      sol_vec = makeDiscVector(disc);
      sol_vec->set(0);

      makeHeatEquationSolar();
    }

    void makeHeatEquationSolar()
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
      //Real wall_temp = 320;
      //Real t_start = 0;
      Heat::SolarPositionCalculator solar_position_calc(0, 0, 0, 0);

      
      auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
      auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
      auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
      auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
      auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

      auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
        air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp, hvac_restore_time);
      heat_eqn = std::make_shared<Heat::HeatEquationSolar>(disc, solar_position_calc, env_interface, air_updator);;

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
      heat_eqn->addNeumannBC(wall_bc, false);
      heat_eqn->getAuxEquations()->getBlockSolution(1)[0] = initial_air_temp;

      heat_eqn->addVolumeGroupParams({1, 1,1 });
    }

  protected:
    std::shared_ptr<Heat::HeatEquationSolar> heat_eqn;
    DiscVectorPtr sol_vec;
};
}

TEST_F(HeatEquationSolarTester, dRdTair)
{
  //TODO: set rhs to non-zero value
  DiscVectorPtr rhs0 = makeDiscVector(disc); rhs0->set(666);
  DiscVectorPtr rhs1 = makeDiscVector(disc); rhs1->set(555);
  DiscVectorPtr dRdTair_fd = makeDiscVector(disc); 
  ArrayType<Real, 1> dRdTair(boost::extents[disc->getDofNumbering()->getNumLocalDofs()]);
  for (size_t i=0; i < dRdTair.shape()[0]; ++i)
    dRdTair[i] = 666;

  Real eps = 1e-7;
  Real t_air = 298;
  sol_vec->set(320);
  heat_eqn->getAuxEquations()->getBlockSolution(1)[0] = t_air;

  heat_eqn->computedRdTinterior_airProduct(sol_vec, 0, 1, dRdTair);

  heat_eqn->computeRhs(sol_vec, 0.0, rhs0);
  heat_eqn->getAuxEquations()->getBlockSolution(1)[0] = t_air + eps;
  heat_eqn->computeRhs(sol_vec, 0.0, rhs1);

  if (!rhs0->isVectorCurrent())
    rhs0->syncArrayToVector();

  if (!rhs1->isVectorCurrent())
    rhs1->syncArrayToVector();

  auto& rhs0_vec = rhs0->getVector();
  auto& rhs1_vec = rhs1->getVector();
  auto& dRdTair_fd_vec = dRdTair_fd->getVector();

  for (size_t i=0; i < dRdTair_fd_vec.shape()[0]; ++i)
    dRdTair_fd_vec[i] = (rhs1_vec[i] - rhs0_vec[i])/eps;

  for (size_t i=0; i < dRdTair_fd_vec.shape()[0]; ++i)
    EXPECT_NEAR(dRdTair[i], dRdTair_fd_vec[i], 1e-5);
}