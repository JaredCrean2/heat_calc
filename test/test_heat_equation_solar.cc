#include "gtest/gtest.h"
#include "bounding_box.h"
#include "discretization/disc_vector.h"
#include "mesh_helper.h"
#include "physics/AuxiliaryEquations.h"
#include "physics/heat/HeatEquation.h"
#include "physics/heat/HeatEquationSolar.h"
#include "physics/heat/air_leakage.h"
#include "physics/heat/bc_defs.h"
#include "physics/heat/interior_temperature_update.h"
#include "physics/heat/environment_interface.h"
#include "physics/heat/source_terms_def.h"
#include "physics/heat/temperature_controller.h"
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
      Heat::EnvironmentData edata{298, 0, {1, 0, 0}, 1000, 100, 100};
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
      auto solar_position_calc = std::make_shared<Heat::SolarPositionCalculatorNaval>(Date{1, 1, 2000}, 0, 0, 0);

      
      auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
      auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
      auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
      auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
      auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);
      auto hvac_model = std::make_shared<Heat::HVACModelSwitch>(min_temp, max_temp, rho*cp, air_volume, hvac_restore_time);

      auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(rho*cp, air_volume, 
        air_leakage, ventilation, interior_loads, window_conduction, hvac_model);
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
      air_updator->setBCs({wall_bc});

      Real collector_area = 2;
      Real collector_efficiency = 0.8;
      Real collector_emissivity = 0.7;
      const std::array<Real, 3>& collector_unit_normal = {0, 0, 1};
      Real output_area= 4;
      utils::BoundingBox box;
      auto controller = std::make_shared<Heat::TemperatureControllerHeatQuadratic>(298, 305);      
      auto nonlinear_source_term = std::make_shared<Heat::SourceTermSolarHeating>(disc->getVolDisc(0), collector_area,
                                                    collector_efficiency, collector_emissivity, collector_unit_normal, output_area,
                                                    box, controller);
      heat_eqn->addSourceTerm(0, nonlinear_source_term);

      sol_aux_vec = makeAuxiliaryEquationsStorage(heat_eqn->getAuxEquations());
      sol_aux_vec->getVector(1)[0] = initial_air_temp;

      heat_eqn->addVolumeGroupParams({1, 1, 1});
    }

  protected:
    std::shared_ptr<Heat::HeatEquationSolar> heat_eqn;
    DiscVectorPtr sol_vec;
    AuxiliaryEquationsStoragePtr sol_aux_vec;
};
}

TEST_F(HeatEquationSolarTester, dRdTair)
{
  DiscVectorPtr rhs0 = makeDiscVector(disc); rhs0->set(666);
  DiscVectorPtr rhs1 = makeDiscVector(disc); rhs1->set(555);
  DiscVectorPtr dRdTair_fd = makeDiscVector(disc); 
  ArrayType<Real, 1> dRdTair(boost::extents[disc->getDofNumbering()->getNumLocalDofs()]);
  for (size_t i=0; i < dRdTair.shape()[0]; ++i)
    dRdTair[i] = 666;

  Real t = 12*60*60;
  Real eps = 1e-7;
  Real t_air = 298;
  sol_vec->set(320);
  sol_aux_vec->getVector(1)[0] = t_air;

  heat_eqn->computedRdTinterior_airProduct(sol_vec, t_air, t, 1, dRdTair);

  heat_eqn->computeRhs(sol_vec, sol_aux_vec, t, rhs0);
  sol_aux_vec->getVector(1)[0] = t_air + eps;

  heat_eqn->computeRhs(sol_vec, sol_aux_vec, t, rhs1);

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