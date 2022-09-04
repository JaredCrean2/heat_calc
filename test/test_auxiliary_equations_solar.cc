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

class AuxiliaryEquationsSolarTester : public StandardDiscSetup, public ::testing::Test
{
  protected:
    AuxiliaryEquationsSolarTester()
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
      Heat::EnvironmentData edata{exterior_air_temp, 0, {1, 0, 0}, 0, 0, 0};
      Real hvac_restore_time = 5;
      Real min_temp = -10000;
      Real max_temp = 10000;
      Real ach50 = 0;
      Real expected_pressure = 5;
      Real interior_load = 0;
      Heat::SolarPositionCalculator solar_position_calc(0, 0, 0, 0);

      
      auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
      auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
      auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
      auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
      auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

      auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
        air_leakage, ventilation, interior_loads, window_conduction, interior_air_temp, hvac_restore_time);
      heat_eqn = std::make_shared<Heat::HeatEquationSolar>(disc, solar_position_calc, env_interface, air_updator);
      heat_eqn->initialize();
      aux_eqn = heat_eqn->getAuxEquations();
    }

    void makeHeatEquationSolarBC()
    {
      Heat::EnvironmentData edata{exterior_air_temp, 0, {1, 0, 0}, 0, 0, 0};
      Real hvac_restore_time = 5;
      Real min_temp = -10000;
      Real max_temp = 10000;
      Real ach50 = 0;
      Real expected_pressure = 5;
      Real interior_load = 0;
      Heat::SolarPositionCalculator solar_position_calc(0, 0, 0, 0);
      window_area = 0;

      auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
      auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
      auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50, expected_pressure, air_volume, cp, rho);
      auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
      auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

      air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
        air_leakage, ventilation, interior_loads, window_conduction, interior_air_temp, hvac_restore_time);
      heat_eqn = std::make_shared<Heat::HeatEquationSolar>(disc, solar_position_calc, env_interface, air_updator);

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

      heat_eqn->initialize();
      aux_eqn = heat_eqn->getAuxEquations();
    }

  protected:
    std::shared_ptr<Heat::HeatEquationSolar> heat_eqn;
    std::shared_ptr<Heat::InteriorAirTemperatureUpdator> air_updator;
    std::shared_ptr<AuxiliaryEquations> aux_eqn;
    DiscVectorPtr sol_vec;

    Real rho = 2;
    Real cp = 3;
    Real air_volume = 400;
    Real exterior_air_temp = 298;
    Real r_val = 7;
    Real window_area = 4;
    Real interior_air_temp = 320;
};
}

TEST_F(AuxiliaryEquationsSolarTester, Sizes)
{
  EXPECT_EQ(aux_eqn->getNumBlocks(), 2);
  EXPECT_EQ(aux_eqn->getBlockSize(0), heat_eqn->getDiscretization()->getDofNumbering()->getNumLocalDofs());
  EXPECT_EQ(aux_eqn->getBlockSize(1), 1);

  auto jac = aux_eqn->getJacobians()->getMatrix(1);
  EXPECT_EQ(jac->getMLocal(), aux_eqn->getBlockSize(1));
  EXPECT_EQ(jac->getNLocal(), aux_eqn->getBlockSize(1));
}


TEST_F(AuxiliaryEquationsSolarTester, Jacobian)
{
  int block_size = aux_eqn->getBlockSize(1);
  ArrayType<Real, 1> rhs0(boost::extents[block_size]), rhs1(boost::extents[block_size]);
  ArrayType<Real, 1> ones(boost::extents[block_size]), jac_product(boost::extents[block_size]);
  ones[0] = 1;
  sol_vec->set(500);

  Real eps=1e-7;
  aux_eqn->getBlockSolution(1)[0] = interior_air_temp;
  aux_eqn->computeRhs(1, sol_vec, 0.0, rhs0);
  aux_eqn->getBlockSolution(1)[0] = interior_air_temp + eps;
  aux_eqn->computeRhs(1, sol_vec, 0.0, rhs1);
  aux_eqn->getBlockSolution(1)[0] = interior_air_temp;

  Real deriv_fd = (rhs1[0] - rhs0[0])/eps;

  auto jac = aux_eqn->getJacobians()->getMatrix(1);
  auto assembler = std::make_shared<linear_system::SimpleAssembler>(jac);
  aux_eqn->computeJacobian(1, sol_vec, 0, assembler);
  jac->finishMatrixAssembly();
  jac->matVec(ones, jac_product);

  EXPECT_NEAR(jac_product[0], deriv_fd, 1e-5);
}


TEST_F(AuxiliaryEquationsSolarTester, MassMatrix)
{
  ArrayType<Real, 1> ones(boost::extents[1]), product(boost::extents[1]);
  ones[0] = 1;

  aux_eqn->multiplyMassMatrix(1, 0, ones, product);
  EXPECT_NEAR(product[0], rho*cp*air_volume, 1e-13);

  auto jac = aux_eqn->getJacobians()->getMatrix(1);
  auto assembler = std::make_shared<linear_system::SimpleAssembler>(jac);
  aux_eqn->computeMassMatrix(1, 0.0, assembler);

  jac->finishMatrixAssembly();
  jac->matVec(ones, product);
  EXPECT_NEAR(product[0], rho*cp*air_volume, 1e-13);
}

TEST_F(AuxiliaryEquationsSolarTester, dRdTair)
{
  makeHeatEquationSolarBC();
  sol_vec->set(320);

  ArrayType<Real, 1> rhs_aux(boost::extents[aux_eqn->getBlockSize(0)]),
                     rhs_direct(boost::extents[aux_eqn->getBlockSize(0)]);
  ArrayType<Real, 1> ones(boost::extents[aux_eqn->getBlockSize(1)]);
  ones[0] = 1;


  aux_eqn->multiplyOffDiagonal(0, 1, sol_vec, 0.0, ones, rhs_aux);
  heat_eqn->computedRdTinterior_airProduct(sol_vec, 0.0, 1, rhs_direct);

  for (int i=0; i < aux_eqn->getBlockSize(0); ++i)
    EXPECT_NEAR(rhs_aux[i], rhs_direct[i], 1e-13);
}


TEST_F(AuxiliaryEquationsSolarTester, dTairdU)
{
  makeHeatEquationSolarBC();
  sol_vec->set(320);

  ArrayType<Real, 1> rhs_aux(boost::extents[aux_eqn->getBlockSize(1)]);
  ArrayType<Real, 1> ones(boost::extents[aux_eqn->getBlockSize(0)]);
  for (int i=0; i < aux_eqn->getBlockSize(0); ++i)
    ones[i] = 1;

  aux_eqn->multiplyOffDiagonal(1, 0, sol_vec, 0.0, ones, rhs_aux);

  auto sol_vec_bar = makeDiscVector(disc);
  sol_vec_bar->set(0);
  air_updator->computeNetFlux_rev(sol_vec, interior_air_temp, 0.0, sol_vec_bar, 1);

  if (!sol_vec_bar->isVectorCurrent())
    sol_vec_bar->syncArrayToVector();

  Real rhs_direct = 0;
  auto& sol_vec_bar_vec = sol_vec_bar->getVector();
  for (int i=0; i < aux_eqn->getBlockSize(0); ++i)
    rhs_direct += sol_vec_bar_vec[i];

  EXPECT_NEAR(rhs_aux[0], rhs_direct, 1e-13);
}
