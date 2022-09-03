#include "gtest/gtest.h"
#include "ProjectDefs.h"
#include "mesh_helper.h"
#include "physics/heat/HeatEquation.h"
#include "physics/heat/HeatEquationSolar.h"
#include "physics/heat/interior_loads.h"
#include "physics/heat/window_conduction_model.h"
#include "time_solver/crank_nicolson.h"
#include "physics/heat/interior_temperature_update.h"
#include "discretization/DirichletBC_defs.h"
#include "discretization/NeumannBC_defs.h"



namespace {

Real degree_space = 5;
Real degree_time = 4;

Real ex_sol(Real x, Real y, Real z, Real t)
{ 
  return std::pow(t, degree_time) + std::pow(x, degree_space) + std::pow(y, degree_space) + std::pow(z, degree_space);
};

std::array<Real, 3> ex_sol_deriv(Real x, Real y, Real z, Real t)
{ 
  std::array<Real, 3> derivs{0, 0, 0};
  if (degree_space > 0)
  {
    derivs[0] = degree_space * std::pow(x, degree_space - 1);
    derivs[1] = degree_space * std::pow(y, degree_space - 1);
    derivs[2] = degree_space * std::pow(z, degree_space - 1);
  }

  return derivs;
}

Real src_func(Real x, Real y, Real z, Real t)
{ 
  Real spacial_term = 0;
  if (degree_space >= 2)
    spacial_term = -(degree_space * (degree_space - 1) * std::pow(x, degree_space-2));

  return spacial_term;
};

class CNDependentTester : public StandardDiscSetup,
                 public testing::Test
{
  protected:
    using HeatPtr = std::shared_ptr<Heat::HeatEquation>;

    CNDependentTester()
    {
      Mesh::MeshSpec spec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, 3, 3, 3);
      setup(3, 1, spec, {false, false, false, false, false, false});
      //setup(3, 1, spec, {true, true, true, true, true, true});

    }

    template <typename Tex, typename Tderiv, typename Tsrc>
    void setSolution(Tex ex_sol, Tderiv deriv, Tsrc src,  
                     std::shared_ptr<Heat::EnvironmentInterface> environment_interface,
                     std::shared_ptr<Heat::InteriorAirTemperatureUpdator> air_temp_updator) 
    {
      const Heat::VolumeGroupParams params{1, 1, 1};
      heat_model = std::make_shared<Heat::HeatEquationSolar>(disc, solar_position_calc, environment_interface, air_temp_updator);
      u_vec   = makeDiscVector(disc);
      //res_vec = makeDiscVector(disc);

      for (int i=0; i < disc->getNumVolDiscs(); ++i)
      {
        heat_model->addSourceTerm(makeSourcetermMMS(disc->getVolDisc(i), src));
        heat_model->addVolumeGroupParams(params);
      }

      for (int i=0; i < disc->getNumBCSurfDiscs(); ++i)
      {
        auto surf = disc->getBCSurfDisc(i);
        if (surf->getIsDirichlet())
          heat_model->addDirichletBC(makeDirichletBCMMS(surf, ex_sol, &(impl::zeroFunc)));
        else
          heat_model->addNeumannBC(makeNeumannBCMMS(surf, deriv), true);
      }


      //res_vec->set(0);
      auto f = [&](Real x, Real y, Real z)
                  { return ex_sol(x, y, z, 0); };
      u_vec->setFunc(f);

      heat_model->initialize(u_vec, 0);

    }

    std::shared_ptr<Heat::HeatEquationSolar> heat_model;
    Heat::SolarPositionCalculator solar_position_calc{0, 0, 0, 0};
    DiscVectorPtr u_vec;
};

}

/*

TEST_F(CNDependentTester, InteriorLoad)
{
  Heat::EnvironmentData edata{298, 0, {1, 0, 0}, 0, 0, 0};
  Real min_temp = 0;
  Real max_temp = 10000;
  Real rho = 2;
  Real cp = 3;
  Real air_volume = 4;
  Real ach50_air = 0;
  Real ach50_vent = 0;
  Real expected_pressure = 5;
  Real interior_load = 6;
  Real r_val = 7;
  Real window_area = 0;
  Real initial_air_temp = 298;
  Real delta_t = 1.0/60;  // time is in hours
  
  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_air, expected_pressure, air_volume, cp, rho);
  auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_vent, expected_pressure, air_volume, cp, rho);
  auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
  auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

  auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
    air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp);

  Real energy_change = interior_load * delta_t * 3600;
  Real temperature_change = energy_change/(rho*cp*air_volume);

  setSolution(ex_sol, ex_sol_deriv, src_func,  env_interface, air_updator);

  timesolvers::TimeStepperOpts opts;
  opts.t_start = 0.0;
  opts.t_end   = 2.5*delta_t;
  opts.delta_t = delta_t;
  opts.mat_type = linear_system::LargeMatrixType::Dense;
  opts.matrix_opts = std::make_shared<linear_system::LargeMatrixOpts>();
  opts.nonlinear_abs_tol = 1e-12;
  opts.nonlinear_rel_tol = 1e-12;
  opts.nonlinear_itermax = 5;  //TODO: test 1

  timesolvers::CrankNicolson crank(heat_model, u_vec, opts);
  crank.solve();

  // the -0.5 is because the previous flux is zero at the initial condition, and the trapizoid
  // rule gives 1/2 of the first timestep flux
  Real timesteps = opts.t_end/delta_t - 0.5;
  Real total_temperature_change = temperature_change * timesteps;
  EXPECT_NEAR(air_updator->getTemperature(), initial_air_temp + total_temperature_change, 1e-10);
}

*/