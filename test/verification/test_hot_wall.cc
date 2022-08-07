#include "gtest/gtest.h"

#include "discretization/DirichletBC_defs.h"
#include "discretization/NeumannBC_defs.h"
#include "discretization/disc_vector.h"
#include "discretization/source_term.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix.h"
#include "physics/heat/HeatEquation.h"
#include "physics/heat/basis_vals.h"
#include "physics/heat/interior_loads.h"
#include "physics/heat/window_conduction_model.h"
#include "physics/heat/interior_temperature_update.h"
#include "linear_system/large_matrix_petsc.h"
#include "linear_system/sparsity_pattern_mesh.h"
#include "time_solver/crank_nicolson.h"
#include "mesh_helper.h"

namespace {

  class HotWallTester : public StandardDiscSetup,
                                           public testing::Test
  {
    protected:
      HotWallTester()
      {
        setup();
      }

      using StandardDiscSetup::setup;

      void setup_final(const int quad_degree, int sol_degree, const Mesh::MeshSpec& spec)
      {
        const std::vector<bool> is_surf_dirichlet = {true, false, false, false, false, false};
        StandardDiscSetup::setup(quad_degree, sol_degree, spec, is_surf_dirichlet);
      }

      void setSolution(std::shared_ptr<Heat::EnvironmentInterface> environment_interface,
                       std::shared_ptr<Heat::InteriorAirTemperatureUpdator> air_temp_updator)
      {

        Heat::SolarPositionCalculator solar_position_calc{0, 0, 0, 0};

        Heat::VolumeGroupParams params = Heat::VolumeGroupParams{0.04, 45, 2020};
        heat        = std::make_shared<Heat::HeatEquationSolar>(disc, solar_position_calc, environment_interface, air_temp_updator);
        u_vec       = makeDiscVector(disc);
        Real surface_temp = 320;
        Real initial_wall_temp = 310;

        for (int i=0; i < disc->getNumVolDiscs(); ++i)
          heat->addVolumeGroupParams(params);

        for (int i=0; i < disc->getNumBCSurfDiscs(); ++i)
        {
          auto surf = disc->getBCSurfDisc(i);

          if (i == 0)
          {

            // y minus surface is held at constant temperature
            assertAlways(surf->getIsDirichlet(), "first surface must be dirichlet");
            auto f = [=](Real x, Real y, Real z, Real t) { return surface_temp; };
            heat->addDirichletBC(makeDirichletBCMMS(surf, f, &(impl::zeroFunc)));
          } else if (i == 2)
          {
            auto bc = std::make_shared<Heat::TarpBC>(surf, 80.0*144/(39*39), 36.0*12/39, 0, std::array<Real, 3>{0, 0, 1});
            heat->addNeumannBC(bc, false);
          } else
          {
            auto f = [](Real x, Real y, Real z, Real t) { return std::array<Real, 3>{0, 0, 0}; };
            heat->addNeumannBC(makeNeumannBCMMS(surf, f), false);
          }
        }

        auto f = [&](Real x, Real y, Real z)
                    { return initial_wall_temp; };
        u_vec->setFunc(f);
        
        heat->initialize(u_vec, 0);

      }

      std::shared_ptr<Heat::HeatEquationSolar> heat;
      DiscVectorPtr u_vec;
  };  

  linear_system::LargeMatrixOptsPetsc get_options()
  {
    linear_system::LargeMatrixOptsPetsc opts;
    opts.is_structurally_symmetric = false;
    opts.is_value_symmetric        = false;
    opts.factor_in_place           = false;
    opts.petsc_opts["ksp_atol"] = "1e-12";
    opts.petsc_opts["ksp_rtol"] = "1e-50";
    opts.petsc_opts["ksp_monitor"] = "";

    return opts;
  }
}



TEST_F(HotWallTester, CaseOne)
{
  //Real kappa = 2;
  //Heat::VolumeGroupParams params{kappa, 3, 4};

  Real kappa = 1;
  Heat::VolumeGroupParams params{kappa, 1, 1};


  Heat::EnvironmentData edata{298, 0, {1, 0, 0}, 0, 0, 0};
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
  Real window_area = 0;
  Real initial_air_temp = 298;
  
  auto env_interface = std::make_shared<Heat::EnvironmentInterfaceConstant>(edata);
  auto air_leakage   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_air, expected_pressure, air_volume, cp, rho);
  auto ventilation   = std::make_shared<Heat::AirLeakageModelPressure>(ach50_vent, expected_pressure, air_volume, cp, rho);
  auto interior_loads = std::make_shared<Heat::InteriorLoadsConstant>(interior_load);
  auto window_conduction = std::make_shared<Heat::WindowConductionModel>(r_val, window_area);

  auto air_updator = std::make_shared<Heat::InteriorAirTemperatureUpdator>(min_temp, max_temp, rho*cp, air_volume, 
    air_leakage, ventilation, interior_loads, window_conduction, initial_air_temp);


  timesolvers::TimeStepperOpts opts;
  opts.t_start = 0.0;
  opts.delta_t = 1.0/60; 
  opts.t_end   = 10*opts.delta_t;
  opts.mat_type = linear_system::LargeMatrixType::Petsc;
  opts.matrix_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>(get_options());
  opts.nonlinear_abs_tol = 1e-12;
  opts.nonlinear_rel_tol = 1e-12;
  opts.nonlinear_itermax = 5;  //TODO: test 1

  int sol_degree = 1;
  int nelem = 4;
  auto meshspec = Mesh::getMeshSpec(0, 1, 0, 1, 0, 1, nelem, nelem, nelem);

  HotWallTester::setup_final(2*sol_degree, sol_degree, meshspec);
  setSolution(env_interface, air_updator);

  //mesh->getFieldDataManager().attachVector(u_vec, "solution");
  //mesh->writeVtkFiles(std::string("mesh") + std::to_string(i));

  timesolvers::CrankNicolson crank(heat, u_vec, opts);
  crank.solve();
}