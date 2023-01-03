#include "gtest/gtest.h"
#include "ProjectDefs.h"
#include "physics/heat/hvac_model.h"

namespace {

void test_derivative_interior_temp(Heat::HVACModel& model, Real interior_temp, Real load_flux, int sign)
{
  const Real eps = sign*1e-7;
  Real flux = model.enforceTemperatureLimit(interior_temp, load_flux);
  Real flux2 = model.enforceTemperatureLimit(interior_temp + eps, load_flux);
  Real flux_deriv = model.enforceTemperatureLimit_dot(interior_temp, 1, load_flux, 0);
  Real flux_fd = (flux2 - flux)/eps;
  std::cout << "flux1 = " << flux << ", flux2 = " << flux2 << ", flux_fd = " << flux_fd << ", flux_deriv = " << flux_deriv << std::endl;
  EXPECT_NEAR(flux_deriv, flux_fd, 1e-6); 
}

void test_derivative_load_flux(Heat::HVACModel& model, Real interior_temp, Real load_flux, int sign)
{
  const Real eps = sign * 1e-7;
  Real flux = model.enforceTemperatureLimit(interior_temp, load_flux);
  Real flux2 = model.enforceTemperatureLimit(interior_temp, load_flux + eps);
  Real flux_deriv = model.enforceTemperatureLimit_dot(interior_temp, 0, load_flux, 1);
  Real flux_fd = (flux2 - flux)/eps;
  std::cout << "flux1 = " << flux << ", flux2 = " << flux2 << ", flux_fd = " << flux_fd << ", flux_deriv = " << flux_deriv << std::endl;

  EXPECT_NEAR(flux_deriv, flux_fd, 1e-6); 
}

void test_derivative_combined(Heat::HVACModel& model, Real interior_temp, Real load_flux, int sign)
{
  const Real eps = sign*1e-7;
  Real interior_temp_dot = 2;  //TODO: make this 2
  Real load_flux_dot = 3;

  Real flux = model.enforceTemperatureLimit(interior_temp, load_flux);
  Real flux2 = model.enforceTemperatureLimit(interior_temp + interior_temp_dot*eps, load_flux + load_flux_dot*eps);
  Real flux_deriv = model.enforceTemperatureLimit_dot(interior_temp, interior_temp_dot, load_flux, load_flux_dot);
  Real flux_fd = (flux2 - flux)/eps;
  std::cout << "flux1 = " << flux << ", flux2 = " << flux2 << ", flux_fd = " << flux_fd << ", flux_deriv = " << flux_deriv << std::endl;

  EXPECT_NEAR(flux_deriv, flux_fd, 1e-6); 
}

void test_reverse_mode(Heat::HVACModel& model, Real interior_temp, Real load_flux)
{
  std::cout << "\nTesting reverse mode" << std::endl;

  Real hvac_flux_bar = 2;
  Real interior_temp_dot = hvac_flux_bar;
  Real load_flux_dot     = hvac_flux_bar;

  Real dflux_d_interior_temp_forward = model.enforceTemperatureLimit_dot(interior_temp, interior_temp_dot, load_flux, 0);
  Real dflux_d_load_flux_forward = model.enforceTemperatureLimit_dot(interior_temp, 0, load_flux, load_flux_dot);

  Real interior_temp_bar = 0;
  Real load_flux_bar = 0;
  model.enforceTemperatureLimit_rev(interior_temp, interior_temp_bar, load_flux, load_flux_bar, hvac_flux_bar);

  EXPECT_NEAR(dflux_d_interior_temp_forward, interior_temp_bar, 1e-13);
  EXPECT_NEAR(dflux_d_load_flux_forward, load_flux_bar, 1e-13);
}


void test_derivative(Heat::HVACModel& model, Real interior_temp, Real load_flux)
{
  test_derivative_interior_temp(model, interior_temp, load_flux,  1);
  test_derivative_interior_temp(model, interior_temp, load_flux, -1);

  test_derivative_load_flux(model, interior_temp, load_flux,  1);
  test_derivative_load_flux(model, interior_temp, load_flux, -1);  

  test_derivative_combined(model, interior_temp, load_flux,  1);
  test_derivative_combined(model, interior_temp, load_flux, -1);  

  test_reverse_mode(model, interior_temp, load_flux);
}

}


TEST(HVACModelSwitch, All)
{
  Real min_temp = 295;
  Real max_temp = 305;
  Real rho_cp = 2;
  Real air_volume = 3;
  Real hvac_restore_time = 4;

  Heat::HVACModelSwitch model(min_temp, max_temp, rho_cp, air_volume, hvac_restore_time);

  Real load_flux = 5;
  EXPECT_EQ(model.enforceTemperatureLimit(296, load_flux), 0);
  EXPECT_EQ(model.enforceTemperatureLimit(300, load_flux), 0);
  EXPECT_EQ(model.enforceTemperatureLimit(304, load_flux), 0);

  EXPECT_NEAR(model.enforceTemperatureLimit(307, load_flux), -rho_cp*air_volume*2/hvac_restore_time - load_flux, 1e-13);
  EXPECT_NEAR(model.enforceTemperatureLimit(293, load_flux), rho_cp*air_volume*2/hvac_restore_time - load_flux, 1e-13);

  test_derivative(model, 296, load_flux);
  test_derivative(model, 300, load_flux);
  test_derivative(model, 304, load_flux);

  test_derivative(model, 293, load_flux);
  test_derivative(model, 307, load_flux);
}

/*
TEST(HVACModelSpline, OutsideSplineRange)
{

  Real min_temp = 295;
  Real max_temp = 305;
  Real rho_cp = 2;
  Real air_volume = 3;
  Real hvac_restore_time = 4;
  Real load_flux = 5;

  Heat::HVACModelSpline model(min_temp, max_temp, rho_cp, air_volume, hvac_restore_time);

  EXPECT_NEAR(model.enforceTemperatureLimit(400, load_flux), -rho_cp*air_volume*95/hvac_restore_time - load_flux, 1e-13);
  EXPECT_NEAR(model.enforceTemperatureLimit(200, load_flux),  rho_cp*air_volume*95/hvac_restore_time - load_flux, 1e-13);

  test_derivative(model, 400, load_flux);
  test_derivative(model, 200, load_flux);

  test_derivative(model, 303, load_flux);
}
*/