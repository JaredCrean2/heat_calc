#include "gtest/gtest.h"
#include "physics/heat/air_leakage.h"

namespace {

void test_finite_difference(Heat::AirLeakageModel& model, Real t_interior)
{
  Real eps = 1e-7;
  Real flux1 = model.computeAirLeakagePower(t_interior);
  Real flux2 = model.computeAirLeakagePower(t_interior + eps);
  Real deriv_fd = (flux2 - flux1)/eps;

  Real deriv_ad;
  Real flux3 = model.computeAirLeakagePowerDot(t_interior, deriv_ad);

  EXPECT_NEAR(flux1, flux3, 1e-13);
  EXPECT_NEAR(deriv_fd, deriv_ad, 1e-6);
}
}

TEST(AirLeakageNatural, Value)
{
  Real ach50 = 2;
  Real expected_pressure = 3;
  Real volume = 4;
  Real cp = 5;
  Real rho = 6;
  Real t_interior = 7;
  Real t_exterior = 10;

  Heat::AirLeakageModelPressure model(ach50, expected_pressure, volume, cp, rho);
  model.setExteriorTemperature(t_exterior);

  EXPECT_NEAR(model.computeAirLeakagePower(t_interior), (t_interior - t_exterior) * (expected_pressure/50) * ach50 * volume * rho * cp / 3600, 1e-13);
  EXPECT_LT(model.computeAirLeakagePower(t_interior), 0);
  test_finite_difference(model, t_interior);
}


TEST(AirLeakagePressure, Value)
{
  Real flow_rate = 3;
  Real efficiency = 0.9;
  Real cp = 5;
  Real rho = 6;
  Real t_interior = 7;
  Real t_exterior = 10;

  Heat::HRVModel model(flow_rate, efficiency, cp, rho);
  model.setExteriorTemperature(t_exterior);
  
  EXPECT_NEAR(model.computeAirLeakagePower(t_interior), (1 - efficiency) * (t_interior - t_exterior) * flow_rate * rho * cp, 1e-13);
  EXPECT_LT(model.computeAirLeakagePower(t_interior), 0);
  
  test_finite_difference(model, t_interior);
}

