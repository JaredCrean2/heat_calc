#include "gtest/gtest.h"
#include "physics/heat/air_leakage.h"

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

  EXPECT_EQ(model.computeAirLeakagePower(t_interior, t_exterior), (t_interior - t_exterior) * (expected_pressure/50) * ach50 * volume * rho * cp / 3600);
}