#include "gtest/gtest.h"
#include "physics/heat/wind_velocity.h"


TEST(AshraeWindVelocity, All)
{
  Real v_met = 2;
  Real z_met = 500;
  Real z_local = 200;
  std::vector<Real> alphas{0.14, 0.22, 0.33, 0.10, 0.22};
  std::vector<Real> deltas{270, 370, 460, 210, 370};

  for (size_t met_type=0; met_type < alphas.size(); ++met_type)
    for (size_t local_type=0; local_type < alphas.size(); ++local_type)
    {
      Heat::AshraeWindVelocity calc(met_type, z_met, local_type);
      Real expected_val = v_met * std::pow(deltas[met_type] / z_met, alphas[met_type]) * std::pow(z_local / deltas[local_type], alphas[local_type]);
      EXPECT_NEAR(calc.computeLocalVelocity(v_met, z_local), expected_val, 1e-13);
    }
}