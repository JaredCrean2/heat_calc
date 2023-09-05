#include "gtest/gtest.h"
#include "physics/heat/temperature_controller.h"

namespace {

void test_finite_difference(Heat::TemperatureController& controller, Real temp)
{
  Real eps =  1e-7;
  Real val1 = controller.get_value(temp);
  Real val2 = controller.get_value(temp + eps);
  Real val_fd = (val2 - val1)/eps;

  Real val_dot = controller.get_value_dot(temp, 2)/2;
  EXPECT_NEAR(val_fd, val_dot, 1e-6);
}
}

TEST(TemperatureControllerQuadratic, Values)
{
  Heat::TemperatureControllerHeatQuadratic controller(290, 300);

  EXPECT_FLOAT_EQ(controller.get_value(270), 1.0);
  EXPECT_FLOAT_EQ(controller.get_value(290), 0.25);
  EXPECT_FLOAT_EQ(controller.get_value(295), 0.0);
}

TEST(TemperatureControllerQuadratic, Derivative)
{
  Heat::TemperatureControllerHeatQuadratic controller(290, 300);

  test_finite_difference(controller, 270);
  test_finite_difference(controller, 290);
  test_finite_difference(controller, 296);
}