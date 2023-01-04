#include "gtest/gtest.h"
#include "physics/heat/window_conduction_model.h"

namespace {

void test_finite_difference(Heat::WindowConductionModel& model, Real t_interior)
{
  Real eps = 1e-7;
  Real flux1 = model.computeConductionPower(t_interior);
  Real flux2 = model.computeConductionPower(t_interior + eps);
  Real deriv_fd = (flux2 - flux1)/eps;

  Real deriv_ad;
  Real flux3 = model.computeConductionPowerDot(t_interior, deriv_ad);

  EXPECT_NEAR(flux1, flux3, 1e-13);
  EXPECT_NEAR(deriv_fd, deriv_ad, 1e-6);
}
}

TEST(WindowConductionModel, Value)
{
  Real r_val = 2;
  Real area = 3;
  Real t_interior = 5;
  Real t_exterior = 10;

  Heat::WindowConductionModel model(r_val, area);

  model.setExteriorTemperature(t_exterior);
  EXPECT_EQ(model.computeConductionPower(t_interior), (t_exterior - t_interior)*area/r_val);
  test_finite_difference(model, t_interior);
}