#include "gtest/gtest.h"
#include "physics/heat/window_conduction_model.h"

TEST(WindowConductionModel, Value)
{
  Real r_val = 2;
  Real area = 3;
  Real t_interior = 5;
  Real t_exterior = 10;

  Heat::WindowConductionModel model(r_val, area);

  EXPECT_EQ(model.computeConductionPower(t_interior, t_exterior), (t_exterior - t_interior)*area/r_val);
}