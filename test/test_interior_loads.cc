#include "gtest/gtest.h"
#include "physics/heat/interior_loads.h"

TEST(InteriorLoadsConstant, Value)
{
  Real load = 100;
  Heat::InteriorLoadsConstant model(load);

  EXPECT_EQ(model.computeLoadPower(), load);
}