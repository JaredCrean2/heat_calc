#include "gtest/gtest.h"
#include "mesh_helper.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/source_terms_def.h"

namespace {
class HeatSourceTermTester : public StandardDiscSetup,
                             public ::testing::Test
{

};
}

TEST_F(HeatSourceTermTester, Values)
{
  setup(3, 1);
  Real collector_area = 2.0, efficiency = 0.8, output_area = 10.0;
  Heat::SourceTermSolarHeating heating(disc->getVolDisc(0), collector_area, efficiency, {0, 0, 1}, output_area);

  Real direct_normal_radiation = 3;
  Real diffuse_radiation = 2.0;
  Heat::DirectionCosines solar_position = {std::sqrt(1 - 1.0/16), 0, 0.25};
  heating.setDirectNormalRadiation(direct_normal_radiation);
  heating.setDiffuseRadiation(diffuse_radiation);
  heating.setSolarDirection(solar_position);

  Real direct_term = solar_position.cs3 * direct_normal_radiation;

  std::vector<Real> vals(disc->getVolDisc(0)->getNumQuadPtsPerElement());
  heating.getValues(0, 0.0, vals.data());
  
  for (int i=0; i < disc->getVolDisc(0)->getNumQuadPtsPerElement(); ++i)
    EXPECT_NEAR(vals[i], collector_area * efficiency * (direct_term + diffuse_radiation) / output_area, 1e-13);
}