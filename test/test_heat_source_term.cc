#include "gtest/gtest.h"
#include "bounding_box.h"
#include "mesh_helper.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/source_terms_def.h"
#include "physics/heat/temperature_controller.h"

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
  Heat::SourceTermSolarHeating heating(disc->getVolDisc(0), collector_area, efficiency, 0, {0, 0, 1}, output_area);

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

TEST_F(HeatSourceTermTester, Derivative)
{
  setup(3, 1);
  Real t_lower = 290, t_upper = 300;
  Real collector_area = 2.0, efficiency = 0.8, output_area = 10.0;
  utils::BoundingBox box;
  auto controller = std::make_shared<Heat::TemperatureControllerHeatQuadratic>(t_lower, t_upper);
  Heat::SourceTermSolarHeating heating(disc->getVolDisc(0), collector_area, efficiency, 0, {0, 0, 1}, output_area,
                                       box, controller);

  Real direct_normal_radiation = 3;
  Real diffuse_radiation = 2.0;
  Real t_air = 291;
  Heat::DirectionCosines solar_position = {std::sqrt(1 - 1.0/16), 0, 0.25};
  heating.setInteriorAirTemperature(t_air);
  heating.setDirectNormalRadiation(direct_normal_radiation);
  heating.setDiffuseRadiation(diffuse_radiation);
  heating.setSolarDirection(solar_position);


  Real eps = 1e-7;
  int npts = disc->getVolDisc(0)->getNumQuadPtsPerElement();
  std::vector<Real> vals1(npts), vals2(npts), vals_fd(npts), vals_ad(npts);

  heating.getValues(0, 0.0, vals1.data());
  heating.setInteriorAirTemperature(t_air + eps);
  heating.getValues(0, 0.0, vals2.data());
  heating.setInteriorAirTemperature(t_air);

  for (int i=0; i < npts; ++i)
    vals_fd[i] = (vals2[i] - vals1[i])/eps;

  heating.getValues_dTair(0, 0.0, 1, vals_ad.data());

  for (int i=0; i < npts; ++i)
    EXPECT_NEAR(vals_fd[i], vals_ad[i], 1e-6);
}