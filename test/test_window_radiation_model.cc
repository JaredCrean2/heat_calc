#include "gtest/gtest.h"
#include "physics/heat/window_radiation_model.h"

TEST(WindowRadiationModel, NoShadow)
{
  Real window_area = 2;
  Real floor_area = 3;
  Real solar_heat_gain_coefficient = 0.9;
  std::array<Real, 3> window_normal = {1, 0, 0};
  Real horizontal_ir_flux = 4;
  Real direct_normal_flux = 5;
  Real diffuse_flux       = 6;
  Heat::DirectionCosines solar_direction{std::sqrt(3)/2, 1.0/2.0, 0};

  Heat::WindowRadiationModel model(window_area, floor_area, solar_heat_gain_coefficient, window_normal);
  model.setDirectNormalRadiation(direct_normal_flux);
  model.setDiffuseRadiation(diffuse_flux);
  model.setIRHorizontalRadiation(horizontal_ir_flux);
  model.setSolarDirection(solar_direction);

  Real total_flux = solar_heat_gain_coefficient * (horizontal_ir_flux + direct_normal_flux * std::sqrt(3)/2 + diffuse_flux);
  EXPECT_NEAR(model.computeFlux(), total_flux * window_area / floor_area, 1e-13);
}


TEST(WindowRadiationModel, BuildingShadow)
{
  Real window_area = 2;
  Real floor_area = 3;
  Real solar_heat_gain_coefficient = 0.9;
  std::array<Real, 3> window_normal = {-1, 0, 0};
  Real horizontal_ir_flux = 4;
  Real direct_normal_flux = 5;
  Real diffuse_flux       = 6;
  Heat::DirectionCosines solar_direction{std::sqrt(3)/2, 1.0/2.0, 0};

  Heat::WindowRadiationModel model(window_area, floor_area, solar_heat_gain_coefficient, window_normal);
  model.setDirectNormalRadiation(direct_normal_flux);
  model.setDiffuseRadiation(diffuse_flux);
  model.setIRHorizontalRadiation(horizontal_ir_flux);
  model.setSolarDirection(solar_direction);

  Real total_flux = solar_heat_gain_coefficient * (horizontal_ir_flux + direct_normal_flux * 0 + diffuse_flux);
  EXPECT_NEAR(model.computeFlux(), total_flux * window_area / floor_area, 1e-13);
}


TEST(WindowRadiationModelNoShading, OverhangShadow)
{
  Real window_area = 2;
  Real floor_area = 3;
  Real solar_heat_gain_coefficient = 0.9;
  std::array<Real, 3> window_normal = {1, 0, 0};
  Real horizontal_ir_flux = 4;
  Real direct_normal_flux = 5;
  Real diffuse_flux       = 6;
  Heat::DirectionCosines solar_direction{std::sqrt(3)/2, 0, 1.0/2.0};
  Real shading_cos_zenith = 1.0/2.0 - 0.1;

  Heat::WindowRadiationModel model(window_area, floor_area, solar_heat_gain_coefficient, window_normal, shading_cos_zenith);
  model.setDirectNormalRadiation(direct_normal_flux);
  model.setDiffuseRadiation(diffuse_flux);
  model.setIRHorizontalRadiation(horizontal_ir_flux);
  model.setSolarDirection(solar_direction);

  Real total_flux = solar_heat_gain_coefficient * (horizontal_ir_flux + direct_normal_flux * 0 + diffuse_flux);
  EXPECT_NEAR(model.computeFlux(), total_flux * window_area / floor_area, 1e-13);
}