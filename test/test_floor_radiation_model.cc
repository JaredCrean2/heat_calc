#include "gtest/gtest.h"
#include "physics/heat/floor_radiation_model.h"

TEST(FloorRadiationModel, Example)
{
  Real window_area = 2;
  std::array<Real, 3> window_normal = {1, 0, 0};
  Real shgc = 0.9;
  Real floor_area = 4;
  Real floor_absorptivity = 0.8;

  Real direct_normal_flux = 100;
  Real diffuse_flux = 50;
  Heat::DirectionCosines cosines = {1, 0, 0};

  Heat::FloorRadiationModel model(window_area, window_normal, shgc, floor_area, floor_absorptivity);
  model.setDirectNormalRadiation(direct_normal_flux);
  model.setDiffuseRadiation(diffuse_flux);
  model.setSolarDirection(cosines);

  EXPECT_NEAR(model.computeFlux(), floor_absorptivity * shgc * window_area * (direct_normal_flux + diffuse_flux)/floor_area, 1e-13);
}

TEST(FloorRadiationModel, Shading)
{
  Real window_area = 2;
  std::array<Real, 3> window_normal = {1, 0, 0};
  Real shgc = 0.9;
  Real floor_area = 4;
  Real floor_absorptivity = 0.8;

  Real direct_normal_flux = 100;
  Real diffuse_flux = 50;
  Heat::DirectionCosines cosines = {-1, 0, 0};

  Heat::FloorRadiationModel model(window_area, window_normal, shgc, floor_area, floor_absorptivity);
  model.setDirectNormalRadiation(direct_normal_flux);
  model.setDiffuseRadiation(diffuse_flux);
  model.setSolarDirection(cosines);

  EXPECT_NEAR(model.computeFlux(), floor_absorptivity * shgc * window_area * diffuse_flux/floor_area, 1e-13);
}