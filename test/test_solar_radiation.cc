#include "gtest/gtest.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/solar_radiation.h"

TEST(SolarRadiation, Value)
{
  Heat::DirectionCosines cosines{std::sqrt(3.0)/2, 1.0/2, 0};
  Real direct_flux = 2;
  Real diffuse_flux = 3;
  Real absorbtivity = 0.7;
  std::array<Real, 3> unit_normal{1, 0, 0};

  Heat::SolarRadiationModel model(absorbtivity);

  model.setDirectNormalRadiation(direct_flux);
  model.setDiffuseRadiation(diffuse_flux);
  model.setSolarDirection(cosines);

  Real flux = model.computeFlux(unit_normal);

  EXPECT_NEAR(flux, absorbtivity*(diffuse_flux + direct_flux*cosines.cs1), 1e-13);

}

TEST(SolarRadiation, ValueBelowHorizon)
{
  Heat::DirectionCosines cosines{-std::sqrt(3.0)/2, 1.0/2, 0};
  Real direct_flux = 2;
  Real diffuse_flux = 3;
  Real absorbtivity = 0.7;
  std::array<Real, 3> unit_normal{1, 0, 0};

  Heat::SolarRadiationModel model(absorbtivity);

  model.setDirectNormalRadiation(direct_flux);
  model.setDiffuseRadiation(diffuse_flux);
  model.setSolarDirection(cosines);

  Real flux = model.computeFlux(unit_normal);

  EXPECT_NEAR(flux, absorbtivity*(diffuse_flux + direct_flux*0), 1e-13);

}