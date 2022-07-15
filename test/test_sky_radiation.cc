#include "gtest/gtest.h"
#include "physics/heat/sky_radiation.h"

TEST(SkyRadiation, Value)
{
  Real emissivity = 0.9;
  Real sigma = 5.6697e-8;
  std::array<Real, 3> vertical{0, 0, 1};
  std::array<Real, 3> unit_normal{0, 1.0/2, std::sqrt(3)/2};
  Real cos_eps = std::sqrt(3)/2;
  Real t_air = 200;
  Real t_ground = t_air;
  Real t_sky = 300;
  Real t_surf = 400;
  Real ir_radiation = sigma * std::pow(t_sky, 4);

  Heat::SkyRadiationModel model(emissivity, vertical);
  model.setIRHorizontalRadiation(ir_radiation);
  model.setAirTemperature(t_air);

  Real beta = std::sqrt(0.5*(1 + cos_eps));
  Real f_ground = 0.5*(1 - cos_eps);
  Real f_sky = beta * 0.5 * (1 + cos_eps);
  Real f_air = (1 - beta) * 0.5 * (1 + cos_eps);

  Real flux_ex = emissivity * sigma * (f_ground * (std::pow(t_ground, 4) - std::pow(t_surf, 4)) +
                                       f_sky    * (std::pow(t_sky, 4)    - std::pow(t_surf, 4)) +
                                       f_air    * (std::pow(t_air, 4)    - std::pow(t_surf, 4)));

  Real flux = model.computeFlux(t_surf, unit_normal);
  EXPECT_NEAR(flux, flux_ex, 1e-13);
}