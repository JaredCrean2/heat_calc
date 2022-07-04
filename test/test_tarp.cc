#include "gtest/gtest.h"
#include "physics/heat/tarp.h"

TEST(Tarp, ForcedConvection)
{
  Real temp = 20;    // the Tarp model has natural convection
                     // heat transfer coefficient independent of the
                     // delta_t, so set this to zero to test only
                     // forced convection term
  Real surface_area   = 4;
  Real perimeter      = 8;
  std::array<Real, 3> vertical_vector{0, 1, 0};
  std::array<Real, 3> pt_at_zero_altitude{1, 1, 0};
  std::array<Real, 3> air_direction{0, 1, 0};
  int met_terrain_index   = 0;
  Real met_altitude       = 10;
  int local_terrain_index = 1;
  Real air_speed          = 2;
  Real local_air_speed    = air_speed * std::pow(270 / met_altitude, 0.14) * std::pow(1.0 / 370, 0.22);
  std::array<Real, 3> pt{2, 2, 0};
  std::array<Real, 3> normal{0, 1, 0};
  Real Wf = 0.5;

  std::vector<Real> roughness_factors{2.17, 1.67, 1.52, 1.13, 1.11, 1.00};

  for (size_t roughness_type=0; roughness_type < roughness_factors.size(); ++roughness_type)
  {
    Heat::TarpModel model(surface_area, perimeter, roughness_type, vertical_vector,
                          pt_at_zero_altitude, met_terrain_index, met_altitude,
                          local_terrain_index);

    model.setAirTemperature(temp);
    model.setAirSpeed(air_speed);
    model.setAirDirection(air_direction);

    Real expected_val = 2.537 * Wf * roughness_factors[roughness_type] * std::sqrt(perimeter * local_air_speed / surface_area);

    EXPECT_NEAR(model.computeHeatTransferCoeff(temp, pt, normal), expected_val, 1e-13);
  }
}


TEST(Tarp, NaturalConvection)
{
  Real delta_t        = 10;
  Real air_temp       = 20;    
  Real surface_area   = 4;
  Real perimeter      = 8;
  std::array<Real, 3> vertical_vector{0, 1, 0};
  std::array<Real, 3> pt_at_zero_altitude{1, 1, 0};
  std::array<Real, 3> air_direction{0, 1, 0};
  int met_terrain_index   = 0;
  Real met_altitude       = 10;
  int local_terrain_index = 1;
  Real air_speed          = 0;
  std::array<Real, 3> pt{2, 2, 0};
  std::array<Real, 3> normal_upward{0, std::sqrt(3)/2, 0.5};
  std::array<Real, 3> normal_downward{0, -std::sqrt(3)/2, 0.5};
  std::array<Real, 3> normal_horizontal{0, 0, 1};
  Real cos_epsilon = std::sqrt(3)/2;
  int roughness_type = 0;

  Heat::TarpModel model(surface_area, perimeter, roughness_type, vertical_vector,
                        pt_at_zero_altitude, met_terrain_index, met_altitude,
                        local_terrain_index);

  model.setAirTemperature(air_temp);
  model.setAirSpeed(air_speed);
  model.setAirDirection(air_direction);

  Real expected_val_upward = 1.810 * std::pow(delta_t, 1.0/3.0) / (1.382 + cos_epsilon);
  Real expected_val_downward = 9.482 * std::pow(delta_t, 1.0/3.0) / (7.283 - cos_epsilon);
  Real expected_val_horizontal = 1.31 * std::pow(delta_t, 1.0/3.0);


  EXPECT_NEAR(model.computeHeatTransferCoeff(air_temp - delta_t, pt, normal_upward),   expected_val_upward, 1e-13);
  EXPECT_NEAR(model.computeHeatTransferCoeff(air_temp + delta_t, pt, normal_downward), expected_val_upward, 1e-13);

  EXPECT_NEAR(model.computeHeatTransferCoeff(air_temp - delta_t, pt, normal_downward), expected_val_downward, 1e-13);
  EXPECT_NEAR(model.computeHeatTransferCoeff(air_temp + delta_t, pt, normal_upward),   expected_val_downward, 1e-13);

  EXPECT_NEAR(model.computeHeatTransferCoeff(air_temp + delta_t, pt, normal_horizontal), expected_val_horizontal, 1e-2);
  EXPECT_NEAR(model.computeHeatTransferCoeff(air_temp - delta_t, pt, normal_horizontal), expected_val_horizontal, 1e-2);

}

