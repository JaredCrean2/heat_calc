#include "gtest/gtest.h"
#include "math.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/dates.h"

/*
TEST(SolarPosition, Noon)
{
  Heat::AzimuthZenith az_prev{0, 0};
  for (int i=0; i < 24; ++i)
  {
    auto az = Heat::solar::computeAzimuthZenith(computeJulianDate({1, 1, 2000}), i, 0, 0, 0);
    if (i > 0)
    {
      if (i <= 12)
        EXPECT_LE(std::acos(az.cos_zenith), std::acos(az_prev.cos_zenith));
      else
        EXPECT_GE(std::acos(az.cos_zenith), std::acos(az_prev.cos_zenith));

//      if (i <= 6)
//        EXPECT_LE(std::acos(az.cos_azimuth),  std::acos(az_prev.cos_azimuth));
//      else
//        EXPECT_GE(std::acos(az.cos_azimuth),  std::acos(az_prev.cos_azimuth));
//
    }

    az_prev = az;

    std::cout << "at hour " << i << ", azimuth = " << radiansToDegrees(std::acos(az.cos_azimuth)) << ", zenith = " << radiansToDegrees(std::acos(az.cos_zenith)) << std::endl;

    //std::cout << "at hour " << i << ", azimuth = " << radiansToDegrees(std::acos(az.cos_azimuth)) << ", zenith = " << radiansToDegrees(std::acos(az.cos_zenith)) << std::endl;
  }
}
*/



// Test against onlin calculator
// from https://gml.noaa.gov/grad/solcalc/azel.html



TEST(SolarPosition, Denver6am)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 6;
  auto az = Heat::solar::computeAzimuthZenith({1, 1, 2000}, time, time_zone, longitude, latitude);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -15.06)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(107.63)), 1e-1);
}

TEST(SolarPosition, DenverNoon)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 12;
  auto az = Heat::solar::computeAzimuthZenith({1, 1, 2000}, time, time_zone, longitude, latitude);
  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 27.28)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(179.12)), 1e-1);

}


TEST(SolarPosition, Denver3pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 15;
  auto az = Heat::solar::computeAzimuthZenith({1, 1, 2000}, time, time_zone, longitude, latitude);


  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 15.02)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(221.57)), 1e-1);
}

TEST(SolarPosition, Denver6pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 18;
  auto az = Heat::solar::computeAzimuthZenith({1, 1, 2000}, time, time_zone, longitude, latitude);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -13.79)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(251.42)), 1e-1);
}


TEST(SolarPositionLowPrecision, Denver6am)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 6;
  auto az = Heat::solar::computeAzimuthZenith_lowprecision({1, 1, 2000}, time, time_zone, longitude, latitude);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -15.06)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(107.63)), 1e-1);
}

TEST(SolarPositionLowPrecision, DenverNoon)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99); 
  int time_zone = 7;
  int time = 12;
  auto az = Heat::solar::computeAzimuthZenith_lowprecision({1, 1, 2000}, time, time_zone, longitude, latitude);
  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 27.28)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(179.12)), 1e-1);

}


TEST(SolarPositionLowPrecision, Denver3pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99); 
  int time_zone = 7;
  int time = 15;
  auto az = Heat::solar::computeAzimuthZenith_lowprecision({1, 1, 2000}, time, time_zone, longitude, latitude);


  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 15.02)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(221.57)), 1e-1);
}

TEST(SolarPositionLowPrecision, Denver6pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 18;
  auto az = Heat::solar::computeAzimuthZenith_lowprecision({1, 1, 2000}, time, time_zone, longitude, latitude);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -13.79)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(251.42)), 1e-1);
}



