#include "gtest/gtest.h"
#include "math.h"
#include "physics/heat/solar_position.h"
#include "physics/heat/solar_position_calculator.h"
#include "physics/heat/dates.h"


// Test against online calculator
// from https://gml.noaa.gov/grad/solcalc/azel.html

TEST(SolarPosition, Denver6am)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 6;
  auto az = Heat::solar::computeAzimuthZenith({1, 1, 2000}, time, time_zone, latitude, longitude);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -15.06)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(107.63)), 1e-1);
}

TEST(SolarPositionCalculator, Denver6am)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 6;
  Heat::SolarPositionCalculatorNaval calc(Date{1, 1, 2000}, time_zone, latitude, longitude);
  auto cosines = calc.computePosition(Real(time)*3600);
  auto az = Heat::solar::computeAzimuthZenith(cosines);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -15.06)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(107.63)), 1e-1);
  EXPECT_GE(cosines.cs1, 0);
}

TEST(SolarPosition, DenverNoon)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 12;
  auto az = Heat::solar::computeAzimuthZenith({1, 1, 2000}, time, time_zone, latitude, longitude);
  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 27.28)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(179.12)), 1e-1);
}

TEST(SolarPositionCalculator, DenverNoon)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 12;
  Heat::SolarPositionCalculatorNaval calc({1, 1, 2000}, time_zone, latitude, longitude);

  auto az_low = Heat::solar::computeAzimuthZenith_lowprecision({1, 1, 2000}, time, time_zone, latitude, longitude);

  auto cosines = calc.computePosition(Real(time)*3600);
  auto az = Heat::solar::computeAzimuthZenith(cosines);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 27.28)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(179.12)), 1e-1);

  EXPECT_NEAR(az_low.cos_zenith, std::cos(degreesToRadians(90 - 27.28)), 1e-1);
  EXPECT_NEAR(az_low.cos_azimuth, std::cos(degreesToRadians(179.12)), 1e-1);  
}

TEST(SolarPosition, Denver3pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 15;
  auto az = Heat::solar::computeAzimuthZenith({1, 1, 2000}, time, time_zone, latitude, longitude);


  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 15.02)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(221.57)), 1e-1);
}

TEST(SolarPositionCalculator, Denver3pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 15;
  Heat::SolarPositionCalculatorNaval calc({1, 1, 2000}, time_zone, latitude, longitude);
  auto cosines = calc.computePosition(Real(time)*3600);
  auto az = Heat::solar::computeAzimuthZenith(cosines);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 15.02)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(221.57)), 1e-1);
}

TEST(SolarPosition, Denver6pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 18;
  auto az = Heat::solar::computeAzimuthZenith({1, 1, 2000}, time, time_zone, latitude, longitude);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -13.79)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(251.42)), 1e-1);
}

TEST(SolarPositionCalculator, Denver6pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 18;
  Heat::SolarPositionCalculatorNaval calc({1, 1, 2000}, time_zone, latitude, longitude);
  auto cosines = calc.computePosition(Real(time)*3600);
  auto az = Heat::solar::computeAzimuthZenith(cosines);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -13.79)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(251.42)), 1e-1);
  EXPECT_LE(cosines.cs1, 0);

}


TEST(SolarPositionLowPrecision, Denver6am)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 6;
  auto az = Heat::solar::computeAzimuthZenith_lowprecision({1, 1, 2000}, time, time_zone, latitude, longitude);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -15.06)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(107.63)), 1e-1);
}

TEST(SolarPositionLowPrecision, DenverNoon)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99); 
  int time_zone = 7;
  int time = 12;
  auto az = Heat::solar::computeAzimuthZenith_lowprecision({1, 1, 2000}, time, time_zone, latitude, longitude);
  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 27.28)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(179.12)), 1e-1);

}

TEST(SolarPositionLowPrecision, Denver3pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99); 
  int time_zone = 7;
  int time = 15;
  auto az = Heat::solar::computeAzimuthZenith_lowprecision({1, 1, 2000}, time, time_zone, latitude, longitude);


  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - 15.02)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(221.57)), 1e-1);
}

TEST(SolarPositionLowPrecision, Denver6pm)
{
  Real latitude = degreesToRadians(39.74);
  Real longitude = degreesToRadians(104.99);
  int time_zone = 7;
  int time = 18;
  auto az = Heat::solar::computeAzimuthZenith_lowprecision({1, 1, 2000}, time, time_zone, latitude, longitude);

  EXPECT_NEAR(az.cos_zenith, std::cos(degreesToRadians(90 - -13.79)), 1e-1);
  EXPECT_NEAR(az.cos_azimuth, std::cos(degreesToRadians(251.42)), 1e-1);
}



