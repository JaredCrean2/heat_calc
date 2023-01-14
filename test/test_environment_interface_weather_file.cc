#include "gtest/gtest.h"
#include "physics/heat/environment_interface.h"
#include "physics/heat/environment_interface_weather_file.h"

TEST(EnvironmentInterfaceWeatherFile, ExactValues)
{
  EPWDataPoint pt1;
  pt1.year = 2000;
  pt1.month = 2;
  pt1.day = 15;
  pt1.hour = 2;
  pt1.minute = 3;

  pt1.temperature = 10;
  pt1.horizontal_ir_radiation = 2;
  pt1.direct_normal_radiation = 3;
  pt1.diffuse_radiation = 4;
  pt1.wind_speed = 5;
  pt1.wind_direction = 6;

  EPWDataPoint pt2 = pt1;
  pt2.hour += 1;
  pt2.temperature += 1;


  std::vector<EPWDataPoint> data = {pt1, pt2};

  Heat::EnvironmentInterfaceWeatherFile interface(data);

  Heat::EnvironmentData data1 = interface.getEnvironmentData(0);
  EXPECT_EQ(data1.air_temp, pt1.temperature);

  Heat::EnvironmentData data2 = interface.getEnvironmentData(60*60);
  EXPECT_EQ(data2.air_temp, pt2.temperature);
}

TEST(EnvironmentInterfaceWeatherFile, InterpolatedValues)
{
  EPWDataPoint pt1;
  pt1.year = 2000;
  pt1.month = 2;
  pt1.day = 15;
  pt1.hour = 2;
  pt1.minute = 3;

  pt1.temperature = 10;
  pt1.horizontal_ir_radiation = 2;
  pt1.direct_normal_radiation = 3;
  pt1.diffuse_radiation = 4;
  pt1.wind_speed = 5;
  pt1.wind_direction = 6;

  EPWDataPoint pt2 = pt1;
  pt2.hour += 1;
  pt2.temperature += 1;


  std::vector<EPWDataPoint> data = {pt1, pt2};

  Heat::EnvironmentInterfaceWeatherFile interface(data);

  Heat::EnvironmentData data1 = interface.getEnvironmentData(60*60.0/4);
  EXPECT_NEAR(data1.air_temp, pt1.temperature + 1.0/4, 1e-8);

  Heat::EnvironmentData data2 = interface.getEnvironmentData(60*60.0/2);
  EXPECT_NEAR(data2.air_temp, pt1.temperature + 1.0/2, 1e-8);

    Heat::EnvironmentData data3 = interface.getEnvironmentData(60*60.0*3/4);
  EXPECT_NEAR(data3.air_temp, pt1.temperature + 3.0/4, 1e-8);

}