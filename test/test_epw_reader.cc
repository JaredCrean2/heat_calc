#include "gtest/gtest.h"
#include "file/EpwReader.h"
#include "physics/heat/environment_interface.h"

TEST(EPWReader, Location)
{
  EPWReader reader("./data/USA_AZ_Phoenix-Sky.Harbor.Intl.AP.722780_TMY3.epw");

  EPWLocation location = reader.getLocation();
  EXPECT_EQ(location.local_name, "Phoenix Sky Harbor Intl Ap");
  EXPECT_EQ(location.state, "AZ");
  EXPECT_EQ(location.country, "USA");
  EXPECT_EQ(location.latitude, 33.45);
  EXPECT_EQ(location.longitude, -111.98);
  EXPECT_EQ(location.time_zone, -7);
}

TEST(EPWReader, RunPeriods)
{
  EPWReader reader("./data/USA_AZ_Phoenix-Sky.Harbor.Intl.AP.722780_TMY3.epw");

  EPWDataPeriods data_periods = reader.getDataPeriods();

  EXPECT_EQ(data_periods.num_periods, 1);
  EXPECT_EQ(data_periods.num_records_per_hour, 1);
}


TEST(EPWReader, FirstDataPoint)
{
  EPWReader reader("./data/USA_AZ_Phoenix-Sky.Harbor.Intl.AP.722780_TMY3.epw");

  auto& data = reader.getData();

  EPWDataPoint pt1 = data[0];
  EXPECT_EQ(pt1.year, 2002);
  EXPECT_EQ(pt1.month, 1);
  EXPECT_EQ(pt1.day, 1);
  EXPECT_EQ(pt1.hour, 1);
  EXPECT_EQ(pt1.minute, 0);

  EXPECT_EQ(pt1.temperature, 10 + 273.15);
  EXPECT_EQ(pt1.horizontal_ir_radiation, 292);
  EXPECT_EQ(pt1.direct_normal_radiation, 0);
  EXPECT_EQ(pt1.diffuse_radiation, 0);
  EXPECT_EQ(pt1.wind_direction, 110);
  EXPECT_EQ(pt1.wind_speed, 1.5);
}

TEST(EPWReader, TenthDataPoint)
{
  EPWReader reader("./data/USA_AZ_Phoenix-Sky.Harbor.Intl.AP.722780_TMY3.epw");

  auto& data = reader.getData();

  EPWDataPoint pt = data[9];
  EXPECT_EQ(pt.year, 2002);
  EXPECT_EQ(pt.month, 1);
  EXPECT_EQ(pt.day, 1);
  EXPECT_EQ(pt.hour, 10);
  EXPECT_EQ(pt.minute, 0);

  EXPECT_EQ(pt.temperature, 11.7 + 273.15);
  EXPECT_EQ(pt.horizontal_ir_radiation, 298);
  EXPECT_EQ(pt.direct_normal_radiation, 782);
  EXPECT_EQ(pt.diffuse_radiation, 47);
  EXPECT_EQ(pt.wind_direction, 40);
  EXPECT_EQ(pt.wind_speed, 1.5);
}

TEST(EPWReader, ConvertToEnvData)
{
  EPWDataPoint pt;
  pt.year = 2000;
  pt.month = 1;
  pt.day = 2;
  pt.hour = 3;
  pt. minute = 4;
  pt.temperature = 300;
  pt.direct_normal_radiation = 5;
  pt.horizontal_ir_radiation = 6;
  pt.diffuse_radiation = 7;
  pt.wind_speed = 8;
  pt.wind_direction = 0;

  Heat::EnvironmentData data = convertToEnvData(pt);

  EXPECT_EQ(data.air_temp, pt.temperature);
  EXPECT_EQ(data.direct_normal_radiation, pt.direct_normal_radiation);
  EXPECT_EQ(data.ir_horizontal_radiation, pt.horizontal_ir_radiation);
  EXPECT_EQ(data.diffuse_radiation, pt.diffuse_radiation);
  EXPECT_EQ(data.air_speed, pt.wind_speed);
  EXPECT_NEAR(data.air_direction[0], 0, 1e-13);
  EXPECT_NEAR(data.air_direction[1], 1, 1e-13);
  EXPECT_NEAR(data.air_direction[2], 0, 1e-13);


  pt.wind_direction = 90;
  data = convertToEnvData(pt);  
  EXPECT_NEAR(data.air_direction[0], 1, 1e-13);
  EXPECT_NEAR(data.air_direction[1], 0, 1e-13);
  EXPECT_NEAR(data.air_direction[2], 0, 1e-13);

  pt.wind_direction = 180;
  data = convertToEnvData(pt);  
  EXPECT_NEAR(data.air_direction[0], 0, 1e-13);
  EXPECT_NEAR(data.air_direction[1], -1, 1e-13);
  EXPECT_NEAR(data.air_direction[2], 0, 1e-13); 

  pt.wind_direction = 270;
  data = convertToEnvData(pt);  
  EXPECT_NEAR(data.air_direction[0], -1, 1e-13);
  EXPECT_NEAR(data.air_direction[1], 0, 1e-13);
  EXPECT_NEAR(data.air_direction[2], 0, 1e-13); 
}