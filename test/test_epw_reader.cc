#include "gtest/gtest.h"
#include "EpwReader.h"

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