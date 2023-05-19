#include "gtest/gtest.h"
#include "file/WeatherFileReader.h"
#include "file/WeatherFileWriter.h"

namespace {
void compareEPW(const EPWDataPoint& pt1, const EPWDataPoint& pt2)
{
  EXPECT_EQ(pt1.year, pt2.year);
  EXPECT_EQ(pt1.month, pt2.month);
  EXPECT_EQ(pt1.day, pt2.day);
  EXPECT_EQ(pt1.hour, pt2.hour);
  EXPECT_EQ(pt1.minute, pt2.minute);

  double tol=1e-13;
  EXPECT_NEAR(pt1.temperature, pt2.temperature, tol);
  EXPECT_NEAR(pt1.horizontal_ir_radiation, pt2.horizontal_ir_radiation, tol);
  EXPECT_NEAR(pt1.direct_normal_radiation, pt2.direct_normal_radiation, tol);
  EXPECT_NEAR(pt1.diffuse_radiation, pt2.diffuse_radiation, tol);
  EXPECT_NEAR(pt1.wind_speed, pt2.wind_speed, tol);
  EXPECT_NEAR(pt1.wind_direction, pt2.wind_direction, tol);
}
}

TEST(WeatherFileIO, WriteRead)
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

  EPWDataPoint pt2;
  pt2.year = 2000 + 1;
  pt2.month = 2 + 1;
  pt2.day = 15 + 1;
  pt2.hour = 2 + 1;
  pt2.minute = 3 + 1;

  pt2.temperature = 10 + 1;
  pt2.horizontal_ir_radiation = 2 + 1;
  pt2.direct_normal_radiation = 3 + 1;
  pt2.diffuse_radiation = 4 + 1;
  pt2.wind_speed = 5 + 1;
  pt2.wind_direction = 6 + 1;  

  std::vector<EPWDataPoint> data;
  data.push_back(pt1);
  data.push_back(pt2);

  WeatherFileWriter writer("tmp.wea");
  writer.write(data);

  WeatherFileReader reader("tmp.wea");
  auto data2 = reader.read();

  EXPECT_EQ(data2.size(), 2u);
  compareEPW(data[0], data2[0]);
  compareEPW(data[1], data2[1]);
}