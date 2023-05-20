#include "gtest/gtest.h"
#include "file/WeatherCat.h"
#include "file/WeatherFileReader.h"
#include "file/WeatherFileWriter.h"
#include "arguments_setup.h"

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


TEST(WeatherCat, CommandLineArgs)
{
  CommandLineArguments args({"outputfile", "--file", "file1.wea"});

  WeatherCatParsedData data = parseWeatherCatData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 1u);
  EXPECT_EQ(data.date_ranges.size(), 1u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));
  
  args = CommandLineArguments({"outputfile", "--file", "file1.wea", "1/2/2020", "3/2/2020"});
  data = parseWeatherCatData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 1u);
  EXPECT_EQ(data.date_ranges.size(), 1u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{0, 0}}));  
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{3, 3, 2020}, Time{0, 0}}));  

  args = CommandLineArguments({"outputfile", "--file", "file1.wea", "1/2/2020"});
  EXPECT_ANY_THROW(parseWeatherCatData(args.getArgc(), args.getArgv()));

  args = CommandLineArguments({"outputfile", "--file", "file1.wea", "1/2/2020", "3/2/2020", "--file", "file2.wea"});
  data = parseWeatherCatData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 2u);
  EXPECT_EQ(data.date_ranges.size(), 2u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_EQ(data.filenames[1], "file2.wea");
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{0, 0}}));  
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{3, 3, 2020}, Time{0, 0}}));  
  EXPECT_TRUE(data.date_ranges[1].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));


  args = CommandLineArguments({"outputfile", "--file", "file1.wea", "1/2/2020", "3/2/2020", "--file", "file2.wea", "2/1/2021", "2/4/2021"});
  data = parseWeatherCatData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 2u);
  EXPECT_EQ(data.date_ranges.size(), 2u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_EQ(data.filenames[1], "file2.wea");
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{1, 1, 1}, Time{0, 0}}));
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{0, 0}}));  
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{3, 3, 2020}, Time{0, 0}}));  

  EXPECT_FALSE(data.date_ranges[1].is_in_range(DateTime{Date{31, 1, 2021}, Time{0, 0}}));
  EXPECT_TRUE(data.date_ranges[1].is_in_range(DateTime{Date{2, 2, 2021}, Time{0, 0}}));  
  EXPECT_FALSE(data.date_ranges[1].is_in_range(DateTime{Date{5, 2, 2021}, Time{0, 0}})); 

  args = CommandLineArguments({"outputfile", "--file", "file1.wea", "1/2/2020", "3/2/2020", "--file", "file2.wea", "2/1/2021"});
  EXPECT_ANY_THROW(parseWeatherCatData(args.getArgc(), args.getArgv()));  

  args = CommandLineArguments({"outputfile", "--file", "file1.wea", "1/2/2020-5:00", "3/2/2020-13:00"});
  data = parseWeatherCatData(args.getArgc(), args.getArgv());

  EXPECT_EQ(data.filenames.size(), 1u);
  EXPECT_EQ(data.date_ranges.size(), 1u);
  EXPECT_EQ(data.filenames[0], "file1.wea");
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{4, 0}}));
  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 1, 2020}, Time{5, 0}}));

  EXPECT_TRUE(data.date_ranges[0].is_in_range(DateTime{Date{2, 3, 2020}, Time{13, 0}}));  
  EXPECT_FALSE(data.date_ranges[0].is_in_range(DateTime{Date{2, 3, 2020}, Time{14, 0}}));  
}

TEST(WeatherCat, CatFiles)
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
  pt2.year = 2000;
  pt2.month = 2;
  pt2.day = 15;
  pt2.hour = 2;
  pt2.minute = 3 + 1;

  pt2.temperature = 10 + 1;
  pt2.horizontal_ir_radiation = 2 + 1;
  pt2.direct_normal_radiation = 3 + 1;
  pt2.diffuse_radiation = 4 + 1;
  pt2.wind_speed = 5 + 1;
  pt2.wind_direction = 6 + 1;  

  EPWDataPoint pt3;
  pt3.year = 2000;
  pt3.month = 2;
  pt3.day = 15;
  pt3.hour = 2;
  pt3.minute = 3 + 2;

  pt3.temperature = 10 + 2;
  pt3.horizontal_ir_radiation = 2 + 2;
  pt3.direct_normal_radiation = 3 + 2;
  pt3.diffuse_radiation = 4 + 2;
  pt3.wind_speed = 5 + 2;
  pt3.wind_direction = 6 + 2;

  EPWDataPoint pt4;
  pt4.year = 2000;
  pt4.month = 2;
  pt4.day = 15;
  pt4.hour = 2;
  pt4.minute = 3 + 3;

  pt4.temperature = 10 + 3;
  pt4.horizontal_ir_radiation = 2 + 3;
  pt4.direct_normal_radiation = 3 + 3;
  pt4.diffuse_radiation = 4 + 3;
  pt4.wind_speed = 5 + 3;
  pt4.wind_direction = 6 + 3;

  std::vector<EPWDataPoint> data1;
  data1.push_back(pt1);
  data1.push_back(pt2);

  std::vector<EPWDataPoint> data2;
  data2.push_back(pt3);
  data2.push_back(pt4);

  WeatherFileWriter writer1("tmp1.wea");
  writer1.write(data1);  

  WeatherFileWriter writer2("tmp2.wea");
  writer2.write(data2);

  CommandLineArguments args({"outputfile.wea", "--file", "tmp1.wea", "--file", "tmp2.wea"});
  WeatherCatParsedData parsedData = parseWeatherCatData(args.getArgc(), args.getArgv());
  WeatherCat cat(parsedData);
  cat.catFiles();

  WeatherFileReader reader("outputfile.wea");
  auto data = reader.read();

  compareEPW(data[0], pt1);
  compareEPW(data[1], pt2);
  compareEPW(data[2], pt3);
  compareEPW(data[3], pt4);
}


TEST(WeatherCat, CatFilesWithJump)
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
  pt2.year = 2000;
  pt2.month = 2;
  pt2.day = 15;
  pt2.hour = 2;
  pt2.minute = 3 + 1;

  pt2.temperature = 10 + 1;
  pt2.horizontal_ir_radiation = 2 + 1;
  pt2.direct_normal_radiation = 3 + 1;
  pt2.diffuse_radiation = 4 + 1;
  pt2.wind_speed = 5 + 1;
  pt2.wind_direction = 6 + 1;  

  EPWDataPoint pt3;
  pt3.year = 2000;
  pt3.month = 2;
  pt3.day = 15;
  pt3.hour = 2;
  pt3.minute = 3 + 10;

  pt3.temperature = 10 + 2;
  pt3.horizontal_ir_radiation = 2 + 2;
  pt3.direct_normal_radiation = 3 + 2;
  pt3.diffuse_radiation = 4 + 2;
  pt3.wind_speed = 5 + 2;
  pt3.wind_direction = 6 + 2;

  EPWDataPoint pt4;
  pt4.year = 2000;
  pt4.month = 2;
  pt4.day = 15;
  pt4.hour = 2;
  pt4.minute = 3 + 11;

  pt4.temperature = 10 + 3;
  pt4.horizontal_ir_radiation = 2 + 3;
  pt4.direct_normal_radiation = 3 + 3;
  pt4.diffuse_radiation = 4 + 3;
  pt4.wind_speed = 5 + 3;
  pt4.wind_direction = 6 + 3;

  std::vector<EPWDataPoint> data1;
  data1.push_back(pt1);
  data1.push_back(pt2);

  std::vector<EPWDataPoint> data2;
  data2.push_back(pt3);
  data2.push_back(pt4);

  WeatherFileWriter writer1("tmp1.wea");
  writer1.write(data1);  

  WeatherFileWriter writer2("tmp2.wea");
  writer2.write(data2);

  CommandLineArguments args({"outputfile.wea", "--file", "tmp1.wea", "--file", "tmp2.wea"});
  WeatherCatParsedData parsedData = parseWeatherCatData(args.getArgc(), args.getArgv());
  WeatherCat cat(parsedData);
  cat.catFiles();

  WeatherFileReader reader("outputfile.wea");
  auto data = reader.read();

  pt3.minute = 5;
  pt4.minute = 6;

  compareEPW(data[0], pt1);
  compareEPW(data[1], pt2);
  compareEPW(data[2], pt3);
  compareEPW(data[3], pt4);
}



TEST(WeatherCat, CatFilesWithDateRange)
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
  pt2.year = 2000;
  pt2.month = 2;
  pt2.day = 15;
  pt2.hour = 2;
  pt2.minute = 3 + 1;

  pt2.temperature = 10 + 1;
  pt2.horizontal_ir_radiation = 2 + 1;
  pt2.direct_normal_radiation = 3 + 1;
  pt2.diffuse_radiation = 4 + 1;
  pt2.wind_speed = 5 + 1;
  pt2.wind_direction = 6 + 1;  

  EPWDataPoint pt3;
  pt3.year = 2000;
  pt3.month = 2;
  pt3.day = 15;
  pt3.hour = 2;
  pt3.minute = 3 + 2;

  pt3.temperature = 10 + 2;
  pt3.horizontal_ir_radiation = 2 + 2;
  pt3.direct_normal_radiation = 3 + 2;
  pt3.diffuse_radiation = 4 + 2;
  pt3.wind_speed = 5 + 2;
  pt3.wind_direction = 6 + 2;

  EPWDataPoint pt4;
  pt4.year = 2000;
  pt4.month = 2;
  pt4.day = 15;
  pt4.hour = 2;
  pt4.minute = 3 + 3;

  pt4.temperature = 10 + 3;
  pt4.horizontal_ir_radiation = 2 + 3;
  pt4.direct_normal_radiation = 3 + 3;
  pt4.diffuse_radiation = 4 + 3;
  pt4.wind_speed = 5 + 3;
  pt4.wind_direction = 6 + 3;

  EPWDataPoint pt5;
  pt5.year = 2000;
  pt5.month = 2;
  pt5.day = 15;
  pt5.hour = 2;
  pt5.minute = 3 + 4;

  pt5.temperature = 10 + 3;
  pt5.horizontal_ir_radiation = 2 + 3;
  pt5.direct_normal_radiation = 3 + 3;
  pt5.diffuse_radiation = 4 + 3;
  pt5.wind_speed = 5 + 3;
  pt5.wind_direction = 6 + 3;

  EPWDataPoint pt6;
  pt6.year = 2000;
  pt6.month = 2;
  pt6.day = 15;
  pt6.hour = 2;
  pt6.minute = 3 + 5;

  pt6.temperature = 10 + 3;
  pt6.horizontal_ir_radiation = 2 + 3;
  pt6.direct_normal_radiation = 3 + 3;
  pt6.diffuse_radiation = 4 + 3;
  pt6.wind_speed = 5 + 3;
  pt6.wind_direction = 6 + 3;  

  std::vector<EPWDataPoint> data1;
  data1.push_back(pt1);
  data1.push_back(pt2);

  std::vector<EPWDataPoint> data2;
  data2.push_back(pt3);
  data2.push_back(pt4);
  data2.push_back(pt5);
  data2.push_back(pt6);

  WeatherFileWriter writer1("tmp1.wea");
  writer1.write(data1);  

  WeatherFileWriter writer2("tmp2.wea");
  writer2.write(data2);

  CommandLineArguments args({"outputfile.wea", "--file", "tmp1.wea", "--file", "tmp2.wea", "2/15/2000-2:06", "2/15/2000-2:07"});
  WeatherCatParsedData parsedData = parseWeatherCatData(args.getArgc(), args.getArgv());
  WeatherCat cat(parsedData);
  cat.catFiles();

  WeatherFileReader reader("outputfile.wea");
  auto data = reader.read();

  pt4.minute = 5;
  pt5.minute = 6;

  compareEPW(data[0], pt1);
  compareEPW(data[1], pt2);
  compareEPW(data[2], pt4);

  compareEPW(data[3], pt5);
}