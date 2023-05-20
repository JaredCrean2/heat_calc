#include "gtest/gtest.h"
#include "arguments_setup.h"
#include "file/DataExtractor.h"
#include "file/EpwReader.h"
#include "file/WeatherFileWriter.h"

TEST(DataExtractor, File)
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

  std::vector<EPWDataPoint> data = {pt1, pt2};
  WeatherFileWriter writer("tmp.wea");
  writer.write(data);

  std::ofstream of("tmp.txt");
  of << "header line" << std::endl;
  of << "0 60" << std::endl;
  of << "1 120" << std::endl;
  of << "2 180" << std::endl;
  of << "3 240" << std::endl;
  of.close();

  CommandLineArguments args({"tmp.txt", "tmp.wea", "2/15/2000-2:05", "2/15/2000-2:06", "outfile.txt"});
  DataExtractorParsedData parsed_data = parseDataExtractor(args.getArgc(), args.getArgv());
  DataExtractor extractor(parsed_data);
  extractor.extract();

  std::ifstream infile("outfile.txt");
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(infile, line))
  {
    lines.push_back(line);
  }

  EXPECT_EQ(lines.size(), 3);
  EXPECT_EQ(lines[0], "header line");
  EXPECT_EQ(lines[1], "1 120");
  EXPECT_EQ(lines[2], "2 180");
}