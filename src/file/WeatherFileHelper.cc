#include "file/WeatherFileHelper.h"

std::string WeatherFileHelper::getHeaders()
{
  return "year,month,day,hour,minute,temperature,horizontal_ir_radiation,direct_normal_radiation,diffuse_radiation,wind_direction,wind_speed";
}

std::string WeatherFileHelper::formatOutputPoint(const EPWDataPoint& pt)
{
  std::stringstream ss;
  ss << std::setprecision(16) << std::scientific;
  ss << pt.year << "," << pt.month << "," << pt.day << "," << pt.hour << "," << pt.minute;
  ss << "," << pt.temperature << "," << pt.horizontal_ir_radiation << ", " << pt.direct_normal_radiation;
  ss << "," << pt.diffuse_radiation << "," << pt.wind_direction << "," << pt.wind_speed;

  return ss.str();
}

EPWDataPoint WeatherFileHelper::parseDataLine(const std::string& line)
{
  std::vector<std::string> line_split = splitLine(line, ",");

  assertAlways(line_split.size() == 11, "Number of fields is incorrect");

  Parser parser;
  EPWDataPoint data;
  data.year   = parser.get<int>(line_split[0]);
  data.month  = parser.get<int>(line_split[1]);
  data.day    = parser.get<int>(line_split[2]);
  data.hour   = parser.get<int>(line_split[3]);
  data.minute = parser.get<int>(line_split[4]);

  data.temperature             = parser.get<Real>(line_split[5]);
  data.horizontal_ir_radiation = parser.get<Real>(line_split[6]);
  data.direct_normal_radiation = parser.get<Real>(line_split[7]);
  data.diffuse_radiation       = parser.get<Real>(line_split[8]);
  data.wind_direction          = parser.get<Real>(line_split[9]);
  data.wind_speed              = parser.get<Real>(line_split[10]);

  return data;
}