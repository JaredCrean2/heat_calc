#include "WeatherFileReader.h"

std::vector<EPWDataPoint> WeatherFileReader::read()
{
  std::string line;
  std::getline(m_file, line);  // headers

  WeatherFileHelper helper;

  while (std::getline(m_file, line))
  {
    m_data.push_back(helper.parseDataLine(line));
  }

  return m_data;
}