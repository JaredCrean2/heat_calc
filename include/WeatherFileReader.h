#ifndef WEATHER_FILE_READER_H
#define WEATHER_FILE_READER_H

#include <string>
#include <fstream>
#include "EpwReader.h"
#include "WeatherFileHelper.h"

class WeatherFileReader
{
  public:
    explicit WeatherFileReader(const std::string& fname) :
      m_file(fname)
    {}

    std::vector<EPWDataPoint> read()
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
      

  private:
    std::ifstream m_file;
    std::vector<EPWDataPoint> m_data;
};

#endif