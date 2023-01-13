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

    std::vector<EPWDataPoint> read();
      

  private:
    std::ifstream m_file;
    std::vector<EPWDataPoint> m_data;
};

#endif