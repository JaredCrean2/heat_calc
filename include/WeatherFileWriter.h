#ifndef WEATHER_FILE_WRITER_H
#define WEATHER_FILE_WRITER_H

#include "EpwReader.h"
#include "WeatherFileHelper.h"
#include <string>
#include <fstream>

// Write EPW data points to the file format used by
// this code
class WeatherFileWriter
{
  public:
    explicit WeatherFileWriter(const std::string& fname) :
      m_file(fname, std::ios::out)
    {}

    void write(const std::vector<EPWDataPoint>& data);

  private:
    std::ofstream m_file;
};


#endif