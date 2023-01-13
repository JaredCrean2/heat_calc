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

    void write(const std::vector<EPWDataPoint>& data)
    {
      WeatherFileHelper helper;
      m_file << helper.getHeaders() << "\n";

      for (auto& pt : data)
        m_file << helper.formatOutputPoint(pt) << "\n";

      m_file.close();
    }

  private:
    std::ofstream m_file;
};


#endif