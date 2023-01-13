#include "WeatherFileWriter.h"

void WeatherFileWriter::write(const std::vector<EPWDataPoint>& data)
{
  WeatherFileHelper helper;
  m_file << helper.getHeaders() << "\n";

  for (auto& pt : data)
    m_file << helper.formatOutputPoint(pt) << "\n";

  m_file.close();
}