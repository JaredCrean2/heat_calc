#ifndef WEATHER_FILE_HELPER_H
#define WEATHER_FILE_HELPER_H

#include "ProjectDefs.h"
#include "EpwReader.h"
#include "error_handling.h"
#include <string>
#include <iomanip>


// Class for helping write the weather file format used by this code
// The file uses MKS unit system
class WeatherFileHelper
{
  public:

    std::string getHeaders();

    std::string formatOutputPoint(const EPWDataPoint& pt);

    EPWDataPoint parseDataLine(const std::string& line);
};

#endif