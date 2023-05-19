#ifndef ENVRIONMENT_INTERFACE_WEATHER_FILE_H
#define ENVRIONMENT_INTERFACE_WEATHER_FILE_H

#include "file/EpwReader.h"
#include "math.h"
#include "physics/heat/environment_interface.h"
#include "file/WeatherFileReader.h"

#include <iostream>

namespace Heat
{

class EnvironmentInterfaceWeatherFile : public EnvironmentInterface
{
  public:
    EnvironmentInterfaceWeatherFile(const std::vector<EPWDataPoint>& data) :
      m_data(data)
    {
      computeSpacing();
    }

    EnvironmentInterfaceWeatherFile(const std::string& fname) :
      EnvironmentInterfaceWeatherFile(WeatherFileReader(fname).read())
    {}

    
    // t is in seconds since the starting time
    EnvironmentData getEnvironmentData(Real t) override;

    Real getJulianDateStart() const { return getJulianDate(m_data[0]); }

  private:
    void computeSpacing();

    void checkUniformSpacing();

    int getIdx(Real t);

    EnvironmentData interpolateData(Real t);

    Real getTDays(Real t);
    
    std::vector<EPWDataPoint> m_data;
    double m_data_spacing;
};

}

#endif