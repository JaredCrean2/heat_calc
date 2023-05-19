#ifndef EPW_READER_H
#define EPW_READER_H

#include "ProjectDefs.h"
#include "physics/heat/environment_interface.h"
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cassert>
#include <array>
#include <ostream>
#include "utils/string_utils.h"


struct EPWDataPoint
{
  int day;
  int month;
  int year;
  int hour;  // hour ranges from 1 to 24
  int minute;
  //int minute;
  double temperature;  // dry bulb tempeature (Kelvin)
  double horizontal_ir_radiation; // W h /m^2
  double direct_normal_radiation;  // W h/m^2
  double diffuse_radiation;        // W h/m^2
  double wind_direction;           // degrees, North is 0 degrees, measured clockwise
  double wind_speed;               // m/s
};

Real getJulianDate(const EPWDataPoint& data);

Heat::EnvironmentData convertToEnvData(const EPWDataPoint& pt);

struct EPWLocation
{
  std::string local_name;
  std::string state;
  std::string country;
  double latitude;
  double longitude;
  int time_zone;
};

struct EPWDataPeriods
{
  int num_periods;
  int num_records_per_hour;
  // There is other data in the DataPeriods line, but it
  // doesn't seem to accurately describe the data in the file
};

std::ostream& operator<<(std::ostream& os, const EPWLocation& location);


class EPWReader
{
  public:
    EPWReader(const std::string& fname) :
      m_file(fname)
    {
      if (!m_file)
        throw std::runtime_error(std::string("could not open file ") + fname);

      readFile();
    }

    using TData = EPWDataPoint;

    const EPWLocation& getLocation() const { return m_location; }

    const EPWDataPeriods& getDataPeriods() const { return m_data_periods; }

    const std::vector<TData>& getData() const { return m_data; }

  private:
    void readFile();

    EPWLocation readLocation();

    void skipRows(int nrows);

    EPWDataPeriods readDataPeriods();

    void readData();

    void validateData(const EPWDataPoint& data);

    TData parseLine(std::string& line);

    std::string parseField(std::string& line, const std::string& delim);

    EPWDataPeriods m_data_periods;
    EPWLocation m_location;
    std::vector<TData> m_data;
    Parser m_parser;
    std::ifstream m_file;
    std::string m_delim = ",";
};

#endif

