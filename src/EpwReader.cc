#include "EpwReader.h"
#include <iostream>
#include <stdexcept>
#include "physics/heat/dates.h"
#include "utils/math.h"


Real getJulianDate(const EPWDataPoint& data)
{
  Date date{data.day, data.month, data.year};
  Time time{data.hour, data.minute};

  // for comparison purposes, the time zone doesn't matter
  return computeJulianDate(date, time, 0);
}

Heat::EnvironmentData convertToEnvData(const EPWDataPoint& pt)
{
  Heat::EnvironmentData data;
  data.air_temp = pt.temperature;
  data.air_speed = pt.wind_speed;
  data.ir_horizontal_radiation = pt.horizontal_ir_radiation;
  data.direct_normal_radiation = pt.direct_normal_radiation;
  data.diffuse_radiation = pt.diffuse_radiation;

  Real theta = degreesToRadians(pt.wind_direction);
  // theta is measured clockwise from north
  data.air_direction = {std::sin(theta), std::cos(theta), 0};

  return data;
}

std::vector<std::string> splitLine(const std::string& line, const std::string& delim)
{
  std::vector<std::string> words;
  if (line.size() == 0)
    return words;

  std::string::size_type prev_pos = 0, pos = 0;
  while (pos < line.length())
  {
    prev_pos = pos;
    pos = line.find(delim, prev_pos);


    if (pos != std::string::npos)
    {
      words.push_back(line.substr(prev_pos, pos - prev_pos));
    } else
      break;


    pos = pos + 1;
  }

  if (pos == std::string::npos)
  {
    std::string::size_type start = prev_pos;
    if (line.substr(prev_pos, delim.size()) == delim)
      start++;
        
    words.push_back(line.substr(start));
  } else if (line.substr(pos-1) == delim)
  {
    words.push_back("");
  }

  return words;
}

std::ostream& operator<<(std::ostream& os, const EPWLocation& location)
{
  constexpr int buf_size = 64;
  std::array<char, buf_size> longitude, latitude;

  char latitude_n_or_s = location.latitude >= 0 ? 'N' : 'S';
  std::snprintf(latitude.data(), buf_size, "%3.2f degrees %c", location.latitude, latitude_n_or_s);

  char longitude_e_or_w = location.longitude >= 0 ? 'E' : 'W';
  std::snprintf(latitude.data(), buf_size, "%3.2f degrees %c", location.longitude, longitude_e_or_w);

  char time_zone_sign = location.time_zone >= 0 ? '+' : '-';
  os << location.local_name << ", " << location.state << ", " << location.country
      << ", " << latitude.data() << ", " << longitude.data() 
      << ", time zone GMT " << time_zone_sign << " " << std::abs(location.time_zone);

  return os;
}

void EPWReader::readFile()
{
  m_location = readLocation();

  skipRows(6);

  m_data_periods = readDataPeriods();

  readData();

  m_file.close();
}

EPWLocation EPWReader::readLocation()
{
  std::string line;
  std::getline(m_file, line);
  auto words = splitLine(line, m_delim);

  if (words.size() != 10)
    throw std::runtime_error("failed to parse LOCATION line");

  if (words[0] != "LOCATION")
    throw std::runtime_error("failed to parse LOCATION line");

  EPWLocation location;
  location.local_name = words[1];
  location.state      = words[2];
  location.country    = words[3];
  location.latitude   = m_parser.get<double>(words[6]);
  location.longitude  = m_parser.get<double>(words[7]);
  location.time_zone  = m_parser.get<int>(words[8]);

  return location;
}

void EPWReader::skipRows(int nrows)
{
  std::string line;
  for (int i=0; i < nrows; ++i)
    std::getline(m_file, line);
}

EPWDataPeriods EPWReader::readDataPeriods()
{
  std::string line;
  std::getline(m_file, line);

  auto words = splitLine(line, m_delim);
  if (words.size() < 6)
    throw std::runtime_error("failed to parse DATA PERIODS line");

  if (words[0] != "DATA PERIODS")
    throw std::runtime_error("failed to parse DATA PERIODS line");

  EPWDataPeriods data_periods;
  data_periods.num_periods = m_parser.get<int>(words[1]);
  data_periods.num_records_per_hour = m_parser.get<int>(words[2]);

  if (data_periods.num_periods != 1)
    throw std::runtime_error("more than one data period not supported");

  return data_periods;
}


void EPWReader::readData()
{
  std::string line;

  while (std::getline(m_file, line))
  {
    auto datapoint = parseLine(line);
    m_data.push_back(datapoint);


  }
}


EPWReader::TData EPWReader::parseLine(std::string& line)
{
  auto fields = splitLine(line, m_delim);
  //auto year_s = parseField(line, delim);  year = m_parser.get<int>(year_s);
  //auto month_s = parseField(line, delim); month = m_parser.get<int>(month_s);
  //auto day_s = parseField(line, delim);   day = m_parser.get<int>(day_s);
  //auto hour_s = parseField(line, delim);  hour = m_parser.get<int>(hour_s);
  //parseField(line, delim); // time interval?
  //parseField(line, delim); // uncertainty data
  //auto temp_s = parseField(line, delim);  temp = m_parser.get<double>(temp_s);

  EPWDataPoint data;
  data.year   = m_parser.get<int>(fields[0]);
  data.month  = m_parser.get<int>(fields[1]);
  data.day    = m_parser.get<int>(fields[2]);
  data.hour   = m_parser.get<int>(fields[3]);
  data.minute = m_parser.get<int>(fields[4]);
  //parseField(line, delim); // time interval?
  //parseField(line, delim); // uncertainty data
  data.temperature = m_parser.get<double>(fields[6]);
  data.horizontal_ir_radiation = m_parser.get<double>(fields[12]);
  data.direct_normal_radiation = m_parser.get<double>(fields[14]);
  data.diffuse_radiation = m_parser.get<double>(fields[15]);
  data.wind_direction = m_parser.get<double>(fields[20]);
  data.wind_speed = m_parser.get<double>(fields[21]);

  validateData(data);

  data.temperature = data.temperature + 273.15;

  return data;
}

void EPWReader::validateData(const EPWDataPoint& data)
{
  if (data.year < 1500)
    throw std::runtime_error("parsed year is too small");

  if (data.month < 1 || data.month > 13)
    throw std::runtime_error("parsed month is out of range");

  if (data.day < 1 || data.day > 31)
    throw std::runtime_error("parsed day is out of range");

  if (data.hour < 1 || data.hour > 24)
    throw std::runtime_error("parsed hour is out of range");

  if (data.minute < 0 || data.minute > 60)
    throw std::runtime_error("parsed minute is out of range"); 

  if (data.temperature < -70 || data.temperature > 70)
    throw std::runtime_error("parsed temperature is out of range or is missing"); 

  if (data.horizontal_ir_radiation < 0 || data.horizontal_ir_radiation == 9999)
    throw std::runtime_error("parsed horizontal IR radiation is out of range or is missing"); 

  if (data.direct_normal_radiation < 0 || data.direct_normal_radiation == 9999)
    throw std::runtime_error("parsed direct normal radiation is out of range or is missing"); 

  if (data.diffuse_radiation < 0 || data.diffuse_radiation == 9999)
    throw std::runtime_error("parsed diffuse radiation is out of range or is missing"); 

  if (data.wind_direction < 0 || data.wind_direction > 360)
    throw std::runtime_error("parsed wind direction is out of range or is missing"); 

  if (data.wind_speed < 0 || data.wind_speed > 40)
    throw std::runtime_error("parsed wind speed is out of range or is missing");     
}



std::string EPWReader::parseField(std::string& line, const std::string& delim)
{
  std::string::size_type pos = line.find(delim);
  std::string field = line.substr(0, pos);
  line.erase(0, pos + delim.length());

  return field;
}
