#include "EpwReader.h"
#include <iostream>

void EPWReader::readFile(const std::string& fname)
{
  std::ifstream f(fname);
  std::string line;

  // 8 lines of header
  for (int i=0; i < 8; ++i)
    std::getline(f, line);


  int linenum = 0;
  while (std::getline(f, line))
  {
    std::cout << "parsing line " << linenum << std::endl;
    auto datapoint = parseLine(line);
    m_data.push_back(datapoint);

    // we require hourly data
    if (m_data.size() > 1 && datapoint.hour != 1)
      assert(datapoint.hour - m_data[m_data.size() - 2].hour == 1);

    linenum++;
  }
}


EPWReader::TData EPWReader::parseLine(std::string& line)
{
  int day, month, year, hour, minute;
  double temp;
  const char* delim = ",";

  auto year_s = parseField(line, delim);  year = m_parser.get<int>(year_s);
  auto month_s = parseField(line, delim); month = m_parser.get<int>(month_s);
  auto day_s = parseField(line, delim);   day = m_parser.get<int>(day_s);
  auto hour_s = parseField(line, delim);  hour = m_parser.get<int>(hour_s);
  parseField(line, delim); // time interval?
  parseField(line, delim); // uncertainty data
  auto temp_s = parseField(line, delim);  temp = m_parser.get<double>(temp_s);

  std::cout << "\nParsed line:" << std::endl;
  std::cout << "day = " << day << std::endl;
  std::cout << "month = " << month << std::endl;
  std::cout << "year = " << year << std::endl;
  std::cout << "hour = " << hour << std::endl;
  std::cout << "temp = " << temp << std::endl;
  return EPWDataPoint{day, month, year, hour, temp};
}


std::string EPWReader::parseField(std::string& line, const std::string& delim)
{
  std::string::size_type pos = line.find(delim);
  std::string field = line.substr(0, pos);
  line.erase(0, pos + delim.length());

  return field;
}
