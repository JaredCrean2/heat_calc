#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cassert>
#include <array>
#include <ostream>

std::vector<std::string> splitLine(const std::string& line, const std::string& delim);


class Parser
{
  public:
    template<typename T>
    T get(const std::string& key)
    {
      T val;

      // use the machinery in std::stringstream to do the parsing
      ss.clear();
      ss.str(key);
      ss.seekg(0);
      ss >> val;

      return val;
    }

  private:
    std::stringstream ss;
};

struct EPWDataPoint
{
  int day;
  int month;
  int year;
  int hour;
  //int minute;
  double temperature;  // dry bulb tempeature
};

struct EPWLocation
{
  std::string local_name;
  std::string state;
  std::string country;
  double latitude;
  double longitude;
  int time_zone;
};

std::ostream& operator<<(std::ostream& os, const EPWLocation& location);


class EPWReader
{
  public:
    EPWReader(const std::string& fname) { readFile(fname);}

    using TData = EPWDataPoint;

    const std::vector<TData> get() const
    {
      return m_data;
    }

  private:

    void readFile(const std::string& fname);

    TData parseLine(std::string& line);

    std::string parseField(std::string& line, const std::string& delim);

    std::vector<TData> m_data;
    Parser m_parser;
};



