#ifndef HEAT_CALC_FILE_WATHER_CAT
#define HEAT_CALC_FILE_WATHER_CAT

#include <string>
#include <vector>
#include "physics/heat/dates.h"

#include "file/WeatherFileReader.h"
#include "file/WeatherFileWriter.h"

class DateRange
{
  public:
    DateRange() :
      m_range_start{Date{1, 1, 1}, Time{0, 0}},
      m_range_end({Date{1, 1, 1}, Time{0, 0}}),
      m_have_dates(false)
    {}

    DateRange(const DateTime& range_start, const DateTime& range_end) :
      m_range_start(range_start),
      m_range_end(range_end),
      m_have_dates(true)
    {}

    bool is_in_range(const DateTime& datetime) const
    {
      if (!m_have_dates)
        return true;
      else
      {
        return datetime >= m_range_start && datetime <= m_range_end;
      }
    }

  private:
    DateTime m_range_start;
    DateTime m_range_end;
    bool m_have_dates;

  friend std::ostream& operator<<(std::ostream& os, const DateRange& range);
};

inline std::ostream& operator<<(std::ostream& os, const DateRange& range)
{
  if (range.m_have_dates)
    os << range.m_range_start << " to " << range.m_range_end;
  else
    os << "infinite date range";

  return os;

}

struct WeatherCatParsedData
{
  std::string output_filename;
  std::vector<std::string> filenames;
  std::vector<DateRange> date_ranges;
};

WeatherCatParsedData parseData(int argc, char* argv[]);

class WeatherCat
{
  public:
    explicit WeatherCat(const WeatherCatParsedData& parsed_data) :
      m_parsed_data(parsed_data)
    {}

    void catFiles();

  private:
    void appendSegment(const std::string& fname, const DateRange& date_range);

    void filterData(const DateRange& date_range, std::vector<EPWDataPoint>& data);

    void rewriteDates(Real jd_start, std::vector<EPWDataPoint>& input_data);

    Real getJulianDate(const EPWDataPoint& data);

    void writeJulianDate(Real jd, EPWDataPoint& data);

    WeatherCatParsedData m_parsed_data;
    std::vector<EPWDataPoint> m_output_data;
};

#endif