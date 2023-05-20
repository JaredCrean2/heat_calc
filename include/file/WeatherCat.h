#ifndef HEAT_CALC_FILE_WATHER_CAT
#define HEAT_CALC_FILE_WATHER_CAT

#include <string>
#include <vector>
#include "physics/heat/dates.h"

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

    bool is_in_range(const DateTime& datetime)
    {
      if (!m_have_dates)
        return true;
      else
        return datetime >= m_range_start && datetime <= m_range_end;
    }

  private:
    DateTime m_range_start;
    DateTime m_range_end;
    bool m_have_dates;
};


struct WeatherCatParsedData
{
  std::vector<std::string> filenames;
  std::vector<DateRange> date_ranges;
};



WeatherCatParsedData parseData(int argc, char* argv[]);

#endif