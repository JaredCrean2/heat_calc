#ifndef HEAT_CALC_FILE_DATE_RANGE
#define HEAT_CALC_FILE_DATE_RANGE

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

#endif