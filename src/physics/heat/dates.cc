#include "physics/heat/dates.h"
#include <iostream>

bool isLeapYear(int year)
{
  if ( !(year % 4 == 0) || ( year % 100 == 0 && (year % 400 != 0)))
    return false;
  else
    return true;
}

// computes number of days since the beginning of the year
int computeDayOfYear(const Date& date)
{
  int leap_year_offset = isLeapYear(date.year) ? 1 : 0;
  std::array<int, 12> days_per_month{31, 28 + leap_year_offset, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  int dy = 0;
  for (int i=0; i < date.month-1; ++i)
    dy += days_per_month[i];

  return dy + date.day;  
}

int computeJulianDate(const Date& date)
{
  int L = (date.month - 14) / 12;
  return date.day - 32075 + (1461 * (date.year + 4800 + L)) / 4
         + (367 * (date.month -2 - 12*L)) / 12
         - (3 * ( (date.year + 4900 + L) / 100 )) /4;
}

Real computeJulianDate(const Date& date, const Time& time, int time_zone)
{
  // time in universal time zone is time in local zone  + time_zone
  Real hours = time.hour + time.minute / 60.0;
  
  int julian_day_whole = computeJulianDate(date);
  return julian_day_whole + (hours - 12)/24.0 + time_zone;
}
