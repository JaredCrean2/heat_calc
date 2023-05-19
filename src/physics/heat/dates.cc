#include "physics/heat/dates.h"
#include <iostream>
#include "utils/string_utils.h"

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


Date parseDate(const std::string& date_str)
{
  std::vector<std::string> words = splitLine(date_str, "/");
  if (words.size() != 3)
    throw std::runtime_error("cannot parse " + date_str + " as a date, must be format month/day/year");

  Parser parser;
  int month = parser.get<int>(words[0]);
  int day   = parser.get<int>(words[1]);
  int year  = parser.get<int>(words[2]);

  Date date{day, month, year};

  validateDate(date);

  return date;
}

void validateDate(const Date& date)
{
  if (date.month < 1 || date.month > 12)
    throw std::runtime_error("month is out of range");

  if (date.day < 1)
    throw std::runtime_error("day is out of range");

  if (date.year < 0)
    throw std::runtime_error("year is out of range");
}

Time parseTime(const std::string time_str)
{
  std::vector<std::string> words = splitLine(time_str, ":");

  if (words.size() != 2)
    throw std::runtime_error("cannot parse " + time_str + " as a time, must be format HH:MM");

  Parser parser;
  int hour = parser.get<int>(words[0]);
  int minute = parser.get<int>(words[1]);

  Time time{hour, minute};

  validateTime(time);

  return time;
}

void validateTime(const Time& time)
{
  if (time.hour < 0 || time.hour > 23)
    throw std::runtime_error("hour is out of range");

  if (time.minute < 0 || time.minute > 59)
    throw std::runtime_error("minute is out of range");
}



DateTime parseDateTime(const std::string& datetime, const Time& time_default)
{
  size_t dash_pos         = datetime.find('-');
  bool have_time          = dash_pos != std::string::npos;
  std::string date_string = datetime.substr(0, dash_pos);

  Date date = parseDate(date_string);
  Time time;
  if (have_time)
  {
    if (datetime.size() <= dash_pos + 1)
      throw std::runtime_error("found the - that separates the date and time, but the time string is malformed");

    time = parseTime(datetime.substr(dash_pos+1, datetime.size()));
  } else
  {
    time = time_default;
  }

  return DateTime{date, time};
}
