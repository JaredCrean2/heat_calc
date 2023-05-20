#include "physics/heat/dates.h"
#include <iostream>
#include "utils/string_utils.h"

// note: month is in the range [1, 12]
int getDaysPerMonth(int month, int year)
{
  int leap_year_offset = isLeapYear(year) ? 1 : 0;
  std::array<int, 12> days_per_month{31, 28 + leap_year_offset, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};  

  return days_per_month[month-1];
}

Date incrementDate(const Date& date)
{

  if (date.day < getDaysPerMonth(date.month, date.year))
    return Date{date.day+1, date.month, date.year};
  else
  {
    Date date_output{1, date.month+1, date.year};
    if (date_output.month > 12)
    {
      date_output.year++;
      date_output.month = 1;
    }

    return date_output;
  }
}

Date decrementDate(const Date& date)
{
  if (date.day > 1)
  {
    return Date{date.day-1, date.month, date.year};
  } else
  {
    int year = date.year;
    int month = date.month - 1;
    if (month < 1)
    {
      month = 12;
      year--;
    }
    int day = getDaysPerMonth(month, year);
    return Date{day, month, year};
  }
}

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
  int dy = 0;
  for (int i=1; i < date.month; ++i)
    dy += getDaysPerMonth(i, date.year); // days_per_month[i];

  return dy + date.day;  
}

int computeJulianDate(const Date& date)
{
  int L = (date.month - 14) / 12;
  return date.day - 32075 + (1461 * (date.year + 4800 + L)) / 4
         + (367 * (date.month -2 - 12*L)) / 12
         - (3 * ( (date.year + 4900 + L) / 100 )) /4;
}

Date computeDateFromJulian(int julian_date)
{
  int J = julian_date + 68569;
  int N = 4*J/146097;
  int K = J - (146097*N + 3) / 4;
  int year_prime = 4000*(K+1)/1461001;
  int L = K - 1461*year_prime/4 + 31;
  int month_prime = 80*L/2447;
  int DM = L - 2447*month_prime/80;
  int M = month_prime/11;
  int MO = month_prime + 2 - 12*M;
  int year = 100*(N-49) + year_prime + M;

  return Date{DM, MO, year};
}


Real computeJulianDate(const Date& date, const Time& time, int time_zone)
{
  // time in universal time zone is time in local zone  + time_zone
  Real hours = time.hour + time.minute / 60.0;
  
  int julian_day_whole = computeJulianDate(date);
  return julian_day_whole + (hours - 12)/24.0 + time_zone;
}

DateTime computeDateTime(Real julian_date, int time_zone)
{
  int julian_day_whole = std::floor(julian_date);
  Date date = computeDateFromJulian(julian_day_whole);

  Real day_remainder = julian_date - julian_day_whole;
  Real hours_remainder = day_remainder * 24;
  int hours_since_noon = std::floor(hours_remainder);
  int minutes = std::round(60*(hours_remainder - hours_since_noon));
  int hours = hours_since_noon + 12;

  if (minutes >= 60)
  {
    hours += minutes / 60;
    minutes = minutes % 60;
  }

  if (hours > 23)
  {
    date = incrementDate(date);
    hours -= 24;
  }

  hours += time_zone;
  if (hours > 23)
  {
    date = incrementDate(date);
    hours -= 24;
  } else if (hours < 0)
  {
    date = decrementDate(date);
    hours += 24;
  }

  return DateTime{date, Time{hours, minutes}};
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
