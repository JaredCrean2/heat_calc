#ifndef PHYSICS_HEAT_DATE_H
#define PHYSICS_HEAT_DATE_H

#include <array>
#include "ProjectDefs.h"

#include <iostream>

// note: month is in the range [1, 12]
int getDaysPerMonth(int month, int year);

// Note: 1 based values
struct Date
{
  int day;
  int month;
  int year;
};

Date parseDate(const std::string& date);

inline bool operator==(const Date& lhs, const Date& rhs)
{
  return lhs.day   == rhs.day &&
         lhs.month == rhs.month &&
         lhs.year  == rhs.year;
}

inline bool operator!=(const Date& lhs, const Date& rhs)
{
  return !(lhs == rhs);
}

inline bool operator<(const Date& lhs, const Date& rhs)
{
  if (lhs.year != rhs.year)
    return lhs.year < rhs.year;
  else if (lhs.month != rhs.month)
    return lhs.month < rhs.month;
  else if (lhs.day != rhs.day)
    return lhs.day < rhs.day;
  else
    return false;
}

inline bool operator>(const Date& lhs, const Date& rhs)
{
  return rhs < lhs;
}

inline bool operator<=(const Date& lhs, const Date& rhs)
{
  return !(lhs > rhs);
}

inline bool operator>=(const Date& lhs, const Date& rhs)
{
  return !(lhs < rhs);
}

inline std::ostream& operator<<(std::ostream& os, const Date& date)
{
  os << date.month << "/" << date.day << "/" << date.year;
  return os;
}

void validateDate(const Date& date);

Date incrementDate(const Date& date);

Date decrementDate(const Date& date);


struct Time
{
  int hour;
  int minute;
};

inline bool operator==(const Time& lhs, const Time& rhs)
{
  return lhs.hour   == rhs.hour &&
         lhs.minute == rhs.minute;
}

inline bool operator!=(const Time& lhs, const Time& rhs)
{
  return !(lhs == rhs);
}

inline std::ostream& operator<<(std::ostream& os, const Time& date)
{
  std::string minute_str;
  if (date.minute < 10)
    minute_str = "0" + std::to_string(date.minute);
  else
    minute_str = std::to_string(date.minute);

  return os << date.hour << ":" << minute_str;
}

inline bool operator<(const Time& lhs, const Time& rhs)
{
  if (lhs.hour != rhs.hour)
    return lhs.hour < rhs.hour;
  else if (lhs.minute != rhs.minute)
    return lhs.minute < rhs.minute;
  else
    return false;
}

inline bool operator>(const Time& lhs, const Time& rhs)
{
  return rhs < lhs;
}

inline bool operator<=(const Time& lhs, const Time& rhs)
{
  return !(lhs > rhs);
}

inline bool operator>=(const Time& lhs, const Time& rhs)
{
  return !(lhs < rhs);
}


Time parseTime(const std::string time);

void validateTime(const Time& time);

bool isLeapYear(int year);

// computes number of days since the beginning of the year
int computeDayOfYear(const Date& date);

int computeJulianDate(const Date& date);

Date computeDateFromJulian(int julian_date);

// time_zone ex. Mountain time is GMT - 7, time zone is 7
Real computeJulianDate(const Date& date, const Time& time, int time_zone);


struct DateTime
{
  Date date;
  Time time;
};

inline bool operator==(const DateTime& lhs, const DateTime& rhs)
{
  return lhs.date == rhs.date &&
         lhs.time == rhs.time;
}

inline bool operator!=(const DateTime& lhs, const DateTime& rhs)
{
  return !(lhs == rhs);
}

inline std::ostream& operator<<(std::ostream& os, const DateTime& date)
{
  os << date.date << " " << date.time;

  return os;
}

inline bool operator<(const DateTime& lhs, const DateTime& rhs)
{
  if (lhs.date != rhs.date)
  {
    return lhs.date < rhs.date;
  } else
  {
    return lhs.time < rhs.time;
  }
}

inline bool operator>(const DateTime& lhs, const DateTime& rhs)
{
  return rhs < lhs;
}

inline bool operator<=(const DateTime& lhs, const DateTime& rhs)
{
  return !(lhs > rhs);
}

inline bool operator>=(const DateTime& lhs, const DateTime& rhs)
{
  return !(lhs < rhs);
}


// parses a string month/year/date[-[H]H:[M]M]
// time_default is used if the hour and minutes part of the string is absent
DateTime parseDateTime(const std::string& datetime, const Time& time_default);

DateTime computeDateTime(Real julian_date, int time_zone);

#endif