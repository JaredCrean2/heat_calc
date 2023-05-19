#ifndef PHYSICS_HEAT_DATE_H
#define PHYSICS_HEAT_DATE_H

#include <array>
#include "ProjectDefs.h"

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

inline std::ostream& operator<<(std::ostream& os, const Date& date)
{
  os << date.month << "/" << date.day << "/" << date.year;
  return os;
}

void validateDate(const Date& date);


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


Time parseTime(const std::string time);

void validateTime(const Time& time);

bool isLeapYear(int year);

// computes number of days since the beginning of the year
int computeDayOfYear(const Date& date);

int computeJulianDate(const Date& date);

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

// parses a string month/year/date[-[H]H:[M]M]
// time_default is used if the hour and minutes part of the string is absent
DateTime parseDateTime(const std::string& datetime, const Time& time_default);

#endif