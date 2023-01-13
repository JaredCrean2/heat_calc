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

struct Time
{
  int hour;
  int minute;
};

bool isLeapYear(int year);

// computes number of days since the beginning of the year
int computeDayOfYear(const Date& date);

int computeJulianDate(const Date& date);

// time_zone ex. Mountain time is GMT - 7, time zone is 7
Real computeJulianDate(const Date& date, const Time& time, int time_zone);

#endif