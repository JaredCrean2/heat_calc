#ifndef PHYSICS_HEAT_DATE_H
#define PHYSICS_HEAT_DATE_H

#include <array>

struct Date
{
  int day;
  int month;
  int year;
};

bool isLeapYear(int year);

// computes number of days since the beginning of the year
int computeDayOfYear(const Date& date);

int computeJulianDate(const Date& date);

#endif