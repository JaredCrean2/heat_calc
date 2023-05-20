#include "gtest/gtest.h"
#include "physics/heat/dates.h"

TEST(Dates, LeapYear)
{
  int idx=0;
  for (int year=1800; year <= 2200; ++year)
  {
    bool is_leap_year = isLeapYear(year);

    if ( (year % 100 == 0) && (year % 400 != 0) )
      EXPECT_FALSE(is_leap_year);
    else if (idx % 4 == 0)
      EXPECT_TRUE(is_leap_year);
    else
      EXPECT_FALSE(isLeapYear(year));
    idx++;
  }
}

TEST(Dates, DayOfYear)
{
  EXPECT_EQ(computeDayOfYear({1, 1, 1971}), 1);
  EXPECT_EQ(computeDayOfYear({5, 1, 1971}), 5);
  EXPECT_EQ(computeDayOfYear({1, 2, 1971}), 32);
  EXPECT_EQ(computeDayOfYear({1, 3, 1971}), 60);
  EXPECT_EQ(computeDayOfYear({31, 12, 1971}), 365);
}


TEST(Dates, computeJulianDate)
{
  int day0 = computeJulianDate({1, 1, 1971});
  
  {
    int day1 = computeJulianDate({2, 1, 1971});
    EXPECT_EQ(day1 - day0, 1);
  }

  {
    int day1 = computeJulianDate({1, 2, 1971});
    EXPECT_EQ(day1 - day0, 31);
  }

  {
    int day1 = computeJulianDate({31, 12, 1971});
    EXPECT_EQ(day1 - day0, 364);
  }
}


TEST(Dates, computeJulianDayExact)
{
  EXPECT_EQ(computeJulianDate({1, 1, 2000}), 2451545);
}


TEST(Dates, computeJulianDateLeapYear)
{
  int day0 = computeJulianDate({28, 2, 1988});
  int day1 = computeJulianDate({29, 2, 1988});
  int day2 = computeJulianDate({1, 3, 1988});
  EXPECT_EQ(day0 + 1, day1);
  EXPECT_EQ(day1 + 1, day2);
}

TEST(Dates, FractionalDayUniversalTime)
{
  Date date{1, 5, 2000};

  EXPECT_EQ(computeJulianDate(date, {12, 00}, 0), computeJulianDate(date));
  EXPECT_EQ(computeJulianDate(date, {13, 00}, 0), computeJulianDate(date) + 1.0/24);
  EXPECT_EQ(computeJulianDate(date, {11, 00}, 0), computeJulianDate(date) - 1.0/24);

  EXPECT_EQ(computeJulianDate(date, {13, 30}, 0), computeJulianDate(date) + 1.5/24);
  EXPECT_EQ(computeJulianDate(date, {10, 30}, 0), computeJulianDate(date) - 1.5/24);

  EXPECT_EQ(computeJulianDate(date, {24, 30}, 0), computeJulianDate(date) + 12.5/24);
}

TEST(Dates, FractionalDayTimeZone)
{
  Date date{1, 5, 2000};

  EXPECT_EQ(computeJulianDate(date, {12, 00}, 2), computeJulianDate(date) + 2);
  EXPECT_EQ(computeJulianDate(date, {13, 00}, 2), computeJulianDate(date) + 1.0/24 + 2);
  EXPECT_EQ(computeJulianDate(date, {11, 00}, 2), computeJulianDate(date) - 1.0/24 + 2);

  EXPECT_EQ(computeJulianDate(date, {13, 30}, 2), computeJulianDate(date) + 1.5/24 + 2);
  EXPECT_EQ(computeJulianDate(date, {10, 30}, 2), computeJulianDate(date) - 1.5/24 + 2);

  EXPECT_EQ(computeJulianDate(date, {24, 30}, 2), computeJulianDate(date) + 12.5/24 + 2);
}

TEST(Dates, CalendarToJulian)
{
  
  {
    DateTime datetime_ex{Date{1, 1, 2020}, Time{12, 00}};
    Real julian_date = computeJulianDate(datetime_ex.date, datetime_ex.time, 0);
    DateTime datetime = computeDateTime(julian_date, 0);
    EXPECT_EQ(datetime_ex, datetime);
  }

  {
    DateTime datetime_ex{Date{1, 1, 2020}, Time{23, 00}};
    Real julian_date = computeJulianDate(datetime_ex.date, datetime_ex.time, 0);
    DateTime datetime = computeDateTime(julian_date, 0);
    EXPECT_EQ(datetime_ex, datetime);
  } 

  {
    DateTime datetime_ex{Date{1, 1, 2020}, Time{1, 00}};
    Real julian_date = computeJulianDate(datetime_ex.date, datetime_ex.time, 0);
    DateTime datetime = computeDateTime(julian_date, 0);
    EXPECT_EQ(datetime_ex, datetime);
  }  

  {
    DateTime datetime_ex{Date{1, 1, 2020}, Time{23, 05}};
    Real julian_date = computeJulianDate(datetime_ex.date, datetime_ex.time, 0);
    DateTime datetime = computeDateTime(julian_date, 0);
    EXPECT_EQ(datetime_ex, datetime);
  }    
}
