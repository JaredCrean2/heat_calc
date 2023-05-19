#include "gtest/gtest.h"
#include "physics/heat/dates.h"

TEST(Date, OperatorEquality)
{
  Date lhs{1, 2, 3}, rhs {1, 2, 3};
  EXPECT_TRUE(lhs == rhs);
  EXPECT_FALSE(lhs != rhs);
}

TEST(Date, OperatorInquality)
{
  Date lhs{1, 2, 3}, rhs{5, 2, 3};
  EXPECT_FALSE(lhs == rhs);
  EXPECT_TRUE(lhs != rhs);

  rhs = {1, 5, 3};
  EXPECT_FALSE(lhs == rhs);
  EXPECT_TRUE(lhs != rhs);

  rhs = {1, 2, 5};
  EXPECT_FALSE(lhs == rhs);
  EXPECT_TRUE(lhs != rhs);  
}

TEST(Date, ComparisonOperators)
{
  Date lhs{1, 2, 3}, rhs{1, 2, 4};
  EXPECT_TRUE(lhs <  rhs);   EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <= rhs);   EXPECT_FALSE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_TRUE(rhs > lhs);
  EXPECT_FALSE(lhs >= rhs);  EXPECT_TRUE(rhs >= lhs);

  rhs = {1, 4, 3};
  EXPECT_TRUE(lhs <  rhs);   EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <= rhs);   EXPECT_FALSE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_TRUE(rhs > lhs);
  EXPECT_FALSE(lhs >= rhs);  EXPECT_TRUE(rhs >= lhs);

  rhs = {4, 2, 3};
  EXPECT_TRUE(lhs <  rhs);   EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <= rhs);   EXPECT_FALSE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_TRUE(rhs > lhs);
  EXPECT_FALSE(lhs >= rhs);  EXPECT_TRUE(rhs >= lhs); 

  rhs = {1, 2, 3};
  EXPECT_FALSE(lhs <  rhs);  EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <=  rhs);  EXPECT_TRUE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_FALSE(rhs > lhs);
  EXPECT_TRUE(lhs >=  rhs);  EXPECT_TRUE(rhs >= lhs);
}

TEST(Date, Parsing)
{
  Date date_ex{1, 2, 2020};
  Date date = parseDate("2/1/2020");
  EXPECT_EQ(date, date_ex);

  EXPECT_ANY_THROW(parseDate("h/2/2020"));

  EXPECT_ANY_THROW(parseDate("13/2/2020"));
  EXPECT_ANY_THROW(parseDate("1/-1/2020"));
  EXPECT_ANY_THROW(parseDate("1/2/-1"));
}


TEST(Time, OperatorEquality)
{
  Time lhs{2, 0}, rhs{2, 0};
  EXPECT_TRUE(lhs == rhs);
}

TEST(Time, OperatorInequality)
{
  Time lhs{2, 0}, rhs{2, 1};
  EXPECT_FALSE(lhs == rhs);
  EXPECT_TRUE(lhs != rhs);

  lhs = {3, 0};
  rhs = {2, 0};
  EXPECT_FALSE(lhs == rhs);
  EXPECT_TRUE(lhs != rhs);  
}


TEST(Time, ComparisonOperators)
{
  Time lhs = {2, 15}, rhs = {2, 16};
  EXPECT_TRUE(lhs <  rhs);   EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <= rhs);   EXPECT_FALSE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_TRUE(rhs > lhs);
  EXPECT_FALSE(lhs >= rhs);  EXPECT_TRUE(rhs >= lhs); 

  rhs = {3, 15};
  EXPECT_TRUE(lhs <  rhs);   EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <= rhs);   EXPECT_FALSE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_TRUE(rhs > lhs);
  EXPECT_FALSE(lhs >= rhs);  EXPECT_TRUE(rhs >= lhs);

  rhs = {2, 15};
  EXPECT_FALSE(lhs <  rhs);  EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <=  rhs);  EXPECT_TRUE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_FALSE(rhs > lhs);
  EXPECT_TRUE(lhs >=  rhs);  EXPECT_TRUE(rhs >= lhs);   
}


TEST(Time, Parsing)
{
  Time time_ex{2, 0};
  Time time = parseTime("2:00");
  EXPECT_EQ(time, time_ex);

  EXPECT_ANY_THROW(parseTime("-2:00"));
  EXPECT_ANY_THROW(parseTime("24:00"));
  EXPECT_ANY_THROW(parseTime("2:-1"));
  EXPECT_ANY_THROW(parseTime("-2:60"));
}

TEST(DateTime, Parsing)
{
  DateTime datetime_ex{Date{2, 1, 2020}, Time{0, 0}};
  DateTime datetime = parseDateTime("1/2/2020", Time{0, 0});
  EXPECT_EQ(datetime, datetime_ex);

  datetime_ex = {Date{2, 1, 2020}, Time{2, 30}};
  datetime = parseDateTime("1/2/2020-2:30", Time{0, 0});
  EXPECT_EQ(datetime, datetime_ex);

  Time time_default{0, 0};
  EXPECT_ANY_THROW(parseDateTime("1/2/2020-2", time_default));
  EXPECT_ANY_THROW(parseDateTime("1/2/2020-2222", time_default));
}

TEST(DateTime, ComparisonOperators)
{
  DateTime lhs{Date{1, 2, 3}, Time{2, 15}};
  DateTime rhs{Date{1, 2, 3}, Time{2, 16}};
  EXPECT_TRUE(lhs <  rhs);   EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <= rhs);   EXPECT_FALSE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_TRUE(rhs > lhs);
  EXPECT_FALSE(lhs >= rhs);  EXPECT_TRUE(rhs >= lhs);

  lhs = {Date{1, 2, 3}, Time{2, 15}};
  rhs = {Date{1, 2, 4}, Time{2, 15}};
  EXPECT_TRUE(lhs <  rhs);   EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <= rhs);   EXPECT_FALSE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_TRUE(rhs > lhs);
  EXPECT_FALSE(lhs >= rhs);  EXPECT_TRUE(rhs >= lhs);

  lhs = {Date{1, 2, 3}, Time{2, 15}};
  rhs = {Date{1, 2, 3}, Time{2, 15}};
  EXPECT_FALSE(lhs <  rhs);  EXPECT_FALSE(rhs < lhs);
  EXPECT_TRUE(lhs <=  rhs);  EXPECT_TRUE(rhs <= lhs);
  EXPECT_FALSE(lhs >  rhs);  EXPECT_FALSE(rhs > lhs);
  EXPECT_TRUE(lhs >=  rhs);  EXPECT_TRUE(rhs >= lhs);
}  