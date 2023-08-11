#include "gtest/gtest.h"
#include "utils/string_utils.h"

TEST(Parser, Integer)
{
  Parser parser;

  EXPECT_EQ(parser.get<int>("1"), 1);
  EXPECT_EQ(parser.get<int>("5"), 5);
  EXPECT_EQ(parser.get<int>("-1"), -1);

  EXPECT_ANY_THROW(parser.get<int>("5a"));
}

TEST(Parser, Double)
{
  Parser parser;

  EXPECT_DOUBLE_EQ(parser.get<double>("1.0"), 1.0);
  EXPECT_DOUBLE_EQ(parser.get<double>("5.0"), 5.0);
  EXPECT_DOUBLE_EQ(parser.get<double>("-1.0"), -1.0);

  EXPECT_ANY_THROW(parser.get<double>("1.0a"));
}

TEST(Parser, Boolean)
{
  Parser parser;

  EXPECT_EQ(parser.get<bool>("true"), true);
  EXPECT_EQ(parser.get<bool>("false"), false);
  EXPECT_ANY_THROW(parser.get<double>("blah"));
}
