#include "gtest/gtest.h"
#include "EpwReader.h"

TEST(SplitLine, EmptyString)
{
  std::string line = "";
  auto words = splitLine(line, ",");

  EXPECT_EQ(words.size(), 0u);
}

TEST(SplitLine, NoDelim)
{
  std::string line = "abc";
  auto words = splitLine(line, ",");

  EXPECT_EQ(words.size(), 1u);
  EXPECT_EQ(words[0], "abc");
}

TEST(SplitLine, LeadingDelim)
{
  std::string line = ",abc";
  auto words = splitLine(line, ",");

  EXPECT_EQ(words.size(), 2u);
  EXPECT_EQ(words[0], "");
  EXPECT_EQ(words[1], "abc");
}

TEST(SplitLine, TwoLeadingDelims)
{
  std::string line = ",,abc";
  auto words = splitLine(line, ",");

  EXPECT_EQ(words.size(), 3u);
  EXPECT_EQ(words[0], "");
  EXPECT_EQ(words[1], "");
  EXPECT_EQ(words[2], "abc");
}

TEST(SplitLine, TwoWords)
{
  std::string line = "abc,def";
  auto words = splitLine(line, ",");

  EXPECT_EQ(words.size(), 2u);
  EXPECT_EQ(words[0], "abc");
  EXPECT_EQ(words[1], "def");
}

TEST(SplitLine, ThreeWords)
{
  std::string line = "abc,def,ghi";
  auto words = splitLine(line, ",");

  EXPECT_EQ(words.size(), 3u);
  EXPECT_EQ(words[0], "abc");
  EXPECT_EQ(words[1], "def");
  EXPECT_EQ(words[2], "ghi");
}

TEST(SplitLine, MiddleDoubleDelim)
{
  std::string line = "abc,,def";
  auto words = splitLine(line, ",");

  EXPECT_EQ(words.size(), 3u);
  EXPECT_EQ(words[0], "abc");
  EXPECT_EQ(words[1], "");
  EXPECT_EQ(words[2], "def");
}

TEST(SplitLine, TrailingDelim)
{
  std::string line = "abc,";
  auto words = splitLine(line, ",");

  EXPECT_EQ(words.size(), 2u);
  EXPECT_EQ(words[0], "abc");
  EXPECT_EQ(words[1], "");
}

TEST(SplitLine, TwoTrailingDelim)
{
  std::string line = "abc,,";
  auto words = splitLine(line, ",");

  EXPECT_EQ(words.size(), 3u);
  EXPECT_EQ(words[0], "abc");
  EXPECT_EQ(words[1], "");
  EXPECT_EQ(words[2], "");
}