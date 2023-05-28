#include "gtest/gtest.h"
#include "file/input_parser.h"

TEST(InputParser, TrimWhitespace)
{
  EXPECT_EQ(trimWhiteSpace("abc"), "abc");
  EXPECT_EQ(trimWhiteSpace(" abc"), "abc");
  EXPECT_EQ(trimWhiteSpace("abc "), "abc");

  EXPECT_EQ(trimWhiteSpace("  abc"), "abc");
  EXPECT_EQ(trimWhiteSpace("abc  "), "abc");

  EXPECT_EQ(trimWhiteSpace("  abc def  "), "abc def");
}

TEST(InputParser, ParseFile)
{
  auto ss = std::make_shared<std::stringstream>();
  (*ss) << "key1 : val1" << std::endl;
  (*ss) << "  key2 : val2  " << std::endl;
  (*ss) << "  key3 : val3  #hello: " << std::endl;
  (*ss) << std::endl;
  (*ss) << std::endl;
  (*ss) << "# key4 : val4" << std::endl;

  std::map<std::string, std::string> defaults = { std::make_pair("key5", "val5")};
  InputParser parser(ss);

  std::map<std::string, std::string> inputs = parser.parse(defaults);

  EXPECT_EQ(inputs.size(), 4u);
  EXPECT_EQ(defaults.size(), 1u);
  EXPECT_EQ(inputs.at("key1"), "val1");
  EXPECT_EQ(inputs.at("key2"), "val2");
  EXPECT_EQ(inputs.at("key3"), "val3");
  EXPECT_EQ(inputs.at("key5"), "val5");
}

TEST(InputParser, ParseValues)
{
  ValueParser parser;
  EXPECT_EQ(parser.parseScalar<double>("2.0"), 2.0);
  EXPECT_EQ(parser.parseScalar<double>("2"), 2.0);

  std::vector<double> vals_ex{1, 2};
  EXPECT_EQ(parser.parseArray<double>("[1, 2]"), vals_ex);

  vals_ex = {1};
  EXPECT_EQ(parser.parseArray<double>("[1]"), vals_ex);
}