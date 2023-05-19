#include "gtest/gtest.h"
#include "file/EpwExtractor.h"

TEST(EpwExtractor, Example)
{
  EpwExtractor extractor("./data/epw_shortened.epw");

  std::istringstream is("0,1\n");

  auto data = extractor.extract(is);
  EXPECT_EQ(data.size(), 8u);
}