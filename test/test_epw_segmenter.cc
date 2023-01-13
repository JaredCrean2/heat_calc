#include "gtest/gtest.h"
#include "EpwSegmenter.h"


TEST(EpwSegmenter, Example)
{
  EpwSegmenter segmenter("./data/epw_shortened.epw");

  EXPECT_EQ(segmenter.getNumSegments(), 2);
  EXPECT_EQ(segmenter.getSegments()[0].first, 0);
  EXPECT_EQ(segmenter.getSegments()[0].second, 2);

  EXPECT_EQ(segmenter.getSegments()[1].first, 3);
  EXPECT_EQ(segmenter.getSegments()[1].second, 7);

  auto data1 = segmenter.getSegment(0);
  EXPECT_EQ(data1.size(), 3u);
  EXPECT_EQ(data1.front().day, 1);
  EXPECT_EQ(data1.front().hour, 1);

  EXPECT_EQ(data1.back().day, 1);
  EXPECT_EQ(data1.back().hour, 3);  

  auto data2 = segmenter.getSegment(1);
  EXPECT_EQ(data2.size(), 5u);
  EXPECT_EQ(data2.front().day, 4);
  EXPECT_EQ(data2.front().hour, 3);

  EXPECT_EQ(data2.back().day, 4);
  EXPECT_EQ(data2.back().hour, 7);  
}