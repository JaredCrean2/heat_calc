#ifndef EPW_EXTRACTOR_H
#define EPW_EXTRACTOR_H

#include "EpwSegmenter.h"

#include <fstream>
#include <ios>

// extracts segments of an EPW file
class EpwExtractor
{
  public:
    explicit EpwExtractor(std::shared_ptr<EpwSegmenter> epw_segments) :
      m_epw_segments(epw_segments)
    {}

    explicit EpwExtractor(const std::string& input_fname) :
      EpwExtractor(std::make_shared<EpwSegmenter>(input_fname))
    {}

    std::vector<EPWDataPoint> extract(std::istream& is=std::cin);

  private:
    void printSegments();

    std::vector<int> getOutputSegments(std::istream& is);

    std::array<char, 64> formatDateAndTime(const EPWDataPoint& pt);

    std::vector<EPWDataPoint> getOutputData(const std::vector<int>& segments);

    std::shared_ptr<EpwSegmenter> m_epw_segments;
};


#endif