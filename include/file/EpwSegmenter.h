#ifndef EPWSEGMENTER_H
#define EPWSEGMENTER_H

#include <memory>
#include <algorithm>

#include "EpwReader.h"
#include "physics/heat/dates.h"
#include "utils/error_handling.h"

#include <iostream>



class EpwSegmenter
{
  public:
    explicit EpwSegmenter(std::shared_ptr<EPWReader> epw_file) :
      m_epw_file(epw_file)
    {
      computeSegments();
    }

    explicit EpwSegmenter(const std::string& fname) :
      m_epw_file(std::make_shared<EPWReader>(fname))
    {
      computeSegments();
    }

    const std::vector<EPWDataPoint>& getData() const { return m_epw_data; }

    int getNumSegments() const { return m_segment_idxs.size(); }

    const std::vector<std::pair<int, int>>& getSegments() const { return m_segment_idxs; }

    std::vector<EPWDataPoint> getSegment(int idx) const;

  private:

    void computeSegments();

    void sortData();


    void findSegmentBoundaries();

    std::shared_ptr<EPWReader> m_epw_file;
    std::vector<EPWDataPoint> m_epw_data;
    std::vector<std::pair<int, int>> m_segment_idxs;
};

#endif
