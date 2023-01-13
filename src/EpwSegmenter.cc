#include "EpwSegmenter.h"

std::vector<EPWDataPoint> EpwSegmenter::getSegment(int idx) const
{
  auto range = m_segment_idxs.at(idx);
  auto range_start = std::next(m_epw_data.begin(), range.first);
  auto range_end = std::next(m_epw_data.begin(), range.second + 1);
  std::vector<EPWDataPoint> data(range_start, range_end);

  return data;
}

void EpwSegmenter::computeSegments()
{
  sortData();
  findSegmentBoundaries();
}

void EpwSegmenter::sortData()
{
  m_epw_data = m_epw_file->getData();

  auto comp = [](const EPWDataPoint& lhs, const EPWDataPoint& rhs)
  {
    return getJulianDate(lhs) < getJulianDate(rhs);
  };

  std::sort(m_epw_data.begin(), m_epw_data.end(), comp);
}


void EpwSegmenter::findSegmentBoundaries()
{
  const double tol = 1e-6;
  assertAlways(m_epw_data.size() > 1, "EPW file must have more than one line");
  Real expected_spacing = m_epw_file->getDataPeriods().num_records_per_hour/24.0;

  int idx_start = 0;
  for (size_t i=1; i < m_epw_data.size(); ++i)
  {
    double spacing = getJulianDate(m_epw_data[i]) - getJulianDate(m_epw_data[i-1]);
    if (std::abs(spacing - expected_spacing) > tol)
    {
      m_segment_idxs.push_back(std::make_pair(idx_start, i-1));
      idx_start = i;
    }
  }

  m_segment_idxs.push_back(std::make_pair(idx_start, m_epw_data.size()-1));
}