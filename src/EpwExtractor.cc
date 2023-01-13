#include "EpwExtractor.h"
#include "WeatherFileHelper.h"


std::vector<EPWDataPoint> EpwExtractor::extract(std::istream& is)
{
  printSegments();
  auto segments = getOutputSegments(is);
  return getOutputData(segments);
}

void EpwExtractor::printSegments()
{
  std::cout << "EPW file has " << m_epw_segments->getNumSegments() << " segments.  They are" << std::endl;
  for (int i=0; i < m_epw_segments->getNumSegments(); ++i)
  {
    auto segment_idxs = m_epw_segments->getSegments()[i];
    EPWDataPoint start_pt = m_epw_segments->getData()[segment_idxs.first];
    EPWDataPoint end_pt = m_epw_segments->getData()[segment_idxs.second];

    std::cout << "segment " << i << ": " << formatDateAndTime(start_pt).data() << " to " << formatDateAndTime(end_pt).data() << std::endl;
  }
}

std::vector<int> EpwExtractor::getOutputSegments(std::istream& is)
{
  std::cout << "Select segments to extract to (comma separated values): ";
  std::string vals_string;
  if (is >> vals_string)
  {
    std::vector<std::string> vals_split = splitLine(vals_string, ",");
    Parser parser;
    std::vector<int> segments;
    for (size_t i=0; i < vals_split.size(); ++i)
    {
      segments.push_back(parser.get<int>(vals_split[i]));
    }

    std::cout << std::endl;

    is >> std::ws;
    return segments;
  } else
  {
    throw std::runtime_error("could not read segments to extract");
  }

}

std::array<char, 64> EpwExtractor::formatDateAndTime(const EPWDataPoint& pt)
{
  std::array<char, 64> buf;
  std::snprintf(buf.data(), buf.size(), "%2d/%2d/%4d %2d:%02d", pt.month, pt.day, pt.year, pt.hour, pt.minute);

  return buf;
}

std::vector<EPWDataPoint> EpwExtractor::getOutputData(const std::vector<int>& segments)
{
  if (segments.size() == 0)
    throw std::runtime_error("cannot output zero segments");

  std::vector<EPWDataPoint> data;
  for (auto& segment : segments)
  {
    auto segment_data = m_epw_segments->getSegment(segment);
    for (auto& data_pt : segment_data)
      data.push_back(data_pt);
  }

  return data;
}