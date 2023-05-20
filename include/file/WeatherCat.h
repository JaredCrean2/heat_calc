#ifndef HEAT_CALC_FILE_WATHER_CAT
#define HEAT_CALC_FILE_WATHER_CAT

#include <string>
#include <vector>

#include "file/WeatherFileReader.h"
#include "file/WeatherFileWriter.h"
#include "file/date_range.h"


struct WeatherCatParsedData
{
  std::string output_filename;
  std::vector<std::string> filenames;
  std::vector<DateRange> date_ranges;
};

WeatherCatParsedData parseWeatherCatData(int argc, char* argv[]);

class WeatherCat
{
  public:
    explicit WeatherCat(const WeatherCatParsedData& parsed_data) :
      m_parsed_data(parsed_data)
    {}

    void catFiles();

  private:
    void appendSegment(const std::string& fname, const DateRange& date_range);

    void filterData(const DateRange& date_range, std::vector<EPWDataPoint>& data);

    void rewriteDates(Real jd_start, std::vector<EPWDataPoint>& input_data);

    Real getJulianDate(const EPWDataPoint& data);

    void writeJulianDate(Real jd, EPWDataPoint& data);

    WeatherCatParsedData m_parsed_data;
    std::vector<EPWDataPoint> m_output_data;
};

#endif