#ifndef HEAT_CALC_FILE_WATHER_CAT
#define HEAT_CALC_FILE_WATHER_CAT

#include <string>
#include <vector>
#include "physics/heat/dates.h"

#include "file/WeatherFileReader.h"
#include "file/WeatherFileWriter.h"

class DateRange
{
  public:
    DateRange() :
      m_range_start{Date{1, 1, 1}, Time{0, 0}},
      m_range_end({Date{1, 1, 1}, Time{0, 0}}),
      m_have_dates(false)
    {}

    DateRange(const DateTime& range_start, const DateTime& range_end) :
      m_range_start(range_start),
      m_range_end(range_end),
      m_have_dates(true)
    {}

    bool is_in_range(const DateTime& datetime) const
    {
      if (!m_have_dates)
        return true;
      else
      {
        return datetime >= m_range_start && datetime <= m_range_end;
      }
    }

  private:
    DateTime m_range_start;
    DateTime m_range_end;
    bool m_have_dates;

  friend std::ostream& operator<<(std::ostream& os, const DateRange& range);
};

inline std::ostream& operator<<(std::ostream& os, const DateRange& range)
{
  if (range.m_have_dates)
    os << range.m_range_start << " to " << range.m_range_end;
  else
    os << "infinite date range";

  return os;

}

struct WeatherCatParsedData
{
  std::string output_filename;
  std::vector<std::string> filenames;
  std::vector<DateRange> date_ranges;
};

WeatherCatParsedData parseData(int argc, char* argv[]);

class WeatherCat
{
  public:
    explicit WeatherCat(const WeatherCatParsedData& parsed_data) :
      m_parsed_data(parsed_data)
    {}

    void catFiles()
    {
      for (size_t i=0; i < m_parsed_data.filenames.size(); ++i)
      {
        if (i != 0)
          std::cout << std::endl;

        std::cout << "processing file " << m_parsed_data.filenames[i] << std::endl;
        appendSegment(m_parsed_data.filenames[i], m_parsed_data.date_ranges[i]);
      }

      WeatherFileWriter writer(m_parsed_data.output_filename);
      writer.write(m_output_data);
    }

  private:
    void appendSegment(const std::string& fname, const DateRange& date_range)
    {
      WeatherFileReader reader(fname);
      std::vector<EPWDataPoint> input_data = reader.read();
      filterData(date_range, input_data);

      if (input_data.size() == 0)
      {
        std::cerr << "Warning: the date range for file " << fname << " selects no data" << std::endl;
        return;
      }

      if (m_output_data.size() >= 2)
      {
        Real jd0 = getJulianDate(m_output_data[m_output_data.size()-2]);
        Real jd1 = getJulianDate(m_output_data[m_output_data.size()-1]);
        Real jd_start = jd1 + (jd1 - jd0);
        
        std::cout << "selected data starts at " << computeDateTime(getJulianDate(input_data[0]), 0) 
                  << ", rewriting to " << computeDateTime(jd_start, 0) << std::endl;
        rewriteDates(jd_start, input_data);
      }

      for (auto& data : input_data)
        m_output_data.push_back(data);
    }

    void filterData(const DateRange& date_range, std::vector<EPWDataPoint>& data)
    {
      int start_idx=-1, end_idx=-1;
      for (size_t i=0; i < data.size(); ++i)
      {
        Real jd = getJulianDate(data[i]);
        DateTime datetime = computeDateTime(jd, 0);

        if (start_idx == -1 && date_range.is_in_range(datetime))
        {
          start_idx = i;
        }

        if (start_idx != -1 && end_idx == -1 && !date_range.is_in_range(datetime))
        {
          end_idx = i;
        }
      }

      if (end_idx == -1)
        end_idx = data.size()-1;
       
      data.resize(end_idx+1);
      data.erase(data.begin(), data.begin() += start_idx);
    }

    void rewriteDates(Real jd_start, std::vector<EPWDataPoint>& input_data)
    {
      auto& data = input_data[0];
      Date date{data.day, data.month, data.year};
      Time time{data.hour, data.minute};
      Real jd_orig_prev = computeJulianDate(date, time, 0); 

      Real jd_new_i = jd_start;
      writeJulianDate(jd_new_i, data);

      // write jd_new_i to data
      for (size_t i=1; i < input_data.size(); ++i)
      {
        Real jd_orig_i   = getJulianDate(input_data[i]);
        Real jd_new_prev = getJulianDate(input_data[i-1]);
        Real jd_new_i    = jd_orig_i - jd_orig_prev + jd_new_prev;

        writeJulianDate(jd_new_i, input_data[i]);

        jd_orig_prev = jd_orig_i;
      }
    }

    Real getJulianDate(const EPWDataPoint& data)
    {
      Date date{data.day, data.month, data.year};
      Time time{data.hour, data.minute};
      return computeJulianDate(date, time, 0);
    }

    void writeJulianDate(Real jd, EPWDataPoint& data)
    {
      DateTime datetime = computeDateTime(jd, 0);

      data.day = datetime.date.day;
      data.month = datetime.date.month;
      data.year = datetime.date.year;
      data.hour = datetime.time.hour;
      data.minute = datetime.time.minute;
    }

    WeatherCatParsedData m_parsed_data;
    std::vector<EPWDataPoint> m_output_data;
};

#endif