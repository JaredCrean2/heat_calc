#ifndef HEAT_CALC_FILE_DATA_EXTRACTOR
#define HEAT_CALC_FILE_DATA_EXTRACTOR

#include "file/EpwReader.h"
#include "file/date_range.h"
#include "file/WeatherFileReader.h"
#include "string_utils.h"


// extracts a section of a data file writtern by PostProcessorManager
// Usage: data_extractor simple_house_data.txt weather.wea date_start date_end
// The weather file is needed because it contains absolute calendar dates, while the
// data file contains dates relative to the start of the weather file
// The start_data and end_date of the format month/day/year-[[H]H:[M]M]

struct DataExtractorParsedData
{
  std::string data_filename;
  std::string weather_filename;
  DateRange date_range;
  std::string output_filename;
};

DataExtractorParsedData parseDataExtractor(int argc, char* argv[]);

class DataExtractor
{
  public:
    explicit DataExtractor(const DataExtractorParsedData& parsed_data) :
      m_parsed_data(parsed_data)
    {}

    void extract()
    {
      get_start_date();
      extract_data();
    }

  private:

    void get_start_date()
    {
      WeatherFileReader reader(m_parsed_data.weather_filename);
      EPWDataPoint data0 = reader.read()[0];

      Date date{data0.day, data0.month, data0.year};
      Time time{data0.hour, data0.minute};
      m_date_at_time_zero = DateTime{date, time};

      std::cout << "date at time 0 = " << m_date_at_time_zero << std::endl;
    }

    void extract_data()
    {
      std::ofstream of(m_parsed_data.output_filename);
      std::ifstream infile(m_parsed_data.data_filename);
      if (!infile)
        throw std::runtime_error(std::string("cannot open input data file ") + m_parsed_data.data_filename);

      std::string headerline, line;
      Parser parser;

      std::getline(infile, headerline);
      of << headerline << "\n";

      Real julian_date_start = computeJulianDate(m_date_at_time_zero.date, m_date_at_time_zero.time, 0);
      while (std::getline(infile, line))
      {
        std::vector<std::string> words = splitLine(line, " ");
        double time_seconds = parser.get<double>(words[1]);
        Real julian_date_line = julian_date_start + time_seconds/(60*60*24);
        DateTime date_line = computeDateTime(julian_date_line, 0);
        //std::cout << "date_line = " << date_line << std::endl;

        if (m_parsed_data.date_range.is_in_range(date_line))
        {
          of << line << "\n";
        }
      }
    }


    DataExtractorParsedData m_parsed_data;
    DateTime m_date_at_time_zero;
};

#endif