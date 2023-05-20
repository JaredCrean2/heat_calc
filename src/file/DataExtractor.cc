#include "file/DataExtractor.h"

DataExtractorParsedData parseDataExtractor(int argc, char* argv[])
{
  if (argc != 6)
  {
    std::cerr << "Usage: " << argv[0] << "data_filename weather_filename.wea date_start date_end output_filename" << std::endl;
    throw std::runtime_error("incorrect number of command line arguments");
  }

  DataExtractorParsedData data;
  data.data_filename = argv[1];
  data.weather_filename = argv[2];
  
  DateTime start_date = parseDateTime(argv[3], Time{00, 00});
  DateTime end_date = parseDateTime(argv[4], Time{23, 59});
  data.date_range = DateRange(start_date, end_date);

  data.output_filename = argv[5];

  return data;
}