#include "file/WeatherCat.h"
#include "file/EpwReader.h"
#include "file/WeatherFileReader.h"
#include "file/WeatherFileWriter.h"
#include "physics/heat/dates.h"


WeatherCatParsedData parseData(int argc, char* argv[])
{
  WeatherCatParsedData data;
  int idx=1;
  if (idx >= argc)
    throw std::runtime_error("first command line argument must be output filename");

  data.output_filename = argv[idx];
  idx++;

  while (idx < argc)
  {
    if (std::string(argv[idx]) != "--file")
      throw std::runtime_error("expected --file at index " + std::to_string(idx));

    idx++;
    if (idx >= argc)
      throw std::runtime_error("reached end of arguments while parsing a --file section");

    data.filenames.emplace_back(argv[idx]);

    idx++;
    if (idx >= argc)
    {
      data.date_ranges.emplace_back();
      break;
    } else if (std::string(argv[idx]) != "--file")
    {
      DateTime range_start = parseDateTime(argv[idx], Time{0, 0});
      idx++;
      if (idx >= argc)
        throw std::runtime_error("found only one date in a --file section, must have two");

      DateTime range_end = parseDateTime(argv[idx], Time{23, 59});
      data.date_ranges.emplace_back(DateTime(range_start), DateTime(range_end));

      idx++;
    } else  // nex arg is --file: no date range is given for current file
    {
      data.date_ranges.emplace_back();
    }
  }

  if (data.filenames.size() == 0)
    throw std::runtime_error("command line arguments must contain at least 1 --file section");

  return data;
}
