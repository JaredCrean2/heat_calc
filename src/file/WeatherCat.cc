#include "file/WeatherCat.h"


WeatherCatParsedData parseData(int argc, char* argv[])
{
  WeatherCatParsedData data;
  int idx=1;
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

  return data;
}