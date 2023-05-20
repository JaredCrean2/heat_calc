#include "physics/heat/dates.h"
#include "string_utils.h"
#include "file/WeatherCat.h"

// combines several .wea files into one.  This rewrites the date and time so that the combined
// file forms a single time range
// Usage: wea_cat --file fname1 [date_start date_end] --file fname2 [date_start date_end]
// The --file argument can be repeated as many times as necessary.
// The dates can be specified as month/year/day[-hour:minute].
// The month must be specified as a number.
// The hour:minute are optional.  If not specified, any time on the specified day will be
// included
int main(int argc, char* argv[])
{

  WeatherCatParsedData parsed_data = parseData(argc, argv);
  WeatherCat cat(parsed_data);
  cat.catFiles();

  return 0;
}