#include "file/WeatherCat.h"
#include "file/EpwReader.h"
#include "file/WeatherFileReader.h"
#include "file/WeatherFileWriter.h"
#include "physics/heat/dates.h"


WeatherCatParsedData parseWeatherCatData(int argc, char* argv[])
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


void WeatherCat::catFiles()
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

void WeatherCat::appendSegment(const std::string& fname, const DateRange& date_range)
{
  std::cout << std::setprecision(16);
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
              << " and ends on " << computeDateTime(getJulianDate(input_data.back()), 0) << std::endl;
    std::cout << "existing data end date = " << computeDateTime(jd1, 0) << std::endl;
    std::cout << "jd = " << jd1 << std::endl;
    std::cout << "New start date " << computeDateTime(jd_start, 0) << std::endl;
    rewriteDates(jd_start, input_data);
  }

  for (auto& data : input_data)
    m_output_data.push_back(data);
}

void WeatherCat::filterData(const DateRange& date_range, std::vector<EPWDataPoint>& data)
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

void WeatherCat::rewriteDates(Real jd_start, std::vector<EPWDataPoint>& input_data)
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

Real WeatherCat::getJulianDate(const EPWDataPoint& data)
{
  Date date{data.day, data.month, data.year};
  Time time{data.hour, data.minute};
  return computeJulianDate(date, time, 0);
}

void WeatherCat::writeJulianDate(Real jd, EPWDataPoint& data)
{
  DateTime datetime = computeDateTime(jd, 0);

  data.day = datetime.date.day;
  data.month = datetime.date.month;
  data.year = datetime.date.year;
  data.hour = datetime.time.hour;
  data.minute = datetime.time.minute;
}
