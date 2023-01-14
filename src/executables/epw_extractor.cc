#include "EpwExtractor.h"
#include "WeatherFileWriter.h"

int main(int argc, char* argv[])
{
  std::string input_fname, output_fname;
  if (argc == 3)
  {
    input_fname = argv[1];
    output_fname = argv[2];
  } else
  {
    std::cerr << "usage: " << argv[0] << " input_fname.epw output_fname.wea" << std::endl;
    return 1;
  }

  EpwExtractor extractor(input_fname);
  auto data = extractor.extract(std::cin);

  WeatherFileWriter writer(output_fname);
  writer.write(data);

  return 0;
}