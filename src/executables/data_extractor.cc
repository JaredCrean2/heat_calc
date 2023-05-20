#include "file/DataExtractor.h"

int main(int argc, char* argv[])
{
  DataExtractorParsedData parsed_data = parseDataExtractor(argc, argv);
  DataExtractor extractor(parsed_data);
  extractor.extract();

  return 0;
}