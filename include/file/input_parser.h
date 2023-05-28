#ifndef HEAT_CALC_FILE_INPUT_PARSER_H
#define HEAT_CALC_FILE_INPUT_PARSER_H

#include <memory>
#include <string>
#include <map>
#include <fstream>
#include "utils/string_utils.h"


// parses input files where each line is of the form
// key : value
// keys are always strings.  The values are initially parsed into string, but
// can then be parsed into other concrete datatypes.  Arrays can be parsed if
// the string is of the form "[ val1, val2 ]"
// Comments start with a # character.  The remainder of the line is ignored.
// Blank lines are allowed
class InputParser
{
  public:
    InputParser(std::shared_ptr<std::istream> infile) :
      m_infile(infile)
    {}

    InputParser(std::string& fname) :
      m_infile(std::make_shared<std::ifstream>(fname))
    {}

    std::map<std::string, std::string> parse(std::map<std::string, std::string> defaults = {});

  private:

    bool isBlankLine(const std::string& line);

    std::string trimComment(const std::string& line);

    std::pair<std::string, std::string> parseKeyAndValue(const std::string& line);

    std::shared_ptr<std::istream> m_infile;
    const char m_comment_char = '#';
    const char m_separator = ':';
};

class ValueParser
{
  public:

    template <typename T>
    T parseScalarValue(const std::string& val)
    {
      return m_parser.get<T>(val);
    }

    template <typename T>
    std::vector<T> parseArrayValue(const std::string& val)
    {
      if (val.front() != '[' || val.back() != ']')
        throw std::runtime_error(std::string("unable to parse string as array: ") + val);

      auto start = val.begin(); std::advance(start, 1);
      auto end = val.end(); std::advance(end, -1);
      
      std::string values(start, end);
      std::vector<std::string> words = splitLine(values, ",");

      std::vector<T> vals;
      for (const auto& word : words)
        vals.push_back(m_parser.get<T>(trimWhiteSpace(word)));
        
      return vals;
    }

  private:
    Parser m_parser;
};

#endif