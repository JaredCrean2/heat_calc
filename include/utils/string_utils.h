#ifndef HEAT_CALC_UTILS_STRING_UTILS
#define HEAT_CALC_UTILS_STRING_UTILS


#include <string>
#include <vector>
#include <sstream>

#include <iostream>

std::vector<std::string> splitLine(const std::string& line, const std::string& delim);

// removes the leading and trailing whitespace from string.  Leaves any intermediate whitespace
std::string trimWhiteSpace(const std::string& str);

class Parser
{
  public:
    template<typename T>
    T get(const std::string& key)
    {
      T val;

      // use the machinery in std::stringstream to do the parsing
      ss.clear();
      ss.str(key);
      ss.seekg(0);
      ss >> val;

      char c;      
      if (ss.fail() || ss.get(c))
        throw std::runtime_error("failed to parse string");

      return val;
    }

  private:
    std::stringstream ss;
};


#endif