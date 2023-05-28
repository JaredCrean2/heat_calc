#include "utils/string_utils.h"
#include <algorithm>

std::vector<std::string> splitLine(const std::string& line, const std::string& delim)
{
  std::vector<std::string> words;
  if (line.size() == 0)
    return words;

  std::string::size_type prev_pos = 0, pos = 0;
  while (pos < line.length())
  {
    prev_pos = pos;
    pos = line.find(delim, prev_pos);


    if (pos != std::string::npos)
    {
      words.push_back(line.substr(prev_pos, pos - prev_pos));
    } else
      break;


    pos = pos + 1;
  }

  if (pos == std::string::npos)
  {
    std::string::size_type start = prev_pos;
    if (line.substr(prev_pos, delim.size()) == delim)
      start++;
        
    words.push_back(line.substr(start));
  } else if (line.substr(pos-1) == delim)
  {
    words.push_back("");
  }

  return words;
}

// removes the leading and trailing whitespace from string.  Leaves any intermediate whitespace
std::string trimWhiteSpace(const std::string& str)
{
  auto is_not_whitespace = [&](unsigned char c) { return !std::isspace(c); };
  auto str_start = std::find_if(str.begin(), str.end(), is_not_whitespace);

  if (str_start == str.end())
    return "";

  auto str_end  = std::find_if(str.rbegin(), str.rend(), is_not_whitespace);
  size_t startidx = std::distance(str.begin(), str_start);
  size_t idxend = str.size() - std::distance(str.rbegin(), str_end) - 1;
  size_t len = idxend - startidx + 1;

  return str.substr(startidx, len);
}