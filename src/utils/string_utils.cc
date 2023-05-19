#include "utils/string_utils.h"

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