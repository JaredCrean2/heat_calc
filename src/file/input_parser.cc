#include "file/input_parser.h"
#include <cctype>
#include <algorithm>
#include <set>


std::map<std::string, std::string> InputParser::parse(std::map<std::string, std::string> defaults)
{
  if (!(*m_infile))
    throw std::runtime_error("unable to open file");

  std::string line;

  std::set<std::string> seen_keys;
  while (std::getline(*m_infile, line))
  {
    line = trimComment(line);
    if (isBlankLine(line))
      continue;

    std::pair<std::string, std::string> keyval = parseKeyAndValue(line);
    if (seen_keys.count(keyval.first) > 0)
      throw std::runtime_error(std::string("duplicate key found on line: ") + line);
    else
      seen_keys.insert(keyval.first);

    if (defaults.count(keyval.first) > 0)
      defaults[keyval.first] = keyval.second;
    else
      defaults.insert(keyval);
  }

  return defaults;
}


bool InputParser::isBlankLine(const std::string& line)
{
  return trimWhiteSpace(line).size() == 0;
}

std::string InputParser::trimComment(const std::string& line)
{
  auto it = std::find(line.begin(), line.end(), m_comment_char);
  if (it == line.end())
    return line;
  else
  {
    return std::string(line.begin(), it);
  }
}

std::pair<std::string, std::string> InputParser::parseKeyAndValue(const std::string& line)
{
  auto it = std::find(line.begin(), line.end(), m_separator);

  if (it == line.end())
    throw std::runtime_error(std::string("unable to parse line: ") + line +", line does not contain a separator");

  std::string key(line.begin(), it);
  std::advance(it, 1);
  std::string val(it, line.end());

  key = trimWhiteSpace(key);
  val = trimWhiteSpace(val);

  return std::make_pair(key, val);
}