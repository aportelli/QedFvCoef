#pragma once

#include <OptParser.hpp>
#include <qedfvcoef.hpp>
#include <sstream>
#include <string>

namespace OPT_PARSER_NS
{

template <>
inline qedfv::DVec3 strTo<qedfv::DVec3>(const std::string &str)
{
  qedfv::DVec3 v;
  std::stringstream vstr(str);
  std::string buf;

  std::getline(vstr, buf, ',');
  v[0] = strTo<double>(buf);
  std::getline(vstr, buf, ',');
  v[1] = strTo<double>(buf);
  std::getline(vstr, buf, ',');
  v[2] = strTo<double>(buf);

  return v;
}

} // namespace OPT_PARSER_NS
