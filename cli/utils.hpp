#pragma once

#include <OptParser.hpp>
#include <exception>
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

  for (unsigned int i = 0; i < 3; ++i)
  {
    std::getline(vstr, buf, ',');
    if (vstr.fail())
    {
      throw std::runtime_error("cannot parse vector '" + str + "'");
    }
    v[i] = strTo<double>(buf);
  }

  return v;
}

} // namespace OPT_PARSER_NS
