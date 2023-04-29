/*
Copyright © 2023 Antonin Portelli <antonin.portelli@me.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

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

template <>
inline qedfv::QedFvCoef::Params strTo<qedfv::QedFvCoef::Params>(const std::string &str)
{
  qedfv::QedFvCoef::Params par;
  std::stringstream vstr(str);
  std::string buf;

  std::getline(vstr, buf, ',');
  if (vstr.fail())
  {
    throw std::runtime_error("cannot parse parameters '" + str + "'");
  }
  par.eta = strTo<double>(buf);
  std::getline(vstr, buf, ',');
  if (vstr.fail())
  {
    throw std::runtime_error("cannot parse parameters '" + str + "'");
  }
  par.nmax = strTo<unsigned int>(buf);

  return par;
}

} // namespace OPT_PARSER_NS
