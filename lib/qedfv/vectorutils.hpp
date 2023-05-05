/*
Copyright © 2023
Matteo Di Carlo <matteo.dicarlo93@gmail.com>
Antonin Portelli <antonin.portelli@me.com>

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
#include <array>
#include <cmath>
#include <qedfv/global.hpp>

namespace qedfv
{

// Vector and basic operations /////////////////////////////////////////////////////////////////////
template <typename T>
using Vec3 = std::array<T, 3>;

typedef Vec3<double> DVec3;
typedef Vec3<int> IVec3;

template <typename T, typename U>
inline auto dot(const Vec3<T> &v, const Vec3<U> &w) -> decltype(v[0] * w[0])
{
  return v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
}

template <typename T>
inline T norm2(const Vec3<T> &v)
{
  return dot(v, v);
}

inline DVec3 sphericalToCartesian(const double r, const double theta, const double phi)
{
  return {r * cos(phi) * sin(theta), r * sin(phi) * sin(theta), r * cos(theta)};
}

inline DVec3 cartesianToSpherical(const double x, const double y, const double z)
{
  return {sqrt(x * x + y * y + z * z), atan2(sqrt(x * x + y * y), z), atan2(y, x)};
}

} // namespace qedfv
