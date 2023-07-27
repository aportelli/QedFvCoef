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
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>

namespace qedfv
{

// Debug printf ////////////////////////////////////////////////////////////////////////////////////
inline void dgbPrintf(const bool debug, const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  if (debug)
  {
    printf("[QedFv]: ");
    vprintf(fmt, args);
  }
  va_end(args);
}

// Floating-point equal operator ///////////////////////////////////////////////////////////////////
constexpr double cmpEps = 1.0e-10;

inline bool isEqual(double a, double b)
{
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * cmpEps);
}

// High-resolution timer ///////////////////////////////////////////////////////////////////////////
inline double clockMs(void)
{
  auto t = std::chrono::high_resolution_clock::now().time_since_epoch();
  auto d = std::chrono::duration_cast<std::chrono::nanoseconds>(t);
  return static_cast<double>(d.count()) * 1e-6;
}

// Simple square ///////////////////////////////////////////////////////////////////////////////////
inline double sq(const double x) { return x * x; }

} // namespace qedfv
