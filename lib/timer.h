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

#include <stdint.h>
#if defined(__linux)
#define HAVE_POSIX_TIMER
#include <time.h>
#ifdef CLOCK_MONOTONIC
#define CLOCKID CLOCK_MONOTONIC
#else
#define CLOCKID CLOCK_REALTIME
#endif
#elif defined(__APPLE__)
#define HAVE_MACH_TIMER
#include <mach/mach_time.h>
#elif defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

static double _qedfv_ns()
{
  static uint64_t is_init = 0;
#if defined(__APPLE__)
  static mach_timebase_info_data_t info;
  if (0 == is_init)
  {
    mach_timebase_info(&info);
    is_init = 1;
  }
  uint64_t now;
  now = mach_absolute_time();
  now *= info.numer;
  now /= info.denom;
  return (double)now;
#elif defined(__linux)
  static struct timespec linux_rate;
  if (0 == is_init)
  {
    clock_getres(CLOCKID, &linux_rate);
    is_init = 1;
  }
  double now;
  struct timespec spec;
  clock_gettime(CLOCKID, &spec);
  now = spec.tv_sec * 1.0e9 + spec.tv_nsec;
  return now;
#elif defined(_WIN32)
  static LARGE_INTEGER win_frequency;
  if (0 == is_init)
  {
    QueryPerformanceFrequency(&win_frequency);
    is_init = 1;
  }
  LARGE_INTEGER now;
  QueryPerformanceCounter(&now);
  return (uint64_t)((1.0e9 * now.QuadPart) / win_frequency.QuadPart);
#endif
}

static double _qedfv_ms() { return _qedfv_ns() / 1.0e6; }