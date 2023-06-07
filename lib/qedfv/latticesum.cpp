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

#include <qedfv/latticesum.hpp>

using namespace qedfv;

double ThreadedSum::sum(Summand func, const unsigned int nmax, const bool debug)
{
  double sum = 0., time;
  int n0, n1, n2, inmax = nmax;
  unsigned int nit = (2 * inmax - 1) * (2 * inmax - 1) * (2 * inmax - 1);

  time = -clockMs();
#pragma omp parallel for collapse(3) reduction(+ : sum)
  for (n0 = -inmax; n0 <= inmax; ++(n0))
    for (n1 = -inmax; n1 <= inmax; ++(n1))
      for (n2 = -inmax; n2 <= inmax; ++(n2))
      {
        IVec3 n = {n0, n1, n2};
        sum += func(n);
      }
  time += clockMs();
  dgbPrintf(debug, "threaded sum, result= %.10e, nmax= %d, time= %f ms, perf= %.1e pt/s\n", sum,
            nmax, time, nit / (1.0e-3 * time));

  return sum;
}
