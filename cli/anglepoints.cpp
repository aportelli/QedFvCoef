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

#include "utils.hpp"
#include <OptParser.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <qedfv/latticesum.hpp>

using namespace optp;
using namespace std;
using namespace qedfv;

struct Point
{
  double theta, phi, n;
};

int main(int argc, const char *argv[])
{
  OptParser opt;
  bool parsed, debug;
  unsigned int nPoints;
  string filename;
  double radius, theta;

  opt.addOption("n", "npoints", OptParser::OptType::value, false,
                "number of points for each angle");
  opt.addOption("o", "output", OptParser::OptType::value, false, "output file");
  opt.addOption("e", "error", OptParser::OptType::value, true, "target relative error",
                strFrom(QEDFV_DEFAULT_ERROR));
  opt.addOption("d", "debug", OptParser::OptType::trigger, true, "show debug messages");
  opt.addOption("", "help", OptParser::OptType::trigger, true, "show this help message and exit");
  parsed = opt.parse(argc, argv);
  if (!parsed or (opt.getArgs().size() != 2) or opt.gotOption("help"))
  {
    cerr << "usage: " << argv[0] << " <radius> <theta>" << endl;
    cerr << endl << "Possible options:" << endl << opt << endl;

    return EXIT_FAILURE;
  }
  debug = opt.gotOption("d");
  nPoints = opt.optionValue<unsigned int>("n");
  filename = opt.optionValue("o");
  radius = strTo<double>(opt.getArgs()[0]);
  theta = strTo<double>(opt.getArgs()[1]);

  double dtheta = M_PI / nPoints, dphi = 2. * M_PI / nPoints;
  unsigned int nmax = static_cast<unsigned int>(radius + 1.);
  vector<Point> result(nPoints * nPoints);
  unsigned int i = 0;

  for (unsigned int ith = 0; ith < nPoints; ++ith)
    for (unsigned int iph = 0; iph < nPoints; ++iph)
    {
      result[i].theta = ith * dtheta;
      result[i].phi = -M_PI + iph * dphi;
      i++;
    }
#pragma omp parallel for
  for (unsigned int k = 0; k < result.size(); ++k)
  {
    DVec3 dir = sphericalToCartesian(1., result[k].theta, result[k].phi);
    SolidAngleCharacteristic chiS(dir, radius, theta);

    result[k].n = ThreadedSum::sum(chiS, nmax, debug);
  }

  char buf[256];
  ofstream file(filename);
  for (unsigned int k = 0; k < result.size(); ++k)
  {
    snprintf(buf, 256, "%10f %10f %.15e", result[k].theta, result[k].phi, result[k].n);
    file << buf << endl;
  }
  file.close();

  return 0;
}