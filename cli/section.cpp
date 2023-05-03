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
#include <random>
#include <qedfvcoef.hpp>

using namespace optp;
using namespace std;
using namespace qedfv;

struct Coef
{
  double theta, phi, c;
};

int main(int argc, const char *argv[])
{
  OptParser opt;
  bool parsed, debug, tuned = false;
  double j, error;
  double vn;
  unsigned int nPoints;
  QedFvCoef::Params par;
  string filename;

  opt.addOption("v", "velocity", OptParser::OptType::value, false, "velocity norm");
  opt.addOption("n", "npoints", OptParser::OptType::value, false,
                "number of points for scan in phi");
  opt.addOption("o", "output", OptParser::OptType::value, false, "output file");
  opt.addOption("e", "error", OptParser::OptType::value, true, "target relative error",
                strFrom(QEDFV_DEFAULT_ERROR));
  opt.addOption("p", "parameters", OptParser::OptType::value, true,
                "algorithm parameters as eta,nmax (e.g. 0.5,50), (default: auto-tuned)");
  opt.addOption("d", "debug", OptParser::OptType::trigger, true, "show debug messages");
  opt.addOption("", "help", OptParser::OptType::trigger, true, "show this help message and exit");
  parsed = opt.parse(argc, argv);
  if (!parsed or (opt.getArgs().size() != 1) or opt.gotOption("help"))
  {
    cerr << "usage: " << argv[0] << " <j>" << endl;
    cerr << endl << "Possible options:" << endl << opt << endl;

    return EXIT_FAILURE;
  }
  j = strTo<double>(opt.getArgs()[0]);
  error = opt.optionValue<double>("e");
  debug = opt.gotOption("d");
  vn = opt.optionValue<double>("v");
  nPoints = opt.optionValue<unsigned int>("n");
  filename = opt.optionValue("o");
  if (opt.gotOption("p"))
  {
    par = opt.optionValue<QedFvCoef::Params>("p");
    tuned = true;
  }

  QedFvCoef coef(debug);
  double phiMin = 0., phiMax = 0.25 * M_PI;
  double dphi = (phiMax - phiMin) / nPoints;
  double thetaMin = 0., thetaMax = 0.5 * M_PI;
  vector<Coef> result(nPoints * nPoints);
  double sphi = 0.;
  unsigned int i = 0;

  if (!tuned)
  {
    double c;

    par = coef.tune(j, {vn, 0., 0.}, error);
    c = coef(j, {vn, 0., 0.}, par);
    printf("Converged for eta= %.2f, nmax= %d, error= %.2e\n", par.eta, par.nmax, fabs(error * c));
  }
  else
  {
    coef(j, {vn, 0., 0.}, par);
  }
  for (unsigned int iph = 0; iph < nPoints; ++iph)
  {
    sphi = sin(dphi * iph);
    thetaMin = -2. * atan(sphi - sqrt(1. + pow(sphi,2.)));

    for (double ith = thetaMin; ith <= thetaMax; ith+=dphi)
    {
      result[i].phi = phiMin + iph * dphi;
      result[i].theta = ith;
      i++;
    }
  }
  result.resize(i);

#pragma omp parallel for
  for (unsigned int k = 0; k < result.size(); ++k)
  {
    DVec3 v = sphericalToCartesian(vn, result[k].theta, result[k].phi);
    result[k].c = coef(j, v, par);
  }

  char buf[256];
  ofstream file(filename);
  for (unsigned int k = 0; k < result.size(); ++k)
  {
    snprintf(buf, 256, "%10f %10f %.15e", result[k].theta, result[k].phi, result[k].c);
    file << buf << endl;
  }
  file.close();

  return 0;
}