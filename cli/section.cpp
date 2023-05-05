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

#include "utils.hpp"
#include <OptParser.hpp>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <qedfv/coef.hpp>

using namespace optp;
using namespace std;
using namespace qedfv;

struct CoefPoint
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
  Coef::Params par;
  Coef::Qed qed;
  string filename;

  opt.addOption("v", "velocity", OptParser::OptType::value, false, "velocity norm");
  opt.addOption("n", "npoints", OptParser::OptType::value, false,
                "number of points for scan in phi");
  opt.addOption("o", "output", OptParser::OptType::value, false, "output file");
  opt.addOption("e", "error", OptParser::OptType::value, true, "target relative error",
                strFrom(QEDFV_DEFAULT_ERROR));
  opt.addOption("p", "parameters", OptParser::OptType::value, true,
                "algorithm parameters as eta,nmax (e.g. 0.5,50), (default: auto-tuned)");
  opt.addOption("q", "qed", OptParser::OptType::value, true,
                "QED theory to use, 'x' for QED_x, possible choices are L,r", "L");
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
  qed = opt.optionValue<Coef::Qed>("q");
  vn = opt.optionValue<double>("v");
  nPoints = opt.optionValue<unsigned int>("n");
  filename = opt.optionValue("o");
  if (opt.gotOption("p"))
  {
    par = opt.optionValue<Coef::Params>("p");
    tuned = true;
  }

  Coef coef(qed, debug);
  double phiMin = 0., phiMax = 0.25 * M_PI;
  double dphi = (phiMax - phiMin) / nPoints;
  double thetaMin = 0., thetaMax = 0.5 * M_PI;
  vector<CoefPoint> result(nPoints * nPoints);
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
    thetaMin = -2. * atan(sphi - sqrt(1. + pow(sphi, 2.)));

    for (double ith = thetaMin; ith <= thetaMax; ith += dphi)
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

  vector<CoefPoint> unfold(i * 48);
  unsigned int r = 0;

  for (unsigned int k = 0; k < result.size(); ++k)
  {
    DVec3 v = sphericalToCartesian(vn, result[k].theta, result[k].phi);
    sort(v.begin(), v.end());
    do
    {
      for (unsigned int i = 0; i < 1 << v.size(); ++i)
      {
        DVec3 w = v;
        for (unsigned int j = 0; j < v.size(); ++j)
        {
          if ((i >> j) & 1)
          {
            w[j] = -w[j];
          }
        }
        DVec3 vs = cartesianToSpherical(w[0], w[1], w[2]);
        unfold[r] = {vs[1], vs[2], result[k].c};
        r++;
      }
    } while (next_permutation(v.begin(), v.end()));
  }

  unfold.resize(r);

  char buf[256];
  ofstream file(filename);
  for (unsigned int k = 0; k < result.size(); ++k)
  {
    snprintf(buf, 256, "%10f %10f %.15e", result[k].theta, result[k].phi, result[k].c);
    file << buf << endl;
  }
  file.close();

  ofstream file2("unfold_" + filename);
  for (unsigned int k = 0; k < unfold.size(); ++k)
  {
    snprintf(buf, 256, "%10f %10f %.15e", unfold[k].theta, unfold[k].phi, unfold[k].c);
    file2 << buf << endl;
  }
  file2.close();

  return 0;
}