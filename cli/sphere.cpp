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
                "number of points for each angle");
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
  double dtheta = M_PI / nPoints, dphi = 2. * M_PI / nPoints;
  vector<CoefPoint> result(nPoints * nPoints);
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