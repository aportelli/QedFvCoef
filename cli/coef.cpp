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
#include <iostream>
#include <qedfv/coef.hpp>

using namespace optp;
using namespace std;
using namespace qedfv;

int main(int argc, const char *argv[])
{
  OptParser opt;
  bool parsed, debug, rest = true, tuned = false;
  double j, error;
  DVec3 v;
  Coef::Params par;
  Coef::Qed qed;

  opt.addOption("v", "velocity", OptParser::OptType::value, true,
                "velocity as comma-separated list (e.g. 0.1,0.2,0.3)");
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
  if (opt.gotOption("v"))
  {
    v = opt.optionValue<DVec3>("v");
    rest = false;
  }
  if (opt.gotOption("p"))
  {
    par = opt.optionValue<Coef::Params>("p");
    tuned = true;
  }

  double c;

  Coef coef(qed, debug);

  if (rest)
  {
    if (!tuned)
    {
      par = coef.tune(j, error);
    }
    c = coef(j, par);
  }
  else
  {
    if (!tuned)
    {
      par = coef.tune(j, v, error);
    }
    c = coef(j, v, par);
  }
  if (!tuned)
  {
    printf("Converged for eta= %.2f, nmax= %d, error= %.2e\n", par.eta, par.nmax, fabs(error * c));
  }
  printf("Coefficient: %.15e\n", c);

  return 0;
}
