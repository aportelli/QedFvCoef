#include "utils.hpp"
#include <OptParser.hpp>
#include <cstdio>
#include <iostream>
#include <qedfvcoef.hpp>
#include <sstream>
#include <vector>

using namespace optp;
using namespace std;
using namespace qedfv;

int main(int argc, const char *argv[])
{
  OptParser opt;
  bool parsed, debug, rest = true;
  double j, error;
  DVec3 v;

  opt.addOption("v", "velocity", OptParser::OptType::value, true,
                "velocity as comma-separated list (e.g. 0.1,0.2,0.3)");
  opt.addOption("e", "error", OptParser::OptType::value, true, "target relative error",
                strFrom(QEDFV_EPSILON));
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
  if (opt.gotOption("v"))
  {
    v = opt.optionValue<DVec3>("v");
    rest = false;
  }

  double c;

  QedFvCoef coef(debug);
  QedFvCoef::Params par;

  if (rest)
  {
    par = coef.tune(j, error);
    c = coef(j, par);
  }
  else
  {
    par = coef.tune(j, v, error);
    c = coef(j, v, par);
  }
  printf("Converged for eta= %.2f, nmax= %d, error= %.2e\n", par.eta, par.nmax, fabs(error * c));
  printf("Coefficient: %.15e\n", c);

  return 0;
}
