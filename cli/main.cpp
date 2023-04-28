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
  double j, eta;
  unsigned int nmax;
  DVec3 v;

  opt.addOption("v", "velocity", OptParser::OptType::value, true,
                "velocity as comma-separated list (e.g. 0.1,0.2,0.3)");
  opt.addOption("d", "debug", OptParser::OptType::trigger, true, "show debug messages");
  opt.addOption("", "help", OptParser::OptType::trigger, true, "show this help message and exit");
  parsed = opt.parse(argc, argv);
  if (!parsed or (opt.getArgs().size() != 3) or opt.gotOption("help"))
  {
    cerr << "usage: " << argv[0] << " <j> <eta> <nmax>" << endl;
    cerr << endl << "Possible options:" << endl << opt << endl;

    return EXIT_FAILURE;
  }
  j = strTo<double>(opt.getArgs()[0]);
  eta = strTo<double>(opt.getArgs()[1]);
  nmax = strTo<unsigned int>(opt.getArgs()[2]);
  debug = opt.gotOption("d");
  if (opt.gotOption("v"))
  {
    v = strTo<DVec3>(opt.optionValue("v"));
    rest = false;
  }

  double c;

  QedFvCoef coef(debug);
  if (rest)
  {
    c = coef(j, eta, nmax);
  }
  else
  {
    c = coef(j, v, eta, nmax);
  }
  cout << c << endl;

  return 0;
}
