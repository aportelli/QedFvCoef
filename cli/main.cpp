#include "utils.hpp"
#include <OptParser.hpp>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <vector>
extern "C"
{
#include <qedfvcoef.h>
}

using namespace optp;
using namespace std;

int main(int argc, const char *argv[])
{
  OptParser opt;
  bool parsed, debug, rest = true;
  double j, eta;
  unsigned int nmax;
  dvec3 v;

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
    string buf;
    stringstream vstr(opt.optionValue("v"));
    getline(vstr, buf, ',');
    v[0] = strTo<double>(buf);
    getline(vstr, buf, ',');
    v[1] = strTo<double>(buf);
    getline(vstr, buf, ',');
    v[2] = strTo<double>(buf);
    rest = false;
  }

  qedfv_context *ctx = qedfv_create_context();
  double coef;

  ctx->debug = debug;
  if (rest)
  {
    coef = qedfv_coef_rest(j, eta, nmax, ctx);
  }
  else
  {
    coef = qedfv_coef(j, v, eta, nmax, ctx);
  }
  printf("%f\n", coef);
  qedfv_destroy_context(ctx);

  return 0;
}
