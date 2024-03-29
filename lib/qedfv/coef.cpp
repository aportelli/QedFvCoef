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

#include <chrono>
#include <cmath>
#include <exception>
#include <omp.h>
#include <qedfv/coef.hpp>
#include <qedfv/latticesum.hpp>

using namespace qedfv;

const std::map<std::string, Coef::Qed> Coef::qedMap = {{"L", Qed::L}, {"r", Qed::r}};

// Public interface ////////////////////////////////////////////////////////////////////////////////
Coef::Coef(const Coef::Qed qed, const bool debug)
{
  setDebug(debug);
  setQed(qed);
  jCache_ = std::nan("");
  intWorkspace_ = gsl_integration_workspace_alloc(QEDFV_GSL_INT_LIMIT);
}

Coef::~Coef(void) { gsl_integration_workspace_free(intWorkspace_); }

void Coef::setDebug(const bool debug) { debug_ = debug; }

void Coef::setQed(const Qed qed) { qed_ = qed; }

bool Coef::isDebug(void) const { return debug_; }

Coef::Qed Coef::getQed(void) const { return qed_; }

Coef::Qed Coef::parseQed(const std::string str) { return qedMap.at(str); }

double Coef::operator()(const double j, const double eta, const unsigned int nmax)
{
  double result;

  if (!isEqual(j, 0.))
  {
    result = acceleratedSum(j, eta, nmax);
    if (j < 3.)
    {
      if (!isEqual(j, jCache_))
      {
        jCache_ = j;
        intCache_ = r(j);
      }
      result -= 4 * M_PI * pow(eta, j - 3.) * intCache_;
    }
    else if (j > 3.)
    {
      if (!isEqual(j, jCache_))
      {
        jCache_ = j;
        intCache_ = rBar(j);
      }
      result += 4 * M_PI * pow(eta, j - 3.) * intCache_;
    }
    else
    {
      throw std::logic_error("error: j = 3 is not implemented");
    }
  }
  else
  {
    result = -1.;
  }
  switch (getQed())
  {
  case Qed::r:
    result += 1.;
    break;
  default:
    break;
  }

  return result;
}

double Coef::operator()(const double j, const Params par) { return (*this)(j, par.eta, par.nmax); }

double Coef::operator()(const double j, const DVec3 v, const double eta, const unsigned int nmax)
{
  double result, aVal, vn2 = norm2(v);

  if (vn2 >= 1.)
  {
    throw std::logic_error("velocity has norm larger than 1");
  }
  if (isEqual(vn2, 0.))
  {
    return (*this)(j, eta, nmax);
  }
  result = acceleratedSum(j, v, eta, nmax);
  aVal = a(5. / (j + 2.), v);
  if (j < 3.)
  {
    if (!isEqual(j, jCache_))
    {
      jCache_ = j;
      intCache_ = r(j);
    }
    result -= 4 * M_PI * aVal * pow(eta, j - 3.) * intCache_;
  }
  else if (j > 3.)
  {
    if (!isEqual(j, jCache_))
    {
      jCache_ = j;
      intCache_ = rBar(j);
    }
    result += 4 * M_PI * aVal * pow(eta, j - 3.) * intCache_;
  }
  else
  {
    throw std::logic_error("error: j = 3 is not implemented");
  }
  switch (getQed())
  {
  case Qed::r:
    result += qedrTerm(v);
    break;
  default:
    break;
  }

  return result;
}

double Coef::operator()(const double j, const DVec3 v, const Params par)
{
  return (*this)(j, v, par.eta, par.nmax);
}

double Coef::qedrTerm(const DVec3 v)
{
  double term = 0.;

  for (const double &vj : v)
  {
    term += 1. / (1 + vj) + 1. / (1 - vj);
  }

  return term / 6.;
}

double Coef::a(const double k, const DVec3 &v)
{
  double vn = sqrt(norm2(v));

  if (isEqual(k, 1.) && (vn > 0.))
  {
    return atanh(vn) / vn;
  }
  else if (!isEqual(k, 1.) && (vn > 0.))
  {
    return (pow(1 - vn * vn, -k) * (-pow(1 - vn, k) * (1 + vn) - (-1 + vn) * pow(1 + vn, k))) /
           ((-1 + k) * vn) / 2;
  }
  else
  {
    return 1.;
  }
}

double Coef::r(const double j)
{
  double time, rj;

  Integrand i = [j](double r) { return pow(r, 2. - j) * (1. - pow(tanh(sinh(r)), j + 2.)); };
  time = -clockMs();
  rj = integrate(i);
  time += clockMs();
  dgbPrintf(isDebug(), "computed R_j, j= %f, R_j= %f, error= %e, time= %f ms\n", j, intCache_,
            intError_, time);

  return rj;
}

double Coef::rBar(const double j)
{
  double time, rbarj;

  Integrand i = [j](double r) { return pow(r, 2. - j) * pow(tanh(sinh(r)), j + 2.); };
  time = -clockMs();
  rbarj = integrate(i);
  time += clockMs();
  dgbPrintf(isDebug(), "computed Rbar_j, j= %f, Rbar_j= %f, error= %e, time= %f ms\n", j, intCache_,
            intError_, time);

  return rbarj;
}

Coef::Params Coef::tune(const double j, const double residual, const double eta0,
                        const double etaInvStep, const unsigned int nmax0,
                        const unsigned int nmaxStep)
{
  CoefFunc coef = [this, j](Params par) { return (*this)(j, par); };

  return tune(coef, residual, eta0, etaInvStep, nmax0, nmaxStep);
}

Coef::Params Coef::tune(const double j, const Params par, const double residual,
                        const double etaInvStep, const unsigned int nmaxStep)
{
  CoefFunc coef = [this, j](Params par) { return (*this)(j, par); };

  return tune(coef, residual, par.eta, etaInvStep, par.nmax, nmaxStep);
}

Coef::Params Coef::tune(const double j, const DVec3 v, const double residual, const double eta0,
                        const double etaInvStep, const unsigned int nmax0,
                        const unsigned int nmaxStep)
{
  CoefFunc coef = [this, j, v](Params par) { return (*this)(j, v, par); };

  return tune(coef, residual, eta0, etaInvStep, nmax0, nmaxStep);
}

Coef::Params Coef::tune(const double j, const DVec3 v, const Params par, const double residual,
                        const double etaInvStep, const unsigned int nmaxStep)
{
  CoefFunc coef = [this, j, v](Params par) { return (*this)(j, v, par); };

  return tune(coef, residual, par.eta, etaInvStep, par.nmax, nmaxStep);
}

// Private methods /////////////////////////////////////////////////////////////////////////////////
double Coef::summand(IVec3 n, const double j, const double eta)
{
  double n2 = norm2(n);

  if (n2 == 0)
  {
    return 0.;
  }
  else
  {
    double norm = sqrt(n2);
    return (1. - pow(tanh(sinh(eta * norm)), j + 2.)) / pow(norm, j);
  }
}

double Coef::summandJ0(IVec3 n, const double eta)
{
  double n2 = norm2(n);

  if (n2 == 0)
  {
    return 0.;
  }
  else
  {
    double norm = sqrt(n2);
    return (1. - pow(tanh(sinh(eta * norm)), 2));
  }
}

double Coef::summand(IVec3 n, const double j, const DVec3 v, const double eta)
{
  double n2 = norm2(n);

  if (n2 == 0)
  {
    return 0.;
  }
  else
  {
    double norm = sqrt(n2);
    DVec3 nhat = {n[0] / norm, n[1] / norm, n[2] / norm};
    double d = 1. / (1. - dot(nhat, v));
    return d * (1. - pow(tanh(sinh(eta * norm * pow(d, -1. / (j + 2.)))), j + 2.)) / pow(norm, j);
  }
}

double Coef::summandJ0(IVec3 n, const DVec3 v, const double eta)
{
  double n2 = norm2(n);

  if (n2 == 0)
  {
    return 0.;
  }
  else
  {
    double norm = sqrt(n2);
    DVec3 nhat = {n[0] / norm, n[1] / norm, n[2] / norm};
    double d = 1. / (1. - dot(nhat, v));
    return d * (1. - pow(tanh(sinh(eta * norm / sqrt(d))), 2));
  }
}

double Coef::integrate(Integrand &func)
{
  gsl_function gsl_f;

  double (*wrapper)(double, void *) = [](double x, void *func) -> double
  { return (*static_cast<Integrand *>(func))(x); };

  gsl_f.function = wrapper;
  gsl_f.params = &func;
  gsl_integration_qagiu(&gsl_f, 0., 0., QEDFV_GSL_INT_PREC, QEDFV_GSL_INT_LIMIT, intWorkspace_,
                        &intCache_, &intError_);

  return intCache_;
}

Coef::Params Coef::tune(CoefFunc &coef, const double residual, const double eta0,
                        const double etaInvStep, const unsigned int nmax0,
                        const unsigned int nmaxStep)
{
  Params par;
  double previous, buf, res;

  auto converge = [&coef, nmaxStep](Params &par)
  {
    double previous = coef(par);
    double buf, res;
    do
    {
      par.nmax += nmaxStep;
      buf = coef(par);
      res = fabs(buf - previous) / (0.5 * (fabs(buf) + fabs(previous)));
      previous = buf;
    } while (res > 1.0e-2 * QEDFV_DEFAULT_ERROR);

    return previous;
  };

  par.eta = eta0;
  par.nmax = nmax0;
  previous = converge(par);
  dgbPrintf(isDebug(), "eta= %.4f nmax= %d c= %.15e\n", par.eta, par.nmax, previous);
  do
  {
    par.eta = par.eta / (1 + etaInvStep * par.eta);
    par.nmax = (par.nmax > 10) ? (par.nmax - 10) : par.nmax;
    buf = converge(par);
    res = fabs(buf - previous) / (0.5 * (fabs(buf) + fabs(previous)));
    res /= fabs(etaInvStep / (par.eta * etaInvStep - 1));
    previous = buf;
    dgbPrintf(isDebug(), "eta= %.4f nmax= %d c= %.15e residual= %.2e\n", par.eta, par.nmax, buf,
              res);
  } while (res > residual);

  return par;
}

double Coef::acceleratedSum(const double j, const double eta, const unsigned int nmax)
{
  double sum;
  Summand wrapper;

  if (isEqual(j, 0.))
  {
    wrapper = [this, eta](const IVec3 &n) { return summandJ0(n, eta); };
  }
  else
  {
    wrapper = [this, j, eta](const IVec3 &n) { return summand(n, j, eta); };
  }
  sum = ThreadedSum::sum(wrapper, nmax, isDebug());

  return sum;
}

double Coef::acceleratedSum(const double j, const DVec3 v, const double eta,
                            const unsigned int nmax)
{
  double sum;
  Summand wrapper;

  if (isEqual(j, 0.))
  {
    wrapper = [this, v, eta](const IVec3 &n) { return summandJ0(n, v, eta); };
  }
  else
  {
    wrapper = [this, j, v, eta](const IVec3 &n) { return summand(n, j, v, eta); };
  }
  sum = ThreadedSum::sum(wrapper, nmax, isDebug());

  return sum;
}
