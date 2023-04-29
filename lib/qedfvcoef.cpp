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

#include "qedfvcoef.hpp"
#include <chrono>
#include <cmath>
#include <exception>
#include <omp.h>

using namespace qedfv;

// High-resolution timer ///////////////////////////////////////////////////////////////////////////
double clockMs(void)
{
  auto t = std::chrono::high_resolution_clock::now().time_since_epoch();
  auto d = std::chrono::duration_cast<std::chrono::nanoseconds>(t);
  return static_cast<double>(d.count()) * 1e-6;
}

// Public interface ////////////////////////////////////////////////////////////////////////////////
QedFvCoef::QedFvCoef(const bool debug)
{
  debug_ = debug;
  jCache_ = std::nan("");
  intWorkspace_ = gsl_integration_workspace_alloc(QEDFV_GSL_INT_LIMIT);
}

QedFvCoef::~QedFvCoef(void) { gsl_integration_workspace_free(intWorkspace_); }

double QedFvCoef::operator()(const double j, const double eta, const unsigned int nmax)
{
  double result;

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

  return result;
}

double QedFvCoef::operator()(const double j, const Params par)
{
  return (*this)(j, par.eta, par.nmax);
}

double QedFvCoef::operator()(const double j, const DVec3 v, const double eta,
                             const unsigned int nmax)
{
  double result, aVal;

  if (norm2(v) >= 1.)
  {
    throw std::logic_error("velocity has norm larger than 1");
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

  return result;
}

double QedFvCoef::operator()(const double j, const DVec3 v, const Params par)
{
  return (*this)(j, v, par.eta, par.nmax);
}

double QedFvCoef::a(const double k, const DVec3 &v)
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

double QedFvCoef::r(const double j)
{
  double time, rj;

  Integrand i = [j](double r) { return pow(r, 2. - j) * (1. - pow(tanh(sinh(r)), j + 2.)); };
  time = -clockMs();
  rj = integrate(i);
  time += clockMs();
  if (debug_)
  {
    printf("[QedFv]: computed R_j, j= %f, R_j= %f, error= %e, time= %f ms\n", j, intCache_,
           intError_, time);
  }

  return rj;
}

double QedFvCoef::rBar(const double j)
{
  double time, rbarj;

  Integrand i = [j](double r) { return pow(r, 2. - j) * pow(tanh(sinh(r)), j + 2.); };
  time = -clockMs();
  rbarj = integrate(i);
  time += clockMs();
  if (debug_)
  {
    printf("[QedFv]: computed Rbar_j, j= %f, Rbar_j= %f, error= %e, time= %f ms\n", j, intCache_,
           intError_, time);
  }

  return rbarj;
}

QedFvCoef::Params QedFvCoef::tune(const double j, const DVec3 v, const double residual,
                                  const double eta0, const double etaFactor,
                                  const unsigned int nmax0, const unsigned int nmaxStep)
{
  CoefFunc coef = [this, j, v](Params par) { return (*this)(j, v, par); };

  return tune(coef, residual, eta0, etaFactor, nmax0, nmaxStep);
}

QedFvCoef::Params QedFvCoef::tune(const double j, const double residual, const double eta0,
                                  const double etaFactor, const unsigned int nmax0,
                                  const unsigned int nmaxStep)
{
  CoefFunc coef = [this, j](Params par) { return (*this)(j, par); };

  return tune(coef, residual, eta0, etaFactor, nmax0, nmaxStep);
}

// Private methods /////////////////////////////////////////////////////////////////////////////////
double QedFvCoef::summand(IVec3 n, const double j, const double eta)
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

double QedFvCoef::summandJ0(IVec3 n, const double eta)
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

double QedFvCoef::summand(IVec3 n, const double j, const DVec3 v, const double eta)
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

double QedFvCoef::summandJ0(IVec3 n, const DVec3 v, const double eta)
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

double QedFvCoef::integrate(Integrand &func)
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

QedFvCoef::Params QedFvCoef::tune(CoefFunc &coef, const double residual, const double eta0,
                                  const double etaFactor, const unsigned int nmax0,
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
    } while (res > QEDFV_DEFAULT_ERROR);

    return previous;
  };

  par.eta = eta0;
  par.nmax = nmax0;
  previous = converge(par);
  if (debug_)
  {
    printf("[QedFv]: eta= %.2f, nmax= %d, c= %.15e\n", par.eta, par.nmax, previous);
  }
  do
  {
    par.eta *= etaFactor;
    par.nmax = (par.nmax > 10) ? (par.nmax - 10) : par.nmax;
    buf = converge(par);
    res = fabs(buf - previous) / (0.5 * (fabs(buf) + fabs(previous)));
    previous = buf;
    if (debug_)
    {
      printf("[QedFv]: eta= %.2f, nmax= %d, c= %.15e, residual= %.2e\n", par.eta, par.nmax, buf,
             res);
    }
  } while (res > residual);

  return par;
}

double QedFvCoef::threadedSum(Summand &func, const unsigned int nmax)
{
  double sum = 0., time;
  int n0, n1, n2, inmax = nmax;

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
  if (debug_)
  {
    printf("[QedFv]: threaded sum, result= %f, nmax= %d, time= %f ms\n", sum, nmax, time);
  }

  return sum;
}

double QedFvCoef::acceleratedSum(const double j, const double eta, const unsigned int nmax)
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
  sum = threadedSum(wrapper, nmax);

  return sum;
}

double QedFvCoef::acceleratedSum(const double j, const DVec3 v, const double eta,
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
  sum = threadedSum(wrapper, nmax);

  return sum;
}
