#include "qedfvcoef.hpp"
#include "timer.h"
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <omp.h>

using namespace qedfv;

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
    fprintf(stderr, "error: j = 3 is not implemented");
    exit(EXIT_FAILURE);
  }

  return result;
}

double QedFvCoef::operator()(const double j, const DVec3 v, const double eta,
                             const unsigned int nmax)
{
  double result, aVal;

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
    fprintf(stderr, "error: j = 3 is not implemented");
    exit(EXIT_FAILURE);
  }

  return result;
}

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
  time = -_qedfv_ms();
  rj = integrate(i);
  time += _qedfv_ms();
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
  time = -_qedfv_ms();
  rbarj = integrate(i);
  time += _qedfv_ms();
  if (debug_)
  {
    printf("[QedFv]: computed Rbar_j, j= %f, Rbar_j= %f, error= %e, time= %f ms\n", j, intCache_,
           intError_, time);
  }

  return rbarj;
}

double QedFvCoef::threadedSum(Summand &func, const unsigned int nmax)
{
  double sum = 0., time;
  int n0, n1, n2, inmax = nmax;

  time = -_qedfv_ms();
#pragma omp parallel for collapse(3) reduction(+ : sum)
  for (n0 = -inmax; n0 <= inmax; ++(n0))
    for (n1 = -inmax; n1 <= inmax; ++(n1))
      for (n2 = -inmax; n2 <= inmax; ++(n2))
      {
        IVec3 n = {n0, n1, n2};
        sum += func(n);
      }
  time += _qedfv_ms();
  if (debug_)
  {
    printf("[QedFv]: threaded sum, result= %f, nmax= %d, time= %f ms\n", sum, nmax, time);
  }

  return sum;
}

double QedFvCoef::acceleratedSum(const double j, const double eta, const unsigned int nmax)
{
  double sum;

  Summand wrapper = [this, j, eta](const IVec3 &n) { return summand(n, j, eta); };
  sum = threadedSum(wrapper, nmax);

  return sum;
}

double QedFvCoef::acceleratedSum(const double j, const DVec3 v, const double eta,
                                 const unsigned int nmax)
{
  double sum;

  Summand wrapper = [this, j, v, eta](const IVec3 &n) { return summand(n, j, v, eta); };
  sum = threadedSum(wrapper, nmax);

  return sum;
}
