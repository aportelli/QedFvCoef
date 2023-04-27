#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "qedfvcoef.h"
#include "timer.h"

#ifndef QEDFV_EPSILON
#define QEDFV_EPSILON 1.e-7
#endif

#ifndef QEDFV_GSL_INT_LIMIT
#define QEDFV_GSL_INT_LIMIT 1000
#endif
#ifndef QEDFV_GSL_INT_PREC
#define QEDFV_GSL_INT_PREC QEDFV_EPSILON
#endif

// private functions /////////////////////////////////////////////////////////////////////
void _qedfv_free(void *pt)
{
  if (pt != NULL)
  {
    free(pt);
  }
}

void _qedfv_debug_printf(qedfv_context *ctx, const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  if (ctx->debug)
  {
    printf("[QedFv]: ");
    vprintf(fmt, args);
  }
  va_end(args);
}

void _qedfv_fatal(qedfv_context *ctx)
{
  qedfv_destroy_context(ctx);
  exit(EXIT_FAILURE);
}

int _qedfv_inorm2(const ivec3 n) { return n[0] * n[0] + n[1] * n[1] + n[2] * n[2]; }

double _qedfv_dnorm2(const dvec3 n) { return n[0] * n[0] + n[1] * n[1] + n[2] * n[2]; }

bool _qedfv_is_equal(double a, double b)
{
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * QEDFV_EPSILON);
}

double _qedfv_rest_summand(ivec3 n, const double j, const double eta)
{
  double norm = sqrt(_qedfv_inorm2(n));

  return (1. - pow(tanh(sinh(eta * norm)), j + 2.)) / pow(norm, j);
}

double _qedfv_summand(ivec3 n, const double j, const dvec3 v, const double eta)
{
  double norm = sqrt(_qedfv_inorm2(n));
  dvec3 nhat = {n[0] / norm, n[1] / norm, n[2] / norm};
  double d = 1. / (1. - (nhat[0] * v[0] + nhat[1] * v[1] + nhat[2] * v[2]));

  return d * (1. - pow(tanh(sinh(eta * norm * pow(d, -1. / (j + 2.)))), j + 2.)) /
         pow(norm, j);
}

double _qedfv_r_integrand(double r, void *param)
{
  double j = *(double *)(param);
  return pow(r, 2. - j) * (1. - pow(tanh(sinh(r)), j + 2.));
}

double _qedfv_rbar_integrand(double r, void *param)
{
  double j = *(double *)(param);
  return pow(r, 2. - j) * pow(tanh(sinh(r)), j + 2.);
}

double _qedfv_A(const double k, const dvec3 v)
{
  double vn = sqrt(_qedfv_dnorm2(v));

  if (_qedfv_is_equal(k, 1.) && (vn > 0.))
  {
    return atanh(vn) / vn;
  }
  else if (!_qedfv_is_equal(k, 1.) && (vn > 0.))
  {
    return (pow(1 - vn * vn, -k) *
            (-pow(1 - vn, k) * (1 + vn) - (-1 + vn) * pow(1 + vn, k))) /
           ((-1 + k) * vn) / 2;
  }
  else
  {
    return 1.;
  }
}
void _qedfv_integrate(double (*func)(double, void *), void *params, qedfv_context *ctx)
{
  gsl_function gsl_f;

  gsl_f.function = func;
  gsl_f.params = params;
  gsl_integration_qagiu(&gsl_f, 0., 0., QEDFV_GSL_INT_PREC, QEDFV_GSL_INT_LIMIT,
                        ctx->int_workspace, &ctx->int_cache, &ctx->int_error);
}

double _qedfv_rest_accelerated_sum(const double j, const double eta,
                                   const unsigned int nmax, qedfv_context *ctx)
{
  double sum = 0., time;
  int n0, n1, n2, inmax = nmax;

  time = -_qedfv_ms();
#pragma omp parallel for collapse(3) reduction(+ : sum)
  for (n0 = -inmax; n0 <= inmax; ++(n0))
    for (n1 = -inmax; n1 <= inmax; ++(n1))
      for (n2 = -inmax; n2 <= inmax; ++(n2))
      {
        ivec3 n = {n0, n1, n2};
        sum += (_qedfv_inorm2(n) != 0) ? _qedfv_rest_summand(n, j, eta) : 0.;
      }
  time += _qedfv_ms();
  _qedfv_debug_printf(ctx, "rest_accelerated_sum, result= %f, nmax= %d, time= %f ms\n",
                      sum, nmax, time);

  return sum;
}

double _qedfv_accelerated_sum(const double j, const dvec3 v, const double eta,
                              const unsigned int nmax, qedfv_context *ctx)
{
  double sum = 0., time;
  int n0, n1, n2, inmax = nmax;

  time = -_qedfv_ms();
#pragma omp parallel for collapse(3) reduction(+ : sum)
  for (n0 = -inmax; n0 <= inmax; ++(n0))
    for (n1 = -inmax; n1 <= inmax; ++(n1))
      for (n2 = -inmax; n2 <= inmax; ++(n2))
      {
        ivec3 n = {n0, n1, n2};
        sum += (_qedfv_inorm2(n) != 0) ? _qedfv_summand(n, j, v, eta) : 0.;
      }
  time += _qedfv_ms();
  _qedfv_debug_printf(ctx, "accelerated_sum, result= %f, nmax= %d, time= %f ms\n", sum,
                      nmax, time);

  return sum;
}

// public functions //////////////////////////////////////////////////////////////////////
qedfv_context *qedfv_create_context(void)
{
  qedfv_context *ctx = (qedfv_context *)malloc(sizeof(qedfv_context));
  ctx->debug = false;
  ctx->j_cache = NAN;
  ctx->int_workspace = gsl_integration_workspace_alloc(QEDFV_GSL_INT_LIMIT);

  return ctx;
}

void qedfv_destroy_context(qedfv_context *ctx)
{
  gsl_integration_workspace_free(ctx->int_workspace);
  _qedfv_free(ctx);
}

double qedfv_r(const double j, qedfv_context *ctx)
{
  double time;

  time = -_qedfv_ms();
  _qedfv_integrate(&_qedfv_r_integrand, (void *)(&j), ctx);
  time += _qedfv_ms();
  _qedfv_debug_printf(ctx, "computed R_j, j= %f, R_j= %f, error= %e, time= %f ms\n", j,
                      ctx->int_cache, ctx->int_error, time);

  return ctx->int_cache;
}

double qedfv_rbar(const double j, qedfv_context *ctx)
{
  double time;

  time = -_qedfv_ms();
  _qedfv_integrate(&_qedfv_rbar_integrand, (void *)(&j), ctx);
  time += _qedfv_ms();
  _qedfv_debug_printf(ctx, "computed Rbar_j, j= %f, Rbar_j= %f, error= %e, time= %f ms\n",
                      j, ctx->int_cache, ctx->int_error, time);

  return ctx->int_cache;
}

double qedfv_coef_rest(const double j, const double eta, const unsigned int nmax,
                       qedfv_context *ctx)
{
  double result;

  result = _qedfv_rest_accelerated_sum(j, eta, nmax, ctx);
  if (j < 3.)
  {
    if (!_qedfv_is_equal(j, ctx->j_cache))
    {
      ctx->j_cache = j;
      ctx->int_cache = qedfv_r(j, ctx);
    }
    result -= 4 * M_PI * pow(eta, j - 3.) * ctx->int_cache;
  }
  else if (j > 3.)
  {
    if (!_qedfv_is_equal(j, ctx->j_cache))
    {
      ctx->j_cache = j;
      ctx->int_cache = qedfv_rbar(j, ctx);
    }
    result += 4 * M_PI * pow(eta, j - 3.) * ctx->int_cache;
  }
  else
  {
    fprintf(stderr, "error: j = 3 is not implemented");
    _qedfv_fatal(ctx);
  }

  return result;
}

double qedfv_coef(const double j, const dvec3 v, const double eta,
                  const unsigned int nmax, qedfv_context *ctx)
{
  double result, a;

  result = _qedfv_accelerated_sum(j, v, eta, nmax, ctx);
  a = _qedfv_A(5. / (j + 2.), v);
  if (j < 3.)
  {
    if (!_qedfv_is_equal(j, ctx->j_cache))
    {
      ctx->j_cache = j;
      ctx->int_cache = qedfv_r(j, ctx);
    }
    result -= 4 * M_PI * a * pow(eta, j - 3.) * ctx->int_cache;
  }
  else if (j > 3.)
  {
    if (!_qedfv_is_equal(j, ctx->j_cache))
    {
      ctx->j_cache = j;
      ctx->int_cache = qedfv_rbar(j, ctx);
    }
    result += 4 * M_PI * a * pow(eta, j - 3.) * ctx->int_cache;
  }
  else
  {
    fprintf(stderr, "error: j = 3 is not implemented");
    _qedfv_fatal(ctx);
  }

  return result;
}
