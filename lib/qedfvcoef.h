#pragma once
#include <stdbool.h>
#include <gsl/gsl_integration.h>

typedef double dvec3[3];
typedef int ivec3[3];

typedef struct qedfv_context_t
{
  double j_cache;
  double int_cache;
  double int_error;
  bool debug;
  gsl_integration_workspace *int_workspace;
} qedfv_context;

qedfv_context *qedfv_create_context(void);
void qedfv_destroy_context(qedfv_context *ctx);
double qedfv_coef_rest(const double j, const double eta, const unsigned int nmax,
                       qedfv_context *ctx);
double qedfv_coef(const double j, const dvec3 v, const double eta,
                  const unsigned int nmax, qedfv_context *ctx);
