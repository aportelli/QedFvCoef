#pragma once
#include <array>
#include <functional>
#include <gsl/gsl_integration.h>
#include <memory>

#ifndef QEDFV_EPSILON
#define QEDFV_EPSILON 1.e-8
#endif
#ifndef QEDFV_GSL_INT_LIMIT
#define QEDFV_GSL_INT_LIMIT 1000
#endif
#ifndef QEDFV_GSL_INT_PREC
#define QEDFV_GSL_INT_PREC QEDFV_EPSILON
#endif

namespace qedfv
{

// Vector and basic operations /////////////////////////////////////////////////////////////////////
template <typename T>
using Vec3 = std::array<T, 3>;

typedef Vec3<double> DVec3;
typedef Vec3<int> IVec3;

template <typename T, typename U>
inline auto dot(const Vec3<T> &v, const Vec3<U> &w) -> decltype(v[0] * w[0])
{
  return v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
}

template <typename T>
inline T norm2(const Vec3<T> &v)
{
  return dot(v, v);
}

// Floating-point equal operator ///////////////////////////////////////////////////////////////////
inline bool isEqual(double a, double b)
{
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * QEDFV_EPSILON);
}

// Main class //////////////////////////////////////////////////////////////////////////////////////
class QedFvCoef
{
public:
  struct Params
  {
    double eta;
    unsigned int nmax;
  };

public:
  QedFvCoef(const bool debug = false);
  ~QedFvCoef(void);
  double operator()(const double j, const double eta, const unsigned int nmax);
  double operator()(const double j, const Params par);
  double operator()(const double j, const DVec3 v, const double eta, const unsigned int nmax);
  double operator()(const double j, const DVec3 v, const Params par);
  double a(const double k, const DVec3 &v);
  double r(const double j);
  double rBar(const double j);
  Params tune(const double j, const DVec3 v, const double residual = QEDFV_EPSILON,
              const double eta0 = 1.0, const double etaFactor = 0.98, const unsigned int nmax0 = 5,
              const unsigned int nmaxStep = 5);
  Params tune(const double j, const double residual = QEDFV_EPSILON, const double eta0 = 1.0,
              const double etaFactor = 0.98, const unsigned int nmax0 = 5,
              const unsigned int nmaxStep = 5);

private:
  typedef std::function<double(double)> Integrand;
  typedef std::function<double(const IVec3 &n)> Summand;
  typedef std::function<double(const Params par)> CoefFunc;

private:
  double summand(IVec3 n, const double j, const double eta);
  double summandJ0(IVec3 n, const double eta);
  double summand(IVec3 n, const double j, const DVec3 v, const double eta);
  double summandJ0(IVec3 n, const DVec3 v, const double eta);
  double integrate(Integrand &func);
  Params tune(CoefFunc &coef, const double residual, const double eta0, const double etaFactor,
              const unsigned int nmax0, const unsigned int nmaxStep);
  double threadedSum(Summand &func, const unsigned int nmax);
  double acceleratedSum(const double j, const double eta, const unsigned int nmax);
  double acceleratedSum(const double j, const DVec3 v, const double eta, const unsigned int nmax);

private:
  bool debug_{false};
  double jCache_, intCache_, intError_;
  gsl_integration_workspace *intWorkspace_;
};

} // namespace qedfv
