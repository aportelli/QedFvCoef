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

#pragma once
#include <array>
#include <functional>
#include <gsl/gsl_integration.h>
#include <map>
#include <memory>
#include <string>

#ifndef QEDFV_DEFAULT_ERROR
#define QEDFV_DEFAULT_ERROR 1.e-8
#endif
#ifndef QEDFV_GSL_INT_LIMIT
#define QEDFV_GSL_INT_LIMIT 1000
#endif
#ifndef QEDFV_GSL_INT_PREC
#define QEDFV_GSL_INT_PREC 1.e-10
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

DVec3 sphericalToCartesian(const double r, const double theta, const double phi)
{
  return {r * cos(phi) * sin(theta), r * sin(phi) * sin(theta), r * cos(theta)};
}

DVec3 CartesianToSpherical(const double x, const double y, const double z)
{
  return {sqrt(x*x + y*y + z*z), atan2(sqrt(x*x + y*y), z), atan2(y, x)};
}

// Floating-point equal operator ///////////////////////////////////////////////////////////////////
constexpr double cmpEps = 1.0e-10;

inline bool isEqual(double a, double b)
{
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * cmpEps);
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

  enum class Qed
  {
    L,
    r
  };

public:
  QedFvCoef(const Qed qed = Qed::L, const bool debug = false);
  ~QedFvCoef(void);
  void setDebug(const bool debug = true);
  void setQed(const Qed qed);
  bool isDebug(void) const;
  Qed getQed(void) const;
  static Qed parseQed(const std::string str);
  double operator()(const double j, const double eta, const unsigned int nmax);
  double operator()(const double j, const Params par);
  double operator()(const double j, const DVec3 v, const double eta, const unsigned int nmax);
  double operator()(const double j, const DVec3 v, const Params par);
  double qedrTerm(const DVec3 v);
  double a(const double k, const DVec3 &v);
  double r(const double j);
  double rBar(const double j);
  Params tune(const double j, const double residual = QEDFV_DEFAULT_ERROR, const double eta0 = 1.0,
              const double etaFactor = 0.98, const unsigned int nmax0 = 5,
              const unsigned int nmaxStep = 5);
  Params tune(const double j, const Params par, const double residual = QEDFV_DEFAULT_ERROR,
              const double etaFactor = 0.98, const unsigned int nmaxStep = 5);
  Params tune(const double j, const DVec3 v, const double residual = QEDFV_DEFAULT_ERROR,
              const double eta0 = 1.0, const double etaFactor = 0.98, const unsigned int nmax0 = 5,
              const unsigned int nmaxStep = 5);
  Params tune(const double j, const DVec3 v, const Params par,
              const double residual = QEDFV_DEFAULT_ERROR, const double etaFactor = 0.98,
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
  Qed qed_{Qed::L};
  double jCache_, intCache_, intError_;
  gsl_integration_workspace *intWorkspace_;
  static const std::map<std::string, Qed> qedMap;
};

} // namespace qedfv
