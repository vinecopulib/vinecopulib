// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <algorithm>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <cmath>
#include <utility>
#include <vector>

namespace vinecopulib {

namespace tools_integration {

//! @brief Integrates `f` over the unit interval.
//!
//! Suited to integrands that are singular or steep at `0` and `1`, as the
//! Kendall's tau integrals of the BB families are: tanh-sinh clusters its
//! abscissas at the endpoints. Use `integrate_zero_to_one_split()` when the
//! integrand instead peaks in the interior.
template<typename F>
inline double
integrate_zero_to_one(F&& f)
{
  // The bounds stay clear of 0 and 1: the integrands are singular there.
  const double lb = 1e-12;
  boost::math::quadrature::tanh_sinh<double> integrator;
  return integrator.integrate(f, lb, 1.0 - lb);
}

//! @brief Integrates `f` over the unit interval, splitting it at `split`.
//!
//! A narrow interior peak is invisible to a rule whose abscissas cluster at the
//! endpoints. Splitting there makes the peak an endpoint of both subintervals,
//! which is exactly where tanh-sinh puts its resolution. Needed for Tawn's tau
//! integral, whose peak is both narrow and far from `1/2` when the parameters
//! are asymmetric.
//! @brief Integrates `f` over the unit interval, splitting it at every
//! interior point of `splits`, with Gauss-Legendre panels graded toward the
//! panel ends.
//!
//! For integrands with steep but finite features at known locations, such as
//! the h-functions of a concentrated circular copula. Each panel between two
//! consecutive split points is divided into sub-panels whose widths halve
//! toward both ends, so a feature sitting at a split point is resolved down to
//! \f$ 2^{-9} \f$ of the panel width, and a fixed 15-point rule is applied on
//! every sub-panel. The cost is bounded and independent of the integrand,
//! unlike an adaptive rule, which refines without limit on steep features.
//! Panels narrower than `1e-7` are skipped; for an integrand bounded by one
//! their contribution is below the accuracy of the rule.
template<typename F>
inline double
integrate_zero_to_one_splits(F&& f, std::vector<double> splits)
{
  using rule = boost::math::quadrature::gauss<double, 15>;
  const int levels = 9;
  auto panel = [&](double a, double b) {
    const double half = 0.5 * (b - a);
    const double mid = a + half;
    double total = 0.0;
    // sub-panels [a + half / 2^k, a + half / 2^(k-1)] and their mirror images
    double lo = a;
    for (int k = levels; k >= 1; --k) {
      const double hi = a + half / std::pow(2.0, k - 1);
      total += rule::integrate(f, lo, hi) +
               rule::integrate(f, b - (hi - a), b - (lo - a));
      lo = hi;
    }
    (void)mid;
    return total;
  };

  const double min_width = 1e-7;
  std::sort(splits.begin(), splits.end());
  double total = 0.0;
  double a = 0.0;
  for (double s : splits) {
    if (s - a > min_width && 1.0 - s > min_width) {
      total += panel(a, s);
      a = s;
    }
  }
  return total + panel(a, 1.0);
}

template<typename F>
inline double
integrate_zero_to_one_split(F&& f, const double split)
{
  const double lb = 1e-12;
  const double ub = 1.0 - lb;
  if (!(split > lb) || !(split < ub)) {
    return integrate_zero_to_one(std::forward<F>(f));
  }
  boost::math::quadrature::tanh_sinh<double> integrator;
  return integrator.integrate(f, lb, split) +
         integrator.integrate(f, split, ub);
}
}
}
