// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <complex>
#include <vector>

namespace vinecopulib {

//! @brief Numerical helpers shared by the circular copula families.
namespace tools_circular {

//! \f$ 2\pi \f$.
inline double
two_pi()
{
  return boost::math::constants::two_pi<double>();
}

//! reduces an angle to \f$ [-\pi, \pi) \f$.
inline double
wrap_pi(double theta)
{
  const double t =
    std::fmod(theta + boost::math::constants::pi<double>(), two_pi());
  return (t < 0 ? t + two_pi() : t) - boost::math::constants::pi<double>();
}

//! @brief Solves \f$ f(\theta) = \text{target} \f$ for an increasing `f`
//! with derivative `df` by safeguarded Newton steps inside `[lo, hi]`.
//!
//! The bracket is halved whenever a Newton step leaves it or does not shrink
//! the residual enough, so the iteration converges even where `df` vanishes
//! (the cardioid at its concentration bound).
template<typename F, typename D>
inline double
invert_increasing(const F& f,
                  const D& df,
                  double target,
                  double lo,
                  double hi,
                  double tol = 1e-12,
                  int max_iter = 100)
{
  double x = 0.5 * (lo + hi);
  double fx = f(x) - target;
  for (int it = 0; it < max_iter; ++it) {
    if (std::abs(fx) < tol) {
      return x;
    }
    if (fx > 0) {
      hi = x;
    } else {
      lo = x;
    }
    const double d = df(x);
    double x_new = (d > 0) ? x - fx / d : 0.5 * (lo + hi);
    if (!(x_new > lo && x_new < hi)) {
      x_new = 0.5 * (lo + hi);
    }
    x = x_new;
    fx = f(x) - target;
    if (hi - lo < tol) {
      return x;
    }
  }
  return x;
}

//! @brief The weighted mean of `factor_i * exp(i theta_i)`.
//!
//! @param theta Angles in radians.
//! @param weights Observation weights; empty means equal weights.
//! @param factor A real factor per observation; empty means one. For a plain
//!   resultant, the modulus is the mean resultant length and the argument the
//!   mean direction.
//! @return The weighted mean, or zero when the total weight is not positive.
inline std::complex<double>
weighted_resultant(const Eigen::VectorXd& theta,
                   const Eigen::VectorXd& weights,
                   const Eigen::VectorXd& factor = Eigen::VectorXd())
{
  std::complex<double> z(0.0, 0.0);
  double wsum = 0.0;
  for (Eigen::Index i = 0; i < theta.size(); ++i) {
    const double w = (weights.size() > 0) ? weights(i) : 1.0;
    const double f = (factor.size() > 0) ? factor(i) : 1.0;
    z += w * f * std::polar(1.0, theta(i));
    wsum += w;
  }
  return (wsum > 0.0) ? z / wsum : std::complex<double>(0.0, 0.0);
}

//! @brief The number of terms after which the Fourier series of the von
//! Mises CDF is below `tol`; the coefficients decay like
//! \f$ \exp(-j^2 / 2\kappa) \f$.
inline size_t
von_mises_series_length(double kappa, double tol = 1e-12)
{
  const double n = std::sqrt(2.0 * std::max(kappa, 1.0) * -std::log(tol));
  return static_cast<size_t>(std::ceil(n)) + 10;
}

//! @brief The Bessel ratios \f$ I_j(\kappa) / I_0(\kappa) \f$ for
//! \f$ j = 1, \dots, n \f$.
inline std::vector<double>
bessel_i_ratios(double kappa, size_t n)
{
  std::vector<double> ratios(n);
  if (kappa <= 0.0) {
    std::fill(ratios.begin(), ratios.end(), 0.0);
    return ratios;
  }
  const double i0 = boost::math::cyl_bessel_i(0, kappa);
  for (size_t j = 1; j <= n; ++j) {
    ratios[j - 1] = boost::math::cyl_bessel_i(static_cast<int>(j), kappa) / i0;
  }
  return ratios;
}

//! @brief The mean resultant length \f$ A(\kappa) = I_1(\kappa) / I_0(\kappa)
//! \f$ of the von Mises distribution.
inline double
von_mises_a(double kappa)
{
  if (kappa <= 0.0) {
    return 0.0;
  }
  return boost::math::cyl_bessel_i(1, kappa) /
         boost::math::cyl_bessel_i(0, kappa);
}

//! @brief The inverse of `von_mises_a()` on `[0, upper]`, clipped to it.
inline double
von_mises_a_inverse(double rbar, double upper)
{
  if (rbar <= 0.0) {
    return 0.0;
  }
  if (von_mises_a(upper) <= rbar) {
    return upper;
  }
  double lo = 0.0, hi = upper;
  for (int it = 0; it < 200 && hi - lo > 1e-12 * std::max(1.0, hi); ++it) {
    const double mid = 0.5 * (lo + hi);
    (von_mises_a(mid) < rbar ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

} // namespace tools_circular
} // namespace vinecopulib
