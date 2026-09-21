// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

namespace vinecopulib {

inline CardioidBicop::CardioidBicop()
{
  family_ = BicopFamily::cardioid;
  parameters_ = Eigen::VectorXd(2);
  parameters_lower_bounds_ = Eigen::VectorXd(2);
  parameters_upper_bounds_ = Eigen::VectorXd(2);
  const double inf = std::numeric_limits<double>::infinity();
  parameters_ << 0, 0;
  parameters_lower_bounds_ << 0, -inf;
  parameters_upper_bounds_ << 0.5, inf;
}

inline double
CardioidBicop::g(double theta, double concentration) const
{
  return (1.0 + 2.0 * concentration * std::cos(theta)) /
         tools_circular::two_pi();
}

inline double
CardioidBicop::lifted_cdf(double theta, double concentration) const
{
  return (theta + 2.0 * concentration * std::sin(theta)) /
         tools_circular::two_pi();
}

//! solves the Kepler-type equation \f$ \theta + 2\rho\sin\theta = 2\pi w \f$;
//! the root lies within \f$ 2\rho \f$ of \f$ 2\pi w \f$.
inline double
CardioidBicop::lifted_cdf_inverse(double w, double concentration) const
{
  const double center = tools_circular::two_pi() * w;
  const double half_width = 2.0 * concentration + 1e-9;
  return tools_circular::invert_increasing(
    [&](double t) { return lifted_cdf(t, concentration); },
    [&](double t) { return g(t, concentration); },
    w,
    center - half_width,
    center + half_width);
}

//! the mean resultant length of the cardioid is its concentration.
inline double
CardioidBicop::concentration_from_resultant(double rbar) const
{
  return std::min(std::max(rbar, 0.0), parameters_upper_bounds_(0));
}
}
