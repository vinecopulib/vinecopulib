// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

namespace vinecopulib {

inline WrappedCauchyBicop::WrappedCauchyBicop()
{
  family_ = BicopFamily::wrapped_cauchy;
  parameters_ = Eigen::VectorXd(2);
  parameters_lower_bounds_ = Eigen::VectorXd(2);
  parameters_upper_bounds_ = Eigen::VectorXd(2);
  const double inf = std::numeric_limits<double>::infinity();
  parameters_ << 0, 0;
  parameters_lower_bounds_ << 0, -inf;
  parameters_upper_bounds_ << 0.99, inf;
}

inline double
WrappedCauchyBicop::g(double theta, double concentration) const
{
  const double r = concentration;
  return (1.0 - r * r) /
         (tools_circular::two_pi() * (1.0 + r * r - 2.0 * r * std::cos(theta)));
}

//! the arctangent term is continuous and bounded for \f$ \rho < 1 \f$, so the
//! expression is a valid lift.
inline double
WrappedCauchyBicop::lifted_cdf(double theta, double concentration) const
{
  const double r = concentration;
  return (theta +
          2.0 * std::atan(r * std::sin(theta) / (1.0 - r * std::cos(theta)))) /
         tools_circular::two_pi();
}

//! closed form within one turn, with the number of completed turns restored.
inline double
WrappedCauchyBicop::lifted_cdf_inverse(double w, double concentration) const
{
  const double r = concentration;
  const double k = std::round(w);
  const double pi = boost::math::constants::pi<double>();
  return tools_circular::two_pi() * k +
         2.0 * std::atan((1.0 - r) / (1.0 + r) * std::tan(pi * (w - k)));
}

//! the mean resultant length of the wrapped Cauchy is its concentration.
inline double
WrappedCauchyBicop::concentration_from_resultant(double rbar) const
{
  return std::min(std::max(rbar, 0.0), parameters_upper_bounds_(0));
}
}
