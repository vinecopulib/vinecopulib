// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

namespace vinecopulib {

inline QuadSectionsBicop::QuadSectionsBicop()
{
  family_ = BicopFamily::quad_sections;
  parameters_ = Eigen::VectorXd(2);
  parameters_lower_bounds_ = Eigen::VectorXd(2);
  parameters_upper_bounds_ = Eigen::VectorXd(2);
  const double inf = std::numeric_limits<double>::infinity();
  parameters_ << 0, 0;
  parameters_lower_bounds_ << 0, -inf;
  parameters_upper_bounds_ << 1, inf;
}

//! quadratic sections are the cubic family with equal amplitudes.
inline std::pair<double, double>
QuadSectionsBicop::amplitudes(
  const Eigen::Ref<const Eigen::VectorXd>& parameters) const
{
  return { parameters(0), parameters(0) };
}

//! the two moment amplitudes estimate the same quantity; average them.
inline Eigen::VectorXd
QuadSectionsBicop::parameters_from_moments(double a, double b, double mu) const
{
  Eigen::VectorXd start(2);
  start << std::min(std::max(0.5 * (a + b), 0.0), 1.0), mu;
  return start;
}
}
