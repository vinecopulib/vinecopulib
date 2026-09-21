// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

namespace vinecopulib {

inline CubicSectionsBicop::CubicSectionsBicop()
{
  family_ = BicopFamily::cubic_sections;
  parameters_ = Eigen::VectorXd(3);
  parameters_lower_bounds_ = Eigen::VectorXd(3);
  parameters_upper_bounds_ = Eigen::VectorXd(3);
  const double inf = std::numeric_limits<double>::infinity();
  parameters_ << 0, 0, 0;
  parameters_lower_bounds_ << 0, -1, -inf;
  parameters_upper_bounds_ << 1, 1, inf;
}

inline std::pair<double, double>
CubicSectionsBicop::amplitudes(
  const Eigen::Ref<const Eigen::VectorXd>& parameters) const
{
  return { parameters(0), parameters(1) };
}

inline Eigen::VectorXd
CubicSectionsBicop::parameters_from_moments(double a, double b, double mu) const
{
  Eigen::VectorXd start(3);
  start << std::min(std::max(a, 0.0), 1.0), std::min(std::max(b, -1.0), 1.0),
    mu;
  return start;
}
}
