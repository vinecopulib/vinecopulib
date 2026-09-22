// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <limits>

namespace vinecopulib {

inline Eigen::Index
CircularBicop::phase_index() const
{
  return parameters_.size() - 1;
}

inline void
CircularBicop::fit(const Eigen::MatrixXd& data,
                   std::string method,
                   double,
                   size_t,
                   const Eigen::VectorXd& weights)
{
  check_fit_method(method);
  if (method != "mle") {
    throw std::runtime_error("only method 'mle' is available for the " +
                             get_family_name() + " copula");
  }
  Eigen::VectorXd start = moment_start(data, weights);
  start = start.cwiseMax(parameters_lower_bounds_);
  start = start.cwiseMin(parameters_upper_bounds_);
  fit_mle(data,
          weights,
          method,
          start,
          parameters_lower_bounds_,
          parameters_upper_bounds_);
  parameters_(phase_index()) =
    tools_circular::wrap_pi(parameters_(phase_index()));
}

//! the densities are bounded, so no corner has tail dependence.
inline Eigen::MatrixXd
CircularBicop::parameters_to_taildep(const Eigen::MatrixXd&)
{
  return Eigen::MatrixXd::Zero(2, 2);
}

inline Eigen::MatrixXd
CircularBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}

//! Kendall's tau of a pair with a circular variable depends on where the
//! circle is cut, so it is not reported.
inline double
CircularBicop::parameters_to_tau(const Eigen::MatrixXd&)
{
  return std::numeric_limits<double>::quiet_NaN();
}

//! unused: the fit starts from `moment_start()`.
inline Eigen::VectorXd
CircularBicop::get_start_parameters(const double)
{
  return parameters_;
}
}
