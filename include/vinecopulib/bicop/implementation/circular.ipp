// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_integration.hpp>

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

//! Kendall's tau of the copula with the cut at zero,
//! \f$ \tau = 1 - 4 \int_0^1 \int_0^1 h_1(v | u) h_2(u | v) \, du \, dv \f$,
//! by nested quadrature of the h-function leaves.
inline double
CircularBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  const Eigen::MatrixXd par_row =
    (parameters.cols() == 1) ? parameters.transpose() : parameters;
  Eigen::MatrixXd uv(1, 2);
  auto inner = [&](double u) {
    auto f = [&](double v) {
      uv(0, 0) = u;
      uv(0, 1) = v;
      return hfunc1_raw(uv, par_row)(0) * hfunc2_raw(uv, par_row)(0);
    };
    return tools_integration::integrate_zero_to_one(f);
  };
  return 1.0 - 4.0 * tools_integration::integrate_zero_to_one(inner);
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

//! unused: the fit starts from `moment_start()`.
inline Eigen::VectorXd
CircularBicop::get_start_parameters(const double)
{
  return parameters_;
}
}
