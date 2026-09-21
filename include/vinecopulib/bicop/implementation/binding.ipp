// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_integration.hpp>

namespace vinecopulib {

inline Eigen::VectorXd
BindingBicop::pdf_raw(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters)
{
  const double two_pi = tools_circular::two_pi();
  auto f = [this, two_pi](const double& u1,
                          const double& u2,
                          const Eigen::Ref<const Eigen::VectorXd>& par) {
    return two_pi * g(two_pi * (u2 - u1) - par(1), par(0));
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

//! \f$ C(u, v) = \int_0^u h_1(v | s) \, ds \f$ by quadrature, split where
//! \f$ h_1 \f$ is steep in \f$ s \f$ (the density peak at
//! \f$ v - s \equiv \mu / 2\pi \f$); the CDF has no closed form for every
//! family and is not needed in the fitting or simulation paths.
inline Eigen::VectorXd
BindingBicop::cdf(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  const double two_pi = tools_circular::two_pi();
  auto f = [this, two_pi](const double& u1,
                          const double& u2,
                          const Eigen::Ref<const Eigen::VectorXd>& par) {
    auto h1 = [&](double s) {
      return lifted_cdf(two_pi * (u2 - s) - par(1), par(0)) -
             lifted_cdf(-two_pi * s - par(1), par(0));
    };
    // rescale [0, u1] to the unit interval; h1 is steep in s where either
    // lifted-CDF argument crosses the peak of g, at s = u2 - mu/2pi and at
    // s = -mu/2pi (mod 1)
    auto frac = [](double x) { return x - std::floor(x); };
    const double shift = par(1) / two_pi;
    return u1 * tools_integration::integrate_zero_to_one_splits(
                  [&](double t) { return h1(u1 * t); },
                  { frac(u2 - shift) / u1, frac(-shift) / u1 });
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
BindingBicop::hfunc1_raw(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters)
{
  const double two_pi = tools_circular::two_pi();
  auto f = [this, two_pi](const double& u1,
                          const double& u2,
                          const Eigen::Ref<const Eigen::VectorXd>& par) {
    return lifted_cdf(two_pi * (u2 - u1) - par(1), par(0)) -
           lifted_cdf(-two_pi * u1 - par(1), par(0));
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
BindingBicop::hfunc2_raw(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters)
{
  const double two_pi = tools_circular::two_pi();
  auto f = [this, two_pi](const double& u1,
                          const double& u2,
                          const Eigen::Ref<const Eigen::VectorXd>& par) {
    return lifted_cdf(two_pi * u2 - par(1), par(0)) -
           lifted_cdf(two_pi * (u2 - u1) - par(1), par(0));
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
BindingBicop::hinv1_raw(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters)
{
  const double two_pi = tools_circular::two_pi();
  auto f = [this, two_pi](const double& u1,
                          const double& w,
                          const Eigen::Ref<const Eigen::VectorXd>& par) {
    const double target = w + lifted_cdf(-two_pi * u1 - par(1), par(0));
    const double v =
      u1 + (par(1) + lifted_cdf_inverse(target, par(0))) / two_pi;
    return std::min(std::max(v, 0.0), 1.0);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
BindingBicop::hinv2_raw(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters)
{
  const double two_pi = tools_circular::two_pi();
  auto f = [this, two_pi](const double& w,
                          const double& u2,
                          const Eigen::Ref<const Eigen::VectorXd>& par) {
    const double target = lifted_cdf(two_pi * u2 - par(1), par(0)) - w;
    const double v =
      u2 - (par(1) + lifted_cdf_inverse(target, par(0))) / two_pi;
    return std::min(std::max(v, 0.0), 1.0);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

//! Kendall's tau by nested quadrature, splitting each integral where the
//! h-functions are steep: the density concentrates along
//! \f$ v - u \equiv \mu / 2\pi \f$, and the h-functions jump where their
//! lifted-CDF arguments cross that peak.
inline double
BindingBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  const Eigen::MatrixXd par_row =
    (parameters.cols() == 1) ? parameters.transpose() : parameters;
  const double two_pi = tools_circular::two_pi();
  const double shift = par_row(0, 1) / two_pi; // the peak of g in turns
  auto frac = [](double x) { return x - std::floor(x); };

  Eigen::MatrixXd uv(1, 2);
  auto inner = [&](double u) {
    auto f = [&](double v) {
      uv(0, 0) = u;
      uv(0, 1) = v;
      return hfunc1_raw(uv, par_row)(0) * hfunc2_raw(uv, par_row)(0);
    };
    // steep where v - u or v sits at the peak
    return tools_integration::integrate_zero_to_one_splits(
      f, { frac(u + shift), frac(shift) });
  };
  // steep where -u sits at the peak
  const double integral =
    tools_integration::integrate_zero_to_one_splits(inner, { frac(-shift) });
  return 1.0 - 4.0 * integral;
}

//! exchanging the arguments negates the phase.
inline void
BindingBicop::flip()
{
  parameters_(1) = -parameters_(1);
}

//! the residual angles \f$ 2\pi(v - u) \f$ have density \f$ g(\cdot - \mu) \f$,
//! so the mean resultant gives the phase and the concentration.
inline Eigen::VectorXd
BindingBicop::moment_start(const Eigen::MatrixXd& data,
                           const Eigen::VectorXd& weights)
{
  Eigen::VectorXd theta =
    tools_circular::two_pi() * (data.col(1) - data.col(0));
  auto resultant = tools_circular::mean_resultant(theta, weights);
  Eigen::VectorXd start(2);
  start << concentration_from_resultant(resultant.second), resultant.first;
  return start;
}
}
