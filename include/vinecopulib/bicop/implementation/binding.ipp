// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <tuple>
#include <vinecopulib/misc/tools_eigen.hpp>

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

//! \f$ C(u, v) = \int_0^u h_1(v | s) \, ds \f$ in closed form: with
//! \f$ K \f$ an antiderivative of the lifted CDF,
//! \f$ 2\pi C(u, v) = K(2\pi v - \mu) - K(2\pi (v - u) - \mu) - K(-\mu)
//! + K(-2\pi u - \mu) \f$, and the Fourier series of \f$ g \f$ gives
//! \f$ K(\theta) = \theta^2 / 4\pi - \sum_j \rho_j \cos(j\theta) / (\pi j^2)
//! \f$.
inline Eigen::VectorXd
BindingBicop::cdf(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  const double two_pi = tools_circular::two_pi();
  const double pi = boost::math::constants::pi<double>();
  // the coefficients of the last (family, concentration) seen on this thread
  thread_local std::tuple<BicopFamily, double, std::vector<double>> series{
    BicopFamily::indep, -1.0, {}
  };
  auto f = [&](const double& u1,
               const double& u2,
               const Eigen::Ref<const Eigen::VectorXd>& par) {
    if (std::get<0>(series) != family_ || std::get<1>(series) != par(0)) {
      series = { family_, par(0), fourier_coefficients(par(0)) };
    }
    const std::vector<double>& coefs = std::get<2>(series);
    auto K = [&](double theta) {
      double sum = 0.0;
      for (size_t j = 1; j <= coefs.size(); ++j) {
        const double jj = static_cast<double>(j);
        sum += coefs[j - 1] * std::cos(jj * theta) / (jj * jj);
      }
      return theta * theta / (4.0 * pi) - sum / pi;
    };
    const double mu = par(1);
    const double c = (K(two_pi * u2 - mu) - K(two_pi * (u2 - u1) - mu) -
                      K(-mu) + K(-two_pi * u1 - mu)) /
                     two_pi;
    return std::min(std::max(c, 0.0), 1.0);
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
  const auto resultant = tools_circular::weighted_resultant(theta, weights);
  Eigen::VectorXd start(2);
  start << concentration_from_resultant(std::abs(resultant)),
    std::arg(resultant);
  return start;
}
}
