// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_var_types.hpp>

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

inline double
CubicSectionsBicop::p_of(double a, double b, double v)
{
  return a * (1.0 - v) * (1.0 - 3.0 * v) + b * v * (2.0 - 3.0 * v);
}

inline double
CubicSectionsBicop::P_of(double a, double b, double v)
{
  return a * v * (1.0 - v) * (1.0 - v) + b * v * v * (1.0 - v);
}

inline double
CubicSectionsBicop::pdf_at(double u,
                           double v,
                           const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  return 1.0 + std::cos(tools_circular::two_pi() * u - par(2)) *
                 p_of(par(0), par(1), v);
}

inline double
CubicSectionsBicop::cdf_at(double u,
                           double v,
                           const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  const double two_pi = tools_circular::two_pi();
  return u * v + (std::sin(two_pi * u - par(2)) + std::sin(par(2))) / two_pi *
                   P_of(par(0), par(1), v);
}

inline double
CubicSectionsBicop::hfunc1_at(
  double u,
  double v,
  const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  return v + std::cos(tools_circular::two_pi() * u - par(2)) *
               P_of(par(0), par(1), v);
}

inline double
CubicSectionsBicop::hfunc2_at(
  double u,
  double v,
  const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  const double two_pi = tools_circular::two_pi();
  return u + (std::sin(two_pi * u - par(2)) + std::sin(par(2))) / two_pi *
               p_of(par(0), par(1), v);
}

inline double
CubicSectionsBicop::hinv1_at(double u,
                             double w,
                             const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  return tools_circular::invert_increasing(
    [&](double v) { return hfunc1_at(u, v, par); },
    [&](double v) { return pdf_at(u, v, par); },
    w,
    0.0,
    1.0);
}

inline double
CubicSectionsBicop::hinv2_at(double w,
                             double v,
                             const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  return tools_circular::invert_increasing(
    [&](double u) { return hfunc2_at(u, v, par); },
    [&](double u) { return pdf_at(u, v, par); },
    w,
    0.0,
    1.0);
}

inline Eigen::VectorXd
CubicSectionsBicop::eval(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters,
                         Leaf leaf,
                         Leaf swapped) const
{
  const bool swap = tools_var_types::is_circular(var_types_[1]) &&
                    !tools_var_types::is_circular(var_types_[0]);
  const Leaf fn = swap ? swapped : leaf;
  auto f = [this, fn](const double& x,
                      const double& y,
                      const Eigen::Ref<const Eigen::VectorXd>& par) {
    return (this->*fn)(x, y, par);
  };
  if (!swap) {
    return tools_eigen::binaryExpr_or_nan(u, parameters, f);
  }
  return tools_eigen::binaryExpr_or_nan(
    tools_eigen::swap_cols(u), parameters, f);
}

inline Eigen::VectorXd
CubicSectionsBicop::pdf_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters)
{
  return eval(
    u, parameters, &CubicSectionsBicop::pdf_at, &CubicSectionsBicop::pdf_at);
}

inline Eigen::VectorXd
CubicSectionsBicop::cdf(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters)
{
  return eval(
    u, parameters, &CubicSectionsBicop::cdf_at, &CubicSectionsBicop::cdf_at);
}

inline Eigen::VectorXd
CubicSectionsBicop::hfunc1_raw(const Eigen::MatrixXd& u,
                               const Eigen::MatrixXd& parameters)
{
  return eval(u,
              parameters,
              &CubicSectionsBicop::hfunc1_at,
              &CubicSectionsBicop::hfunc2_at);
}

inline Eigen::VectorXd
CubicSectionsBicop::hfunc2_raw(const Eigen::MatrixXd& u,
                               const Eigen::MatrixXd& parameters)
{
  return eval(u,
              parameters,
              &CubicSectionsBicop::hfunc2_at,
              &CubicSectionsBicop::hfunc1_at);
}

inline Eigen::VectorXd
CubicSectionsBicop::hinv1_raw(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& parameters)
{
  return eval(u,
              parameters,
              &CubicSectionsBicop::hinv1_at,
              &CubicSectionsBicop::hinv2_at);
}

inline Eigen::VectorXd
CubicSectionsBicop::hinv2_raw(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& parameters)
{
  return eval(u,
              parameters,
              &CubicSectionsBicop::hinv2_at,
              &CubicSectionsBicop::hinv1_at);
}

//! the leaves read the argument order from `var_types`, so the parameters
//! are unchanged.
inline void
CubicSectionsBicop::flip()
{
}

//! moment estimates from E[exp(i 2 pi U) p_k(V)] = (a / 15 + b / 60) exp(i mu)
//! and (a / 60 + b / 15) exp(i mu) for p_1(v) = (1 - v)(1 - 3v) and
//! p_2(v) = v(2 - 3v).
inline Eigen::VectorXd
CubicSectionsBicop::moment_start(const Eigen::MatrixXd& data,
                                 const Eigen::VectorXd& weights)
{
  const Eigen::Index ic = tools_var_types::is_circular(var_types_[0]) ? 0 : 1;
  const Eigen::VectorXd theta = tools_circular::two_pi() * data.col(ic);
  const Eigen::ArrayXd v = data.col(1 - ic).array();
  const auto m1 = tools_circular::weighted_resultant(
    theta, weights, ((1.0 - v) * (1.0 - 3.0 * v)).matrix());
  const auto m2 = tools_circular::weighted_resultant(
    theta, weights, (v * (2.0 - 3.0 * v)).matrix());
  double mu = std::arg(m1 + m2);
  const std::complex<double> rot = std::polar(1.0, -mu);
  const double r1 = std::real(m1 * rot);
  const double r2 = std::real(m2 * rot);
  // invert the 2 x 2 system [[1/15, 1/60], [1/60, 1/15]] (a, b)' = (r1, r2)'
  const double det = 1.0 / 225.0 - 1.0 / 3600.0;
  double a = (r1 / 15.0 - r2 / 60.0) / det;
  double b = (r2 / 15.0 - r1 / 60.0) / det;
  if (a < 0) {
    // (a, b, mu) and (-a, -b, mu + pi) define the same copula
    a = -a;
    b = -b;
    mu = tools_circular::wrap_pi(mu + boost::math::constants::pi<double>());
  }
  Eigen::VectorXd start(3);
  start << std::min(std::max(a, 0.0), 1.0), std::min(std::max(b, -1.0), 1.0),
    mu;
  return start;
}
}
