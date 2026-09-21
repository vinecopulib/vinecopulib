// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_var_types.hpp>

namespace vinecopulib {

inline bool
SectionsBicop::circular_second() const
{
  return tools_var_types::is_circular(var_types_[1]) &&
         !tools_var_types::is_circular(var_types_[0]);
}

//! \f$ p(v) = a (1 - v)(1 - 3v) + b\, v (2 - 3v) \f$.
inline double
SectionsBicop::p_of(double a, double b, double v)
{
  return a * (1.0 - v) * (1.0 - 3.0 * v) + b * v * (2.0 - 3.0 * v);
}

//! the antiderivative of `p_of` with \f$ P(0) = P(1) = 0 \f$.
inline double
SectionsBicop::P_of(double a, double b, double v)
{
  return a * v * (1.0 - v) * (1.0 - v) + b * v * v * (1.0 - v);
}

inline double
SectionsBicop::pdf_internal(double u,
                            double v,
                            const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  const auto ab = amplitudes(par);
  const double mu = par(par.size() - 1);
  return 1.0 + std::cos(tools_circular::two_pi() * u - mu) *
                 p_of(ab.first, ab.second, v);
}

inline double
SectionsBicop::cdf_internal(double u,
                            double v,
                            const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  const auto ab = amplitudes(par);
  const double mu = par(par.size() - 1);
  const double two_pi = tools_circular::two_pi();
  return u * v + (std::sin(two_pi * u - mu) + std::sin(mu)) / two_pi *
                   P_of(ab.first, ab.second, v);
}

//! \f$ h_1(v | u) = \partial C / \partial u \f$.
inline double
SectionsBicop::hfunc1_internal(
  double u,
  double v,
  const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  const auto ab = amplitudes(par);
  const double mu = par(par.size() - 1);
  return v + std::cos(tools_circular::two_pi() * u - mu) *
               P_of(ab.first, ab.second, v);
}

//! \f$ h_2(u | v) = \partial C / \partial v \f$.
inline double
SectionsBicop::hfunc2_internal(
  double u,
  double v,
  const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  const auto ab = amplitudes(par);
  const double mu = par(par.size() - 1);
  const double two_pi = tools_circular::two_pi();
  return u + (std::sin(two_pi * u - mu) + std::sin(mu)) / two_pi *
               p_of(ab.first, ab.second, v);
}

//! both inverses are bracketed root solves on the unit interval; the
//! derivative of each h-function in its own argument is the density.
inline double
SectionsBicop::hinv1_internal(
  double u,
  double w,
  const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  return tools_circular::invert_increasing(
    [&](double v) { return hfunc1_internal(u, v, par); },
    [&](double v) { return pdf_internal(u, v, par); },
    w,
    0.0,
    1.0);
}

inline double
SectionsBicop::hinv2_internal(
  double w,
  double v,
  const Eigen::Ref<const Eigen::VectorXd>& par) const
{
  return tools_circular::invert_increasing(
    [&](double u) { return hfunc2_internal(u, v, par); },
    [&](double u) { return pdf_internal(u, v, par); },
    w,
    0.0,
    1.0);
}

inline Eigen::VectorXd
SectionsBicop::pdf_raw(const Eigen::MatrixXd& u,
                       const Eigen::MatrixXd& parameters)
{
  const bool swap = circular_second();
  auto f = [this, swap](const double& u1,
                        const double& u2,
                        const Eigen::Ref<const Eigen::VectorXd>& par) {
    return swap ? pdf_internal(u2, u1, par) : pdf_internal(u1, u2, par);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
SectionsBicop::cdf(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  const bool swap = circular_second();
  auto f = [this, swap](const double& u1,
                        const double& u2,
                        const Eigen::Ref<const Eigen::VectorXd>& par) {
    return swap ? cdf_internal(u2, u1, par) : cdf_internal(u1, u2, par);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
SectionsBicop::hfunc1_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters)
{
  const bool swap = circular_second();
  auto f = [this, swap](const double& u1,
                        const double& u2,
                        const Eigen::Ref<const Eigen::VectorXd>& par) {
    return swap ? hfunc2_internal(u2, u1, par) : hfunc1_internal(u1, u2, par);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
SectionsBicop::hfunc2_raw(const Eigen::MatrixXd& u,
                          const Eigen::MatrixXd& parameters)
{
  const bool swap = circular_second();
  auto f = [this, swap](const double& u1,
                        const double& u2,
                        const Eigen::Ref<const Eigen::VectorXd>& par) {
    return swap ? hfunc1_internal(u2, u1, par) : hfunc2_internal(u1, u2, par);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
SectionsBicop::hinv1_raw(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters)
{
  const bool swap = circular_second();
  auto f = [this, swap](const double& u1,
                        const double& w,
                        const Eigen::Ref<const Eigen::VectorXd>& par) {
    return swap ? hinv2_internal(w, u1, par) : hinv1_internal(u1, w, par);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
SectionsBicop::hinv2_raw(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters)
{
  const bool swap = circular_second();
  auto f = [this, swap](const double& w,
                        const double& u2,
                        const Eigen::Ref<const Eigen::VectorXd>& par) {
    return swap ? hinv1_internal(u2, w, par) : hinv2_internal(w, u2, par);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

//! \f$ \tau = (a + b) \sin\mu / (3\pi) \f$.
inline double
SectionsBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  const Eigen::VectorXd par = (parameters.cols() == 1)
                                ? Eigen::VectorXd(parameters)
                                : Eigen::VectorXd(parameters.row(0));
  const auto ab = amplitudes(par);
  const double mu = par(par.size() - 1);
  return (ab.first + ab.second) * std::sin(mu) /
         (3.0 * boost::math::constants::pi<double>());
}

//! the argument order is read from the variable types, which `Bicop::flip()`
//! exchanges; the parameters do not change.
inline void
SectionsBicop::flip()
{
}

//! moment estimates from \f$ E[e^{i 2\pi U} p_k(V)] \f$ with
//! \f$ p_1(v) = (1 - v)(1 - 3v) \f$ and \f$ p_2(v) = v(2 - 3v) \f$, which equal
//! \f$ (a/15 + b/60) e^{i\mu} \f$ and \f$ (a/60 + b/15) e^{i\mu} \f$.
inline Eigen::VectorXd
SectionsBicop::moment_start(const Eigen::MatrixXd& data,
                            const Eigen::VectorXd& weights)
{
  const Eigen::Index ic = circular_second() ? 1 : 0;
  const Eigen::VectorXd u = data.col(ic);
  const Eigen::VectorXd v = data.col(1 - ic);
  const double two_pi = tools_circular::two_pi();

  std::complex<double> m1(0.0, 0.0), m2(0.0, 0.0);
  double wsum = 0.0;
  for (Eigen::Index i = 0; i < u.size(); ++i) {
    const double w = (weights.size() > 0) ? weights(i) : 1.0;
    const std::complex<double> z(std::cos(two_pi * u(i)),
                                 std::sin(two_pi * u(i)));
    m1 += w * z * ((1.0 - v(i)) * (1.0 - 3.0 * v(i)));
    m2 += w * z * (v(i) * (2.0 - 3.0 * v(i)));
    wsum += w;
  }
  m1 /= wsum;
  m2 /= wsum;
  double mu = std::arg(m1 + m2);
  const std::complex<double> rot(std::cos(-mu), std::sin(-mu));
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
  return parameters_from_moments(a, b, mu);
}
}
