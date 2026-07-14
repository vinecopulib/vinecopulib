// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_integration.hpp>

namespace vinecopulib {
inline Eigen::VectorXd
ExtremeValueBicop::cdf(const tools_eigen::ConstMatRef& u,
                       const tools_eigen::ConstMatRef& parameters)
{
  if (parameters.rows() == 1) {
    Eigen::ArrayXd l1 = u.col(0).array().log();
    Eigen::ArrayXd l2 = u.col(1).array().log();
    Eigen::ArrayXd ls = l1 + l2;
    Eigen::ArrayXd t = l2 / ls;
    Eigen::ArrayXd a = pickands_arr(t, parameters.row(0).transpose());
    return (ls * a).exp().matrix();
  }
  auto f = [this](const double& u1,
                  const double& u2,
                  const Eigen::Ref<const Eigen::VectorXd>& par) {
    double t = std::log(u2) / std::log(u1 * u2);
    t = pickands(t, par);
    t = (std::log(u1) + std::log(u2)) * t;
    return std::exp(t);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

//! default fused Pickands evaluation: the three separate virtual calls
inline void
ExtremeValueBicop::pickands_all(
  const double& t,
  const Eigen::Ref<const Eigen::VectorXd>& parameters,
  double& A,
  double& A1,
  double& A2)
{
  A = pickands(t, parameters);
  A1 = pickands_derivative(t, parameters);
  A2 = pickands_derivative2(t, parameters);
}

//! default array-valued Pickands evaluation: a scalar loop
inline Eigen::ArrayXd
ExtremeValueBicop::pickands_arr(
  const Eigen::ArrayXd& t,
  const Eigen::Ref<const Eigen::VectorXd>& parameters)
{
  Eigen::ArrayXd a(t.size());
  for (Eigen::Index i = 0; i < t.size(); ++i) {
    a(i) = pickands(t(i), parameters);
  }
  return a;
}

//! default array-valued fused Pickands evaluation: a scalar loop; families
//! with vectorizable shared subexpressions override this
inline void
ExtremeValueBicop::pickands_all(
  const Eigen::ArrayXd& t,
  const Eigen::Ref<const Eigen::VectorXd>& parameters,
  Eigen::ArrayXd& A,
  Eigen::ArrayXd& A1,
  Eigen::ArrayXd& A2)
{
  const Eigen::Index n = t.size();
  A.resize(n);
  A1.resize(n);
  A2.resize(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    pickands_all(t(i), parameters, A(i), A1(i), A2(i));
  }
}

inline Eigen::VectorXd
ExtremeValueBicop::pdf_raw(const tools_eigen::ConstMatRef& u,
                           const tools_eigen::ConstMatRef& parameters)
{
  if (parameters.rows() == 1) {
    Eigen::ArrayXd l1 = u.col(0).array().log();
    Eigen::ArrayXd l2 = u.col(1).array().log();
    Eigen::ArrayXd ls = l1 + l2;
    Eigen::ArrayXd t = l2 / ls;
    Eigen::ArrayXd a, a1, a2;
    pickands_all(t, parameters.row(0).transpose(), a, a1, a2);
    Eigen::ArrayXd fac =
      a * a + (1.0 - 2.0 * t) * a1 * a - (1.0 - t) * t * (a1 * a1 + a2 / ls);
    return ((ls * a).exp() * fac / (u.col(0).array() * u.col(1).array()))
      .matrix();
  }
  auto f = [this](const double& u1,
                  const double& u2,
                  const Eigen::Ref<const Eigen::VectorXd>& par) {
    double t = std::log(u2) / std::log(u1 * u2);
    double t2, t3, t4;
    pickands_all(t, par, t2, t3, t4);

    t3 = std::pow(t2, 2) + (1 - 2 * t) * t3 * t2 -
         (1 - t) * t * (std::pow(t3, 2) + t4 / std::log(u1 * u2));
    t2 = (std::log(u1) + std::log(u2)) * t2;

    return std::exp(t2) * t3 / (u1 * u2);
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
ExtremeValueBicop::hfunc1_raw(const tools_eigen::ConstMatRef& u,
                              const tools_eigen::ConstMatRef& parameters)
{
  if (parameters.rows() == 1) {
    Eigen::ArrayXd l1 = u.col(0).array().log();
    Eigen::ArrayXd l2 = u.col(1).array().log();
    Eigen::ArrayXd ls = l1 + l2;
    Eigen::ArrayXd t = l2 / ls;
    Eigen::ArrayXd a, a1, a2;
    pickands_all(t, parameters.row(0).transpose(), a, a1, a2);
    return ((ls * a).exp() * (a - t * a1) / u.col(0).array()).matrix();
  }
  auto f = [this](const double& u1,
                  const double& u2,
                  const Eigen::Ref<const Eigen::VectorXd>& par) {
    double t = std::log(u2) / std::log(u1 * u2);
    double t2 = pickands(t, par);
    double t3 = pickands_derivative(t, par);
    t3 = t2 - t * t3;
    t2 = (std::log(u1) + std::log(u2)) * t2;

    return std::exp(t2) * t3 / u1;
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
ExtremeValueBicop::hfunc2_raw(const tools_eigen::ConstMatRef& u,
                              const tools_eigen::ConstMatRef& parameters)
{
  if (parameters.rows() == 1) {
    Eigen::ArrayXd l1 = u.col(0).array().log();
    Eigen::ArrayXd l2 = u.col(1).array().log();
    Eigen::ArrayXd ls = l1 + l2;
    Eigen::ArrayXd t = l2 / ls;
    Eigen::ArrayXd a, a1, a2;
    pickands_all(t, parameters.row(0).transpose(), a, a1, a2);
    return ((ls * a).exp() * (a + (1.0 - t) * a1) / u.col(1).array()).matrix();
  }
  auto f = [this](const double& u1,
                  const double& u2,
                  const Eigen::Ref<const Eigen::VectorXd>& par) {
    double t = std::log(u2) / std::log(u1 * u2);
    double t2 = pickands(t, par);
    double t3 = pickands_derivative(t, par);
    t3 = t2 + (1 - t) * t3;
    t2 = (std::log(u1) + std::log(u2)) * t2;

    return std::exp(t2) * t3 / u2;
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
ExtremeValueBicop::hinv1_raw(const tools_eigen::ConstMatRef& u,
                             const tools_eigen::ConstMatRef& parameters)
{
  return hinv1_num_raw(u, parameters);
}

inline Eigen::VectorXd
ExtremeValueBicop::hinv2_raw(const tools_eigen::ConstMatRef& u,
                             const tools_eigen::ConstMatRef& parameters)
{
  return hinv2_num_raw(u, parameters);
}

inline double
ExtremeValueBicop::parameters_to_tau(const Eigen::MatrixXd& par)
{
  auto f = [this, &par](const double t) {
    double A, A1, A2;
    pickands_all(t, par.col(0), A, A1, A2);
    return t * (1 - t) * A2 / A;
  };
  return tools_integration::integrate_zero_to_one(f);
}

inline Eigen::MatrixXd
ExtremeValueBicop::parameters_to_taildep(const Eigen::MatrixXd& par)
{
  // extreme-value copulas have no lower tail dependence; the upper tail
  // dependence coefficient is 2 * (1 - A(0.5)), where A is the Pickands
  // dependence function.
  Eigen::MatrixXd taildep = Eigen::MatrixXd::Zero(2, 2);
  taildep(1, 1) = 2 * (1 - pickands(0.5, par.col(0)));
  return taildep;
}
}
