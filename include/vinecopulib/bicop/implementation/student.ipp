// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <algorithm>
#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {
inline StudentBicop::StudentBicop()
{
  family_ = BicopFamily::student;
  parameters_ = Eigen::VectorXd(2);
  parameters_lower_bounds_ = Eigen::VectorXd(2);
  parameters_upper_bounds_ = Eigen::VectorXd(2);
  parameters_ << 0, 50;
  parameters_lower_bounds_ << -1, 2;
  parameters_upper_bounds_ << 1, 50;
}

inline Eigen::VectorXd
StudentBicop::pdf_impl(const Eigen::MatrixXd& u, double rho, double nu)
{
  Eigen::VectorXd f = Eigen::VectorXd::Ones(u.rows());
  Eigen::MatrixXd tmp = tools_stats::qt(u, nu);

  f = tmp.col(0).cwiseAbs2() + tmp.col(1).cwiseAbs2() -
      (2 * rho) * tmp.rowwise().prod();
  f /= nu * (1.0 - pow(rho, 2.0));
  f = f + Eigen::VectorXd::Ones(u.rows());
  f = f.array().pow(-(nu + 2.0) / 2.0);
  f = f.cwiseQuotient(tools_stats::dt(tmp, nu).rowwise().prod());
  f *= boost::math::tgamma_ratio((nu + 2.0) / 2.0, nu / 2.0);
  f /= (nu * constant::pi * sqrt(1.0 - pow(rho, 2.0)));

  return f;
}

inline Eigen::VectorXd
StudentBicop::cdf_impl(const Eigen::MatrixXd& u, double rho, double nu)
{
  using namespace tools_stats;

  // for integer nu, just use pbvt
  // otherwise, interpolate linearly between floor(nu) and ceil(nu)
  if (nu == round(nu)) {
    int inu = static_cast<int>(nu);
    return pbvt(qt(u, inu), inu, rho);
  } else {
    int nu1 = static_cast<int>(std::floor(nu));
    int nu2 = static_cast<int>(std::ceil(nu));
    double weight = (nu - static_cast<double>(nu1)) /
                    (static_cast<double>(nu2) - static_cast<double>(nu1));
    return pbvt(qt(u, nu1), nu1, rho) * (1 - weight) +
           pbvt(qt(u, nu2), nu2, rho) * weight;
  }
}

inline Eigen::VectorXd
StudentBicop::hfunc1_impl(const Eigen::MatrixXd& u, double rho, double nu)
{
  Eigen::VectorXd h = Eigen::VectorXd::Ones(u.rows());
  Eigen::MatrixXd tmp = tools_stats::qt(u, nu);
  h = nu * h + tmp.col(0).cwiseAbs2();
  h *= (1.0 - pow(rho, 2)) / (nu + 1.0);
  h = h.cwiseSqrt().cwiseInverse().cwiseProduct(tmp.col(1) - rho * tmp.col(0));
  h = tools_stats::pt(h, nu + 1.0);

  return h;
}

inline Eigen::VectorXd
StudentBicop::hinv1_impl(const Eigen::MatrixXd& u, double rho, double nu)
{
  Eigen::VectorXd hinv = Eigen::VectorXd::Ones(u.rows());
  Eigen::VectorXd tmp = u.col(1);
  Eigen::VectorXd tmp2 = u.col(0);
  tmp = tools_stats::qt(tmp, nu + 1.0);
  tmp2 = tools_stats::qt(tmp2, nu);

  hinv = nu * hinv + tmp2.cwiseAbs2();
  hinv *= (1.0 - pow(rho, 2)) / (nu + 1.0);
  hinv = hinv.cwiseSqrt().cwiseProduct(tmp) + rho * tmp2;
  hinv = tools_stats::pt(hinv, nu);

  return hinv;
}

inline Eigen::VectorXd
StudentBicop::pdf_raw(const Eigen::MatrixXd& u,
                      const Eigen::MatrixXd& parameters)
{
  if (parameters.rows() == 1) {
    return pdf_impl(u, parameters(0, 0), parameters(0, 1));
  }
  Eigen::VectorXd out(u.rows());
  for (Eigen::Index i = 0; i < u.rows(); ++i) {
    out(i) = pdf_impl(u.row(i), parameters(i, 0), parameters(i, 1))(0);
  }
  return out;
}

inline Eigen::VectorXd
StudentBicop::cdf(const Eigen::MatrixXd& u, const Eigen::MatrixXd& parameters)
{
  if (parameters.rows() == 1) {
    return cdf_impl(u, parameters(0, 0), parameters(0, 1));
  }
  Eigen::VectorXd out(u.rows());
  for (Eigen::Index i = 0; i < u.rows(); ++i) {
    out(i) = cdf_impl(u.row(i), parameters(i, 0), parameters(i, 1))(0);
  }
  return out;
}

inline Eigen::VectorXd
StudentBicop::hfunc1_raw(const Eigen::MatrixXd& u,
                         const Eigen::MatrixXd& parameters)
{
  if (parameters.rows() == 1) {
    return hfunc1_impl(u, parameters(0, 0), parameters(0, 1));
  }
  Eigen::VectorXd out(u.rows());
  for (Eigen::Index i = 0; i < u.rows(); ++i) {
    out(i) = hfunc1_impl(u.row(i), parameters(i, 0), parameters(i, 1))(0);
  }
  return out;
}

inline Eigen::VectorXd
StudentBicop::hinv1_raw(const Eigen::MatrixXd& u,
                        const Eigen::MatrixXd& parameters)
{
  if (parameters.rows() == 1) {
    return hinv1_impl(u, parameters(0, 0), parameters(0, 1));
  }
  Eigen::VectorXd out(u.rows());
  for (Eigen::Index i = 0; i < u.rows(); ++i) {
    out(i) = hinv1_impl(u.row(i), parameters(i, 0), parameters(i, 1))(0);
  }
  return out;
}

inline Eigen::VectorXd
StudentBicop::get_start_parameters(const double tau)
{
  Eigen::VectorXd parameters = get_parameters();
  parameters(0) = std::sin(tau * constant::pi / 2);
  parameters(1) = 5;
  return parameters;
}

inline Eigen::MatrixXd
StudentBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}

inline Eigen::MatrixXd
StudentBicop::parameters_to_taildep(const Eigen::MatrixXd& parameters)
{
  double rho = parameters(0);
  double nu = parameters(1);
  // the t-copula has tail dependence in all four corners; the concordant
  // (lower-lower, upper-upper) corners use rho, the discordant (lower-upper,
  // upper-lower) corners use -rho.
  Eigen::MatrixXd arg(2, 1);
  arg(0) = -std::sqrt((nu + 1.0) * (1.0 - rho) / (1.0 + rho));
  arg(1) = -std::sqrt((nu + 1.0) * (1.0 + rho) / (1.0 - rho));
  Eigen::MatrixXd lambda = 2.0 * tools_stats::pt(arg, nu + 1.0);
  Eigen::MatrixXd taildep(2, 2);
  taildep(0, 0) = lambda(0); // lower-lower
  taildep(1, 1) = lambda(0); // upper-upper
  taildep(0, 1) = lambda(1); // lower-upper
  taildep(1, 0) = lambda(1); // upper-lower
  return taildep;
}

//! @name Analytic derivative leaves
//!
//! Closed forms ported from the VineCopula C sources (tcopuladeriv.c,
//! tcopuladeriv_new.c, logderiv.c). The family is exchangeable, so
//! `u2`-flavored pdf selectors are routed through a column/selector swap.
//! @{

//! applies a scalar kernel to each row of `u` with per-row parameters
inline Eigen::VectorXd
StudentBicop::apply_kernel(const Eigen::MatrixXd& u,
                           const Eigen::MatrixXd& parameters,
                           double (*kernel)(double, double, double, double))
{
  auto f = [kernel](const double& u1,
                    const double& u2,
                    const Eigen::Ref<const Eigen::VectorXd>& par) {
    return kernel(u1, u2, par(0), par(1));
  };
  return tools_eigen::binaryExpr_or_nan(u, parameters, f);
}

inline Eigen::VectorXd
StudentBicop::pdf_deriv_raw(const Eigen::MatrixXd& u,
                            const Eigen::MatrixXd& parameters,
                            const std::string& deriv)
{
  if (tools_deriv::is_u2_only(deriv)) {
    // exchangeability: c(u1, u2) = c(u2, u1)
    return pdf_deriv_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }
  if (deriv == "par1") {
    return apply_kernel(u, parameters, &diff_pdf_rho);
  } else if (deriv == "par2") {
    return apply_kernel(u, parameters, &diff_pdf_nu);
  } else if (deriv == "u1") {
    return apply_kernel(u, parameters, &diff_pdf_u1);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
StudentBicop::pdf_deriv2_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters,
                             const std::string& deriv)
{
  if (tools_deriv::is_u2_only(deriv)) {
    // exchangeability: c(u1, u2) = c(u2, u1)
    return pdf_deriv2_raw(
      tools_eigen::swap_cols(u), parameters, tools_deriv::swap_args(deriv));
  }
  if (deriv == "par1par1") {
    return apply_kernel(u, parameters, &diff2_pdf_rho);
  } else if (deriv == "par1par2") {
    return apply_kernel(u, parameters, &diff2_pdf_rho_nu);
  } else if (deriv == "par2par2") {
    return apply_kernel(u, parameters, &diff2_pdf_nu);
  } else if (deriv == "par1u1") {
    return apply_kernel(u, parameters, &diff2_pdf_rho_u1);
  } else if (deriv == "par2u1") {
    return apply_kernel(u, parameters, &diff2_pdf_nu_u1);
  } else if (deriv == "u1u1") {
    return apply_kernel(u, parameters, &diff2_pdf_u1);
  } else if (deriv == "u1u2") {
    return apply_kernel(u, parameters, &diff2_pdf_u1_u2);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
StudentBicop::hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                               const Eigen::MatrixXd& parameters,
                               const std::string& deriv)
{
  if (deriv == "par1") {
    return apply_kernel(u, parameters, &diff_hfunc1_rho);
  } else if (deriv == "par2") {
    return apply_kernel(u, parameters, &diff_hfunc1_nu);
  } else if (deriv == "u1") {
    return apply_kernel(u, parameters, &diff_hfunc1_u1);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
StudentBicop::hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  if (deriv == "par1par1") {
    return apply_kernel(u, parameters, &diff2_hfunc1_rho);
  } else if (deriv == "par1par2") {
    return apply_kernel(u, parameters, &diff2_hfunc1_rho_nu);
  } else if (deriv == "par2par2") {
    return apply_kernel(u, parameters, &diff2_hfunc1_nu);
  } else if (deriv == "par1u1") {
    return apply_kernel(u, parameters, &diff2_hfunc1_rho_u1);
  } else if (deriv == "par2u1") {
    return apply_kernel(u, parameters, &diff2_hfunc1_nu_u1);
  } else if (deriv == "u1u1") {
    return apply_kernel(u, parameters, &diff2_hfunc1_u1);
  }
  throw std::runtime_error("unexpected derivative selector: " + deriv);
}

inline Eigen::VectorXd
StudentBicop::logpdf_deriv_raw(const Eigen::MatrixXd& u,
                               const Eigen::MatrixXd& parameters,
                               const std::string& deriv)
{
  if (deriv == "par1") {
    return apply_kernel(u, parameters, &diff_lpdf_rho);
  } else if (deriv == "par2") {
    return apply_kernel(u, parameters, &diff_lpdf_nu);
  }
  return AbstractBicop::logpdf_deriv_raw(u, parameters, deriv);
}

inline Eigen::VectorXd
StudentBicop::logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  if (deriv == "par1par1") {
    return apply_kernel(u, parameters, &diff2_lpdf_rho);
  } else if (deriv == "par1par2") {
    return apply_kernel(u, parameters, &diff2_lpdf_rho_nu);
  } else if (deriv == "par2par2") {
    return apply_kernel(u, parameters, &diff2_lpdf_nu);
  }
  return AbstractBicop::logpdf_deriv2_raw(u, parameters, deriv);
}

//! @}

//! @name Scalar derivative kernels
//!
//! Ported from the VineCopula C sources; the Maple-generated temporary
//! variable chains are kept verbatim. The pdf kernels use the same argument
//! order as the C code; the h-function kernels differentiate the C code's
//! `h(u|v)` (conditioning on the second argument), while `hfunc1` conditions
//! on the first, so they are called with their `u := u2`, `v := u1`.
//! @{

//! scalar Student t copula density (mirrors `pdf_impl`)
inline double
StudentBicop::pdf_scalar(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double x1 = boost::math::quantile(dist, u1);
  double x2 = boost::math::quantile(dist, u2);
  double f =
    (x1 * x1 + x2 * x2 - 2.0 * rho * x1 * x2) / (nu * (1.0 - rho * rho)) + 1.0;
  f = std::pow(f, -(nu + 2.0) / 2.0);
  f /= boost::math::pdf(dist, x1) * boost::math::pdf(dist, x2);
  f *= boost::math::tgamma_ratio((nu + 2.0) / 2.0, nu / 2.0);
  f /= nu * constant::pi * std::sqrt(1.0 - rho * rho);
  return f;
}

//! ported from VineCopula logderiv.c difflPDF_rho_tCopula
inline double
StudentBicop::diff_lpdf_rho(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u1);
  double t2 = boost::math::quantile(dist, u2);
  double t3 = -(nu + 2.0) / 2.0;
  double t10 = nu * (1.0 - rho * rho);
  double t4 = -2.0 * t1 * t2 / t10;
  double t11 = t1 * t1 + t2 * t2 - 2.0 * rho * t1 * t2;
  double t5 = 2.0 * t11 * rho / t10 / (1.0 - rho * rho);
  double t6 = 1.0 + (t11 / t10);
  double t7 = rho / (1.0 - rho * rho);
  return t3 * (t4 + t5) / t6 + t7;
}

//! ported from VineCopula logderiv.c difflPDF_nu_tCopula_new
inline double
StudentBicop::diff_lpdf_nu(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::digamma((nu + 1.0) / 2.0);
  double t2 = boost::math::digamma(nu / 2.0);
  double t3 = 0.5 * std::log(1.0 - rho * rho);
  double t4 = (nu - 2.0) / (2.0 * nu);
  double t5 = 0.5 * std::log(nu);
  double t6 = -t1 + t2 + t3 - t4 - t5;
  double t10 = (nu + 2.0) / 2.0;

  double x1 = boost::math::quantile(dist, u1);
  double x2 = boost::math::quantile(dist, u2);
  double out1 = tools_stats::diff_x_nu(x1, nu);
  double out2 = tools_stats::diff_x_nu(x2, nu);
  double t7 = 1.0 + 2.0 * x1 * out1;
  double t8 = 1.0 + 2.0 * x2 * out2;
  double t9 = (nu + 1.0) / 2.0 * (t7 / (nu + x1 * x1) + t8 / (nu + x2 * x2));
  double M = nu * (1.0 - rho * rho) + x1 * x1 + x2 * x2 - 2.0 * rho * x1 * x2;
  double t11 = 1.0 - rho * rho + 2.0 * x1 * out1 + 2.0 * x2 * out2 -
               2.0 * rho * (x1 * out2 + x2 * out1);
  double t12 = 0.5 * std::log((nu + x1 * x1) * (nu + x2 * x2));
  double t13 = 0.5 * std::log(M);

  return t6 + t9 + t12 - t10 * t11 / M - t13;
}

//! ported from VineCopula logderiv.c diff2lPDF_rho_tCopula
inline double
StudentBicop::diff2_lpdf_rho(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u1);
  double t2 = boost::math::quantile(dist, u2);

  double t4 = 1.0 - rho * rho;
  double t3 = -(nu + 1.0) * (1.0 + rho * rho) / t4 / t4;
  double M = nu * t4 + t1 * t1 + t2 * t2 - 2.0 * rho * t1 * t2;
  double t5 = (nu + 2.0) * nu / M;
  double t6 = 2.0 * (nu + 2.0) * std::pow(nu * rho + t1 * t2, 2.0) / M / M;

  return t3 + t5 + t6;
}

//! ported from VineCopula logderiv.c diff2lPDF_nu_tCopula_new
inline double
StudentBicop::diff2_lpdf_nu(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = (nu + 1.0) / 2.0;
  double t2 = nu / 2.0;
  double t3 = 1.0 / (nu * nu);
  double t4 = 1.0 / (2.0 * nu);
  double t5 = 0.5 * boost::math::trigamma(t1);
  double t6 = 1.0 - rho * rho;
  double t9 = 0.5 * boost::math::trigamma(t2);
  double t10 = -t5 + t9 - t3 - t4;

  double x1 = boost::math::quantile(dist, u1);
  double x2 = boost::math::quantile(dist, u2);
  double out1 = tools_stats::diff_x_nu(x1, nu);
  double out2 = tools_stats::diff_x_nu(x2, nu);
  double M = nu * t6 + x1 * x1 + x2 * x2 - 2.0 * rho * x1 * x2;
  double t8 = x1 * out2 + out1 * x2;
  double M_nu = t6 + 2.0 * x1 * out1 + 2.0 * x2 * out2 - 2.0 * rho * t8;

  double t11 = 1.0 + 2.0 * x1 * out1;
  double t12 = nu + x1 * x1;
  double t13 = t11 / t12;

  double t14 = 1.0 + 2.0 * x2 * out2;
  double t15 = nu + x2 * x2;
  double t16 = t14 / t15;

  double out3 = tools_stats::diff2_x_nu(x1, nu);
  double out4 = tools_stats::diff2_x_nu(x2, nu);

  double t17 = 2.0 * out1 * out1 + 2.0 * x1 * out3;
  double t18 = t17 / t12;

  double t19 = 2.0 * out2 * out2 + 2.0 * x2 * out4;
  double t20 = t19 / t15;

  double t21 = t11 * t11 / t12 / t12;
  double t22 = t14 * t14 / t15 / t15;

  double M_nu_nu = 2.0 * out1 * out1 + 2.0 * x1 * out3 + 2.0 * out2 * out2 +
                   2.0 * x2 * out4 - 4.0 * rho * out1 * out2 -
                   2.0 * rho * (x2 * out3 + x1 * out4);

  return t10 + 0.5 * (t13 + t16) + t1 * (t18 - t21 + t20 - t22) + 0.5 * t13 +
         0.5 * t16 - M_nu / M -
         (nu / 2.0 + 1.0) * (M_nu_nu / M - M_nu * M_nu / M / M);
}

//! ported from VineCopula logderiv.c diff2lPDF_rho_nu_tCopula_new
inline double
StudentBicop::diff2_lpdf_rho_nu(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t4 = 1.0 - rho * rho;
  double t3 = rho / t4;
  double t5 = nu + 2.0;

  double x1 = boost::math::quantile(dist, u1);
  double x2 = boost::math::quantile(dist, u2);
  double out1 = tools_stats::diff_x_nu(x1, nu);
  double out2 = tools_stats::diff_x_nu(x2, nu);
  double M = nu * t4 + x1 * x1 + x2 * x2 - 2.0 * rho * x1 * x2;
  double M_rho = -2.0 * (nu * rho + x1 * x2);
  double t8 = x1 * out2 + out1 * x2;
  double M_nu = t4 + 2.0 * x1 * out1 + 2.0 * x2 * out2 - 2.0 * rho * t8;

  return -t3 + t5 / M * (rho + t8 + 0.5 * M_nu * M_rho / M) - 0.5 * M_rho / M;
}

//! ported from VineCopula tcopuladeriv.c diffPDF_rho_tCopula; the factor in
//! parentheses there is exactly difflPDF_rho_tCopula
inline double
StudentBicop::diff_pdf_rho(double u1, double u2, double rho, double nu)
{
  return pdf_scalar(u1, u2, rho, nu) * diff_lpdf_rho(u1, u2, rho, nu);
}

//! ported from VineCopula tcopuladeriv_new.c diffPDF_nu_tCopula_new; the
//! factor in parentheses there is exactly difflPDF_nu_tCopula_new
inline double
StudentBicop::diff_pdf_nu(double u1, double u2, double rho, double nu)
{
  return pdf_scalar(u1, u2, rho, nu) * diff_lpdf_nu(u1, u2, rho, nu);
}

//! ported from VineCopula tcopuladeriv_new.c diffPDF_u_tCopula_new
inline double
StudentBicop::diff_pdf_u1(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double c = pdf_scalar(u1, u2, rho, nu);
  double x1 = boost::math::quantile(dist, u1);
  double x2 = boost::math::quantile(dist, u2);
  double out1 = tools_stats::diff_dt_u(x1, nu);
  double t1 = c / boost::math::pdf(dist, x1);
  double t2 = (nu + 2.0) * (x1 - rho * x2);
  double t6 = rho * rho;
  double t7 = x1 * x1;
  double t8 = x2 * x2;
  double t3 = nu * (1.0 - t6) + t7 + t8 - 2 * rho * x1 * x2;
  double t4 = t2 / t3;
  return -t1 * (t4 + out1);
}

//! ported from VineCopula tcopuladeriv.c diff2PDF_rho_tCopula: the bracket
//! there is diff2lPDF_rho_tCopula plus (diffPDF_rho / c)^2
inline double
StudentBicop::diff2_pdf_rho(double u1, double u2, double rho, double nu)
{
  double lp = diff_lpdf_rho(u1, u2, rho, nu);
  return pdf_scalar(u1, u2, rho, nu) *
         (diff2_lpdf_rho(u1, u2, rho, nu) + lp * lp);
}

//! ported from VineCopula tcopuladeriv_new.c diff2PDF_nu_tCopula_new: the
//! bracket there is diff2lPDF_nu_tCopula_new plus (diffPDF_nu / c)^2
inline double
StudentBicop::diff2_pdf_nu(double u1, double u2, double rho, double nu)
{
  double lp = diff_lpdf_nu(u1, u2, rho, nu);
  return pdf_scalar(u1, u2, rho, nu) *
         (diff2_lpdf_nu(u1, u2, rho, nu) + lp * lp);
}

//! ported from VineCopula tcopuladeriv_new.c diff2PDF_rho_nu_tCopula_new:
//! the bracket there is diff2lPDF_rho_nu_tCopula_new plus
//! (diffPDF_rho / c) * (diffPDF_nu / c)
inline double
StudentBicop::diff2_pdf_rho_nu(double u1, double u2, double rho, double nu)
{
  return pdf_scalar(u1, u2, rho, nu) *
         (diff2_lpdf_rho_nu(u1, u2, rho, nu) +
          diff_lpdf_rho(u1, u2, rho, nu) * diff_lpdf_nu(u1, u2, rho, nu));
}

//! ported from VineCopula tcopuladeriv_new.c diff2PDF_rho_u_tCopula_new
inline double
StudentBicop::diff2_pdf_rho_u1(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = nu + 2.0;
  double t4 = 1.0 - rho * rho;

  double c = pdf_scalar(u1, u2, rho, nu);
  double x1 = boost::math::quantile(dist, u1);
  double x2 = boost::math::quantile(dist, u2);
  double M = nu * t4 + x1 * x1 + x2 * x2 - 2.0 * rho * x1 * x2;
  double t2 = nu * rho + x1 * x2;
  double t3 = x1 - rho * x2;
  double out1 = diff_pdf_rho(u1, u2, rho, nu);
  double out2 = diff_pdf_u1(u1, u2, rho, nu);

  return c * (t1 / M / boost::math::pdf(dist, x1) * (x2 - 2.0 * t2 / M * t3)) +
         out1 / c * out2;
}

//! ported from VineCopula tcopuladeriv_new.c diff2PDF_nu_u_tCopula_new
inline double
StudentBicop::diff2_pdf_nu_u1(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = nu + 2.0;
  double t3 = (nu + 1.0) / nu;
  double t14 = rho * rho;
  double t4 = 1.0 - t14;
  double t17 = nu * nu;

  double c = pdf_scalar(u1, u2, rho, nu);
  double x1 = boost::math::quantile(dist, u1);
  double x2 = boost::math::quantile(dist, u2);
  double t15 = x1 * x1;
  double t16 = x2 * x2;
  double M = nu * t4 + t15 + t16 - 2.0 * rho * x1 * x2;
  double t2 = boost::math::pdf(dist, x1);
  double diff_pdf = diff_pdf_nu(u1, u2, rho, nu);
  double diff_dt = tools_stats::diff_dt_nu(x1, nu);
  double diff_dt2 = tools_stats::diff_dt_u(x1, nu);
  double diff_dt3 = tools_stats::diff_dt_x(x1, nu);
  double out1 = tools_stats::diff_x_nu(x1, nu);
  double out2 = tools_stats::diff_x_nu(x2, nu);
  double t8 = x1 * out2 + out1 * x2;
  double M_nu = t4 + 2.0 * x1 * out1 + 2.0 * x2 * out2 - 2.0 * rho * t8;

  double t6 = 1.0 + t15 / nu;
  double t7 = x1 - rho * x2;
  double t9 = -diff_pdf / t2 + c / t2 / t2 * (diff_dt + diff_dt3 * out1);
  double t10 = t1 * t7 / M + diff_dt2;
  double t11 = c / t2;
  double t12 = t7 / M - t1 * t7 / M / M * M_nu + t1 * (out1 - rho * out2) / M;
  double t13 = -out1 * t3 / t6 + x1 / (t17 + nu * t15) +
               x1 * t3 / t6 / t6 * (2.0 * x1 * out1 / nu - t15 / t17);

  return t9 * t10 - t11 * (t12 + t13);
}

//! ported from VineCopula tcopuladeriv_new.c diff2PDF_u_tCopula_new
inline double
StudentBicop::diff2_pdf_u1(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = nu + 2.0;
  double t14 = rho * rho;
  double t4 = 1.0 - t14;

  double c = pdf_scalar(u1, u2, rho, nu);
  double x1 = boost::math::quantile(dist, u1);
  double x2 = boost::math::quantile(dist, u2);
  double t15 = x1 * x1;
  double t16 = x2 * x2;
  double M = nu * t4 + t15 + t16 - 2.0 * rho * x1 * x2;
  double t2 = boost::math::pdf(dist, x1);
  double diff_pdf = diff_pdf_u1(u1, u2, rho, nu);
  double diff_dt2 = tools_stats::diff_dt_u(x1, nu);

  double t7 = x1 - rho * x2;

  double t5 = -diff_pdf / t2 + diff_dt2 / t2 / t2 * c;
  double t6 = t1 * t7 / M + diff_dt2;

  double t11 = c / t2;
  double t8 = 1.0 / t2;
  double t9 = t1 / M - 2.0 * t1 * t7 * t7 / M / M;
  double t13 = 1.0 + t15 / nu;

  double t12 = t8 * (-(nu + 1.0) / (nu + t15) +
                     2.0 * t15 * (nu + 1.0) / nu / nu / t13 / t13);

  return t5 * t6 - t11 * (t8 * t9 + t12);
}

//! ported from VineCopula tcopuladeriv_new.c diff2PDF_u_v_tCopula_new
inline double
StudentBicop::diff2_pdf_u1_u2(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = nu + 2.0;
  double t4 = 1.0 - rho * rho;

  double c = pdf_scalar(u1, u2, rho, nu);
  double x1 = boost::math::quantile(dist, u1);
  double x2 = boost::math::quantile(dist, u2);
  double M = nu * t4 + x1 * x1 + x2 * x2 - 2.0 * rho * x1 * x2;
  double t2 = boost::math::pdf(dist, x1);
  double t5 = boost::math::pdf(dist, x2);

  double t6 = x1 - rho * x2;
  double t7 = x2 - rho * x1;

  double diff_dt1 = tools_stats::diff_dt_u(x1, nu);
  double diff_dt2 = tools_stats::diff_dt_u(x2, nu);

  double t8 = c / t2 / t5;
  double t9 = t1 * t6 / M + diff_dt1;
  double t10 = t1 * t7 / M + diff_dt2;
  double t11 = t1 * rho / M;
  double t12 = 2.0 * t1 * t6 * t7 / M / M;

  return t8 * (t9 * t10 + t11 + t12);
}

//! ported from VineCopula tcopuladeriv.c diffhfunc_rho_tCopula
inline double
StudentBicop::diff_hfunc1_rho(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u2);
  double t2 = boost::math::quantile(dist, u1);
  double t3 = t1 - rho * t2;
  double t4 = nu + t2 * t2;
  double t5 = 1.0 - rho * rho;
  double t6 = nu + 1.0;
  double t7 = t4 * t5 / t6;
  double t8 = std::sqrt(t7);
  boost::math::students_t dist1(nu + 1.0);
  double t9 = boost::math::pdf(dist1, t3 / t8);
  double t10 = -t2 / t8;
  double t11 = t3 * t4 * rho / t6 / t7 / t8;
  return t9 * (t10 + t11);
}

//! ported from VineCopula tcopuladeriv_new.c diffhfunc_nu_tCopula_new
inline double
StudentBicop::diff_hfunc1_nu(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u2);
  double t2 = boost::math::quantile(dist, u1);
  double t3 = t1 - rho * t2;
  double t4 = nu + t2 * t2;
  double t5 = 1.0 - rho * rho;
  double t6 = nu + 1.0;
  double t7 = t4 * t5 / t6;
  double t8 = std::sqrt(t7);

  double t10 = t3 / t8;
  boost::math::students_t dist1(t6);
  double t9 = boost::math::pdf(dist1, t10);
  double t11 = nu + 1.0;

  double diff_t = tools_stats::diff_t_nu(t10, t11);
  double out1 = tools_stats::diff_x_nu(t1, nu);
  double out2 = tools_stats::diff_x_nu(t2, nu);

  double t12 = out1 - rho * out2;
  double t13 = t12 / t8;
  double t14 = 1.0 + 2.0 * t2 * out2;
  double t15 = t14 / t6;
  double t16 = t4 / t6 / t6;

  return diff_t + t9 * (t13 - 0.5 * t10 / t7 * t5 * (t15 - t16));
}

//! ported from VineCopula tcopuladeriv_new.c diffhfunc_v_tCopula_new
inline double
StudentBicop::diff_hfunc1_u1(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u2);
  double t2 = boost::math::quantile(dist, u1);
  double t3 = t1 - rho * t2;
  double t4 = nu + t2 * t2;
  double t5 = 1.0 - rho * rho;
  double t6 = nu + 1.0;
  double t7 = t4 * t5 / t6;
  double t8 = std::sqrt(t7);

  double t10 = t3 / t8;
  boost::math::students_t dist1(nu + 1.0);
  double t9 = boost::math::pdf(dist1, t10);
  double t11 = nu + 1.0;

  double t12 = 1.0 / boost::math::pdf(dist, t2);
  double t13 = -rho / t8;
  double t14 = t10 / t7;
  double t15 = t5 / t11 * t2;

  return t9 * t12 * (t13 - t14 * t15);
}

//! ported from VineCopula tcopuladeriv_new.c diff2hfunc_rho_tCopula_new
inline double
StudentBicop::diff2_hfunc1_rho(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u2);
  double t2 = boost::math::quantile(dist, u1);
  double t3 = t1 - rho * t2;
  double t4 = nu + t2 * t2;
  double t5 = 1.0 - rho * rho;
  double t6 = nu + 1.0;
  double t7 = t4 * t5 / t6;
  double t8 = std::sqrt(t7);
  double t14 = t4 / t6;

  double t10 = t3 / t8;
  boost::math::students_t dist1(nu + 1.0);
  double t9 = boost::math::pdf(dist1, t10);
  double t11 = nu + 1.0;

  double diff_t = tools_stats::diff_dt_x(t10, t11);

  double t12 = -t2 / t8;
  double t13 = t10 / t7 * rho * t14;

  return diff_t * (t12 + t13) * (t12 + t13) +
         t9 * ((t2 / t7 / t8 - 1.5 * t13 / t7) * (-2.0 * rho) * t14 +
               t10 / t7 * t14);
}

//! ported from VineCopula tcopuladeriv_new.c diff2hfunc_nu_tCopula_new
inline double
StudentBicop::diff2_hfunc1_nu(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u2);
  double t2 = boost::math::quantile(dist, u1);
  double t3 = t1 - rho * t2;
  double t4 = nu + t2 * t2;
  double t5 = 1.0 - rho * rho;
  double t6 = nu + 1.0;
  double t7 = t4 * t5 / t6;
  double t8 = std::sqrt(t7);

  double t10 = t3 / t8;
  boost::math::students_t dist1(t6);
  double t9 = boost::math::pdf(dist1, t10);
  double t11 = nu + 1.0;

  double diff_t = tools_stats::diff_dt_x(t10, t11);
  double diff_t2 = tools_stats::diff_t_nu_nu(t10, t11);
  double diff_t3 = tools_stats::diff_dt_nu(t10, t11);

  double out1 = tools_stats::diff_x_nu(t1, nu);
  double out2 = tools_stats::diff_x_nu(t2, nu);
  double out3 = tools_stats::diff2_x_nu(t1, nu);
  double out4 = tools_stats::diff2_x_nu(t2, nu);

  double t12 = out1 - rho * out2;
  double t13 = t12 / t8;
  double t14 = 1.0 + 2.0 * t2 * out2;
  double t15 = t14 / t6;
  double t16 = t4 / t6 / t6;
  double t17 = t15 - t16;

  double t18 = t13 - 0.5 * t10 / t7 * t5 * t17;
  double t19 = -t13 / t7 + 0.75 * t10 / t7 / t7 * t5 * t17;
  double t20 = (2.0 * out2 * out2 + 2.0 * t2 * out4) / t6 - 2.0 * t15 / t6 +
               2.0 * t16 / t6;
  double t21 = (out3 - rho * out4) / t8;

  return diff_t2 + 2.0 * diff_t3 * t18 + diff_t * t18 * t18 +
         t9 * (t19 * t5 * t17 + t21 - 0.5 * t10 / t7 * t5 * t20);
}

//! ported from VineCopula tcopuladeriv_new.c diff2hfunc_rho_nu_tCopula_new
inline double
StudentBicop::diff2_hfunc1_rho_nu(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u2);
  double t2 = boost::math::quantile(dist, u1);
  double t3 = t1 - rho * t2;
  double t4 = nu + t2 * t2;
  double t5 = 1.0 - rho * rho;
  double t6 = nu + 1.0;
  double t7 = t4 * t5 / t6;
  double t8 = std::sqrt(t7);

  double t10 = t3 / t8;
  boost::math::students_t dist1(nu + 1.0);
  double t9 = boost::math::pdf(dist1, t10);
  double t11 = nu + 1.0;

  double diff_t = tools_stats::diff_dt_x(t10, t11);
  double diff_t3 = tools_stats::diff_dt_nu(t10, t11);

  double out1 = tools_stats::diff_x_nu(t1, nu);
  double out2 = tools_stats::diff_x_nu(t2, nu);

  double t12 = out1 - rho * out2;
  double t13 = t12 / t8;
  double t14 = 1.0 + 2.0 * t2 * out2;
  double t15 = t14 / t6;
  double t16 = t4 / t6 / t6;
  double t17 = t15 - t16;

  double t18 = t13 - 0.5 * t10 / t7 * t5 * t17;
  double t19 = -t2 / t8 + t10 / t7 * rho * t4 / t6;
  double t20 = 0.5 * t2 / t7 / t8 - 1.5 * t10 / t7 / t7 * rho * t4 / t6;
  double t21 = out2 / t8;
  double t22 = t13 / t7 * rho * t4 / t6;
  double t23 = t10 / t7 * rho * (t6 * t14 - t4) / t6 / t6;

  return (diff_t3 + diff_t * t18) * t19 +
         t9 * (t20 * t5 * t17 - t21 + t22 + t23);
}

//! ported from VineCopula tcopuladeriv_new.c diff2hfunc_rho_v_tCopula_new
inline double
StudentBicop::diff2_hfunc1_rho_u1(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u2);
  double t2 = boost::math::quantile(dist, u1);
  double t3 = t1 - rho * t2;
  double t4 = nu + t2 * t2;
  double t5 = 1.0 - rho * rho;
  double t6 = nu + 1.0;
  double t7 = t4 * t5 / t6;
  double t8 = std::sqrt(t7);

  double t10 = t3 / t8;
  boost::math::students_t dist1(t6);
  double t9 = boost::math::pdf(dist1, t10);
  double t11 = nu + 1.0;

  double diff_t = tools_stats::diff_dt_x(t10, t11);

  double t12 = boost::math::pdf(dist, t2);
  double t13 = -t2 / t8 + t10 / t7 * rho * t4 / t6;
  double t14 = -rho / t8 - t10 / t7 * t2 * t5 / t6;
  double t15 = -1.0 / t8;
  double t16 = t2 * t2 / t8 / t7 * t5 / t6;
  double t17 = t10 / t7 * 2.0 * rho * t2 / t6;
  double t18 = 1.5 * t10 / t7 / t7 * t5 / t6 * t2 + 0.5 * rho / t8 / t7;
  double t19 = -2.0 * rho * t4 / t6;

  return diff_t / t12 * t13 * t14 + t9 / t12 * (t15 + t16 + t17 + t18 * t19);
}

//! ported from VineCopula tcopuladeriv_new.c diff2hfunc_nu_v_tCopula_new
inline double
StudentBicop::diff2_hfunc1_nu_u1(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u2);
  double t2 = boost::math::quantile(dist, u1);
  double t3 = t1 - rho * t2;
  double t4 = nu + t2 * t2;
  double t5 = 1.0 - rho * rho;
  double t6 = nu + 1.0;
  double t7 = t4 * t5 / t6;
  double t8 = std::sqrt(t7);

  double t10 = t3 / t8;
  boost::math::students_t dist1(nu + 1.0);
  double t9 = boost::math::pdf(dist1, t10);
  double t11 = nu + 1.0;

  double diff_t = tools_stats::diff_dt_nu(t2, nu);
  double diff_t2 = tools_stats::diff_dt_x(t10, t11);
  double diff_t3 = tools_stats::diff_dt_nu(t10, t11);
  double diff_t4 = tools_stats::diff_dt_x(t2, nu);

  double out1 = tools_stats::diff_x_nu(t1, nu);
  double out2 = tools_stats::diff_x_nu(t2, nu);

  double t12 = boost::math::pdf(dist, t2);

  double t13 = out1 - rho * out2;
  double t14 = t13 / t8;

  double t16 = (1.0 + 2.0 * t2 * out2) / t6 - t4 / t6 / t6;
  double t15 = t14 - 0.5 * t10 / t7 * t5 * t16;

  double t17 = -rho / t8 - t10 / t7 * t5 / t6 * t2;
  double t18 = 0.5 * rho / t8 / t7 * t5 * t16;
  double t19 = t14 / t7 * t5 / t6 * t2;
  double t20 = t10 / t7 * t5 / t6 * out2;
  double t21 = t10 / t7 * t5 / t6 / t6 * t2;
  double t22 = 1.5 * t10 / t7 / t7 * t5 * t5 / t6 * t16 * t2;

  return (diff_t3 / t12 + diff_t2 / t12 * t15 -
          t9 * (diff_t / t12 / t12 + diff_t4 / t12 / t12 * out2)) *
           t17 +
         t9 / t12 * (t18 - t19 - t20 + t21 + t22);
}

//! ported from VineCopula tcopuladeriv_new.c diff2hfunc_v_tCopula_new
inline double
StudentBicop::diff2_hfunc1_u1(double u1, double u2, double rho, double nu)
{
  boost::math::students_t dist(nu);
  double t1 = boost::math::quantile(dist, u2);
  double t2 = boost::math::quantile(dist, u1);
  double t3 = t1 - rho * t2;
  double t4 = nu + t2 * t2;
  double t5 = 1.0 - rho * rho;
  double t6 = nu + 1.0;
  double t7 = t4 * t5 / t6;
  double t8 = std::sqrt(t7);

  double t10 = t3 / t8;
  boost::math::students_t dist1(t6);
  double t9 = boost::math::pdf(dist1, t10);
  double t11 = nu + 1.0;

  double diff_t = tools_stats::diff_dt_x(t10, t11);
  double diff_t2 = tools_stats::diff_dt_x(t2, nu);
  double t12 = boost::math::pdf(dist, t2);

  double t13 = -rho / t8 - t10 / t7 * t5 / t6 * t2;

  double t14 = 0.5 * rho / t8 / t7 + 1.5 * t10 / t7 / t7 * t5 / t6 * t2;
  double t15 = t5 / t6 * 2.0 * t2 / t12;
  double t16 = t5 / t6 / t12 * (t3 - rho * t2) / t7 / t8;

  return (diff_t / t12 / t12 * t13 - t9 * diff_t2 / t12 / t12 / t12) * t13 +
         t9 / t12 * (t14 * t15 - t16);
}

//! @}
}
