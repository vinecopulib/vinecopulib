// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <vinecopulib/bicop/bb1.hpp>
#include <vinecopulib/bicop/bb6.hpp>
#include <vinecopulib/bicop/bb7.hpp>
#include <vinecopulib/bicop/bb8.hpp>
#include <vinecopulib/bicop/clayton.hpp>
#include <vinecopulib/bicop/frank.hpp>
#include <vinecopulib/bicop/gaussian.hpp>
#include <vinecopulib/bicop/gumbel.hpp>
#include <vinecopulib/bicop/indep.hpp>
#include <vinecopulib/bicop/joe.hpp>
#include <vinecopulib/bicop/student.hpp>
#include <vinecopulib/bicop/tawn.hpp>
#include <vinecopulib/bicop/tll.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

//! virtual destructor
inline AbstractBicop::~AbstractBicop() {}

//! Instantiates a bivariate copula using the default contructor
//!
//! @param family The copula family.
//! @param parameters The copula parameters (optional, must be compatible
//!     with family).
//! @return A pointer to an object that inherits from AbstractBicop.
//! @{
inline BicopPtr
AbstractBicop::create(BicopFamily family, const Eigen::MatrixXd& parameters)
{
  BicopPtr new_bicop;
  switch (family) {
    case BicopFamily::indep:
      new_bicop = BicopPtr(new IndepBicop());
      break;
    case BicopFamily::gaussian:
      new_bicop = BicopPtr(new GaussianBicop());
      break;
    case BicopFamily::student:
      new_bicop = BicopPtr(new StudentBicop());
      break;
    case BicopFamily::clayton:
      new_bicop = BicopPtr(new ClaytonBicop());
      break;
    case BicopFamily::gumbel:
      new_bicop = BicopPtr(new GumbelBicop());
      break;
    case BicopFamily::frank:
      new_bicop = BicopPtr(new FrankBicop());
      break;
    case BicopFamily::joe:
      new_bicop = BicopPtr(new JoeBicop());
      break;
    case BicopFamily::bb1:
      new_bicop = BicopPtr(new Bb1Bicop());
      break;
    case BicopFamily::bb6:
      new_bicop = BicopPtr(new Bb6Bicop());
      break;
    case BicopFamily::bb7:
      new_bicop = BicopPtr(new Bb7Bicop());
      break;
    case BicopFamily::bb8:
      new_bicop = BicopPtr(new Bb8Bicop());
      break;
    case BicopFamily::tawn:
      new_bicop = BicopPtr(new TawnBicop());
      break;
    case BicopFamily::tll:
      new_bicop = BicopPtr(new TllBicop());
      break;

    default:
      throw std::runtime_error(std::string("Family not implemented"));
  }

  if (parameters.size() > 0) {
    new_bicop->set_parameters(parameters);
  }

  return new_bicop;
}

//!@}

inline Eigen::MatrixXd
AbstractBicop::no_tau_to_parameters(const double&)
{
  throw std::runtime_error("Method not implemented for this family");
}

//! Default tail dependence: not implemented for this family, so all four
//! corners are reported as NaN. Families with a closed form override this
//! (including those that genuinely have zero tail dependence, e.g. `indep`,
//! `gaussian`, `frank`).
inline Eigen::MatrixXd
AbstractBicop::parameters_to_taildep(const Eigen::MatrixXd&)
{
  return Eigen::MatrixXd::Constant(2, 2, NAN);
}

//! Blomqvist's beta computed generically from the copula cdf as
//! \f$ \beta = 4 C(0.5, 0.5) - 1 \f$. Works for all families (including the
//! nonparametric kernel estimator).
inline double
AbstractBicop::parameters_to_beta(const Eigen::MatrixXd& parameters)
{
  Eigen::MatrixXd u(1, 2);
  u << 0.5, 0.5;
  // the cdf leaf expects an (m x p) parameter matrix (one row per evaluation
  // point); parameters is a (p x 1) column, so transpose it to a single row.
  Eigen::MatrixXd par_row =
    (parameters.cols() == 1) ? parameters.transpose() : parameters;
  return 4.0 * this->cdf(u, par_row)(0) - 1.0;
}

//! Getters and setters.
//! @{
inline BicopFamily
AbstractBicop::get_family() const
{
  return family_;
}

inline std::string
AbstractBicop::get_family_name() const
{
  return vinecopulib::get_family_name(family_);
}

inline double
AbstractBicop::get_loglik() const
{
  return loglik_;
}

inline void
AbstractBicop::set_loglik(const double loglik)
{
  loglik_ = loglik;
}

inline void
AbstractBicop::set_var_types(const std::vector<std::string>& var_types)
{
  if (var_types.size() != 2) {
    throw std::runtime_error("var_types must have size two.");
  }
  var_types_ = var_types;
  all_continuous_ = (var_types_[0] == "c") && (var_types_[1] == "c");
}

inline const Eigen::MatrixXd&
AbstractBicop::parameters_row() const
{
  static thread_local Eigen::MatrixXd storage;
  storage = get_parameters().transpose();
  return storage;
}
//! @}

//! evaluates the pdf, but truncates it's value by DBL_MIN and DBL_MAX.
//! @param u Matrix of evaluation points.
inline Eigen::VectorXd
AbstractBicop::pdf(const Eigen::MatrixXd& u)
{
  if (!all_continuous_) {
    // discrete margins go through the parameter-aware difference quotients
    return pdf(u, parameters_row());
  }
  Eigen::VectorXd pdf = pdf_raw(u.leftCols(2), parameters_row());
  tools_eigen::trim(pdf, DBL_MIN, DBL_MAX);
  return pdf;
}

inline Eigen::VectorXd
AbstractBicop::hfunc1(const Eigen::MatrixXd& u)
{
  if (var_types_[0] == "d") {
    return hfunc1(u, parameters_row());
  }
  return hfunc1_raw(u.leftCols(2), parameters_row());
}

inline Eigen::VectorXd
AbstractBicop::hfunc2(const Eigen::MatrixXd& u)
{
  if (var_types_[1] == "d") {
    return hfunc2(u, parameters_row());
  }
  return hfunc2_raw(u.leftCols(2), parameters_row());
}

inline Eigen::VectorXd
AbstractBicop::hinv1(const Eigen::MatrixXd& u)
{
  if (var_types_[0] == "c") {
    return hinv1_raw(u.leftCols(2), parameters_row());
  } else {
    return hinv1_num(u);
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv2(const Eigen::MatrixXd& u)
{
  if (var_types_[1] == "c") {
    return hinv2_raw(u.leftCols(2), parameters_row());
  } else {
    return hinv2_num(u);
  }
}

//! default log-density: log of the (trimmed) density.
inline Eigen::VectorXd
AbstractBicop::log_pdf_raw(const tools_eigen::ConstMatRef& u,
                           const tools_eigen::ConstMatRef& parameters)
{
  Eigen::VectorXd pdf = pdf_raw(u, parameters);
  tools_eigen::trim(pdf, DBL_MIN, DBL_MAX);
  return pdf.array().log().matrix();
}

//! evaluates the log-likelihood.
//! @param u Data matrix.
//! @param weights Optional weights for each observation.
inline double
AbstractBicop::loglik(const Eigen::MatrixXd& u, const Eigen::VectorXd& weights)
{
  Eigen::MatrixXd log_pdf;
  if (all_continuous_) {
    // same clamp semantics as trimming the density to [DBL_MIN, DBL_MAX]
    log_pdf = log_pdf_raw(u.leftCols(2), parameters_row());
    tools_eigen::trim(log_pdf, std::log(DBL_MIN), std::log(DBL_MAX));
  } else {
    log_pdf = this->pdf(u).array().log();
  }
  if (weights.size() > 0) {
    log_pdf = log_pdf.cwiseProduct(weights);
  }
  tools_eigen::remove_nans(log_pdf);
  return log_pdf.sum();
}

//! Numerical inversion of h-functions
//!
//! These are generic functions to invert the hfunctions numerically.
//! They can be used in derived classes to define \c hinv1 and \c hinv2.
//!
//! @param u \f$m \times 2\f$ matrix of evaluation points.
//! @return The numerical inverse of h-functions.
//! @{
inline Eigen::VectorXd
AbstractBicop::hinv1_num(const Eigen::MatrixXd& u)
{
  return hinv1_num(u, parameters_row());
}

inline Eigen::VectorXd
AbstractBicop::hinv2_num(const Eigen::MatrixXd& u)
{
  return hinv2_num(u, parameters_row());
}
//! @}

//! @name Parameter-aware leaves and dispatchers
//!
//! `parameters` has shape m x p with m in {1, n}: row i holds the p parameters
//! for observation i; a single row is broadcast to all observations. These are
//! the sole evaluation interface; the state-based methods above call them with
//! the stored parameters as a 1 x p broadcast row.
//! @{

inline Eigen::VectorXd
AbstractBicop::pdf(const tools_eigen::ConstMatRef& u,
                   const tools_eigen::ConstMatRef& parameters)
{
  Eigen::VectorXd pdf;
  if (all_continuous_) {
    pdf = pdf_raw(u.leftCols(2), parameters);
  } else if ((var_types_[0] == "d") && (var_types_[1] == "d")) {
    pdf = pdf_d_d(u, parameters);
  } else {
    pdf = pdf_c_d(u, parameters);
  }
  tools_eigen::trim(pdf, DBL_MIN, DBL_MAX);
  return pdf;
}

inline Eigen::VectorXd
AbstractBicop::pdf_c_d(const tools_eigen::ConstMatRef& u,
                       const tools_eigen::ConstMatRef& parameters)
{
  Eigen::VectorXd pdf(u.rows());
  Eigen::MatrixXd umax = u.leftCols(2);
  Eigen::MatrixXd umin = u.rightCols(2);
  Eigen::VectorXd udiff(u.rows());

  if (var_types_[0] != "c") {
    udiff = (u.col(0) - u.col(2)).cwiseAbs();
  } else {
    udiff = (u.col(1) - u.col(3)).cwiseAbs();
  }

  const bool bc = (parameters.rows() == 1);
  Eigen::MatrixXd par_i;
  if (bc) {
    par_i = parameters.row(0);
  }
  for (Eigen::Index i = 0; i < u.rows(); i++) {
    if (!bc) {
      par_i = parameters.row(i);
    }
    if (udiff(i) > 5e-5) {
      if (var_types_[0] != "c") {
        pdf(i) =
          (hfunc2_raw(umax.row(i), par_i) - hfunc2_raw(umin.row(i), par_i))(0);
      } else {
        pdf(i) =
          (hfunc1_raw(umax.row(i), par_i) - hfunc1_raw(umin.row(i), par_i))(0);
      }
      pdf(i) /= udiff(i);
    } else {
      pdf(i) = pdf_raw((umax.row(i) + umin.row(i)) / 2, par_i)(0);
    }
  }
  return pdf.cwiseAbs();
}

inline Eigen::VectorXd
AbstractBicop::pdf_d_d(const tools_eigen::ConstMatRef& u,
                       const tools_eigen::ConstMatRef& parameters)
{
  Eigen::VectorXd pdf(u.rows());
  Eigen::MatrixXd umax = u.leftCols(2);
  Eigen::MatrixXd umin = u.rightCols(2);
  Eigen::MatrixXd udiff = (umax - umin).cwiseAbs();

  const bool bc = (parameters.rows() == 1);
  Eigen::MatrixXd par_i;
  if (bc) {
    par_i = parameters.row(0);
  }
  for (Eigen::Index i = 0; i < u.rows(); i++) {
    if (!bc) {
      par_i = parameters.row(i);
    }
    // the difference quotient can be instable, use derivative if denominator
    // too small
    if (udiff.row(i).maxCoeff() < 5e-5) {
      pdf(i) = pdf_raw((umax.row(i) + umin.row(i)) / 2, par_i)(0);
    } else if (udiff(i, 0) < 5e-5) {
      umax(i, 0) = (umax(i, 0) + umin(i, 0)) / 2;
      umin(i, 0) = (umax(i, 0) + umin(i, 0)) / 2;
      pdf(i) = (hfunc1_raw(umax.row(i), par_i)(0) -
                hfunc1_raw(umin.row(i), par_i)(0)) /
               udiff(i, 1);
    } else if (udiff(i, 1) < 5e-5) {
      umax(i, 1) = (umax(i, 1) + umin(i, 1)) / 2;
      umin(i, 1) = (umax(i, 1) + umin(i, 1)) / 2;
      pdf(i) = (hfunc2_raw(umax.row(i), par_i)(0) -
                hfunc2_raw(umin.row(i), par_i)(0)) /
               udiff(i, 0);
    } else {
      pdf(i) = cdf(umax.row(i), par_i)(0) + cdf(umin.row(i), par_i)(0);
      std::swap(umax(i, 0), umin(i, 0));
      pdf(i) -= cdf(umax.row(i), par_i)(0) + cdf(umin.row(i), par_i)(0);
      pdf(i) /= udiff(i, 0) * udiff(i, 1);
    }
  }

  return pdf.cwiseAbs();
}

inline Eigen::VectorXd
AbstractBicop::hfunc1(const tools_eigen::ConstMatRef& u,
                      const tools_eigen::ConstMatRef& parameters)
{
  if (var_types_[0] == "d") {
    Eigen::MatrixXd uu = u;
    uu.col(3) = uu.col(1);
    auto u1diff = (uu.col(0) - uu.col(2)).cwiseAbs();
    Eigen::VectorXd h(u.rows());

    const bool bc = (parameters.rows() == 1);
    Eigen::MatrixXd par_i;
    if (bc) {
      par_i = parameters.row(0);
    }
    for (Eigen::Index i = 0; i < u.rows(); i++) {
      if (!bc) {
        par_i = parameters.row(i);
      }
      if (std::abs(u1diff(i)) > 5e-5) {
        h(i) = cdf(uu.row(i).leftCols(2), par_i)(0) -
               cdf(uu.row(i).rightCols(2), par_i)(0);
        h(i) /= u1diff(i);
      } else {
        uu(i, 0) = (uu(i, 0) + uu(i, 2)) / 2;
        h(i) = hfunc1_raw(uu.row(i).leftCols(2), par_i)(0);
      }
    }
    return h.cwiseAbs();
  } else {
    return hfunc1_raw(u.leftCols(2), parameters);
  }
}

inline Eigen::VectorXd
AbstractBicop::hfunc2(const tools_eigen::ConstMatRef& u,
                      const tools_eigen::ConstMatRef& parameters)
{
  if (var_types_[1] == "d") {
    Eigen::MatrixXd uu = u;
    uu.col(2) = uu.col(0);
    auto u2diff = (uu.col(1) - uu.col(3)).cwiseAbs();
    Eigen::VectorXd h(u.rows());

    const bool bc = (parameters.rows() == 1);
    Eigen::MatrixXd par_i;
    if (bc) {
      par_i = parameters.row(0);
    }
    for (Eigen::Index i = 0; i < u.rows(); i++) {
      if (!bc) {
        par_i = parameters.row(i);
      }
      if (u2diff(i) > 5e-5) {
        h(i) = cdf(uu.row(i).leftCols(2), par_i)(0) -
               cdf(uu.row(i).rightCols(2), par_i)(0);
        h(i) /= u2diff(i);
      } else {
        uu(i, 1) = (uu(i, 1) + uu(i, 3)) / 2;
        h(i) = hfunc2_raw(uu.row(i).leftCols(2), par_i)(0);
      }
    }
    return h.cwiseAbs();
  } else {
    return hfunc2_raw(u.leftCols(2), parameters);
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv1(const tools_eigen::ConstMatRef& u,
                     const tools_eigen::ConstMatRef& parameters)
{
  if (var_types_[0] == "c") {
    return hinv1_raw(u.leftCols(2), parameters);
  } else {
    return hinv1_num(u, parameters);
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv2(const tools_eigen::ConstMatRef& u,
                     const tools_eigen::ConstMatRef& parameters)
{
  if (var_types_[1] == "c") {
    return hinv2_raw(u.leftCols(2), parameters);
  } else {
    return hinv2_num(u, parameters);
  }
}

inline Eigen::VectorXd
AbstractBicop::hinv1_num(const tools_eigen::ConstMatRef& u,
                         const tools_eigen::ConstMatRef& parameters)
{
  Eigen::MatrixXd u_new = u;
  auto h1 = [&](const Eigen::VectorXd& v) {
    u_new.col(1) = v;
    return hfunc1(u_new, parameters);
  };

  return tools_eigen::invert_f(u.col(1), h1);
}

inline Eigen::VectorXd
AbstractBicop::hinv2_num(const tools_eigen::ConstMatRef& u,
                         const tools_eigen::ConstMatRef& parameters)
{
  Eigen::MatrixXd u_new = u;
  auto h1 = [&](const Eigen::VectorXd& x) {
    u_new.col(0) = x;
    return hfunc2(u_new, parameters);
  };

  return tools_eigen::invert_f(u.col(0), h1);
}

//! builds the active-row evaluator for the safeguarded Newton h-inverse.
//! `cond_col` is the fixed conditioning column and `solve_col` is filled with
//! the current iterates; the returned callback evaluates the h-function and
//! density (its derivative) on just the unconverged rows, gathering the
//! matching per-row parameters when they vary by row.
inline tools_eigen::NewtonEval
AbstractBicop::make_hinv_eval(const tools_eigen::ConstMatRef& u,
                              const tools_eigen::ConstMatRef& parameters,
                              int cond_col,
                              int solve_col,
                              HinvHFunc hfunc)
{
  const bool per_row = parameters.rows() > 1;
  const Eigen::VectorXd ucond = u.col(cond_col);
  // Capture an OWNING copy of the parameters: the returned callback outlives
  // this call, and in the broadcast case a captured `Ref` would alias
  // `parameters_row()`'s thread_local buffer (clobbered by any state-based
  // evaluation). The copy is 1 x p in the common broadcast case.
  const Eigen::MatrixXd params = parameters;
  return [this, ucond, params, cond_col, solve_col, per_row, hfunc](
           const std::vector<Eigen::Index>& active,
           const Eigen::VectorXd& v_active,
           Eigen::VectorXd& f_out,
           Eigen::VectorXd& fprime_out) {
    const Eigen::Index m = static_cast<Eigen::Index>(active.size());
    Eigen::MatrixXd uu(m, 2);
    for (Eigen::Index k = 0; k < m; ++k) {
      uu(k, cond_col) = ucond(active[k]);
      uu(k, solve_col) = v_active(k);
    }
    if (per_row) {
      Eigen::MatrixXd par_sub(m, params.cols());
      for (Eigen::Index k = 0; k < m; ++k) {
        par_sub.row(k) = params.row(active[k]);
      }
      f_out = (this->*hfunc)(uu, par_sub);
      fprime_out = pdf_raw(uu, par_sub);
    } else {
      f_out = (this->*hfunc)(uu, params);
      fprime_out = pdf_raw(uu, params);
    }
  };
}

inline Eigen::VectorXd
AbstractBicop::hinv1_num_raw(const tools_eigen::ConstMatRef& u,
                             const tools_eigen::ConstMatRef& parameters)
{
  // invert h1(u2 | u1) over u2; d h1 / d u2 = c(u1, u2) = pdf, so the
  // density is the derivative that drives the safeguarded Newton step
  auto eval = make_hinv_eval(u, parameters, 0, 1, &AbstractBicop::hfunc1_raw);
  return tools_eigen::invert_f_newton(u.col(1), eval);
}

inline Eigen::VectorXd
AbstractBicop::hinv2_num_raw(const tools_eigen::ConstMatRef& u,
                             const tools_eigen::ConstMatRef& parameters)
{
  // invert h2(u1 | u2) over u1; d h2 / d u1 = c(u1, u2) = pdf
  auto eval = make_hinv_eval(u, parameters, 1, 0, &AbstractBicop::hfunc2_raw);
  return tools_eigen::invert_f_newton(u.col(0), eval);
}
//! @}

//! @name Derivative leaves (defaults)
//!
//! The defaults difference the value leaves by central finite differences, so
//! every family with parameter bounds and value leaves supports the derivative
//! interface out of the box; the second-order versions difference the
//! (possibly analytic) first-derivative leaves. Analytic families override
//! these with closed forms; nonparametric families override where finite
//! differences are meaningless. Steps are clipped to the parameter bounds /
//! the unit interval, with the effective step kept in the denominator.
//! @{

//! differentiates `f` w.r.t. component `comp` (0-based parameter index, `-1`
//! for the first argument, `-2` for the second) by central differences.
inline Eigen::VectorXd
AbstractBicop::fd_deriv(
  const std::function<Eigen::VectorXd(const Eigen::MatrixXd&,
                                      const Eigen::MatrixXd&)>& f,
  const Eigen::MatrixXd& u,
  const Eigen::MatrixXd& parameters,
  int comp)
{
  if (comp >= 0) {
    Eigen::MatrixXd par_plus = parameters;
    Eigen::MatrixXd par_minus = parameters;
    // each family controls its own clamping through the bound getters
    Eigen::MatrixXd lb_bounds = get_parameters_lower_bounds();
    Eigen::MatrixXd ub_bounds = get_parameters_upper_bounds();
    double lb = lb_bounds(comp);
    double ub = ub_bounds(comp);
    for (Eigen::Index i = 0; i < parameters.rows(); ++i) {
      double par = parameters(i, comp);
      double eps = 1e-4 * std::max(1.0, std::fabs(par));
      par_plus(i, comp) = std::min(par + eps, ub);
      par_minus(i, comp) = std::max(par - eps, lb);
    }
    Eigen::VectorXd diff = f(u, par_plus) - f(u, par_minus);
    if (parameters.rows() == 1) {
      return diff / (par_plus(0, comp) - par_minus(0, comp));
    }
    Eigen::ArrayXd eps = (par_plus.col(comp) - par_minus.col(comp)).array();
    return (diff.array() / eps).matrix();
  }

  // argument derivative; stay strictly inside the unit interval
  Eigen::Index col = (comp == -1) ? 0 : 1;
  Eigen::MatrixXd u_plus = u;
  Eigen::MatrixXd u_minus = u;
  u_plus.col(col) = (u.col(col).array() + 1e-5).min(1 - 1e-10);
  u_minus.col(col) = (u.col(col).array() - 1e-5).max(1e-10);
  Eigen::ArrayXd eps = u_plus.col(col).array() - u_minus.col(col).array();
  return ((f(u_plus, parameters) - f(u_minus, parameters)).array() / eps)
    .matrix();
}

inline Eigen::VectorXd
AbstractBicop::pdf_deriv_raw(const Eigen::MatrixXd& u,
                             const Eigen::MatrixXd& parameters,
                             const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto f = [this](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return pdf_raw(uu, pp);
  };
  return fd_deriv(f, u, parameters, comps[0]);
}

inline Eigen::VectorXd
AbstractBicop::pdf_deriv2_raw(const Eigen::MatrixXd& u,
                              const Eigen::MatrixXd& parameters,
                              const std::string& deriv)
{
  // difference the first-derivative leaf w.r.t. the second component; this
  // uses the analytic first derivative when the family provides one
  auto comps = tools_deriv::parse_components(deriv);
  auto first = tools_deriv::comp_to_string(comps[0]);
  auto f = [this, first](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return pdf_deriv_raw(uu, pp, first);
  };
  return fd_deriv(f, u, parameters, comps[1]);
}

inline Eigen::VectorXd
AbstractBicop::hfunc1_deriv_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto f = [this](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return hfunc1_raw(uu, pp);
  };
  return fd_deriv(f, u, parameters, comps[0]);
}

inline Eigen::VectorXd
AbstractBicop::hfunc1_deriv2_raw(const Eigen::MatrixXd& u,
                                 const Eigen::MatrixXd& parameters,
                                 const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto first = tools_deriv::comp_to_string(comps[0]);
  auto f = [this, first](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return hfunc1_deriv_raw(uu, pp, first);
  };
  return fd_deriv(f, u, parameters, comps[1]);
}

inline Eigen::VectorXd
AbstractBicop::hfunc2_deriv_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto f = [this](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return hfunc2_raw(uu, pp);
  };
  return fd_deriv(f, u, parameters, comps[0]);
}

inline Eigen::VectorXd
AbstractBicop::hfunc2_deriv2_raw(const Eigen::MatrixXd& u,
                                 const Eigen::MatrixXd& parameters,
                                 const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  auto first = tools_deriv::comp_to_string(comps[0]);
  auto f = [this, first](const Eigen::MatrixXd& uu, const Eigen::MatrixXd& pp) {
    return hfunc2_deriv_raw(uu, pp, first);
  };
  return fd_deriv(f, u, parameters, comps[1]);
}

inline Eigen::VectorXd
AbstractBicop::logpdf_deriv_raw(const Eigen::MatrixXd& u,
                                const Eigen::MatrixXd& parameters,
                                const std::string& deriv)
{
  Eigen::ArrayXd c = pdf_raw(u, parameters).array().max(DBL_MIN);
  return (pdf_deriv_raw(u, parameters, deriv).array() / c).matrix();
}

inline Eigen::VectorXd
AbstractBicop::logpdf_deriv2_raw(const Eigen::MatrixXd& u,
                                 const Eigen::MatrixXd& parameters,
                                 const std::string& deriv)
{
  auto comps = tools_deriv::parse_components(deriv);
  Eigen::ArrayXd c = pdf_raw(u, parameters).array().max(DBL_MIN);
  Eigen::ArrayXd c_xy = pdf_deriv2_raw(u, parameters, deriv).array();
  Eigen::ArrayXd c_x =
    pdf_deriv_raw(u, parameters, tools_deriv::comp_to_string(comps[0])).array();
  Eigen::ArrayXd c_y =
    (comps[0] == comps[1])
      ? c_x
      : pdf_deriv_raw(u, parameters, tools_deriv::comp_to_string(comps[1]))
          .array();
  return (c_xy / c - (c_x / c) * (c_y / c)).matrix();
}
//! @}

namespace tools_deriv {

//! splits a derivative selector into components; the encoding is the 0-based
//! parameter index for `"par<k>"`, `-1` for `"u1"`, and `-2` for `"u2"`.
inline std::vector<int>
parse_components(const std::string& deriv)
{
  std::vector<int> comps;
  size_t pos = 0;
  while (pos < deriv.size()) {
    if (deriv.compare(pos, 1, "u") == 0) {
      if ((pos + 1 >= deriv.size()) ||
          ((deriv[pos + 1] != '1') && (deriv[pos + 1] != '2'))) {
        throw std::runtime_error("invalid derivative selector: '" + deriv +
                                 "'");
      }
      comps.push_back((deriv[pos + 1] == '1') ? -1 : -2);
      pos += 2;
    } else if (deriv.compare(pos, 3, "par") == 0) {
      pos += 3;
      size_t k = 0;
      size_t num_digits = 0;
      while ((pos < deriv.size()) && (deriv[pos] >= '0') &&
             (deriv[pos] <= '9')) {
        k = 10 * k + static_cast<size_t>(deriv[pos] - '0');
        pos++;
        num_digits++;
      }
      if (num_digits == 0) {
        k = 1; // "par" is short for "par1"
      }
      if ((k < 1) || (num_digits > 3)) {
        throw std::runtime_error("invalid derivative selector: '" + deriv +
                                 "'");
      }
      comps.push_back(static_cast<int>(k) - 1);
    } else {
      throw std::runtime_error("invalid derivative selector: '" + deriv + "'");
    }
  }
  if (comps.empty()) {
    throw std::runtime_error("derivative selector cannot be empty");
  }
  return comps;
}

inline std::string
comp_to_string(int comp)
{
  if (comp == -1) {
    return "u1";
  } else if (comp == -2) {
    return "u2";
  }
  return "par" + std::to_string(comp + 1);
}

//! sorts components canonically (parameters by index, then u1, then u2) and
//! concatenates them.
inline std::string
components_to_string(std::vector<int> comps)
{
  auto rank = [](int comp) { return (comp >= 0) ? comp : 1000 - comp; };
  std::sort(comps.begin(), comps.end(), [&](int lhs, int rhs) {
    return rank(lhs) < rank(rhs);
  });
  std::string str;
  for (auto comp : comps) {
    str += comp_to_string(comp);
  }
  return str;
}

//! validates a user-facing selector against the derivative order and the
//! number of parameters and returns its canonical form; a single component
//! in a second-order selector means differentiating twice w.r.t. it.
inline std::string
canonicalize(const std::string& deriv, size_t order, size_t npars)
{
  auto comps = parse_components(deriv);
  if ((comps.size() == 1) && (order == 2)) {
    comps.push_back(comps[0]);
  }
  if (comps.size() != order) {
    throw std::runtime_error("derivative selector '" + deriv + "' has " +
                             std::to_string(comps.size()) +
                             " components; expected " + std::to_string(order));
  }
  for (auto comp : comps) {
    if (comp >= static_cast<int>(npars)) {
      throw std::runtime_error(
        "derivative selector '" + deriv + "' refers to parameter " +
        std::to_string(comp + 1) + ", but the family has " +
        std::to_string(npars) + " parameter(s)");
    }
  }
  return components_to_string(comps);
}

//! swaps `"u1"` and `"u2"` in a selector (for exchangeable families).
inline std::string
swap_args(const std::string& deriv)
{
  auto comps = parse_components(deriv);
  for (auto& comp : comps) {
    if (comp == -1) {
      comp = -2;
    } else if (comp == -2) {
      comp = -1;
    }
  }
  return components_to_string(comps);
}

//! whether a selector involves `"u2"` but not `"u1"` (an exchangeable family
//! can then route it to its `"u1"`-flavored leaf via `swap_args`).
inline bool
is_u2_only(const std::string& deriv)
{
  auto comps = parse_components(deriv);
  bool has_u1 = std::find(comps.begin(), comps.end(), -1) != comps.end();
  bool has_u2 = std::find(comps.begin(), comps.end(), -2) != comps.end();
  return has_u2 && !has_u1;
}
}
}
