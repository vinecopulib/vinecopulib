// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_var_types.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {
inline KernelBicop::KernelBicop()
{
  // construct default grid (equally spaced on Gaussian scale)
  size_t grid_size = 30;
  auto grid_points = this->make_normal_grid(grid_size);

  interp_grid_ = std::make_shared<tools_interpolation::InterpolationGrid>(
    grid_points,
    Eigen::MatrixXd::Constant(grid_size, grid_size, 1.0) // independence
  );
  npars_ = 0.0;
}

//! the knots follow the geometry of the axes: when a variable changes between
//! linear and circular, the grid is rebuilt on the new default knots with the
//! density values kept.
inline void
KernelBicop::set_var_types(const std::vector<std::string>& var_types)
{
  const auto old_types = var_types_;
  AbstractBicop::set_var_types(var_types);
  bool same_geometry = true;
  for (size_t a = 0; a < 2; ++a) {
    same_geometry =
      same_geometry && (tools_var_types::is_circular(old_types[a]) ==
                        tools_var_types::is_circular(var_types_[a]));
  }
  if (!same_geometry) {
    const Eigen::MatrixXd values = interp_grid_->get_values();
    interp_grid_ = std::make_shared<tools_interpolation::InterpolationGrid>(
      make_grid_points(var_types_[0], values.rows()),
      make_grid_points(var_types_[1], values.cols()),
      values,
      0);
  }
}

inline std::vector<Eigen::VectorXd>
KernelBicop::get_grid_knots() const
{
  return { interp_grid_->get_grid_points(0), interp_grid_->get_grid_points(1) };
}

//! replaces the grid by one with explicit knots; the values are taken as they
//! are (no margin normalization).
inline void
KernelBicop::set_grid(const std::vector<Eigen::VectorXd>& knots,
                      const Eigen::MatrixXd& values)
{
  if (knots.size() != 2) {
    throw std::runtime_error("grid knots must be given for both axes.");
  }
  for (size_t a = 0; a < 2; ++a) {
    const Eigen::VectorXd& g = knots[a];
    if (g.size() < 3) {
      throw std::runtime_error("each axis needs at least 3 grid knots.");
    }
    for (Eigen::Index i = 0; i < g.size(); ++i) {
      const bool ordered = (i == 0) || (g(i) > g(i - 1));
      if (!ordered || g(i) < 0.0 || g(i) > 1.0) {
        throw std::runtime_error(
          "grid knots must be increasing and lie in [0, 1].");
      }
    }
  }
  if (values.minCoeff() < 0) {
    throw std::runtime_error("density should be larger than 0. ");
  }
  interp_grid_ = std::make_shared<tools_interpolation::InterpolationGrid>(
    knots[0], knots[1], values, 0);
}

inline Eigen::VectorXd
KernelBicop::pdf_raw(const Eigen::MatrixXd& u, const Eigen::MatrixXd&)
{
  auto pdf = interp_grid_->interpolate(u);
  tools_eigen::trim(pdf, 1e-20, DBL_MAX);
  return pdf;
}

inline Eigen::VectorXd
KernelBicop::cdf(const Eigen::MatrixXd& u, const Eigen::MatrixXd&)
{
  return interp_grid_->integrate_2d(u);
}

inline double
KernelBicop::rect_prob(double a1,
                       double b1,
                       double a2,
                       double b2,
                       const Eigen::MatrixXd&)
{
  return interp_grid_->rect_mass(a1, b1, a2, b2);
}

inline double
KernelBicop::cond_interval_prob(double u_cond,
                                double lo,
                                double hi,
                                size_t cond_var,
                                const Eigen::MatrixXd&)
{
  return interp_grid_->cond_interval_mass(u_cond, lo, hi, cond_var);
}

inline Eigen::VectorXd
KernelBicop::hfunc1_raw(const Eigen::MatrixXd& u, const Eigen::MatrixXd&)
{
  return interp_grid_->integrate_1d(u, 1);
}

inline Eigen::VectorXd
KernelBicop::hfunc2_raw(const Eigen::MatrixXd& u, const Eigen::MatrixXd&)
{
  return interp_grid_->integrate_1d(u, 2);
}

inline Eigen::VectorXd
KernelBicop::hinv1_raw(const Eigen::MatrixXd& u, const Eigen::MatrixXd&)
{
  // direct inversion of the interpolated conditional cdf; replaces the
  // generic bisection (which re-integrated the grid 35 times per point)
  return interp_grid_->inverse_integrate_1d(u, 1);
}

inline Eigen::VectorXd
KernelBicop::hinv2_raw(const Eigen::MatrixXd& u, const Eigen::MatrixXd&)
{
  return interp_grid_->inverse_integrate_1d(u, 2);
}

inline double
KernelBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  auto oldpars = this->get_parameters();
  auto old_types = var_types_;
  this->set_parameters(parameters);
  var_types_ = tools_var_types::as_continuous(var_types_);

  std::vector<int> seeds = {
    204967043, 733593603, 184618802, 399707801, 290266245
  };
  auto u = tools_stats::ghalton(1000, 2, seeds);
  u.col(1) = hinv1_raw(u, Eigen::MatrixXd());

  this->set_parameters(oldpars);
  var_types_ = old_types;
  return wdm::wdm(u, "tau")(0, 1);
}

inline double
KernelBicop::get_npars() const
{
  return npars_;
}

inline void
KernelBicop::set_npars(const double& npars)
{
  if (npars < 0) {
    throw std::runtime_error("npars must be positive.");
  }
  npars_ = npars;
}

inline Eigen::MatrixXd
KernelBicop::get_parameters() const
{
  return interp_grid_->get_values();
}

inline Eigen::MatrixXd
KernelBicop::get_parameters_lower_bounds() const
{
  const auto& values = interp_grid_->get_values();
  return Eigen::MatrixXd::Constant(values.rows(), values.cols(), 0.0);
}

inline Eigen::MatrixXd
KernelBicop::get_parameters_upper_bounds() const
{
  const auto& values = interp_grid_->get_values();
  return Eigen::MatrixXd::Constant(values.rows(), values.cols(), 1e4);
}

//! the values are placed on the default knots of the current variable types;
//! a matrix of the current shape keeps the current knots.
inline void
KernelBicop::set_parameters(const Eigen::MatrixXd& parameters)
{
  Eigen::Index rows = parameters.rows();
  Eigen::Index cols = parameters.cols();
  if (rows < 3 || cols < 3) {
    std::stringstream message;
    message << "parameters must be a matrix with at least 3 rows and 3 "
            << "columns, got " << rows << " rows and " << cols << " columns.";
    throw std::runtime_error(message.str().c_str());
  }
  if (parameters.minCoeff() < 0) {
    std::stringstream message;
    message << "density should be larger than 0. ";
    throw std::runtime_error(message.str().c_str());
  }
  const auto& values = interp_grid_->get_values();
  if (rows == values.rows() && cols == values.cols()) {
    // don't normalize again!
    interp_grid_->set_values(parameters, 0);
  } else {
    // create new interpolation grid with new size
    interp_grid_ = std::make_shared<tools_interpolation::InterpolationGrid>(
      make_grid_points(var_types_[0], rows),
      make_grid_points(var_types_[1], cols),
      parameters,
      0);
  }
}

inline void
KernelBicop::flip()
{
  interp_grid_->flip();
}

inline Eigen::MatrixXd
KernelBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}

inline Eigen::VectorXd
KernelBicop::make_grid_points(const std::string& var_type, size_t m)
{
  if (tools_var_types::is_circular(var_type)) {
    return Eigen::VectorXd::LinSpaced(static_cast<Eigen::Index>(m), 0.0, 1.0);
  }
  Eigen::VectorXd grid_points(m);
  for (size_t i = 0; i < m; ++i)
    grid_points(i) =
      -3.25 + static_cast<double>(i) * (6.5 / static_cast<double>(m - 1));
  return tools_stats::pnorm(grid_points);
}

// construct default grid (equally spaced on Gaussian scale)
inline Eigen::VectorXd
KernelBicop::make_normal_grid(size_t m)
{
  return make_grid_points(tools_var_types::continuous(), m);
}
}
