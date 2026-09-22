// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <array>
#include <utility>
#include <vector>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_interpolation {
//! A class for bilinear interpolation of bivariate copulas
//!
//! The class is used for implementing kernel estimators. It makes storing the
//! observations obsolete and allows for fast numerical integration. The
//! interpolant is piecewise bilinear, so its mass over a grid-aligned
//! rectangle is available in closed form; see `rect_mass()`. Each axis has
//! its own knot vector, so a circular axis can carry uniform knots while a
//! linear axis concentrates them in the tails.
class InterpolationGrid
{
public:
  InterpolationGrid() = default;

  InterpolationGrid(const Eigen::VectorXd& grid_points,
                    const Eigen::MatrixXd& values,
                    int norm_maxiter = 25);

  InterpolationGrid(const Eigen::VectorXd& grid_points1,
                    const Eigen::VectorXd& grid_points2,
                    const Eigen::MatrixXd& values,
                    int norm_maxiter = 25);

  Eigen::MatrixXd get_values() const;

  Eigen::VectorXd get_grid_points(size_t axis) const;

  void set_values(const Eigen::MatrixXd& values, int norm_maxiter = 25);

  void flip();

  Eigen::VectorXd interpolate(const tools_eigen::ConstMatRef& x);

  Eigen::VectorXd integrate_1d(const tools_eigen::ConstMatRef& u,
                               size_t cond_var);

  //! solves `integrate_1d(...) = p` for the non-conditioning coordinate;
  //! direct inversion of the piecewise-quadratic conditional cdf (used for
  //! the h-function inverses of kernel copulas).
  Eigen::VectorXd inverse_integrate_1d(const tools_eigen::ConstMatRef& u,
                                       size_t cond_var);

  Eigen::VectorXd integrate_2d(const tools_eigen::ConstMatRef& u);

  //! @brief probability of one rectangle, without the cancellation a
  //! difference of four `integrate_2d()` values carries.
  double rect_mass(double a1, double b1, double a2, double b2) const;

  //! @brief probability of one interval in the free coordinate, given the
  //! other.
  double cond_interval_mass(double u_cond,
                            double lo,
                            double hi,
                            size_t cond_var) const;

private:
  // the grid line at a fixed conditioning coordinate; `cond_knot` evaluates a
  // knot of it on demand, so no caller has to materialize the line
  //! @brief A grid line at a fixed conditioning coordinate.
  struct CondLine
  {
    ptrdiff_t cell;
    double x2x, xx1, x2x1;
    size_t cond_var;
  };
  CondLine cond_line(double u_cond, size_t cond_var) const;
  double cond_knot(const CondLine& line, ptrdiff_t j) const;
  // the axis along which a conditional line runs: the one not held fixed
  static size_t free_axis(size_t cond_var);

  // the weights of the two nodes of a cell `[g0, g1]`, integrating the linear
  // basis over the cell's overlap with `[a, b]`; the quadrature every integral
  // here is built from
  static std::pair<double, double> cell_weights(double g0,
                                                double g1,
                                                double a,
                                                double b);
  ptrdiff_t interval_weights(size_t axis,
                             double lo,
                             double hi,
                             Eigen::VectorXd& w) const;
  void row_integrals(double u, Eigen::VectorXd& out) const;
  // normalizes the grid margins; internal only (callers must refresh the
  // cached integrals afterwards, as the ctor and set_values do)
  void normalize_margins(int max_iter);
  void init(const Eigen::VectorXd& grid_points1,
            const Eigen::VectorXd& grid_points2,
            const Eigen::MatrixXd& values,
            int norm_maxiter);
  void update_weights(size_t axis);
  ptrdiff_t binary_search(size_t axis, double x) const;
  ptrdiff_t find_cell(size_t axis, double x) const;
  void update_cell_lookup(size_t axis);
  void update_cached_integrals();
  double cond_quantile(double u_cond,
                       double p,
                       size_t cond_var,
                       Eigen::VectorXd& knots) const;
  double int_on_grid(double upr, const Eigen::VectorXd& vals) const;

  // knots of the first (row) and second (column) axis
  std::array<Eigen::VectorXd, 2> grid_points_;
  Eigen::MatrixXd values_;
  // bucket acceleration tables for cell searches, one per axis; built once
  // (the knots are immutable after construction)
  std::array<std::vector<ptrdiff_t>, 2> cell_lookup_;
  // trapezoid weights of each axis' knots, so that `weights_[a].dot(v)`
  // integrates the piecewise linear function through `(grid_points_[a], v)`
  // over [0, 1]; built once alongside `cell_lookup_`
  std::array<Eigen::VectorXd, 2> weights_;
  // cumulative row integrals R(k, j) = int_0^{grid2_j} values_(k, .);
  // refreshed eagerly whenever values_ changes (lazy caching would race
  // when a shared grid is evaluated from multiple threads)
  Eigen::MatrixXd row_cum_int_;
};
}
}

#include <vinecopulib/misc/implementation/tools_interpolation.ipp>
