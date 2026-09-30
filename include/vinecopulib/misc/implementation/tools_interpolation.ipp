// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <stdexcept>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_interpolation {
//! Constructor
//!
//! @param grid_points An ascending sequence of grid_points; used in both
//! dimensions.
//! @param values A dxd matrix of copula density values evaluated at
//! (grid_points_i, grid_points_j).
//! @param norm_maxiter Maximum number of margin-rescaling passes; `0` leaves
//! the values untouched.
inline InterpolationGrid::InterpolationGrid(const Eigen::VectorXd& grid_points,
                                            const Eigen::MatrixXd& values,
                                            int norm_maxiter)
{
  if (values.cols() != values.rows()) {
    throw std::runtime_error("values must be a quadratic matrix");
  }
  if (grid_points.size() != values.rows()) {
    throw std::runtime_error(
      "number of grid_points must equal dimension of values");
  }

  grid_points_ = grid_points;
  values_ = values;

  // move boundary points to 0/1, so we don't have to extrapolate
  grid_points_(0) = 0.0;
  grid_points_(grid_points.size() - 1) = 1.0;

  update_cell_lookup();
  update_weights();
  normalize_margins(norm_maxiter);
  update_cached_integrals();
}

//! builds the bucket acceleration table for cell searches; the grid is
//! immutable after construction, so this runs exactly once.
inline void
InterpolationGrid::update_cell_lookup()
{
  const ptrdiff_t n_buckets = 1024;
  cell_lookup_.resize(n_buckets);
  for (ptrdiff_t k = 0; k < n_buckets; ++k) {
    cell_lookup_[k] =
      binary_search(static_cast<double>(k) / static_cast<double>(n_buckets));
  }
}

//! cumulative trapezoidal integrals of each row of the values and of their
//! transpose (used by `integrate_2d` and `rect_mass`).
inline void
InterpolationGrid::update_cached_integrals()
{
  const ptrdiff_t m = grid_points_.size();
  values_t_ = values_.transpose();
  row_cum_int_ = cumulative_row_integrals(values_);
  row_cum_int_t_ = cumulative_row_integrals(values_t_);
  line_totals_ = row_cum_int_.col(m - 1);
  line_totals_t_ = row_cum_int_t_.col(m - 1);
  margin_cum_ = cumulative_integral(line_totals_);
  margin_cum_t_ = cumulative_integral(line_totals_t_);
  total_mass_ =
    0.5 * (weights_.dot(line_totals_) + weights_.dot(line_totals_t_));
  update_orientation();
}

//! integrates along the transpose if it is the smaller of the two matrices in
//! column-major lexicographic order, so that a grid and its flipped
//! counterpart run the same arithmetic with the roles of the arguments
//! swapped; a symmetric grid gives the same values either way
inline void
InterpolationGrid::update_orientation()
{
  transposed_ = false;
  for (Eigen::Index k = 0; k < values_.size(); ++k) {
    if (values_(k) != values_t_(k)) {
      transposed_ = values_t_(k) < values_(k);
      return;
    }
  }
}

inline Eigen::MatrixXd
InterpolationGrid::cumulative_row_integrals(const Eigen::MatrixXd& values) const
{
  const ptrdiff_t m = grid_points_.size();
  Eigen::MatrixXd cum(m, m);
  for (ptrdiff_t k = 0; k < m; ++k) {
    double c = 0.0;
    cum(k, 0) = 0.0;
    for (ptrdiff_t j = 0; j < m - 1; ++j) {
      c += (values(k, j + 1) + values(k, j)) *
           (grid_points_(j + 1) - grid_points_(j)) / 2.0;
      cum(k, j + 1) = c;
    }
  }
  return cum;
}

//! O(1) cell search: bucket lookup plus a guarded advance (exact for any
//! ascending grid; equivalent to `binary_search`).
inline ptrdiff_t
InterpolationGrid::find_cell(double x) const
{
  const ptrdiff_t n_buckets = static_cast<ptrdiff_t>(cell_lookup_.size());
  const ptrdiff_t m = grid_points_.size();
  ptrdiff_t b = static_cast<ptrdiff_t>(x * static_cast<double>(n_buckets));
  b = std::min(std::max(b, static_cast<ptrdiff_t>(0)), n_buckets - 1);
  ptrdiff_t i = cell_lookup_[b];
  while ((i < m - 2) && (grid_points_(i + 1) <= x)) {
    ++i;
  }
  return i;
}

inline Eigen::MatrixXd
InterpolationGrid::get_values() const
{
  return values_;
}

inline void
InterpolationGrid::set_values(const Eigen::MatrixXd& values, int norm_maxiter)
{
  if (values.size() != values_.size()) {
    if (values.rows() != values_.rows()) {
      std::stringstream message;
      message << "values have has wrong number of rows; "
              << "expected: " << values_.rows() << ", "
              << "actual: " << values.rows() << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
    if (values.cols() != values_.cols()) {
      std::stringstream message;
      message << "values have wrong number of columns; "
              << "expected: " << values_.cols() << ", "
              << "actual: " << values.cols() << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
  }

  values_ = values;
  normalize_margins(norm_maxiter);
  update_cached_integrals();
}

inline void
InterpolationGrid::flip()
{
  values_.swap(values_t_);
  row_cum_int_.swap(row_cum_int_t_);
  line_totals_.swap(line_totals_t_);
  margin_cum_.swap(margin_cum_t_);
  update_orientation();
}

//! trapezoid weights, so that `weights_.dot(v)` is the integral over [0, 1]
//! of the piecewise linear function through `(grid_points_, v)`. They sum to
//! 1, because the sum telescopes to `grid_points_(m - 1) - grid_points_(0)`.
//! The grid is immutable after construction, so this runs exactly once.
inline void
InterpolationGrid::update_weights()
{
  const ptrdiff_t m = grid_points_.size();
  if (m < 2) {
    weights_.resize(0);
    return;
  }
  weights_.resize(m);
  weights_(0) = (grid_points_(1) - grid_points_(0)) / 2.0;
  weights_.segment(1, m - 2) =
    (grid_points_.tail(m - 2) - grid_points_.head(m - 2)) / 2.0;
  weights_(m - 1) = (grid_points_(m - 1) - grid_points_(m - 2)) / 2.0;
}

//! renormalizes the estimate to uniform margins
//!
//! @details Each pass is the elementwise geometric mean of the two ways of
//! rescaling the grid: rows then columns, and columns then rows. Averaging
//! them leaves the two margins equally close to uniform and makes the pass
//! commute with transposition exactly, so a grid and its flipped counterpart
//! normalize to flipped counterparts whether or not the iteration has
//! converged.
//!
//! @param max_iter Maximum number of rescaling passes; `0` leaves the values
//! untouched. Rescaling also stops as soon as both margins integrate to 1
//! within `1e-10`.
inline void
InterpolationGrid::normalize_margins(int max_iter)
{
  const ptrdiff_t m = grid_points_.size();
  if ((max_iter < 1) || (m < 2)) {
    return;
  }

  const double tol = 1e-10;
  const double min_mass = 1e-20; // prevent 0/0
  const Eigen::VectorXd& w = weights_;
  Eigen::MatrixXd vt(m, m);

  for (int k = 0; k < max_iter; ++k) {
    // the transpose is materialized rather than left as an expression, so
    // that both margins are the same product on a column-major matrix and
    // transposing the grid swaps them bit for bit
    vt = values_.transpose();
    const Eigen::VectorXd r = (values_ * w).cwiseMax(min_mass);
    const Eigen::VectorXd c = (vt * w).cwiseMax(min_mass);
    const double err = std::max((r.array() - 1.0).abs().maxCoeff(),
                                (c.array() - 1.0).abs().maxCoeff());
    if (err < tol) {
      break;
    }

    // Both orders are rank-one rescalings of the same values, so the second
    // margin of each is an integral against reweighted grid weights and no
    // intermediate grid is needed: rows-then-columns divides by
    // `r_i * c2_j`, columns-then-rows by `c_j * r2_i`.
    const Eigen::VectorXd r2 =
      (values_ * w.cwiseQuotient(c)).cwiseMax(min_mass);
    const Eigen::VectorXd c2 = (vt * w.cwiseQuotient(r)).cwiseMax(min_mass);
    const Eigen::VectorXd sr = r.cwiseProduct(r2).cwiseSqrt().cwiseInverse();
    const Eigen::VectorXd sc = c.cwiseProduct(c2).cwiseSqrt().cwiseInverse();

    // one fused pass: two successive rank-one scalings would round the two
    // orders differently and lose the equivariance
    for (ptrdiff_t j = 0; j < m; ++j) {
      for (ptrdiff_t i = 0; i < m; ++i) {
        values_(i, j) *= sr(i) * sc(j);
      }
    }
  }
}

inline ptrdiff_t
InterpolationGrid::binary_search(double x)
{
  ptrdiff_t low = 0;
  ptrdiff_t high = grid_points_.size() - 2; // there's one cell less than points
  ptrdiff_t mid;

  while (low < high) {
    mid = (low + high + 1) / 2; // Use upper midpoint
    if (grid_points_(mid) <= x) {
      low = mid; // Move lower bound up
    } else {
      high = mid - 1; // Move upper bound down
    }
  }

  return low;
}

//! @brief Bilinear interpolation in two dimensions.
//!
//! @param x Mx2 matrix of evaluation points.
//! @return a vector of resulting interpolated values
inline Eigen::VectorXd
InterpolationGrid::interpolate(const tools_eigen::ConstMatRef& x)
{
  auto f = [this](double u1, double u2) {
    const ptrdiff_t i = find_cell(u1);
    const ptrdiff_t j = find_cell(u2);
    const double x2x = grid_points_(i + 1) - u1;
    const double xx1 = u1 - grid_points_(i);
    const double y2y = grid_points_(j + 1) - u2;
    const double yy1 = u2 - grid_points_(j);
    return (values_(i, j) * x2x * y2y + values_(i + 1, j) * xx1 * y2y +
            values_(i, j + 1) * x2x * yy1 + values_(i + 1, j + 1) * xx1 * yy1) /
           ((grid_points_(i + 1) - grid_points_(i)) *
            (grid_points_(j + 1) - grid_points_(j)));
  };

  return tools_eigen::binaryExpr_or_nan(x, f);
}

//! the bracketing cell and interpolation distances of the grid line at
//! `u_cond`, for `cond_knot()` to evaluate a knot from
inline InterpolationGrid::CondLine
InterpolationGrid::cond_line(double u_cond, size_t cond_var) const
{
  const ptrdiff_t i = find_cell(u_cond);
  return { i,
           grid_points_(i + 1) - u_cond,
           u_cond - grid_points_(i),
           grid_points_(i + 1) - grid_points_(i),
           cond_var };
}

//! knot `j` of the line; the floor only absorbs rounding, as interpolating a
//! nonnegative grid is nonnegative
inline double
InterpolationGrid::cond_knot(const CondLine& line, ptrdiff_t j) const
{
  const ptrdiff_t i = line.cell;
  const double v =
    (line.cond_var == 1)
      ? (values_(i, j) * line.x2x + values_(i + 1, j) * line.xx1) / line.x2x1
      : (values_(j, i) * line.x2x + values_(j, i + 1) * line.xx1) / line.x2x1;
  return std::max(v, 0.0);
}

//! the weights of the nodes at `g0` and `g1` integrating the linear basis over
//! the cell's overlap with `[a, b]`; zero for a cell outside it. Callers take
//! wholly covered cells themselves, where the trapezoid needs no division.
inline std::pair<double, double>
InterpolationGrid::cell_weights(double g0, double g1, double a, double b)
{
  const double off = std::max(a - g0, 0.0);
  const double len = std::max(std::min(b, g1) - std::max(a, g0), 0.0);
  const double upper = 0.5 * len * (2.0 * off + len) / (g1 - g0);
  return { len - upper, upper };
}

//! inverts `integrate_1d` in its second argument: the conditional cdf is
//! piecewise quadratic and nondecreasing, so the quantile has a closed form
//! within the bracketing cell. Where the density vanishes the cdf is flat and
//! the inverse is not unique; the smallest quantile is returned.
inline double
InterpolationGrid::cond_quantile(double u_cond,
                                 double p,
                                 size_t cond_var,
                                 Eigen::VectorXd& knots) const
{
  const ptrdiff_t m = grid_points_.size();

  // the grid line is walked twice below, so interpolate it once into the
  // caller's buffer
  const CondLine line = cond_line(u_cond, cond_var);
  for (ptrdiff_t j = 0; j < m; ++j) {
    knots(j) = cond_knot(line, j);
  }

  // total mass (normalization of the conditional cdf)
  const double int1 = weights_.dot(knots);
  const double target =
    std::min(std::max(p, 1e-10), 1 - 1e-10) * std::max(int1, 1e-20);

  // locate the bracketing cell and solve the quadratic within it
  double cum = 0.0;
  double v_k = knots(0);
  for (ptrdiff_t k = 0; k < m - 1; ++k) {
    const double v_k1 = knots(k + 1);
    const double g_k = grid_points_(k);
    const double dg = grid_points_(k + 1) - g_k;
    const double cell = (v_k1 + v_k) * dg / 2.0;
    if ((cum + cell >= target) || (k == m - 2)) {
      // target = cum + v_k s + (v_k1 - v_k) / (2 dg) s^2 with s in [0, dg];
      // stable quadratic root (b > 0 always holds, -c >= 0 within the cell)
      const double a = (v_k1 - v_k) / (2.0 * dg);
      const double b = v_k;
      const double c = cum - target;
      const double disc = std::max(b * b - 4.0 * a * c, 0.0);
      const double denom = b + std::sqrt(disc);
      double s;
      if (denom <= 0.0) {
        // the cell carries no mass: the cdf is flat across it, so every point
        // in it is a quantile and the left endpoint is the smallest
        s = 0.0;
      } else if (std::fabs(a) < 1e-300) {
        s = -c / b;
      } else {
        s = 2.0 * (-c) / denom;
      }
      s = std::min(std::max(s, 0.0), dg);
      return g_k + s;
    }
    cum += cell;
    v_k = v_k1;
  }
  return 1.0; // unreachable: the last cell always catches the target
}

//! Integrate the grid along one axis
//!
//! @param u Mx2 matrix of evaluation points
//! @param cond_var Either 1 or 2; the axis considered fixed.
//! @return a vector of resulting integral values
inline Eigen::VectorXd
InterpolationGrid::integrate_1d(const tools_eigen::ConstMatRef& u,
                                size_t cond_var)
{
  auto f = [this, cond_var](double u1, double u2) {
    const double p = (cond_var == 1) ? cond_interval_mass(u1, 0.0, u2, 1)
                                     : cond_interval_mass(u2, 0.0, u1, 2);
    // clipped here rather than in `cond_interval_mass`, which is a mass
    return std::min(std::max(p, 1e-10), 1 - 1e-10);
  };

  return tools_eigen::binaryExpr_or_nan(u, f);
}

//! Inverse of `integrate_1d` w.r.t. the non-conditioning coordinate.
//!
//! @param u Mx2 matrix of evaluation points; for `cond_var == 1` the first
//!   column holds the conditioning coordinate and the second the probability
//!   level (and vice versa for `cond_var == 2`).
//! @param cond_var Either 1 or 2; the axis considered fixed.
inline Eigen::VectorXd
InterpolationGrid::inverse_integrate_1d(const tools_eigen::ConstMatRef& u,
                                        size_t cond_var)
{
  Eigen::VectorXd knots(grid_points_.size());
  auto f = [this, cond_var, &knots](double u1, double u2) {
    return (cond_var == 1) ? cond_quantile(u1, u2, 1, knots)
                           : cond_quantile(u2, u1, 2, knots);
  };

  return tools_eigen::binaryExpr_or_nan(u, f);
}

//! Integrate the grid along the two axis
//!
//! @details The integral is rescaled along each coordinate so that both
//! margins are exactly uniform: `C(x, y) = M(x, y) f(x) g(y) M(1, 1)`, where
//! `M` is the grid's mass, `f(x) = x / M(x, 1)` and `g(y) = y / M(1, y)`. That
//! form is symmetric in the two arguments, and flipping the grid and swapping
//! the coordinates gives exactly the same values.
//!
//! @param u Mx2 matrix of evaluation points
//! @return a vector of resulting integral values
inline Eigen::VectorXd
InterpolationGrid::integrate_2d(const tools_eigen::ConstMatRef& u)
{
  Eigen::VectorXd lines;

  auto f = [this, &lines](double u1, double u2) {
    const double x = std::min(std::max(u1, 0.0), 1.0);
    const double y = std::min(std::max(u2, 0.0), 1.0);
    double p = 0.0;
    if ((x > 0.0) && (y > 0.0)) {
      if (transposed_) {
        row_integrals(values_t_, row_cum_int_t_, x, lines);
      } else {
        row_integrals(values_, row_cum_int_, y, lines);
      }
      const double mass = int_on_grid(transposed_ ? y : x, lines);
      const double fx =
        x / std::max(margin_integral(line_totals_, margin_cum_, x), 1e-20);
      const double gy =
        y / std::max(margin_integral(line_totals_t_, margin_cum_t_, y), 1e-20);
      p = total_mass_ * ((fx * gy) * mass);
    }
    return std::min(std::max(p, 1e-10), 1 - 1e-10);
  };

  return tools_eigen::binaryExpr_or_nan(u, f);
}

//! @brief Nonnegative quadrature weights for the integral over `[lo, hi]`.
//!
//! @param lo,hi Interval bounds, clamped to `[0, 1]`; `hi <= lo` gives zero.
//! @param w Filled with the weights of the nodes the interval covers.
//! @return The index of the first of those nodes, so that
//!   `w.dot(v.segment(first, w.size()))` integrates the piecewise linear
//!   function through `(grid_points_, v)`.
inline ptrdiff_t
InterpolationGrid::interval_weights(double lo,
                                    double hi,
                                    Eigen::VectorXd& w) const
{
  const double a = std::min(std::max(lo, 0.0), 1.0);
  const double b = std::min(std::max(hi, a), 1.0);
  const ptrdiff_t ka = find_cell(a);
  const ptrdiff_t kb = find_cell(b);
  w.setZero(kb - ka + 2);

  for (ptrdiff_t k = ka; k <= kb; ++k) {
    const auto [w0, w1] =
      cell_weights(grid_points_(k), grid_points_(k + 1), a, b);
    w(k - ka) += w0;
    w(k - ka + 1) += w1;
  }
  return ka;
}

//! @brief Partial integrals of every grid line over `[0, u]`.
//!
//! @param u Upper limit, clamped to `[0, 1]`.
//! @param out Filled with one integral per grid line.
inline void
InterpolationGrid::row_integrals(const Eigen::MatrixXd& values,
                                 const Eigen::MatrixXd& cum,
                                 double u,
                                 Eigen::VectorXd& out) const
{
  const ptrdiff_t m = grid_points_.size();
  const double y = std::min(std::max(u, 0.0), 1.0);
  const ptrdiff_t j = find_cell(y);
  const double dg = grid_points_(j + 1) - grid_points_(j);
  const double s = y - grid_points_(j);
  out.resize(m);
  for (ptrdiff_t k = 0; k < m; ++k) {
    out(k) =
      cum(k, j) +
      (2 * values(k, j) + (values(k, j + 1) - values(k, j)) * s / dg) * s / 2.0;
  }
}

//! @brief Probability of the rectangle `(a1, b1] x (a2, b2]`.
//!
//! @details Not clipped, unlike `integrate_2d()`: an empty rectangle is
//! exactly `0`. Symmetric in the two arguments, as `integrate_2d()` is.
//!
//! @param a1,b1 Bounds in the first argument, in either order.
//! @param a2,b2 Bounds in the second argument, in either order.
//! @return The probability, or `0` for an empty rectangle.
inline double
InterpolationGrid::rect_mass(double a1, double b1, double a2, double b2) const
{
  // rotating the data can leave a left limit above its own value
  const double x0 = std::min(a1, b1);
  const double x1 = std::max(a1, b1);
  const double y0 = std::min(a2, b2);
  const double y1 = std::max(a2, b2);
  if (!(x1 > x0) || !(y1 > y0)) {
    return 0.0;
  }

  // `integrate_2d()`'s rescaling, expanded over the four blocks so that the
  // only differences taken are the rescalings' own increments, each over its
  // common denominator; every term is grouped so that swapping the two
  // coordinates swaps operands of a commutative operation
  Blocks b;
  if (transposed_) {
    const Blocks q = blocks_along(values_t_, row_cum_int_t_, y0, y1, x0, x1);
    b = { q.below_left, q.left, q.below, q.inside };
  } else {
    b = blocks_along(values_, row_cum_int_, x0, x1, y0, y1);
  }
  const Rescaling f = rescaling(line_totals_, margin_cum_, x0, x1);
  const Rescaling g = rescaling(line_totals_t_, margin_cum_t_, y0, y1);
  return total_mass_ * (((f.upper * g.upper) * b.inside +
                         b.below_left * (f.increment * g.increment)) +
                        (b.below * (f.upper * g.increment) +
                         b.left * (g.upper * f.increment)));
}

//! the integral over `[lo, hi]` of the piecewise linear function through
//! `(grid_points_, v)`, with nonnegative weights
inline double
InterpolationGrid::interval_integral(double lo,
                                     double hi,
                                     const Eigen::VectorXd& v) const
{
  const double a = std::min(std::max(lo, 0.0), 1.0);
  const double b = std::min(std::max(hi, a), 1.0);
  double total = 0.0;
  for (ptrdiff_t k = find_cell(a), kb = find_cell(b); k <= kb; ++k) {
    const auto [w0, w1] =
      cell_weights(grid_points_(k), grid_points_(k + 1), a, b);
    total += w0 * v(k) + w1 * v(k + 1);
  }
  return total;
}

//! the cumulative trapezoidal integrals of `v` at the nodes
inline Eigen::VectorXd
InterpolationGrid::cumulative_integral(const Eigen::VectorXd& v) const
{
  const ptrdiff_t m = grid_points_.size();
  Eigen::VectorXd cum(m);
  cum(0) = 0.0;
  for (ptrdiff_t k = 0; k < m - 1; ++k) {
    cum(k + 1) = cum(k) + (v(k + 1) + v(k)) *
                            (grid_points_(k + 1) - grid_points_(k)) / 2.0;
  }
  return cum;
}

//! the integral of `v` over `[0, x]`, from its cumulative integrals `cum`
inline double
InterpolationGrid::margin_integral(const Eigen::VectorXd& v,
                                   const Eigen::VectorXd& cum,
                                   double x) const
{
  const double b = std::min(std::max(x, 0.0), 1.0);
  const ptrdiff_t j = find_cell(b);
  const auto [w0, w1] =
    cell_weights(grid_points_(j), grid_points_(j + 1), grid_points_(j), b);
  return cum(j) + (w0 * v(j) + w1 * v(j + 1));
}

//! the four blocks, with each grid line of `values` integrated first
inline InterpolationGrid::Blocks
InterpolationGrid::blocks_along(const Eigen::MatrixXd& values,
                                const Eigen::MatrixXd& cum,
                                double x0,
                                double x1,
                                double y0,
                                double y1) const
{
  // the mass of each grid line below the rectangle and over its own strip;
  // per thread and reused, since a discrete pair asks this once per row
  thread_local Eigen::VectorXd below, strip;
  row_integrals(values, cum, y0, below);
  if (y0 > 0.0) {
    const ptrdiff_t m = grid_points_.size();
    const double a = std::min(std::max(y0, 0.0), 1.0);
    const double b = std::min(std::max(y1, a), 1.0);
    strip.setZero(m);
    for (ptrdiff_t k = find_cell(a), kb = find_cell(b); k <= kb; ++k) {
      const auto [w0, w1] =
        cell_weights(grid_points_(k), grid_points_(k + 1), a, b);
      strip += w0 * values.col(k) + w1 * values.col(k + 1);
    }
  } else {
    // nothing below, so the cached integrals are the strip itself
    row_integrals(values, cum, y1, strip);
  }
  return { int_on_grid(x0, below),
           interval_integral(x0, x1, below),
           int_on_grid(x0, strip),
           interval_integral(x0, x1, strip) };
}

inline InterpolationGrid::Rescaling
InterpolationGrid::rescaling(const Eigen::VectorXd& totals,
                             const Eigen::VectorXd& cum,
                             double x0,
                             double x1) const
{
  const double m0 = margin_integral(totals, cum, x0);
  const double strip = interval_integral(x0, x1, totals);
  const double m1 = std::max(m0 + strip, 1e-20);
  if (!(x0 > 0.0)) {
    // the blocks the increment multiplies are empty
    return { x1 / m1, 0.0 };
  }
  // `x1 / m1 - x0 / m0` over the common denominator, where it cancels against
  // `x1 - x0` rather than against one
  return { x1 / m1, ((x1 - x0) * m0 - x0 * strip) / (m1 * m0) };
}

//! @brief Probability that the free coordinate falls in `(lo, hi]`, given the
//! other.
//!
//! @details Not clipped, unlike `integrate_1d()`: an empty interval is exactly
//! `0` and the whole line exactly `1`, so the masses of a partition sum to one.
//!
//! @param u_cond The coordinate held fixed.
//! @param lo,hi Bounds in the free coordinate, in either order.
//! @param cond_var Either 1 or 2; the axis held fixed, as for `integrate_1d()`.
//! @return The conditional probability.
inline double
InterpolationGrid::cond_interval_mass(double u_cond,
                                      double lo,
                                      double hi,
                                      size_t cond_var) const
{
  const ptrdiff_t m = grid_points_.size();
  const double a = std::min(std::max(std::min(lo, hi), 0.0), 1.0);
  const double b = std::min(std::max(std::max(lo, hi), a), 1.0);

  const CondLine line = cond_line(u_cond, cond_var);
  double mass = 0.0, total = 0.0;
  double v_k = cond_knot(line, 0);
  double g_k = grid_points_(0);
  for (ptrdiff_t k = 0; k < m - 1; ++k) {
    const double v_k1 = cond_knot(line, k + 1);
    const double g_k1 = grid_points_(k + 1);
    total += weights_(k) * v_k;
    if (g_k < b && g_k1 > a) {
      if (a <= g_k && b >= g_k1) {
        mass += 0.5 * (g_k1 - g_k) * (v_k + v_k1);
      } else {
        const auto [w0, w1] = cell_weights(g_k, g_k1, a, b);
        mass += w0 * v_k + w1 * v_k1;
      }
    }
    v_k = v_k1;
    g_k = g_k1;
  }
  total += weights_(m - 1) * v_k;

  return mass / std::max(total, 1e-20);
}

// ---------------- Utility functions for integration ----------------

//! @brief Integral of the piecewise linear function through
//! `(grid_points_, vals)` over `[0, upr]`.
//!
//! @param upr Upper limit, clamped to `[0, 1]`.
//! @param vals One value per grid point.
inline double
InterpolationGrid::int_on_grid(double upr, const Eigen::VectorXd& vals) const
{
  const double b = std::min(std::max(upr, 0.0), 1.0);
  double total = 0.0;
  double g_k = grid_points_(0);
  for (ptrdiff_t k = 0; g_k < b; ++k) {
    const double g_k1 = grid_points_(k + 1);
    if (b < g_k1) {
      const auto [w0, w1] = cell_weights(g_k, g_k1, 0.0, b);
      return total + w0 * vals(k) + w1 * vals(k + 1);
    }
    // a whole cell, where the trapezoid factors a multiply cheaper than the
    // two weights `cell_weights()` returns
    total += (vals(k + 1) + vals(k)) * (g_k1 - g_k) / 2.0;
    g_k = g_k1;
  }
  return total;
}
}
}
