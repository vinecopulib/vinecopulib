// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <stdexcept>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_interpolation {
//! Constructor with the same knots on both axes
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
  init(grid_points, grid_points, values, norm_maxiter);
}

//! Constructor with one knot vector per axis
//!
//! @param grid_points1 An ascending sequence of knots for the first axis
//! (the rows of `values`).
//! @param grid_points2 An ascending sequence of knots for the second axis
//! (the columns of `values`).
//! @param values A matrix of copula density values evaluated at
//! (grid_points1_i, grid_points2_j).
//! @param norm_maxiter Maximum number of margin-rescaling passes; `0` leaves
//! the values untouched.
inline InterpolationGrid::InterpolationGrid(const Eigen::VectorXd& grid_points1,
                                            const Eigen::VectorXd& grid_points2,
                                            const Eigen::MatrixXd& values,
                                            int norm_maxiter)
{
  init(grid_points1, grid_points2, values, norm_maxiter);
}

inline void
InterpolationGrid::init(const Eigen::VectorXd& grid_points1,
                        const Eigen::VectorXd& grid_points2,
                        const Eigen::MatrixXd& values,
                        int norm_maxiter)
{
  if (grid_points1.size() != values.rows()) {
    throw std::runtime_error(
      "number of grid_points must equal dimension of values");
  }
  if (grid_points2.size() != values.cols()) {
    throw std::runtime_error(
      "number of grid_points must equal dimension of values");
  }

  grid_points_[0] = grid_points1;
  grid_points_[1] = grid_points2;
  values_ = values;

  for (size_t axis = 0; axis < 2; ++axis) {
    // move boundary points to 0/1, so we don't have to extrapolate
    grid_points_[axis](0) = 0.0;
    grid_points_[axis](grid_points_[axis].size() - 1) = 1.0;
    update_cell_lookup(axis);
    update_weights(axis);
  }
  normalize_margins(norm_maxiter);
  update_cached_integrals();
}

//! builds the bucket acceleration table for cell searches on one axis; the
//! knots are immutable after construction, so this runs exactly once.
inline void
InterpolationGrid::update_cell_lookup(size_t axis)
{
  const ptrdiff_t n_buckets = 1024;
  cell_lookup_[axis].resize(n_buckets);
  for (ptrdiff_t k = 0; k < n_buckets; ++k) {
    cell_lookup_[axis][k] = binary_search(
      axis, static_cast<double>(k) / static_cast<double>(n_buckets));
  }
}

//! cumulative trapezoidal integrals of each row (used by `integrate_2d`).
inline void
InterpolationGrid::update_cached_integrals()
{
  const ptrdiff_t m1 = grid_points_[0].size();
  const ptrdiff_t m2 = grid_points_[1].size();
  const Eigen::VectorXd& g2 = grid_points_[1];
  row_cum_int_.resize(m1, m2);
  for (ptrdiff_t k = 0; k < m1; ++k) {
    double cum = 0.0;
    row_cum_int_(k, 0) = 0.0;
    for (ptrdiff_t j = 0; j < m2 - 1; ++j) {
      cum += (values_(k, j + 1) + values_(k, j)) * (g2(j + 1) - g2(j)) / 2.0;
      row_cum_int_(k, j + 1) = cum;
    }
  }
}

//! O(1) cell search: bucket lookup plus a guarded advance (exact for any
//! ascending grid; equivalent to `binary_search`).
inline ptrdiff_t
InterpolationGrid::find_cell(size_t axis, double x) const
{
  const Eigen::VectorXd& g = grid_points_[axis];
  const std::vector<ptrdiff_t>& lookup = cell_lookup_[axis];
  const ptrdiff_t n_buckets = static_cast<ptrdiff_t>(lookup.size());
  const ptrdiff_t m = g.size();
  ptrdiff_t b = static_cast<ptrdiff_t>(x * static_cast<double>(n_buckets));
  b = std::min(std::max(b, static_cast<ptrdiff_t>(0)), n_buckets - 1);
  ptrdiff_t i = lookup[b];
  while ((i < m - 2) && (g(i + 1) <= x)) {
    ++i;
  }
  return i;
}

inline Eigen::MatrixXd
InterpolationGrid::get_values() const
{
  return values_;
}

//! @brief The knots of one axis.
//! @param axis `0` for the first (row) axis, `1` for the second.
inline Eigen::VectorXd
InterpolationGrid::get_grid_points(size_t axis) const
{
  return grid_points_[axis];
}

inline void
InterpolationGrid::set_values(const Eigen::MatrixXd& values, int norm_maxiter)
{
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

  values_ = values;
  normalize_margins(norm_maxiter);
  update_cached_integrals();
}

//! exchanges the two axes: the values are transposed and the knots follow.
inline void
InterpolationGrid::flip()
{
  values_.transposeInPlace();
  std::swap(grid_points_[0], grid_points_[1]);
  std::swap(weights_[0], weights_[1]);
  std::swap(cell_lookup_[0], cell_lookup_[1]);
  update_cached_integrals();
}

//! trapezoid weights of one axis, so that `weights_[axis].dot(v)` is the
//! integral over [0, 1] of the piecewise linear function through
//! `(grid_points_[axis], v)`. They sum to 1, because the sum telescopes to
//! `g(m - 1) - g(0)`. The knots are immutable after construction, so this
//! runs exactly once.
inline void
InterpolationGrid::update_weights(size_t axis)
{
  const Eigen::VectorXd& g = grid_points_[axis];
  Eigen::VectorXd& w = weights_[axis];
  const ptrdiff_t m = g.size();
  if (m < 2) {
    w.resize(0);
    return;
  }
  w.resize(m);
  w(0) = (g(1) - g(0)) / 2.0;
  w.segment(1, m - 2) = (g.tail(m - 2) - g.head(m - 2)) / 2.0;
  w(m - 1) = (g(m - 1) - g(m - 2)) / 2.0;
}

//! renormalizes the estimate to uniform margins
//!
//! @details Each pass is the elementwise geometric mean of the two ways of
//! rescaling the grid: rows then columns, and columns then rows. Averaging
//! them leaves the two margins equally close to uniform and makes the pass
//! commute with transposition exactly, so a grid and its flipped counterpart
//! normalize to flipped counterparts whether or not the iteration has
//! converged. Rows that hold the same values (the two ends of a circular
//! axis) receive the same factor, so periodic endpoint equality survives.
//!
//! @param max_iter Maximum number of rescaling passes; `0` leaves the values
//! untouched. Rescaling also stops as soon as both margins integrate to 1
//! within `1e-10`.
inline void
InterpolationGrid::normalize_margins(int max_iter)
{
  const ptrdiff_t m1 = grid_points_[0].size();
  const ptrdiff_t m2 = grid_points_[1].size();
  if ((max_iter < 1) || (m1 < 2) || (m2 < 2)) {
    return;
  }

  const double tol = 1e-10;
  const double min_mass = 1e-20;           // prevent 0/0
  const Eigen::VectorXd& w1 = weights_[0]; // integrates along the rows' axis
  const Eigen::VectorXd& w2 = weights_[1]; // integrates along the columns'
  Eigen::MatrixXd vt(m2, m1);

  for (int k = 0; k < max_iter; ++k) {
    // the transpose is materialized rather than left as an expression, so
    // that both margins are the same product on a column-major matrix and
    // transposing the grid swaps them bit for bit
    vt = values_.transpose();
    const Eigen::VectorXd r = (values_ * w2).cwiseMax(min_mass);
    const Eigen::VectorXd c = (vt * w1).cwiseMax(min_mass);
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
      (values_ * w2.cwiseQuotient(c)).cwiseMax(min_mass);
    const Eigen::VectorXd c2 = (vt * w1.cwiseQuotient(r)).cwiseMax(min_mass);
    const Eigen::VectorXd sr = r.cwiseProduct(r2).cwiseSqrt().cwiseInverse();
    const Eigen::VectorXd sc = c.cwiseProduct(c2).cwiseSqrt().cwiseInverse();

    // one fused pass: two successive rank-one scalings would round the two
    // orders differently and lose the equivariance
    for (ptrdiff_t j = 0; j < m2; ++j) {
      for (ptrdiff_t i = 0; i < m1; ++i) {
        values_(i, j) *= sr(i) * sc(j);
      }
    }
  }
}

inline ptrdiff_t
InterpolationGrid::binary_search(size_t axis, double x) const
{
  const Eigen::VectorXd& g = grid_points_[axis];
  ptrdiff_t low = 0;
  ptrdiff_t high = g.size() - 2; // there's one cell less than points
  ptrdiff_t mid;

  while (low < high) {
    mid = (low + high + 1) / 2; // Use upper midpoint
    if (g(mid) <= x) {
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
  const Eigen::VectorXd& g1 = grid_points_[0];
  const Eigen::VectorXd& g2 = grid_points_[1];
  auto f = [this, &g1, &g2](double u1, double u2) {
    const ptrdiff_t i = find_cell(0, u1);
    const ptrdiff_t j = find_cell(1, u2);
    const double x2x = g1(i + 1) - u1;
    const double xx1 = u1 - g1(i);
    const double y2y = g2(j + 1) - u2;
    const double yy1 = u2 - g2(j);
    return (values_(i, j) * x2x * y2y + values_(i + 1, j) * xx1 * y2y +
            values_(i, j + 1) * x2x * yy1 + values_(i + 1, j + 1) * xx1 * yy1) /
           ((g1(i + 1) - g1(i)) * (g2(j + 1) - g2(j)));
  };

  return tools_eigen::binaryExpr_or_nan(x, f);
}

inline size_t
InterpolationGrid::free_axis(size_t cond_var)
{
  return (cond_var == 1) ? 1 : 0;
}

//! the bracketing cell and interpolation distances of the grid line at
//! `u_cond`, for `cond_knot()` to evaluate a knot from
inline InterpolationGrid::CondLine
InterpolationGrid::cond_line(double u_cond, size_t cond_var) const
{
  const size_t axis = (cond_var == 1) ? 0 : 1;
  const Eigen::VectorXd& g = grid_points_[axis];
  const ptrdiff_t i = find_cell(axis, u_cond);
  return { i, g(i + 1) - u_cond, u_cond - g(i), g(i + 1) - g(i), cond_var };
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
  const size_t axis = free_axis(cond_var);
  const Eigen::VectorXd& g = grid_points_[axis];
  const ptrdiff_t m = g.size();

  // the grid line is walked twice below, so interpolate it once into the
  // caller's buffer
  const CondLine line = cond_line(u_cond, cond_var);
  for (ptrdiff_t j = 0; j < m; ++j) {
    knots(j) = cond_knot(line, j);
  }

  // total mass (normalization of the conditional cdf)
  const double int1 = weights_[axis].dot(knots);
  const double target =
    std::min(std::max(p, 1e-10), 1 - 1e-10) * std::max(int1, 1e-20);

  // locate the bracketing cell and solve the quadratic within it
  double cum = 0.0;
  double v_k = knots(0);
  for (ptrdiff_t k = 0; k < m - 1; ++k) {
    const double v_k1 = knots(k + 1);
    const double g_k = g(k);
    const double dg = g(k + 1) - g_k;
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
  Eigen::VectorXd knots(grid_points_[free_axis(cond_var)].size());
  auto f = [this, cond_var, &knots](double u1, double u2) {
    return (cond_var == 1) ? cond_quantile(u1, u2, 1, knots)
                           : cond_quantile(u2, u1, 2, knots);
  };

  return tools_eigen::binaryExpr_or_nan(u, f);
}

//! Integrate the grid along the two axis
//!
//! @param u Mx2 matrix of evaluation points
//! @return a vector of resulting integral values
inline Eigen::VectorXd
InterpolationGrid::integrate_2d(const tools_eigen::ConstMatRef& u)
{
  Eigen::VectorXd tmpvals2;

  auto f = [this, &tmpvals2](double u1, double u2) {
    row_integrals(u2, tmpvals2);
    double tmpint = int_on_grid(u1, tmpvals2);
    double tmpint1 = weights_[0].dot(tmpvals2);
    return std::min(std::max(tmpint * u2 / tmpint1, 1e-10), 1 - 1e-10);
  };

  return tools_eigen::binaryExpr_or_nan(u, f);
}

//! @brief Nonnegative quadrature weights for the integral over `[lo, hi]`
//! along one axis.
//!
//! @param axis The axis whose knots are integrated over.
//! @param lo,hi Interval bounds, clamped to `[0, 1]`; `hi <= lo` gives zero.
//! @param w Filled with the weights of the nodes the interval covers.
//! @return The index of the first of those nodes, so that
//!   `w.dot(v.segment(first, w.size()))` integrates the piecewise linear
//!   function through `(grid_points_[axis], v)`.
inline ptrdiff_t
InterpolationGrid::interval_weights(size_t axis,
                                    double lo,
                                    double hi,
                                    Eigen::VectorXd& w) const
{
  const Eigen::VectorXd& g = grid_points_[axis];
  const double a = std::min(std::max(lo, 0.0), 1.0);
  const double b = std::min(std::max(hi, a), 1.0);
  const ptrdiff_t ka = find_cell(axis, a);
  const ptrdiff_t kb = find_cell(axis, b);
  w.setZero(kb - ka + 2);

  for (ptrdiff_t k = ka; k <= kb; ++k) {
    const auto [w0, w1] = cell_weights(g(k), g(k + 1), a, b);
    w(k - ka) += w0;
    w(k - ka + 1) += w1;
  }
  return ka;
}

//! @brief Partial integrals of every grid line over `[0, u]` along the second
//! axis.
//!
//! @param u Upper limit, clamped to `[0, 1]`.
//! @param out Filled with one integral per grid line.
inline void
InterpolationGrid::row_integrals(double u, Eigen::VectorXd& out) const
{
  const Eigen::VectorXd& g2 = grid_points_[1];
  const ptrdiff_t m1 = grid_points_[0].size();
  const double y = std::min(std::max(u, 0.0), 1.0);
  const ptrdiff_t j = find_cell(1, y);
  const double dg = g2(j + 1) - g2(j);
  const double s = y - g2(j);
  out.resize(m1);
  for (ptrdiff_t k = 0; k < m1; ++k) {
    out(k) =
      row_cum_int_(k, j) +
      (2 * values_(k, j) + (values_(k, j + 1) - values_(k, j)) * s / dg) * s /
        2.0;
  }
}

//! @brief Probability of the rectangle `(a1, b1] x (a2, b2]`.
//!
//! @details Not clipped, unlike `integrate_2d()`: an empty rectangle is
//! exactly `0`.
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

  // the mass of each grid line over the rectangle's own strip, and below it
  Eigen::VectorXd wx, wy, strip, below;
  const ptrdiff_t i0 = interval_weights(0, x0, x1, wx);
  row_integrals(y0, below);
  if (y0 > 0.0) {
    const ptrdiff_t j0 = interval_weights(1, y0, y1, wy);
    strip = values_.middleCols(j0, wy.size()) * wy;
  } else {
    // nothing below to subtract, so the cached integrals are the strip itself
    row_integrals(y1, strip);
  }

  const Eigen::VectorXd& w1 = weights_[0];
  const double m_strip = w1.dot(strip);
  const double m_below = std::max(w1.dot(below), 1e-20);
  const double total = std::max(m_below + m_strip, 1e-20);

  // `integrate_2d()` rescales each grid line by `lambda(y) = y / M(1, y)`, so
  // the probability is `lambda(y1) strip + (lambda(y1) - lambda(y0)) below`.
  // The lambda difference is the only part that cancels; expanded over the
  // common denominator it cancels against `y1 - y0` rather than against one.
  const double dlambda =
    ((y1 - y0) * m_below - y0 * m_strip) / (total * m_below);

  return y1 * wx.dot(strip.segment(i0, wx.size())) / total +
         dlambda * wx.dot(below.segment(i0, wx.size()));
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
  const size_t axis = free_axis(cond_var);
  const Eigen::VectorXd& g = grid_points_[axis];
  const Eigen::VectorXd& w = weights_[axis];
  const ptrdiff_t m = g.size();
  const double a = std::min(std::max(std::min(lo, hi), 0.0), 1.0);
  const double b = std::min(std::max(std::max(lo, hi), a), 1.0);

  const CondLine line = cond_line(u_cond, cond_var);
  double mass = 0.0, total = 0.0;
  double v_k = cond_knot(line, 0);
  double g_k = g(0);
  for (ptrdiff_t k = 0; k < m - 1; ++k) {
    const double v_k1 = cond_knot(line, k + 1);
    const double g_k1 = g(k + 1);
    total += w(k) * v_k;
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
  total += w(m - 1) * v_k;

  return mass / std::max(total, 1e-20);
}

// ---------------- Utility functions for integration ----------------

//! @brief Integral of the piecewise linear function through
//! `(grid_points_[0], vals)` over `[0, upr]`.
//!
//! @param upr Upper limit, clamped to `[0, 1]`.
//! @param vals One value per knot of the first axis.
inline double
InterpolationGrid::int_on_grid(double upr, const Eigen::VectorXd& vals) const
{
  const Eigen::VectorXd& g = grid_points_[0];
  const double b = std::min(std::max(upr, 0.0), 1.0);
  double total = 0.0;
  double g_k = g(0);
  for (ptrdiff_t k = 0; g_k < b; ++k) {
    const double g_k1 = g(k + 1);
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
