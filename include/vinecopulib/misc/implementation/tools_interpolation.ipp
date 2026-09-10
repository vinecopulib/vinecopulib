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

//! cumulative trapezoidal integrals of each row (used by `integrate_2d`).
inline void
InterpolationGrid::update_cached_integrals()
{
  const ptrdiff_t m = grid_points_.size();
  row_cum_int_.resize(m, m);
  for (ptrdiff_t k = 0; k < m; ++k) {
    double cum = 0.0;
    row_cum_int_(k, 0) = 0.0;
    for (ptrdiff_t j = 0; j < m - 1; ++j) {
      cum += (values_(k, j + 1) + values_(k, j)) *
             (grid_points_(j + 1) - grid_points_(j)) / 2.0;
      row_cum_int_(k, j + 1) = cum;
    }
  }
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
  values_.transposeInPlace();
  update_cached_integrals();
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

inline Eigen::Matrix<ptrdiff_t, 1, 2>
InterpolationGrid::get_indices(double x0, double x1)
{
  Eigen::Matrix<ptrdiff_t, 1, 2> out;
  out(0) = this->find_cell(x0);
  out(1) = this->find_cell(x1);
  return out;
}

//! Interpolate linearly in two dimensions
//!
//! @param z11 Value corresponding to (x1, y1)
//! @param z12 Value corresponding to (x1, y2)
//! @param z21 Value corresponding to (x2, y1)
//! @param z22 Value corresponding to (x2, y2)
//! @param x1 First cell value for the first dimension
//! @param x2 Second cell value for the first dimension
//! @param y1 First cell value for the second dimension
//! @param y2 Second cell value for the second dimension
//! @param x Evaluation point for the first dimension
//! @param y Evaluation point for the second dimension
inline double
InterpolationGrid::bilinear_interpolation(double z11,
                                          double z12,
                                          double z21,
                                          double z22,
                                          double x1,
                                          double x2,
                                          double y1,
                                          double y2,
                                          double x,
                                          double y)
{
  double x2x1, y2y1, x2x, y2y, yy1, xx1;
  x2x1 = x2 - x1;
  y2y1 = y2 - y1;
  x2x = x2 - x;
  y2y = y2 - y;
  yy1 = y - y1;
  xx1 = x - x1;
  return (z11 * x2x * y2y + z21 * xx1 * y2y + z12 * x2x * yy1 +
          z22 * xx1 * yy1) /
         (x2x1 * y2y1);
}

//! Interpolation in two dimensions
//!
//! @param x Mx2 matrix of evaluation points.
//! @return a vector of resulting interpolated values
inline Eigen::VectorXd
InterpolationGrid::interpolate(const tools_eigen::ConstMatRef& x)
{

  auto f = [this](double x0, double x1) {
    auto indices = this->get_indices(x0, x1);
    return bilinear_interpolation(this->values_(indices(0), indices(1)),
                                  this->values_(indices(0), indices(1) + 1),
                                  this->values_(indices(0) + 1, indices(1)),
                                  this->values_(indices(0) + 1, indices(1) + 1),
                                  this->grid_points_(indices(0)),
                                  this->grid_points_(indices(0) + 1),
                                  this->grid_points_(indices(1)),
                                  this->grid_points_(indices(1) + 1),
                                  x0,
                                  x1);
  };

  return tools_eigen::binaryExpr_or_nan(x, f);
}

//! conditional cdf along one axis, fused: one cell search for the
//! conditioning coordinate, the interpolated knot values on the fly, and a
//! single pass accumulating both the partial and the full integral (no
//! allocation per query).
//!
//! @param u_cond The conditioning coordinate.
//! @param u The coordinate up to which the density is integrated.
//! @param cond_var Either 1 or 2; the axis considered fixed.
inline double
InterpolationGrid::cond_cdf(double u_cond, double u, size_t cond_var) const
{
  const ptrdiff_t m = grid_points_.size();
  const ptrdiff_t i = find_cell(u_cond);
  const double x1 = grid_points_(i);
  const double x2 = grid_points_(i + 1);
  const double x2x = x2 - u_cond;
  const double xx1 = u_cond - x1;
  const double x2x1 = x2 - x1;

  // interpolated grid-line value at knot j; bilinear interpolation of a
  // nonnegative grid is nonnegative, so the guard only absorbs rounding
  auto knot = [&](ptrdiff_t j) {
    double v;
    if (cond_var == 1) {
      v = (values_(i, j) * x2x + values_(i + 1, j) * xx1) / x2x1;
    } else {
      v = (values_(j, i) * x2x + values_(j, i + 1) * xx1) / x2x1;
    }
    return std::max(v, 0.0);
  };

  double tmpint = 0.0, int1 = 0.0;
  double v_k = knot(0);
  const bool do_partial = (u > grid_points_(0));
  for (ptrdiff_t k = 0; k < m - 1; ++k) {
    const double v_k1 = knot(k + 1);
    const double g_k = grid_points_(k);
    const double g_k1 = grid_points_(k + 1);
    int1 += weights_(k) * v_k;
    if (do_partial && !(u < g_k)) {
      if (u < g_k1) {
        tmpint +=
          (2 * v_k + (v_k1 - v_k) * (u - g_k) / (g_k1 - g_k)) * (u - g_k) / 2.0;
      } else {
        tmpint += (v_k1 + v_k) * (g_k1 - g_k) / 2.0;
      }
    }
    v_k = v_k1;
  }
  int1 += weights_(m - 1) * v_k;

  return std::min(std::max(tmpint / std::max(int1, 1e-20), 1e-10), 1 - 1e-10);
}

//! inverts `cond_cdf` in its second argument: the conditional cdf is
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
  const ptrdiff_t i = find_cell(u_cond);
  const double x1 = grid_points_(i);
  const double x2 = grid_points_(i + 1);
  const double x2x = x2 - u_cond;
  const double xx1 = u_cond - x1;
  const double x2x1 = x2 - x1;

  // the grid line is walked twice below, so interpolate it once into the
  // caller's buffer
  for (ptrdiff_t j = 0; j < m; ++j) {
    const double v = (cond_var == 1)
                       ? (values_(i, j) * x2x + values_(i + 1, j) * xx1) / x2x1
                       : (values_(j, i) * x2x + values_(j, i + 1) * xx1) / x2x1;
    knots(j) = std::max(v, 0.0);
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
    return (cond_var == 1) ? cond_cdf(u1, u2, 1) : cond_cdf(u2, u1, 2);
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
//! @param u Mx2 matrix of evaluation points
//! @return a vector of resulting integral values
inline Eigen::VectorXd
InterpolationGrid::integrate_2d(const tools_eigen::ConstMatRef& u)
{
  ptrdiff_t m = grid_points_.size();
  Eigen::VectorXd tmpvals2(m);

  auto f = [this, m, &tmpvals2](double u1, double u2) {
    // partial row integrals up to u2, from the cached cumulative integrals
    // plus the remaining partial cell (one O(m) pass instead of m
    // interpolation sweeps per query)
    const ptrdiff_t j = find_cell(u2);
    const double g_j = grid_points_(j);
    const double dg = grid_points_(j + 1) - g_j;
    const double s = u2 - g_j;
    for (ptrdiff_t k = 0; k < m; ++k) {
      double partial = 0.0;
      if (u2 > grid_points_(0)) {
        partial =
          (2 * values_(k, j) + (values_(k, j + 1) - values_(k, j)) * s / dg) *
          s / 2.0;
      }
      tmpvals2(k) = row_cum_int_(k, j) + partial;
    }
    double tmpint = int_on_grid(u1, tmpvals2, grid_points_);
    double tmpint1 = weights_.dot(tmpvals2);
    return std::min(std::max(tmpint * u2 / tmpint1, 1e-10), 1 - 1e-10);
  };

  return tools_eigen::binaryExpr_or_nan(u, f);
}

//! nonnegative quadrature weights for the integral over `[lo, hi]`.
//!
//! @details Fills `w` with the weights of the nodes the interval covers and
//! returns the index of the first of them, so that
//! `w.dot(v.segment(first, w.size()))` is exactly the integral over `[lo, hi]`
//! of the piecewise linear function through `(grid_points_, v)`. Every weight
//! is nonnegative, because the hat functions are, which is what makes a mass
//! built from them free of cancellation however narrow the interval. Endpoints
//! are clamped to `[0, 1]`, which the constructor makes the grid's own range,
//! and an empty or inverted interval gives zero weights.
inline ptrdiff_t
InterpolationGrid::interval_weights(double lo,
                                    double hi,
                                    Eigen::VectorXd& w) const
{
  const ptrdiff_t m = grid_points_.size();
  const double a = std::min(std::max(lo, grid_points_(0)), grid_points_(m - 1));
  const double b = std::min(std::max(hi, a), grid_points_(m - 1));
  const ptrdiff_t ka = find_cell(a);
  const ptrdiff_t kb = find_cell(b);
  w.setZero(kb - ka + 2);

  // a partial cell is parameterized by its width rather than by its upper end,
  // so a narrow interval is never formed as a difference of two numbers of
  // order one -- which would put the 1 / width amplification back into the
  // weights the exact route exists to avoid
  auto add_cell = [&](ptrdiff_t k, double s0, double d) {
    const double h = grid_points_(k + 1) - grid_points_(k);
    const double q = 0.5 * d * (2.0 * s0 + d);
    w(k - ka) += h * (d - q);
    w(k - ka + 1) += h * q;
  };

  if (ka == kb) {
    const double h = grid_points_(ka + 1) - grid_points_(ka);
    add_cell(ka, (a - grid_points_(ka)) / h, (b - a) / h);
    return ka;
  }

  const double ha = grid_points_(ka + 1) - grid_points_(ka);
  const double s0 = (a - grid_points_(ka)) / ha;
  add_cell(ka, s0, 1.0 - s0);
  const double hb = grid_points_(kb + 1) - grid_points_(kb);
  add_cell(kb, 0.0, (b - grid_points_(kb)) / hb);
  // the trapezoid rule is exact for a piecewise linear integrand, so each whole
  // cell in between contributes half its width to both of its nodes
  for (ptrdiff_t k = ka + 1; k < kb; ++k) {
    const double h = grid_points_(k + 1) - grid_points_(k);
    w(k - ka) += h / 2.0;
    w(k - ka + 1) += h / 2.0;
  }
  return ka;
}

//! partial integrals of every grid line up to `u`, from the cached cumulative
//! integrals plus the remaining partial cell (one pass over the grid lines
//! instead of one interpolation sweep each). `u` is clamped as in
//! `interval_weights()`.
inline void
InterpolationGrid::row_integrals(double u, Eigen::VectorXd& out) const
{
  const ptrdiff_t m = grid_points_.size();
  const double y = std::min(std::max(u, grid_points_(0)), grid_points_(m - 1));
  const ptrdiff_t j = find_cell(y);
  const double g_j = grid_points_(j);
  const double dg = grid_points_(j + 1) - g_j;
  const double s = y - g_j;
  out.resize(m);
  for (ptrdiff_t k = 0; k < m; ++k) {
    out(k) =
      row_cum_int_(k, j) +
      (2 * values_(k, j) + (values_(k, j + 1) - values_(k, j)) * s / dg) * s /
        2.0;
  }
}

//! probability of `(a1, b1] x (a2, b2]`, the value a difference of four
//! `integrate_2d()` values defines, arranged so that almost none of it cancels.
//!
//! @details Differencing four corners turns an absolute error `eps` into
//! `~4 eps / (w1 w2)` in the rectangle's widths, which a mixed-discrete density
//! then divides by `w1 w2` again. Writing `M` for the interpolant's own mass
//! and `lambda(y) = y / M(1, y)` for the rescaling `integrate_2d()` applies to
//! each grid line, that difference is identically `lambda(b2) R + (lambda(b2) -
//! lambda(a2)) S`, with `R` the rectangle's mass and `S` the mass of `(a1, b1]
//! x (0, a2]`. Both are sums of nonnegative terms against a nonnegative grid
//! and do not cancel at all; only the two `lambda` cancel, and they multiply a
//! term of order `w1` rather than one of order one. The amplification is
//! therefore `1 / w2` alone, one power instead of two.
//!
//! @param a1,b1 Bounds in the first argument, in either order.
//! @param a2,b2 Bounds in the second argument, in either order.
//! @return The probability, or `0` for an empty rectangle.
inline double
InterpolationGrid::rect_mass(double a1, double b1, double a2, double b2) const
{
  // rotating the data can leave a left limit above its own value, so the
  // rectangle is oriented here rather than by taking an absolute value of the
  // result, which would hide a sign error instead of preventing one
  const double x0 = std::min(a1, b1);
  const double x1 = std::max(a1, b1);
  const double y0 = std::min(a2, b2);
  const double y1 = std::max(a2, b2);
  if (!(x1 > x0) || !(y1 > y0)) {
    return 0.0;
  }

  Eigen::VectorXd wx, wy, rows;
  const ptrdiff_t i0 = interval_weights(x0, x1, wx);
  const ptrdiff_t j0 = interval_weights(y0, y1, wy);
  const auto strip = values_.middleCols(j0, wy.size());

  // the rectangle's own mass, and the mass of the whole strip it sits in
  const double mass = wx.dot(strip.middleRows(i0, wx.size()) * wy);
  const double strip_mass = weights_.dot(strip * wy);

  if (!(y0 > grid_points_(0))) {
    // no lower corner to subtract: the probability is the mass rescaled so
    // that the strip carries exactly its own marginal probability
    return y1 * mass / std::max(strip_mass, 1e-20);
  }

  // Below the rectangle, the mass of `(x0, x1] x (0, y0]` and its marginal
  // counterpart. `lambda(y) = y / M(1, y)` is the rescaling the distribution
  // function applies to each grid line, and the probability is
  // `lambda(y1) mass + (lambda(y1) - lambda(y0)) below`. Only the difference of
  // the two lambdas cancels, and writing it over the strip's own mass keeps
  // that cancellation proportional to the strip's width rather than to one: the
  // two terms of the numerator agree to the extent that the fitted margin is
  // uniform, and a rectangle the model gives almost no mass to would otherwise
  // come out with the wrong sign.
  row_integrals(y0, rows);
  const double below_marginal = std::max(weights_.dot(rows), 1e-20);
  const double below = wx.dot(rows.segment(i0, wx.size()));
  const double marginal = std::max(below_marginal + strip_mass, 1e-20);
  const double dlambda = ((y1 - y0) * below_marginal - y0 * strip_mass) /
                         (marginal * below_marginal);

  return y1 * mass / marginal + dlambda * below;
}

//! probability that the free coordinate falls in `(lo, hi]`, given the other.
//!
//! @details The value a difference of two `integrate_1d()` values defines. A
//! conditional distribution function is the interpolated grid line divided by
//! its own total, so this is a ratio of two nonnegative sums and does not
//! cancel at all -- a partition of the free axis therefore sums to exactly `1`.
//!
//! @param u_cond The coordinate held fixed.
//! @param lo,hi Bounds in the free coordinate, in either order.
//! @param cond_var Either 1 or 2; the axis considered fixed, as for
//!   `integrate_1d()`.
inline double
InterpolationGrid::cond_interval_mass(double u_cond,
                                      double lo,
                                      double hi,
                                      size_t cond_var) const
{
  const ptrdiff_t m = grid_points_.size();
  const ptrdiff_t i = find_cell(u_cond);
  const double x1 = grid_points_(i);
  const double x2 = grid_points_(i + 1);
  const double x2x = x2 - u_cond;
  const double xx1 = u_cond - x1;
  const double x2x1 = x2 - x1;

  // the grid line at the conditioning coordinate; bilinear interpolation of a
  // nonnegative grid is nonnegative, so the guard only absorbs rounding
  Eigen::VectorXd knots(m);
  for (ptrdiff_t j = 0; j < m; ++j) {
    const double v = (cond_var == 1)
                       ? (values_(i, j) * x2x + values_(i + 1, j) * xx1) / x2x1
                       : (values_(j, i) * x2x + values_(j, i + 1) * xx1) / x2x1;
    knots(j) = std::max(v, 0.0);
  }

  Eigen::VectorXd w;
  const ptrdiff_t j0 = interval_weights(std::min(lo, hi), std::max(lo, hi), w);
  return w.dot(knots.segment(j0, w.size())) /
         std::max(weights_.dot(knots), 1e-20);
}

// ---------------- Utility functions for integration ----------------

//! Integrate using a trapezoid rule
//!
//! @param upr Upper limit of integration (lower is 0).
//! @param vals Vector of values to be interpolated and integrated.
//! @param grid Vector of grid points on which vals has been computed.
//!
//! @return integral of a piecewise linear function defined by (grid_i, vals_i).
inline double
InterpolationGrid::int_on_grid(const double& upr,
                               const Eigen::VectorXd& vals,
                               const Eigen::VectorXd& grid)
{
  double tmpint = 0.0;

  if (upr > grid(0)) {
    // go up the grid and integrate
    for (ptrdiff_t k = 0; k < (grid.size() - 1); ++k) {
      // stop loop if fully integrated
      if (upr < grid(k))
        break;

      // don't integrate over full cell if upr is in interior
      if (upr < grid(k + 1)) {
        tmpint += (2 * vals(k) + (vals(k + 1) - vals(k)) * (upr - grid(k)) /
                                   (grid(k + 1) - grid(k))) *
                  (upr - grid(k)) / 2.0;
      } else {
        tmpint += (vals(k + 1) + vals(k)) * (grid(k + 1) - grid(k)) / 2.0;
      }
    }
  }

  return tmpint;
}
}
}
