// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>
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
//! the values untouched. The passes stop at convergence, well before the
//! default bound unless the dependence is extreme.
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
  // A grid's transpose, and its lines' integrals both ways, all built by the
  // same code from `values_` and from its transpose, so a transposed grid
  // holds the same arrays with the roles swapped, bit for bit.
  values_t_ = values_.transpose();
  cumulative_lines(values_, row_cum_int_);
  cumulative_lines(values_t_, col_cum_int_);
}

//! @brief Cumulative integrals of every row of `v` along its columns,
//! `cum(k, j) = int_0^{grid_j} v(k, .)`.
inline void
InterpolationGrid::cumulative_lines(const Eigen::MatrixXd& v,
                                    Eigen::MatrixXd& cum) const
{
  const ptrdiff_t m = grid_points_.size();
  cum.resize(m, m);
  for (ptrdiff_t k = 0; k < m; ++k) {
    double total = 0.0;
    cum(k, 0) = 0.0;
    for (ptrdiff_t j = 0; j < m - 1; ++j) {
      total +=
        (v(k, j + 1) + v(k, j)) * (grid_points_(j + 1) - grid_points_(j)) / 2.0;
      cum(k, j + 1) = total;
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
//! The normalization runs to convergence: it stops once the margins' residual,
//! already at the level of rounding, no longer shrinks. A grid left short of
//! uniform margins is not a copula density, and its distribution function,
//! which rescales one argument's margin only, then depends on which argument
//! is first by as much as the residual. Passes converge slowly under strong
//! dependence, so after the first few dozen `newton_margins()` finishes, in a
//! handful of steps.
//!
//! @param max_iter Maximum number of rescaling passes; `0` leaves the values
//! untouched.
inline void
InterpolationGrid::normalize_margins(int max_iter)
{
  const ptrdiff_t m = grid_points_.size();
  if ((max_iter < 1) || (m < 2)) {
    return;
  }

  // at machine precision there is nothing left to do; below `rounding`, a
  // residual that stops shrinking has reached the floor its sums round to
  const double exact = 8 * std::numeric_limits<double>::epsilon();
  const double rounding = 1e-12;
  const double min_mass = 1e-20; // prevent 0/0
  double previous = std::numeric_limits<double>::infinity();
  const Eigen::VectorXd& w = weights_;
  Eigen::MatrixXd vt(m, m);

  // the passes that converge most grids; past them, a pass gains little
  const int newton_after = 25;
  const int newton_steps = 50;

  for (int k = 0; k < max_iter; ++k) {
    // the transpose is materialized rather than left as an expression, so
    // that both margins are the same product on a column-major matrix and
    // transposing the grid swaps them bit for bit
    vt = values_.transpose();
    const Eigen::VectorXd r = (values_ * w).cwiseMax(min_mass);
    const Eigen::VectorXd c = (vt * w).cwiseMax(min_mass);
    const double err = std::max((r.array() - 1.0).abs().maxCoeff(),
                                (c.array() - 1.0).abs().maxCoeff());
    if ((err <= exact) || ((err < rounding) && (err >= previous))) {
      break;
    }
    previous = err;
    if (k == newton_after) {
      if (newton_margins(newton_steps)) {
        break;
      }
      // the passes resume from wherever the steps left the grid
      previous = std::numeric_limits<double>::infinity();
      continue;
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

//! Newton's method for the row and column scalings that make both margins
//! uniform, on their logarithms
//!
//! @details Scaling row \f$ i \f$ by \f$ e^{a_i} \f$ and column \f$ j \f$
//! by \f$ e^{b_j} \f$ moves the log margins by \f$ a + P b \f$ and
//! \f$ b + Q a \f$ to first order, where \f$ P = \mathrm{diag}(1 / r) V W
//! \f$ and \f$ Q = \mathrm{diag}(1 / c) V^\top W \f$ are row stochastic.
//! Eliminating either unknown leaves an \f$ m \times m \f$ system, singular
//! along the scalings of rows against columns that leave the grid as it is:
//! one per block of its support, which a term per block pins. Both
//! eliminations are solved and their steps averaged, so that a step commutes
//! with transposition exactly, as a pass does. Each step is halved until it
//! reduces the residual.
//!
//! @param max_steps Maximum number of steps.
//! @return Whether the margins converged. A step that fails to reduce the
//! residual stops the method, leaving the values of the last one that did.
inline bool
InterpolationGrid::newton_margins(int max_steps)
{
  const ptrdiff_t m = grid_points_.size();
  const double exact = 8 * std::numeric_limits<double>::epsilon();
  const double rounding = 1e-12;
  const double min_mass = 1e-20;
  const Eigen::VectorXd& w = weights_;
  const Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(m, m);

  // both margins as the same product on a column-major matrix, as in a pass
  Eigen::MatrixXd vt = values_.transpose();
  auto residual = [&](const Eigen::MatrixXd& v,
                      const Eigen::MatrixXd& v_t,
                      Eigen::VectorXd& r,
                      Eigen::VectorXd& c) {
    r = (v * w).cwiseMax(min_mass);
    c = (v_t * w).cwiseMax(min_mass);
    if (!r.allFinite() || !c.allFinite()) {
      // a step that overflowed is no improvement, whatever `maxCoeff` would
      // make of the NaN it leaves
      return std::numeric_limits<double>::infinity();
    }
    return std::max((r.array() - 1.0).abs().maxCoeff(),
                    (c.array() - 1.0).abs().maxCoeff());
  };
  // `diag(1 / margin) v W`, row stochastic. Entries too small to move a step
  // are dropped, so that the products of two never fall below the smallest
  // normal number, where arithmetic is orders of magnitude slower; the grids
  // of strongly dependent pairs hold values near it.
  auto stochastic = [&](const Eigen::VectorXd& margin,
                        const Eigen::MatrixXd& v) {
    const Eigen::MatrixXd s =
      margin.cwiseInverse().asDiagonal() * v * w.asDiagonal();
    return Eigen::MatrixXd((s.array() < 1e-150).select(0.0, s));
  };
  // The blocks of the grid's support: row `i` and column `j` are linked when
  // `p` or `q` holds an entry between them, and a block is what the links
  // join. A block can scale its rows against its columns without changing the
  // grid, which leaves each eliminated system singular along that scaling. The
  // pins fix each one with `1 / |B|` between two rows, or two columns, of a
  // block `B`: on a connected grid, both are the rank-one `1 / m`.
  const Eigen::MatrixXd pin =
    Eigen::MatrixXd::Constant(m, m, 1.0 / static_cast<double>(m));
  Eigen::MatrixXd pin_rows(m, m), pin_cols(m, m);
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> linked(m, m);
  std::vector<ptrdiff_t> block(2 * m), rows, cols, todo;
  // whether the grid is connected; if not, sets `pin_rows` and `pin_cols`
  auto connected = [&](const Eigen::MatrixXd& p, const Eigen::MatrixXd& q) {
    linked = (p.array() > 0.0) || (q.transpose().array() > 0.0);
    // nodes `0, ..., m - 1` are the rows and `m, ..., 2m - 1` the columns
    std::fill(block.begin(), block.end(), -1);
    rows.clear();
    cols.clear();
    for (ptrdiff_t start = 0; start < 2 * m; ++start) {
      if (block[start] >= 0) {
        continue;
      }
      const auto label = static_cast<ptrdiff_t>(rows.size());
      rows.push_back(0);
      cols.push_back(0);
      block[start] = label;
      todo.push_back(start);
      while (!todo.empty()) {
        const ptrdiff_t node = todo.back();
        todo.pop_back();
        const bool is_row = node < m;
        ++(is_row ? rows : cols)[label];
        for (ptrdiff_t k = 0; k < m; ++k) {
          const ptrdiff_t other = is_row ? m + k : k;
          if ((block[other] < 0) &&
              (is_row ? linked(node, k) : linked(k, node - m))) {
            block[other] = label;
            todo.push_back(other);
          }
        }
      }
    }
    if (rows.size() == 1) {
      return true;
    }
    for (ptrdiff_t j = 0; j < m; ++j) {
      for (ptrdiff_t i = 0; i < m; ++i) {
        pin_rows(i, j) = (block[i] == block[j])
                           ? 1.0 / static_cast<double>(rows[block[i]])
                           : 0.0;
        pin_cols(i, j) = (block[m + i] == block[m + j])
                           ? 1.0 / static_cast<double>(cols[block[m + i]])
                           : 0.0;
      }
    }
    return false;
  };

  Eigen::VectorXd r(m), c(m), r_try(m), c_try(m);
  double err = residual(values_, vt, r, c);
  double previous = std::numeric_limits<double>::infinity();
  Eigen::MatrixXd trial(m, m), trial_t(m, m);
  for (int step = 0; step < max_steps; ++step) {
    if ((err <= exact) || ((err < rounding) && (err >= previous))) {
      return true;
    }
    const Eigen::VectorXd lr = r.array().log();
    const Eigen::VectorXd lc = c.array().log();
    const Eigen::MatrixXd p = stochastic(r, values_);
    const Eigen::MatrixXd q = stochastic(c, vt);
    const bool whole = connected(p, q);
    const Eigen::VectorXd b1 =
      Eigen::MatrixXd(eye - q * p + (whole ? pin : pin_cols))
        .partialPivLu()
        .solve(q * lr - lc);
    const Eigen::VectorXd a1 = -lr - p * b1;
    const Eigen::VectorXd a2 =
      Eigen::MatrixXd(eye - p * q + (whole ? pin : pin_rows))
        .partialPivLu()
        .solve(p * lc - lr);
    const Eigen::VectorXd b2 = -lc - q * a2;
    const Eigen::VectorXd a = (a1 + a2) / 2.0;
    const Eigen::VectorXd b = (b1 + b2) / 2.0;

    bool improved = false;
    double t = 1.0;
    for (int halving = 0; halving < 30; ++halving, t /= 2.0) {
      const Eigen::VectorXd sr = (t * a).array().exp();
      const Eigen::VectorXd sc = (t * b).array().exp();
      for (ptrdiff_t j = 0; j < m; ++j) {
        for (ptrdiff_t i = 0; i < m; ++i) {
          trial(i, j) = values_(i, j) * (sr(i) * sc(j));
        }
      }
      trial_t = trial.transpose();
      const double err_try = residual(trial, trial_t, r_try, c_try);
      if (err_try < err) {
        values_.swap(trial);
        vt.swap(trial_t);
        r.swap(r_try);
        c.swap(c_try);
        previous = err;
        err = err_try;
        improved = true;
        break;
      }
    }
    if (!improved) {
      // at the floor of rounding, no step can reduce the residual further
      return err < rounding;
    }
  }
  return (err <= exact) || ((err < rounding) && (err >= previous));
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
  // the rows of the transpose are the columns, read by the same expression
  const Eigen::MatrixXd& v = (line.cond_var == 1) ? values_ : values_t_;
  const double knot = (v(i, j) * line.x2x + v(i + 1, j) * line.xx1) / line.x2x1;
  return std::max(knot, 0.0);
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
//! @param u Mx2 matrix of evaluation points
//! @return a vector of resulting integral values
inline Eigen::VectorXd
InterpolationGrid::integrate_2d(const tools_eigen::ConstMatRef& u)
{
  auto f = [this](double u1, double u2) {
    const double a = std::min(std::max(u1, 0.0), 1.0);
    const double b = std::min(std::max(u2, 0.0), 1.0);
    const ptrdiff_t ia = find_cell(a);
    const ptrdiff_t jb = find_cell(b);
    // Rows to `b` swept across to `a`, or columns to `a` across to `b`: the
    // same mass either way round. Which one is a rule that swaps with the
    // arguments, so a grid and its transpose run the same code on the same
    // arrays and agree bit for bit; on the diagonal, both, averaged.
    const auto rows = [&] {
      return sweep(values_, row_cum_int_, a, ia, b, jb);
    };
    const auto cols = [&] {
      return sweep(values_t_, col_cum_int_, b, jb, a, ia);
    };
    const double c = (a < b)   ? rows()
                     : (b < a) ? cols()
                               : 0.5 * (rows() + cols());
    return std::min(std::max(c, 1e-10), 1 - 1e-10);
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
    const double g0 = grid_points_(k);
    const double g1 = grid_points_(k + 1);
    if ((a <= g0) && (b >= g1)) {
      // a whole cell: the trapezoid, with no division
      const double half = 0.5 * (g1 - g0);
      w(k - ka) += half;
      w(k - ka + 1) += half;
    } else {
      const auto [w0, w1] = cell_weights(g0, g1, a, b);
      w(k - ka) += w0;
      w(k - ka + 1) += w1;
    }
  }
  return ka;
}

//! @brief Mass over `[0, a] x [0, b]` of the rows of `v`: each row's
//! integral up to `b`, integrated across the rows up to `a`.
//!
//! @param v The grid, or its transpose.
//! @param cum The cumulative integrals of its rows, from `cumulative_lines()`.
//! @param a,ia Limit across the rows, and its cell from `find_cell()`.
//! @param b,jb Limit along the rows, and its cell from `find_cell()`.
inline double
InterpolationGrid::sweep(const Eigen::MatrixXd& v,
                         const Eigen::MatrixXd& cum,
                         double a,
                         ptrdiff_t ia,
                         double b,
                         ptrdiff_t jb) const
{
  double total = 0.0;
  double l_k = line_integral(v, cum, 0, jb, b);
  for (ptrdiff_t k = 0; k < ia; ++k) {
    const double l_k1 = line_integral(v, cum, k + 1, jb, b);
    total += (l_k1 + l_k) * (grid_points_(k + 1) - grid_points_(k)) / 2.0;
    l_k = l_k1;
  }
  const auto [w0, w1] =
    cell_weights(grid_points_(ia), grid_points_(ia + 1), 0.0, a);
  return total + w0 * l_k + w1 * line_integral(v, cum, ia + 1, jb, b);
}

//! @brief Integral over `[0, upr]` of row `k` of `v`.
//!
//! @param v The grid, or its transpose.
//! @param cum The cumulative integrals of its rows.
//! @param k The row.
//! @param j The cell holding `upr`, from `find_cell()`.
//! @param upr Upper limit, in `[0, 1]`.
inline double
InterpolationGrid::line_integral(const Eigen::MatrixXd& v,
                                 const Eigen::MatrixXd& cum,
                                 ptrdiff_t k,
                                 ptrdiff_t j,
                                 double upr) const
{
  const double dg = grid_points_(j + 1) - grid_points_(j);
  const double s = upr - grid_points_(j);
  return cum(k, j) + (2 * v(k, j) + (v(k, j + 1) - v(k, j)) * s / dg) * s / 2.0;
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
  // per-thread buffers: this runs once per observation of a discrete edge
  thread_local Eigen::VectorXd wx, wy;
  const ptrdiff_t i0 = interval_weights(x0, x1, wx);
  const ptrdiff_t j0 = interval_weights(y0, y1, wy);
  // `wx' V wy` over the covered nodes, rows first or columns first: a sum of
  // nonnegative terms either way, so nothing cancels however narrow the
  // rectangle. Which one is a rule that swaps with the arguments, so a grid
  // and its transpose run the same code on the same arrays and agree bit for
  // bit; for a rectangle on the diagonal, both, averaged.
  const auto rows = [&] { return block_mass(values_, i0, wx, j0, wy); };
  const auto cols = [&] { return block_mass(values_t_, j0, wy, i0, wx); };
  if ((x0 < y0) || ((x0 == y0) && (x1 < y1))) {
    return rows();
  }
  if ((y0 < x0) || ((y0 == x0) && (y1 < x1))) {
    return cols();
  }
  return 0.5 * (rows() + cols());
}

//! @brief `wa' v[i0:, j0:] wb`, each row's sum first.
inline double
InterpolationGrid::block_mass(const Eigen::MatrixXd& v,
                              ptrdiff_t i0,
                              const Eigen::VectorXd& wa,
                              ptrdiff_t j0,
                              const Eigen::VectorXd& wb)
{
  double total = 0.0;
  for (ptrdiff_t i = 0; i < wa.size(); ++i) {
    double line = 0.0;
    for (ptrdiff_t j = 0; j < wb.size(); ++j) {
      line += v(i0 + i, j0 + j) * wb(j);
    }
    total += wa(i) * line;
  }
  return total;
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

}
}
