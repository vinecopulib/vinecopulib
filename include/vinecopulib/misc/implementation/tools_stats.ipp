// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <algorithm>
#include <array>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/seed_seq.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <unsupported/Eigen/FFT>
#include <vinecopulib/misc/tools_stats_ghalton.hpp>
#include <vinecopulib/misc/tools_stats_sobol.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <wdm/eigen.hpp>
#include <wdm/ranks.hpp>

namespace vinecopulib {

//! Utilities for statistical analysis
namespace tools_stats {

//! @brief Simulates from the multivariate uniform distribution.
//!
//! If `qrng = TRUE`, generalized Halton sequences (see `ghalton()`) are used
//! for \f$ d \leq 300 \f$ and Sobol sequences otherwise (see `sobol()`).
//!
//! @param n Number of observations.
//! @param d Dimension.
//! @param qrng If true, quasi-numbers are generated.
//! @param seeds Seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//! @return An \f$ n \times d \f$ matrix of independent
//! \f$ \mathrm{U}[0, 1] \f$ random variables.
inline Eigen::MatrixXd
simulate_uniform(const size_t& n,
                 const size_t& d,
                 bool qrng,
                 std::vector<int> seeds)
{
  if (qrng) {
    if (d > 300) {
      return tools_stats::sobol(n, d, seeds);
    } else {
      return tools_stats::ghalton(n, d, seeds);
    }
  }
  if ((n < 1) || (d < 1)) {
    throw std::runtime_error("n and d must be at least 1.");
  }
  if (seeds.size() == 0) {
    // no seeds provided, seed randomly
    std::random_device rd{};
    seeds = std::vector<int>(20);
    std::generate(
      seeds.begin(), seeds.end(), [&]() { return static_cast<int>(rd()); });
  }

  // initialize random engine and uniform distribution
  boost::random::seed_seq seq(seeds.begin(), seeds.end());
  boost::random::mt19937 generator(seq);
  boost::random::uniform_real_distribution<double> distribution(0.0, 1.0);

  // NullaryExpr fills the result directly (column-major, same order as the
  // previous unaryExpr-based version) without a second allocation
  return Eigen::MatrixXd::NullaryExpr(
    n, d, [&]() { return distribution(generator); });
}

//! @brief Simulates from independendent normals.
//!
//! @param n Number of observations.
//! @param d Dimension.
//! @param qrng If true, quasi-numbers are generated.
//! @param seeds Seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//!//!
//! @return An \f$ n \times d \f$ matrix of independent
//! \f$ \mathrm{N}(0, 1) \f$ random variables.
inline Eigen::MatrixXd
simulate_normal(const size_t& n,
                const size_t& d,
                bool qrng,
                std::vector<int> seeds)
{
  return qnorm(tools_stats::simulate_uniform(n, d, qrng, seeds));
}

// (internal) 1-d worker with pre-converted weights; moves `xvec` into
// `wdm::impl::rank` to avoid a copy.
inline Eigen::VectorXd
pseudo_obs_1d_impl(std::vector<double>&& xvec,
                   const std::string& ties_method,
                   const std::vector<double>& weights,
                   const std::vector<int>& seeds)
{
  // correction for NaNs (must be counted before the move)
  size_t n = xvec.size();
  for (size_t i = 0; i < xvec.size(); i++) {
    if (std::isnan(xvec[i])) {
      n--;
    }
  }
  auto res = wdm::impl::rank(std::move(xvec), weights, ties_method, seeds);
  return Eigen::Map<Eigen::VectorXd>(res.data(), res.size()) /
         (static_cast<double>(n) + 1.0);
}

//! @brief Applies the empirical probability integral transform to a data
//! matrix.
//!
//! Gives pseudo-observations from the copula by applying the empirical
//! distribution function (scaled by \f$ n + 1 \f$) to each margin/column.
//!
//! @param x A matrix of real numbers.
//! @param ties_method Indicates how to treat ties; same as in R, see
//! https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.
//! @param weights Vector of weights for the observations.
//! @param seeds Seeds for the random number generator, used only when
//! `ties_method = "random"`.
//! @return Pseudo-observations of the copula, i.e. \f$ F_X(x) \f$
//! (column-wise).
inline Eigen::MatrixXd
to_pseudo_obs(Eigen::MatrixXd x,
              const std::string& ties_method,
              const Eigen::VectorXd& weights,
              std::vector<int> seeds)
{
  // convert the weights once instead of once per column
  const auto wvec = wdm::utils::convert_vec(weights);
  const size_t n = x.rows();
  for (int j = 0; j < x.cols(); ++j) {
    std::vector<double> xvec(x.data() + n * j, x.data() + n * (j + 1));
    x.col(j) = pseudo_obs_1d_impl(std::move(xvec), ties_method, wvec, seeds);
  }

  return x;
}

//! @brief Applies the empirical probability integral transform to a data
//! vector.
//!
//! Gives pseudo-observations from the copula by applying the empirical
//! distribution function (scaled by \f$ n + 1 \f$) to each margin/column.
//!
//! @param x A vector of real numbers.
//! @param ties_method Indicates how to treat ties; same as in R, see
//! https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.
//! @param weights Vector of weights for the observations.
//! @param seeds Seeds for the random number generator, used only when
//! `ties_method = "random"`.
//! @return Pseudo-observations of the copula, i.e. \f$ F_X(x) \f$.
inline Eigen::VectorXd
to_pseudo_obs_1d(Eigen::VectorXd x,
                 const std::string& ties_method,
                 const Eigen::VectorXd& weights,
                 std::vector<int> seeds)
{
  return pseudo_obs_1d_impl(wdm::utils::convert_vec(x),
                            ties_method,
                            wdm::utils::convert_vec(weights),
                            seeds);
}

//! @brief The distance below which `soft_pseudo_obs()` does not tell values
//! apart by value alone: the square root of the machine epsilon, the usual
//! scale for what rounding in a computation can move.
inline double
default_soft_scale()
{
  return std::sqrt(std::numeric_limits<double>::epsilon());
}

//! @brief Whether a pair is in the order of its own values.
//!
//! Anything that treats a pair's two arguments by position, a draw or an
//! iteration indexed by column, is made a function of the pair rather than
//! of how it was passed by putting the pair in this order first: at the first
//! row whose two values (then, with four columns, whose two left limits)
//! differ, the smaller first.
//!
//! @param u A pair, `[u1, u2]` or `[u1, u2, u1^-, u2^-]`.
//! @return Whether the columns have to be swapped.
inline bool
swaps_pair(const Eigen::MatrixXd& u)
{
  for (Eigen::Index i = 0; i < u.rows(); ++i) {
    if (u(i, 0) != u(i, 1)) {
      return u(i, 1) < u(i, 0);
    }
    if ((u.cols() == 4) && (u(i, 2) != u(i, 3))) {
      return u(i, 3) < u(i, 2);
    }
  }
  return false;
}

//! @brief Pseudo-observations whose ranks move continuously with the data.
//!
//! Ranks every column as `to_pseudo_obs(x, "random", weights, seeds)` does,
//! with two differences that make the ranks a continuous function of the
//! data. Tied values are ordered by a key per observation, drawn once from
//! `seeds`: a uniformly random order within every tie group, as with
//! `"random"`, but one that depends only on the group's members. And two
//! values closer than a few multiples of `scale` \f$ s \f$ are ranked partly
//! by key and partly by value: observation \f$ i \f$ counts as ranked after
//! \f$ j \f$ with weight
//! \f[
//!   p_{ij} = \kappa_{ij} \, [k_i > k_j] + (1 - \kappa_{ij}) \,
//!            \Phi(g_{ij} / s), \qquad
//!   \kappa_{ij} = e^{-g_{ij}^2 / (2 s^2)}, \quad g_{ij} = x_i - x_j,
//! \f]
//! so that data moving by \f$ \delta \f$ move a rank by \f$ O(\delta / s) \f$
//! rather than by a whole block when near-equal values reorder. Values further
//! apart than about \f$ 9 s \f$ are ranked exactly by value, so data without
//! such near-ties get the ordinary ranks.
//!
//! @param x Data, one variable per column.
//! @param weights Optional weights, one per observation.
//! @param seeds Seeds of the keys.
//! @param scale The distance \f$ s \f$; zero ranks exactly by value, only the
//!   ties being ordered by `ties_method`.
//! @param ties_method `"random"` orders ties by the keys, as above;
//!   `"average"` counts every tie-mate as half ranked before, which gives
//!   tied values their average rank, and uses \f$ [k_i > k_j] = 1/2 \f$ in
//!   \f$ p_{ij} \f$ likewise.
//! @return \f$ r_i / (n + 1) \f$, with \f$ r_i = w_i + \sum_{j \neq i} w_j
//!   p_{ij} \f$ the weighted rank and \f$ n \f$ the number of non-`NaN`
//!   values; `NaN` where `x` is.
inline Eigen::MatrixXd
soft_pseudo_obs(const Eigen::MatrixXd& x,
                const Eigen::VectorXd& weights,
                const std::vector<int>& seeds,
                double scale,
                const std::string& ties_method)
{
  if ((ties_method != "random") && (ties_method != "average")) {
    throw std::runtime_error("ties_method must be 'random' or 'average'.");
  }
  const bool random = (ties_method == "random");
  // a zero scale ranks exactly by value, the ties aside
  const bool soft = (scale > 0.0);
  const Eigen::Index n = x.rows();
  Eigen::MatrixXd out = Eigen::MatrixXd::Constant(
    n, x.cols(), std::numeric_limits<double>::quiet_NaN());
  // beyond `hard`, both terms of p_ij are exactly 0 or 1 in double precision
  const double hard = 9.0;
  const double inv_scale = 1.0 / scale;
  const double inv_sqrt2 = 0.70710678118654752440;

  for (Eigen::Index col = 0; col < x.cols(); ++col) {
    // one key per observation and column, from the engine's raw output, which
    // is specified exactly; a column of its own keeps the tie orders of two
    // columns independent. "average" reads no keys; they are all zero there,
    // which leaves the order to the value and the index.
    std::vector<uint64_t> key(static_cast<size_t>(n), 0);
    if (random) {
      std::vector<int> col_seeds = seeds;
      col_seeds.push_back(static_cast<int>(col));
      boost::random::seed_seq seq(col_seeds.begin(), col_seeds.end());
      boost::random::mt19937 engine(seq);
      for (auto& k : key) {
        const uint64_t high = engine();
        k = (high << 32) | static_cast<uint64_t>(engine());
      }
    }

    std::vector<Eigen::Index> order;
    for (Eigen::Index i = 0; i < n; ++i) {
      if (!std::isnan(x(i, col))) {
        order.push_back(i);
      }
    }
    const double m = static_cast<double>(order.size());
    std::vector<double> w(static_cast<size_t>(n), 1.0);
    if (weights.size() > 0) {
      double total = 0.0;
      for (auto i : order) {
        total += weights(i);
      }
      for (auto i : order) {
        w[i] = weights(i) * m / total;
      }
    }
    // a total order: value, then key, then index
    std::sort(order.begin(), order.end(), [&](Eigen::Index a, Eigen::Index b) {
      if (x(a, col) != x(b, col)) {
        return x(a, col) < x(b, col);
      }
      if (key[a] != key[b]) {
        return key[a] < key[b];
      }
      return a < b;
    });

    // groups of equal values, in order; `prefix` accumulates their weights
    std::vector<size_t> group;
    std::vector<double> prefix(order.size() + 1, 0.0);
    for (size_t k = 0; k < order.size(); ++k) {
      if ((k == 0) || (x(order[k], col) != x(order[k - 1], col))) {
        group.push_back(k);
      }
      prefix[k + 1] = prefix[k] + w[order[k]];
    }
    group.push_back(order.size());

    // `lo` / `hi`: the groups within `hard * scale` of group `a`; every group
    // below `lo` ranks before it and every group from `hi` on after it
    size_t lo = 0, hi = 0;
    for (size_t a = 0; a + 1 < group.size(); ++a) {
      const double va = x(order[group[a]], col);
      if (soft) {
        while ((va - x(order[group[lo]], col)) * inv_scale > hard) {
          ++lo;
        }
        while ((hi + 1 < group.size()) &&
               (x(order[group[hi]], col) - va) * inv_scale <= hard) {
          ++hi;
        }
      } else {
        lo = a;
        hi = a + 1;
      }
      for (size_t k = group[a]; k < group[a + 1]; ++k) {
        const Eigen::Index i = order[k];
        // below the window; then the tie-mates, those with smaller keys
        // under "random", which come first in the group, or half of them all
        // under "average"
        const double mates =
          random ? (prefix[k] - prefix[group[a]])
                 : 0.5 * (prefix[group[a + 1]] - prefix[group[a]] - w[i]);
        double r = prefix[group[lo]] + w[i] + mates;
        for (size_t b = lo; b < hi; ++b) {
          if (b == a) {
            continue;
          }
          const double gs = (va - x(order[group[b]], col)) * inv_scale;
          const double kappa = std::exp(-0.5 * gs * gs);
          const double phi = 0.5 * std::erfc(-gs * inv_sqrt2);
          // group b's weight below i's key: the group is in key order
          auto first = order.begin() + static_cast<std::ptrdiff_t>(group[b]);
          auto last = order.begin() + static_cast<std::ptrdiff_t>(group[b + 1]);
          auto pos = std::lower_bound(
            first, last, key[i], [&](Eigen::Index j, uint64_t kk) {
              return key[j] < kk;
            });
          const double whole = prefix[group[b + 1]] - prefix[group[b]];
          const double below =
            random ? prefix[static_cast<size_t>(pos - order.begin())] -
                       prefix[group[b]]
                   : 0.5 * whole;
          r += kappa * below + (1.0 - kappa) * phi * whole;
        }
        out(i, col) = r / (m + 1.0);
      }
    }
  }
  return out;
}

//! @brief Pseudo-observations of a pair, as a kernel pair copula ranks it.
//!
//! `soft_pseudo_obs()` of the pair's first two columns, drawn in the pair's
//! own order (`swaps_pair()`) and returned in the order passed: each column's
//! keys come from a stream of its own, and this keeps the result a function
//! of the pair rather than of the order of its arguments.
//!
//! @param data The pair, `[u1, u2]` or `[u1, u2, u1^-, u2^-]`.
//! @param weights Optional weights, one per observation.
//! @param scale As for `soft_pseudo_obs()`.
//! @param seeds Seeds of the keys.
//! @return An \f$ n \times 2 \f$ matrix of pseudo-observations.
inline Eigen::MatrixXd
pair_soft_pseudo_obs(const Eigen::MatrixXd& data,
                     const Eigen::VectorXd& weights,
                     double scale,
                     const std::vector<int>& seeds)
{
  const bool swapped = swaps_pair(data);
  Eigen::MatrixXd pair = data.leftCols(2);
  if (swapped) {
    pair.col(0).swap(pair.col(1));
  }
  Eigen::MatrixXd psobs = soft_pseudo_obs(pair, weights, seeds, scale);
  if (swapped) {
    psobs.col(0).swap(psobs.col(1));
  }
  return psobs;
}

//! @brief The number of points close to each point, counted continuously.
//!
//! \f$ m_i = \sum_j e^{-\|x_i - x_j\|^2 / (2 s^2)} \f$: an exactly repeated
//! point counts its copies, a point far from all others counts one, and
//! points moving by \f$ \delta \ll s \f$ change the counts by
//! \f$ O(\delta / s) \f$. Summing \f$ f_i / m_i \f$ over points sums \f$ f \f$
//! over the distinct points, continuously in the data.
//!
//! @param x Points, one per row.
//! @param scale The distance \f$ s \f$.
//! @return The counts, each at least one.
inline Eigen::VectorXd
soft_multiplicity(const Eigen::MatrixXd& x, double scale)
{
  const Eigen::Index n = x.rows();
  std::vector<Eigen::Index> order(static_cast<size_t>(n));
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](Eigen::Index a, Eigen::Index b) {
    if (x(a, 0) != x(b, 0)) {
      return x(a, 0) < x(b, 0);
    }
    if (x(a, 1) != x(b, 1)) {
      return x(a, 1) < x(b, 1);
    }
    return a < b;
  });
  // where the run of points sharing a first coordinate ends
  std::vector<size_t> run_end(order.size());
  for (size_t k = order.size(); k-- > 0;) {
    run_end[k] =
      ((k + 1 < order.size()) && (x(order[k + 1], 0) == x(order[k], 0)))
        ? run_end[k + 1]
        : k + 1;
  }
  const double reach = 9.0 * scale;
  Eigen::VectorXd count = Eigen::VectorXd::Zero(n);
  for (size_t k = 0; k < order.size(); ++k) {
    const Eigen::Index i = order[k];
    for (size_t l = k; l < order.size(); ++l) {
      const Eigen::Index j = order[l];
      if (x(j, 0) - x(i, 0) > reach) {
        break;
      }
      if (x(j, 1) - x(i, 1) > reach) {
        // the rest of this run is further still, in its second coordinate
        l = run_end[l] - 1;
        continue;
      }
      const double d2 = (x.row(i) - x.row(j)).squaredNorm() / (scale * scale);
      const double kernel = std::exp(-0.5 * d2);
      count(i) += kernel;
      if (l != k) {
        count(j) += kernel;
      }
    }
  }
  return count;
}

// Construct a box covering from a matrix of samples.
// @param u A matrix of samples.
// @param K The number of boxes in each dimension.
inline BoxCovering::BoxCovering(const Eigen::MatrixXd& u, uint16_t K)
  : u_(u)
  , K_(K)
{
  boxes_.resize(static_cast<size_t>(K) * K);

  n_ = u.rows();
  for (size_t i = 0; i < n_; i++) {
    boxes_[cell(u(i, 0)) * K_ + cell(u(i, 1))].insert(i);
  }
}

// The cell along one axis containing a coordinate; the unit interval is
// closed at 1, so a coordinate of exactly 1 belongs to the last cell, and
// coordinates outside [0, 1] are assigned to the nearest boundary cell.
inline size_t
BoxCovering::cell(double x) const
{
  const double clamped = std::min(std::max(x, 0.0), 1.0);
  return std::min(static_cast<size_t>(std::floor(clamped * K_)),
                  static_cast<size_t>(K_ - 1));
}

// One past the last cell along one axis that a coordinate can touch, capped
// at the number of cells.
inline size_t
BoxCovering::cell_end(double x) const
{
  const double clamped = std::min(std::max(x, 0.0), 1.0);
  return static_cast<size_t>(std::ceil(clamped * K_));
}

// Get the indices of the samples in a box defined by lower and upper bounds.
// @param lower Lower bounds of the box.
// @param upper Upper bounds of the box.
inline std::vector<size_t>
BoxCovering::get_box_indices(const Eigen::VectorXd& lower,
                             const Eigen::VectorXd& upper) const
{
  std::vector<size_t> indices;
  get_box_indices(lower, upper, indices);
  return indices;
}

// Buffer-reusing variant: clears `indices` and fills it (the capacity
// persists across calls, avoiding an allocation per call).
inline void
BoxCovering::get_box_indices(const Eigen::VectorXd& lower,
                             const Eigen::VectorXd& upper,
                             std::vector<size_t>& indices) const
{
  indices.clear();
  const size_t l0 = cell(lower(0));
  const size_t l1 = cell(lower(1));
  const size_t u0 = cell_end(upper(0));
  const size_t u1 = cell_end(upper(1));

  for (size_t k = l0; k < u0; k++) {
    for (size_t j = l1; j < u1; j++) {
      for (auto& i : boxes_[k * K_ + j]) {
        if ((k == l0) || (k == u0 - 1)) {
          if ((u_(i, 0) < lower(0)) || (u_(i, 0) > upper(0)))
            continue;
        }
        if ((j == l1) || (j == u1 - 1)) {
          if ((u_(i, 1) < lower(1)) || (u_(i, 1) > upper(1)))
            continue;
        }
        indices.push_back(i);
      }
    }
  }
}

// Swap a sample in the box covering.
// @param i Index of the sample to swap.
inline void
BoxCovering::swap_sample(size_t i, const Eigen::VectorXd& new_sample)
{
  boxes_[cell(u_(i, 0)) * K_ + cell(u_(i, 1))].erase(i);
  u_.row(i) = new_sample;
  boxes_[cell(new_sample(0)) * K_ + cell(new_sample(1))].insert(i);
}

// SplitMix64's finalizer: a bijection of 64-bit integers whose outputs are
// statistically uniform, computed in integer arithmetic so that every build
// draws the same keys.
inline uint64_t
mix64(uint64_t z)
{
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

//! @brief Recovers a continuous latent sample from a sample of a discrete
//! copula.
//!
//! Treats the discrete sample as an interval-censored density estimation
//! problem: observation \f$ i \f$ is only known to lie in the rectangle
//! \f$ (u_{i1}^-, u_{i1}] \times (u_{i2}^-, u_{i2}] \f$. A latent point is
//! initialized inside that rectangle and then refined by `niter` sweeps: for
//! each observation, one of the observations whose current latent point falls
//! in its own rectangle (widened by \f$ b \f$ on the normal scale) is picked
//! uniformly at random, and the latent point is set to that neighbor's value
//! plus Gaussian noise of scale \f$ b \f$ --- a draw from a Gaussian kernel
//! density estimate restricted to the compatible observations. The refined
//! points are not confined to the original rectangles.
//!
//! The draw is deterministic: every random component is drawn from a
//! fixed-seed generator, so repeated calls on the same input agree bit for
//! bit. It is also independent of argument order: recovering the latent sample
//! of \f$ (u_1, u_2) \f$ and of \f$ (u_2, u_1) \f$ gives the same pair of
//! columns, swapped.
//!
//! @param u An \f$ n \times 4 \f$ matrix \f$ (u_1, u_2, u_1^-, u_2^-) \f$
//! holding each observation's distribution function values and their left
//! limits; any other number of columns throws.
//! @param b The bandwidth of the kernel density estimator.
//! @param niter The number of sweeps.
//!
//! @return An \f$ n \times 2 \f$ matrix of pseudo-observations of the latent
//! sample, i.e. on the copula scale rather than the normal scale the sweeps run
//! on.
inline Eigen::MatrixXd
find_latent_sample(const Eigen::MatrixXd& u, double b, size_t niter)
{
  using namespace tools_stats;
  size_t n = u.rows();
  if (u.cols() != 4) {
    throw std::runtime_error("u must have four columns.");
  }

  // The draws below are indexed by column position, so ordering the pair by
  // its own values is what makes the result a function of the observations
  // rather than of how they were passed.
  const bool swapped = swaps_pair(u);
  Eigen::MatrixXd v = u;
  if (swapped) {
    v.col(0).swap(v.col(1));
    v.col(2).swap(v.col(3));
  }

  // Pseudo-random, not quasi-random: a low-discrepancy sequence is a
  // deterministic near-lattice whose coordinates are negatively dependent by
  // construction, which is what makes it integrate well and what makes it
  // unusable as per-observation noise -- it imposes that lattice on the latent
  // points and attenuates their dependence. Fixed seeds keep the draw
  // reproducible.
  auto w = simulate_uniform(n, 2, false, { 5 });
  // exactly the value where there is no atom, so that ties stay exact
  Eigen::MatrixXd uu =
    v.rightCols(2).array() +
    w.array() * (v.leftCols(2).array() - v.rightCols(2).array());

  auto covering = BoxCovering(uu);
  std::vector<size_t> indices;

  Eigen::MatrixXd lb = qnorm(v.rightCols(2));
  Eigen::MatrixXd ub = qnorm(v.leftCols(2));
  lb = pnorm(lb.array() - b);
  ub = pnorm(ub.array() + b);

  Eigen::MatrixXd x(n, 2), norm_sim(n, 2);

  for (size_t it = 0; it < niter; it++) {
    // continuous in the cloud, so that values within rounding of each other
    // cannot reorder a block of it
    uu = soft_pseudo_obs(
      uu, Eigen::VectorXd(), {}, default_soft_scale(), "average");
    x = qnorm(uu);
    // the seed vectors hold `int`, which a `size_t` does not narrow to
    // implicitly inside a braced initializer
    const auto seed = static_cast<int>(it);
    norm_sim = simulate_normal(n, 2, false, { seed, 5 }).array() * b;
    const uint64_t sweep = mix64(static_cast<uint64_t>(it) + 55);

    for (size_t i = 0; i < n; i++) {
      covering.get_box_indices(lb.row(i), ub.row(i), indices);
      if (!indices.empty()) {
        // The neighbor is drawn uniformly as the compatible one holding the
        // smallest of a fixed key per (sweep, observation, neighbor). A
        // neighbor entering or leaving the compatible set changes the draw only
        // if it holds that key, where an index into the list of neighbors would
        // move with any change to the list and redraw every later observation.
        const uint64_t target = mix64(sweep ^ mix64(i));
        size_t j = indices[0];
        uint64_t best = mix64(target ^ j);
        for (size_t k = 1; k < indices.size(); ++k) {
          const uint64_t key = mix64(target ^ indices[k]);
          if ((key < best) || ((key == best) && (indices[k] < j))) {
            best = key;
            j = indices[k];
          }
        }
        x.row(i) = x.row(j) + norm_sim.row(i);
        uu.row(i) = pnorm(x.row(i));
        covering.swap_sample(i, uu.row(i));
      }
    }
  }

  Eigen::MatrixXd latent =
    soft_pseudo_obs(x, Eigen::VectorXd(), {}, default_soft_scale(), "average");
  if (swapped) {
    latent.col(0).swap(latent.col(1));
  }
  return latent;
}

// Utility function to compute the next power of 2.
inline size_t
next_power_of_two(size_t n)
{
  size_t power = 1;
  while (power < n) {
    power *= 2;
  }
  return power;
}

//! reusable state for the FFT window smoother: the plan, the window
//! transform, and the scratch buffers only depend on (fftSize, wl), which is
//! constant across all `win` calls of one `ace` run
struct SmoothingWorkspace
{
  Eigen::FFT<double> fft;
  Eigen::VectorXd xx;
  Eigen::VectorXcd win_fft, tmp1, tmp2;
  size_t fft_size{ 0 };
  size_t wl{ 0 };
};

//! window smoother
inline Eigen::VectorXd
win(const Eigen::VectorXd& x, size_t wl, SmoothingWorkspace& ws)
{
  size_t n = x.size();
  // pad length to powers of 2 to force FFT to use its fastest algorithm
  size_t fftSize = next_power_of_two(n + 2 * wl);

  if ((ws.fft_size != fftSize) || (ws.wl != wl)) {
    // the conjugated window transform is constant; compute it once
    Eigen::VectorXd yy = Eigen::VectorXd::Zero(fftSize);
    yy.head(2 * wl + 1) = Eigen::VectorXd::Ones(2 * wl + 1);
    Eigen::VectorXcd tmp = ws.fft.fwd(yy);
    ws.win_fft = tmp.conjugate();
    ws.xx = Eigen::VectorXd::Zero(fftSize);
    ws.fft_size = fftSize;
    ws.wl = wl;
  }

  ws.xx.setZero();
  ws.xx.segment(2 * wl, n) = x;
  ws.fft.fwd(ws.tmp1, ws.xx);
  ws.tmp1 = ws.tmp1.cwiseProduct(ws.win_fft);
  ws.fft.inv(ws.tmp2, ws.tmp1);

  Eigen::VectorXd result = ws.tmp2.real().segment(wl, n);
  result /= 2.0 * static_cast<double>(wl) + 1.0;
  result.head(wl).setConstant(result(wl));
  result.tail(wl).setConstant(result(n - wl - 1));

  return result;
}

//! helper routine for ace (In R, this would be win(x[ind], wl)[ranks])
inline Eigen::VectorXd
cef(const Eigen::VectorXd& x,
    const Eigen::Matrix<size_t, Eigen::Dynamic, 1>& ind,
    const Eigen::Matrix<size_t, Eigen::Dynamic, 1>& ranks,
    size_t wl,
    SmoothingWorkspace& ws)
{
  Eigen::VectorXd cey = x(ind);
  cey = win(cey, wl, ws);
  return cey(ranks);
}

//! alternating conditional expectation algorithm
inline Eigen::MatrixXd
ace(const Eigen::MatrixXd& data,                        // data
    const Eigen::VectorXd& weights = Eigen::VectorXd(), // weights
    size_t wl = 0,                // window length for the smoother
    size_t outer_iter_max = 100,  // max number of outer iterations
    size_t inner_iter_max = 10,   // max number of inner iterations
    double outer_abs_tol = 2e-15, // outer stopping criterion
    double inner_abs_tol = 1e-4)  // inner stopping criterion
{
  // sample size and memory allocation for the outer/inner loops
  size_t n = data.rows();
  Eigen::VectorXd tmp(n);

  size_t nw = weights.size();
  Eigen::VectorXd w(n);
  if (nw == 0) {
    w = Eigen::VectorXd::Ones(n);
  } else {
    if (nw != n) {
      throw std::runtime_error("weights should have a length equal to "
                               "the number of rows in data");
    }
    w = weights;
  }

  // default window size
  double n_dbl = static_cast<double>(n);
  if (wl == 0) {
    wl = static_cast<size_t>(std::ceil(n_dbl / 5));
  }

  // assign order/ranks to ind/ranks
  Eigen::Matrix<size_t, Eigen::Dynamic, 2> ind(n, 2);
  Eigen::Matrix<size_t, Eigen::Dynamic, 2> ranks(n, 2);
  for (size_t i = 0; i < 2; i++) {
    std::vector<double> xvec(data.data() + n * i, data.data() + n * (i + 1));
    auto order = tools_stl::get_order(xvec);
    for (auto j : order) {
      ind(j, i) = order[j];
      ranks(order[j], i) = j;
    }
  }

  // initialize output
  Eigen::MatrixXd phi = ranks.cast<double>();
  phi.array() -= (n_dbl - 1.0) / 2.0 - 1.0;
  phi /= std::sqrt(n_dbl * (n_dbl - 1.0) / 12.0);
  if (nw > 0) {
    phi.col(0) = phi.col(0).cwiseProduct(w);
    phi.col(1) = phi.col(1).cwiseProduct(w);
  }

  // initialize variables for the outer loop
  size_t outer_iter = 1;
  double outer_eps = 1.0;
  double outer_abs_err = 1.0;
  SmoothingWorkspace ws;

  // outer loop (expectation of the first variable given the second)
  while (outer_iter <= outer_iter_max && outer_abs_err > outer_abs_tol) {
    // initialize variables for the inner loop
    size_t inner_iter = 1;
    double inner_eps = 1.0;
    double inner_abs_err = 1.0;

    // inner loop (expectation of the second variable given the first)
    while (inner_iter <= inner_iter_max && inner_abs_err > inner_abs_tol) {
      // conditional expectation
      phi.col(1) =
        cef(phi.col(0).cwiseProduct(w), ind.col(1), ranks.col(1), wl, ws);

      // center and standardize
      double m1 = phi.col(1).sum() / n_dbl;
      phi.col(1).array() -= m1;
      double s1 = std::sqrt(phi.col(1).cwiseAbs2().sum() / (n_dbl - 1));
      phi.col(1) /= s1;

      // compute error and increase step
      inner_abs_err = inner_eps;
      tmp = (phi.col(1) - phi.col(0));
      inner_eps = tmp.cwiseAbs2().sum() / n_dbl;
      inner_abs_err = std::fabs(inner_abs_err - inner_eps);
      inner_iter = inner_iter + 1;
    }

    // conditional expectation
    phi.col(0) =
      cef(phi.col(1).cwiseProduct(w), ind.col(0), ranks.col(0), wl, ws);

    // center and standardize
    double m0 = phi.col(0).sum() / n_dbl;
    phi.col(0).array() -= m0;
    double s0 = std::sqrt(phi.col(0).cwiseAbs2().sum() / (n_dbl - 1));
    phi.col(0) /= s0;

    // compute error and increase step
    outer_abs_err = outer_eps;
    tmp = (phi.col(1) - phi.col(0));
    outer_eps = tmp.cwiseAbs2().sum() / n_dbl;
    outer_abs_err = std::fabs(outer_abs_err - outer_eps);
    outer_iter = outer_iter + 1;
  }

  // return result
  return phi;
}

//! @name Dependence measures
//! @{

//! calculates the pairwise maximum correlation coefficient.
//!
//! @details Symmetric in the two variables, as the measure is by definition.
inline double
pairwise_mcor(const Eigen::MatrixXd& x, const Eigen::VectorXd& weights)
{
  // ACE updates one variable first, and where the dependence is weak the two
  // orders can stop at different correlations; ordering the pair by its own
  // values, as `find_latent_sample` does, makes the result a function of the
  // pair rather than of how it was passed
  Eigen::MatrixXd v = x.leftCols(2);
  if (swaps_pair(v)) {
    v.col(0).swap(v.col(1));
  }
  Eigen::MatrixXd phi = ace(v, weights);
  return wdm::wdm(phi, "cor", weights)(0, 1);
}

//! @brief calculates the pairwise symmetrized Chatterjee's xi.
//!
//! Chatterjee's xi measures how well one variable is a measurable function of
//! the other and is therefore asymmetric. The symmetrized version is the
//! larger of the two directions, so it detects a functional relationship
//! whichever way it runs.
//!
//! @literature
//! Chatterjee, Sourav. *A New Coefficient of Correlation*. Journal of the
//! American Statistical Association 116(536), 2009-2022, 2021
inline double
pairwise_cxi(const Eigen::MatrixXd& x, const Eigen::VectorXd& weights)
{
  double xi12 = wdm::wdm(x.col(0), x.col(1), "cxi", weights);
  double xi21 = wdm::wdm(x.col(1), x.col(0), "cxi", weights);
  if (std::isnan(xi12) || std::isnan(xi21)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::max(xi12, xi21);
}
//! @}

//! @brief Simulates from the multivariate Generalized Halton Sequence.
//!
//! For more information on Generalized Halton Sequence, see
//! Faure, H., Lemieux, C. (2009). Generalized Halton Sequences in 2008:
//! A Comparative Study. ACM-TOMACS 19(4), Article 15.
//!
//! @param n Number of observations.
//! @param d Dimension.
//! @param seeds Seeds to scramble the quasi-random numbers; if empty
//! (default),
//!   the quasi-random number generator is seeded randomly.
//!
//! @return An \f$ n \times d \f$ matrix of quasi-random
//! \f$ \mathrm{U}[0, 1] \f$ variables.
inline Eigen::MatrixXd
ghalton(const size_t& n, const size_t& d, const std::vector<int>& seeds)
{
  if ((n < 1) || (d < 1)) {
    throw std::runtime_error("n and d must be at least 1.");
  }

  Eigen::MatrixXd res(d, n);

  // Coefficients of the shift
  Eigen::MatrixXi shcoeff(d, 32);
  Eigen::VectorXi base = tools_ghalton::primes.block(0, 0, d, 1);
  Eigen::MatrixXd u = Eigen::VectorXd::Zero(d, 1);
  auto U = simulate_uniform(d, 32, false, seeds);
  for (int k = 31; k >= 0; k--) {
    shcoeff.col(k) =
      (base.cast<double>()).cwiseProduct(U.block(0, k, d, 1)).cast<int>();
    u = (u + shcoeff.col(k).cast<double>()).cwiseQuotient(base.cast<double>());
  }
  res.block(0, 0, d, 1) = u;

  Eigen::VectorXi perm = tools_ghalton::permTN2.block(0, 0, d, 1);
  Eigen::MatrixXi coeff(d, 32);
  Eigen::VectorXi tmp(d);
  auto mod = [](const int& u1, const int& u2) { return u1 % u2; };
  for (size_t i = 1; i < n; i++) {

    // Find i in the prime base
    tmp = Eigen::VectorXi::Constant(d, static_cast<int>(i));
    coeff = Eigen::MatrixXi::Zero(d, 32);
    int k = 0;
    while ((tmp.maxCoeff() > 0) && (k < 32)) {
      coeff.col(k) = tmp.binaryExpr(base, mod);
      tmp = tmp.cwiseQuotient(base);
      k++;
    }

    u = Eigen::VectorXd::Zero(d);
    k = 31;
    while (k >= 0) {
      tmp = perm.cwiseProduct(coeff.col(k)) + shcoeff.col(k);
      u = u + tmp.binaryExpr(base, mod).cast<double>();
      u = u.cwiseQuotient(base.cast<double>());
      k--;
    }
    res.block(0, i, d, 1) = u;
  }

  return res.transpose();
}

//! @brief Simulates from the multivariate Sobol sequence.
//!
//! For more information on the Sobol sequence, see S. Joe and F. Y. Kuo
//! (2008), constructing Sobol  sequences with better two-dimensional
//! projections, SIAM J. Sci. Comput. 30, 2635–2654.
//!
//! @param n Number of observations.
//! @param d Dimension.
//! @param seeds Seeds to scramble the quasi-random numbers; if empty
//! (default),
//!   the quasi-random number generator is seeded randomly.
//!
//! @return An \f$ n \times d \f$ matrix of quasi-random
//! \f$ \mathrm{U}[0, 1] \f$ variables.
inline Eigen::MatrixXd
sobol(const size_t& n, const size_t& d, const std::vector<int>& seeds)
{
  if ((n < 1) || (d < 1)) {
    throw std::runtime_error("n and d must be at least 1.");
  }

  // output matrix
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(n, d);

  // L = max number of bits needed
  size_t L =
    static_cast<size_t>(std::ceil(log(static_cast<double>(n)) / log(2.0)));

  // Vector of scrambling factors
  Eigen::MatrixXd scrambling = simulate_uniform(d, 1, false, seeds);

  // C(i) = index from the right of the first zero bit of i + 1
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> C(n);
  C(0) = 1;
  for (size_t i = 1; i < n; i++) {
    C(i) = 1;
    size_t value = i;
    while (value & 1) {
      value >>= 1;
      C(i)++;
    }
  }

  // Compute the first dimension

  // Compute direction numbers scaled by pow(2,32)
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> V(L);
  for (size_t i = 0; i < L; i++) {
    V(i) = static_cast<size_t>(1) << (32 - (i + 1)); // all m's = 1
  }

  // Evaluate X scaled by pow(2,32)
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> X(n);
  X(0) = static_cast<size_t>(scrambling(0) * 4294967296.0);
  for (size_t i = 1; i < n; i++) {
    X(i) = X(i - 1) ^ V(C(i - 1) - 1);
  }
  output.block(0, 0, n, 1) = X.cast<double>();

  // Compute the remaining dimensions
  for (size_t j = 0; j < d - 1; j++) {

    // Get parameters from static vectors
    size_t a = tools_sobol::a_sobol[j];
    size_t s = tools_sobol::s_sobol[j];

    Eigen::Map<Eigen::Matrix<size_t, Eigen::Dynamic, 1>> m(
      tools_sobol::minit_sobol[j], s);

    // Compute direction numbers scaled by pow(2,32)
    for (size_t i = 0; i < std::min(L, s); i++)
      V(i) = m(i) << (32 - (i + 1));

    if (L > s) {
      for (size_t i = s; i < L; i++) {
        V(i) = V(i - s) ^ (V(i - s) >> s);
        for (size_t k = 0; k < s - 1; k++)
          V(i) ^= (((a >> (s - 2 - k)) & 1) * V(i - k - 1));
      }
    }

    // Evaluate X
    X(0) = static_cast<size_t>(scrambling(j + 1) * 4294967296.0);
    for (size_t i = 1; i < n; i++)
      X(i) = X(i - 1) ^ V(C(i - 1) - 1);
    output.block(0, j + 1, n, 1) = X.cast<double>();
  }

  // Scale output by pow(2,32)
  output /= 4294967296.0;

  return output;
}

//! @brief Computes bivariate t probabilities.
//!
//! Based on the method described by
//! Dunnett, C.W. and M. Sobel, (1954),
//! A bivariate generalization of Student's t-distribution
//! with tables for certain special cases,
//! Biometrika 41, pp. 153-169. Translated from the Fortran routines of
//! Alan Genz (www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f).
//!
//! @param z An \f$ n \times 2 \f$ matrix of evaluation points.
//! @param nu Number of degrees of freedom.
//! @param rho Correlation.
//!
//! @return An \f$ n \times 1 \f$ vector of probabilities.
inline Eigen::VectorXd
pbvt(const Eigen::MatrixXd& z, int nu, double rho)
{
  double snu = sqrt(static_cast<double>(nu));
  double ors = 1 - pow(rho, 2.0);
  // even-nu starting value; depends only on rho, so hoisted out of the
  // per-element lambda
  double bvt0 = atan2(sqrt(ors), -rho) / 6.2831853071795862;

  auto f = [snu, nu, ors, rho, bvt0](double h, double k) {
    double d1, d2, bvt, gmph, gmpk, xnkh, xnhk, btnckh, btnchk, btpdkh, btpdhk;
    int hs, ks;

    double hrk = h - rho * k;
    double krh = k - rho * h;
    if (std::fabs(hrk) + ors > 0.) {
      /* Computing 2nd power */
      d1 = hrk;
      /* Computing 2nd power */
      d2 = hrk;
      /* Computing 2nd power */
      double d3 = k;
      xnhk = d1 * d1 / (d2 * d2 + ors * (nu + d3 * d3));
      /* Computing 2nd power */
      d1 = krh;
      /* Computing 2nd power */
      d2 = krh;
      /* Computing 2nd power */
      d3 = h;
      xnkh = d1 * d1 / (d2 * d2 + ors * (nu + d3 * d3));
    } else {
      xnhk = 0.;
      xnkh = 0.;
    }
    d1 = h - rho * k;
    hs = static_cast<int>(d1 >= 0 ? 1 : -1);
    d1 = k - rho * h;
    ks = static_cast<int>(d1 >= 0 ? 1 : -1);
    if (nu % 2 == 0) {
      bvt = bvt0;
      /* Computing 2nd power */
      d1 = h;
      gmph = h / sqrt((nu + d1 * d1) * 16);
      /* Computing 2nd power */
      d1 = k;
      gmpk = k / sqrt((nu + d1 * d1) * 16);
      btnckh = atan2(sqrt(xnkh), sqrt(1 - xnkh)) * 2 / 3.14159265358979323844;
      btpdkh = sqrt(xnkh * (1 - xnkh)) * 2 / 3.14159265358979323844;
      btnchk = atan2(sqrt(xnhk), sqrt(1 - xnhk)) * 2 / 3.14159265358979323844;
      btpdhk = sqrt(xnhk * (1 - xnhk)) * 2 / 3.14159265358979323844;
      size_t i1 = static_cast<size_t>(nu / 2);
      for (size_t j = 1; j <= i1; ++j) {
        double jj = static_cast<double>(j << 1);
        bvt += gmph * (ks * btnckh + 1);
        bvt += gmpk * (hs * btnchk + 1);
        btnckh += btpdkh;
        btpdkh = jj * btpdkh * (1 - xnkh) / (jj + 1);
        btnchk += btpdhk;
        btpdhk = jj * btpdhk * (1 - xnhk) / (jj + 1);
        /* Computing 2nd power */
        d1 = h;
        gmph = gmph * (jj - 1) / (jj * (d1 * d1 / nu + 1));
        /* Computing 2nd power */
        d1 = k;
        gmpk = gmpk * (jj - 1) / (jj * (d1 * d1 / nu + 1));
      }
    } else {
      /* Computing 2nd power */
      d1 = h;
      /* Computing 2nd power */
      d2 = k;
      double qhrk = sqrt(d1 * d1 + d2 * d2 - rho * 2 * h * k + nu * ors);
      double hkrn = h * k + rho * nu;
      double hkn = h * k - nu;
      double hpk = h + k;
      bvt =
        atan2(-snu * (hkn * qhrk + hpk * hkrn), hkn * hkrn - nu * hpk * qhrk) /
        6.2831853071795862;
      if (bvt < -1e-15) {
        bvt += 1;
      }
      /* Computing 2nd power */
      d1 = h;
      gmph = h / (snu * 6.2831853071795862 * (d1 * d1 / nu + 1));
      /* Computing 2nd power */
      d1 = k;
      gmpk = k / (snu * 6.2831853071795862 * (d1 * d1 / nu + 1));
      btnckh = sqrt(xnkh);
      btpdkh = btnckh;
      btnchk = sqrt(xnhk);
      btpdhk = btnchk;
      size_t i1 = static_cast<size_t>((nu - 1) / 2);
      for (size_t j = 1; j <= i1; ++j) {
        double jj = static_cast<double>(j << 1);
        bvt += gmph * (ks * btnckh + 1);
        bvt += gmpk * (hs * btnchk + 1);
        btpdkh = (jj - 1) * btpdkh * (1 - xnkh) / jj;
        btnckh += btpdkh;
        btpdhk = (jj - 1) * btpdhk * (1 - xnhk) / jj;
        btnchk += btpdhk;
        /* Computing 2nd power */
        d1 = h;
        gmph = jj * gmph / ((jj + 1) * (d1 * d1 / nu + 1));
        /* Computing 2nd power */
        d1 = k;
        gmpk = jj * gmpk / ((jj + 1) * (d1 * d1 / nu + 1));
      }
    }
    return bvt;
  };

  return tools_eigen::binaryExpr_or_nan(z, f);
}

//! @brief Compute bivariate normal probabilities.
//!
//! A function for computing bivariate normal probabilities;
//! developed using Drezner, Z. and Wesolowsky, G. O. (1989),
//! On the Computation of the Bivariate Normal Integral,
//! J. Stat. Comput. Simul.. 35 pp. 101-107.
//! with extensive modifications for double precisions by
//! Alan Genz and Yihong Ge. Translated from the Fortran routines of
//! Alan Genz (www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f).
//!
//! @param z An \f$ n \times 2 \f$ matrix of evaluation points.
//! @param rho Correlation.
//!
//! @return An \f$ n \times 1 \f$ vector of probabilities.
inline Eigen::VectorXd
pbvnorm(const Eigen::MatrixXd& z, double rho)
{

  // normal cdf; direct erfc form of the standard normal cdf (avoids the
  // boost distribution-object dispatch on every evaluation)
  static const double inv_sqrt2 = 0.70710678118654752440;
  auto phi = [](double y) { return 0.5 * std::erfc(-y * inv_sqrt2); };

  // set-up required constants
  size_t lg;
  if (std::fabs(rho) < .3f) {
    lg = 3;
  } else if (std::fabs(rho) < .75f) {
    lg = 6;
  } else {
    lg = 10;
  }
  Eigen::VectorXd w(lg), x(lg);
  if (std::fabs(rho) < .3f) {
    w << 0.1713244923791705, 0.3607615730481384, 0.4679139345726904;
    x << -.9324695142031522, -.6612093864662647, -.238619186083197;
  } else if (std::fabs(rho) < .75f) {
    w << 0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
      0.2031674267230659, 0.2334925365383547, 0.2491470458134029;
    x << -.9815606342467191, -.904117256370475, -.769902674194305,
      -.5873179542866171, -.3678314989981802, -.1252334085114692;
  } else {
    w << 0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
      0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
      0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
      0.1527533871307259;
    x << -.9931285991850949, -.9639719272779138, -.9122344282513259,
      -.8391169718222188, -.7463319064601508, -.636053680726515,
      -.5108670019508271, -.3737060887154196, -.2277858511416451,
      -.07652652113349733;
  }

  // everything that depends only on rho and the quadrature nodes is
  // precomputed here instead of once per evaluation point; fixed-size stack
  // tables and branch-gated setup keep single-point calls (the discrete
  // difference quotients evaluate the cdf row by row) as cheap as the
  // pre-hoisting code
  const double asr = asin(rho);
  std::array<double, 10> sn1{}, sn2{}, dn1{}, dn2{};
  std::array<double, 10> xs1{}, rs1{}, xs2{}, rs2{};
  const double as = (1 - rho) * (rho + 1);
  const double a_full = std::sqrt(as);
  const double a_half = a_full / 2;
  if (std::fabs(rho) < .925f) {
    for (size_t i = 0; i < lg; ++i) {
      double sn = std::sin(asr * (x(i) + 1) / 2);
      sn1[i] = sn;
      dn1[i] = 1 - sn * sn;
      sn = std::sin(asr * (-x(i) + 1) / 2);
      sn2[i] = sn;
      dn2[i] = 1 - sn * sn;
    }
  } else {
    for (size_t i = 0; i < lg; ++i) {
      double d1 = a_half * (x(i) + 1);
      xs1[i] = d1 * d1;
      rs1[i] = std::sqrt(1 - xs1[i]);
      d1 = -x(i) + 1;
      xs2[i] = as * (d1 * d1) / 4;
      rs2[i] = std::sqrt(1 - xs2[i]);
    }
  }

  auto f = [lg,
            rho,
            w,
            phi,
            asr,
            &sn1,
            &sn2,
            &dn1,
            &dn2,
            as,
            a_full,
            a_half,
            &xs1,
            &rs1,
            &xs2,
            &rs2](double h, double k) {
    size_t i1;
    double d1, d2, hk, bvn;
    h = -h;
    k = -k;
    hk = h * k;
    bvn = 0.0;
    if (std::fabs(rho) < .925f) {
      double hs = (h * h + k * k) / 2;
      i1 = lg;
      for (size_t i = 0; i < i1; ++i) {
        bvn += w(i) * std::exp((sn1[i] * hk - hs) / dn1[i]);
        bvn += w(i) * std::exp((sn2[i] * hk - hs) / dn2[i]);
      }
      d1 = -h;
      d2 = -k;
      bvn = bvn * asr / 12.566370614359172 + phi(d1) * phi(d2);
    } else {
      if (rho < 0.) {
        k = -k;
        hk = -hk;
      }
      if (std::fabs(rho) < 1.) {
        /* Computing 2nd power */
        d1 = h - k;
        double bs = d1 * d1;
        double c = (4 - hk) / 8;
        double d = (12 - hk) / 16;
        bvn = a_full * std::exp(-(bs / as + hk) / 2) *
              (1 - c * (bs - as) * (1 - d * bs / 5) / 3 + c * d * as * as / 5);
        if (hk > -160.) {
          double b = std::sqrt(bs);
          d1 = -b / a_full;
          bvn -= std::exp(-hk / 2) * std::sqrt(6.283185307179586) * phi(d1) *
                 b * (1 - c * bs * (1 - d * bs / 5) / 3);
        }
        i1 = lg;
        for (size_t i = 0; i < i1; ++i) {
          bvn += a_half * w(i) *
                 (std::exp(-bs / (xs1[i] * 2) - hk / (rs1[i] + 1)) / rs1[i] -
                  std::exp(-(bs / xs1[i] + hk) / 2) *
                    (c * xs1[i] * (d * xs1[i] + 1) + 1));
          /* Computing 2nd power */
          d1 = rs2[i] + 1;
          bvn += a_half * w(i) * std::exp(-(bs / xs2[i] + hk) / 2) *
                 (std::exp(-hk * xs2[i] / (d1 * d1 * 2)) / rs2[i] -
                  (c * xs2[i] * (d * xs2[i] + 1) + 1));
        }
        bvn = -bvn / 6.283185307179586;
      }
      if (rho > 0.) {
        d1 = -std::max(h, k);
        bvn += phi(d1);
      } else {
        bvn = -bvn;
        if (k > h) {
          if (h < 0.) {
            bvn = bvn + phi(k) - phi(h);
          } else {
            d1 = -h;
            d2 = -k;
            bvn = bvn + phi(d1) - phi(d2);
          }
        }
      }
    }
    return bvn;
  };

  return tools_eigen::binaryExpr_or_nan(z, f);
}

}
}
