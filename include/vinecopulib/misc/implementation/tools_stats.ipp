// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/seed_seq.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <memory>
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

  Eigen::MatrixXd u(n, d);
  return u.unaryExpr([&](double) { return distribution(generator); });
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
//! @return Pseudo-observations of the copula, i.e. \f$ F_X(x) \f$
//! (column-wise).
inline Eigen::MatrixXd
to_pseudo_obs(Eigen::MatrixXd x,
              const std::string& ties_method,
              const Eigen::VectorXd& weights,
              std::vector<int> seeds)
{
  for (int j = 0; j < x.cols(); ++j)
    x.col(j) = to_pseudo_obs_1d(
      static_cast<Eigen::VectorXd>(x.col(j)), ties_method, weights, seeds);

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
//! @return Pseudo-observations of the copula, i.e. \f$ F_X(x) \f$.
inline Eigen::VectorXd
to_pseudo_obs_1d(Eigen::VectorXd x,
                 const std::string& ties_method,
                 const Eigen::VectorXd& weights,
                 std::vector<int> seeds)
{
  size_t n = x.size();
  auto xvec = wdm::utils::convert_vec(x);
  auto res =
    wdm::impl::rank(xvec, wdm::utils::convert_vec(weights), ties_method, seeds);
  x = Eigen::Map<Eigen::VectorXd>(res.data(), res.size());

  // correction for NaNs
  if (wdm::utils::any_nan(xvec)) {
    for (size_t i = 0; i < xvec.size(); i++) {
      if (std::isnan(xvec[i])) {
        n--;
      }
    }
  }

  return x.array() / (static_cast<double>(n) + 1.0);
}

// Construct a box covering from a matrix of samples.
// @param u A matrix of samples.
// @param K The number of boxes in each dimension.
inline BoxCovering::BoxCovering(const Eigen::MatrixXd& u, uint16_t K)
  : u_(u)
  , K_(K)
{
  boxes_.resize(K);
  for (size_t k = 0; k < K; k++) {
    boxes_[k].resize(K);
    for (size_t j = 0; j < K; j++) {
      boxes_[k][j] = std::make_unique<Box>(
        std::vector<double>{ static_cast<double>(k) / K,
                             static_cast<double>(j) / K },
        std::vector<double>{ static_cast<double>(k + 1) / K,
                             static_cast<double>(j + 1) / K });
    }
  }

  n_ = u.rows();
  for (size_t i = 0; i < n_; i++) {
    size_t k = static_cast<size_t>(std::floor(u(i, 0) * K));
    size_t j = static_cast<size_t>(std::floor(u(i, 1) * K));
    boxes_[k][j]->indices_.insert(i);
  }
}

// Get the indices of the samples in a box defined by lower and upper bounds.
// @param lower Lower bounds of the box.
// @param upper Upper bounds of the box.
inline std::vector<size_t>
BoxCovering::get_box_indices(const Eigen::VectorXd& lower,
                             const Eigen::VectorXd& upper) const
{
  std::vector<size_t> indices;
  indices.reserve(n_);
  auto l0 = static_cast<size_t>(std::floor(lower(0) * K_));
  auto l1 = static_cast<size_t>(std::floor(lower(1) * K_));
  auto u0 = static_cast<size_t>(std::ceil(upper(0) * K_));
  auto u1 = static_cast<size_t>(std::ceil(upper(1) * K_));

  for (size_t k = l0; k < u0; k++) {
    for (size_t j = l1; j < u1; j++) {
      for (auto& i : boxes_[k][j]->indices_) {
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

  return indices;
}

// Swap a sample in the box covering.
// @param i Index of the sample to swap.
inline void
BoxCovering::swap_sample(size_t i, const Eigen::VectorXd& new_sample)
{
  auto k = static_cast<size_t>(std::floor(u_(i, 0) * K_));
  auto j = static_cast<size_t>(std::floor(u_(i, 1) * K_));
  boxes_[k][j]->indices_.erase(i);

  u_.row(i) = new_sample;
  k = static_cast<size_t>(std::floor(new_sample(0) * K_));
  j = static_cast<size_t>(std::floor(new_sample(1) * K_));
  boxes_[k][j]->indices_.insert(i);
}

//! Create a single box.
inline BoxCovering::Box::Box(const std::vector<double>& lower,
                             const std::vector<double>& upper)
  : lower_(lower)
  , upper_(upper)
{
}

// Recovers a (continuous) latent sample from a sample of a discrete copula by
// treating it as an interval-censored density estimation problem.
// @param u A matrix of samples.
// @param b The bandwidth of the kernel density estimator.
// @param niter The number of iterations.
inline Eigen::MatrixXd
find_latent_sample(const Eigen::MatrixXd& u, double b, size_t niter)
{
  using namespace tools_stats;
  size_t n = u.rows();
  if (u.cols() != 4) {
    throw std::runtime_error("u must have four columns.");
  }

  auto w = simulate_uniform(n, 2, true, { 5 });
  Eigen::MatrixXd uu = w.array() * u.leftCols(2).array() +
                       (1 - w.array()) * u.rightCols(2).array();

  auto covering = BoxCovering(uu);
  std::vector<size_t> indices;

  Eigen::MatrixXd lb = qnorm(u.rightCols(2));
  Eigen::MatrixXd ub = qnorm(u.leftCols(2));
  lb = pnorm(lb.array() - b);
  ub = pnorm(ub.array() + b);

  Eigen::MatrixXd x(n, 2), norm_sim(n, 2);

  for (uint16_t it = 0; it < niter; it++) {
    uu = to_pseudo_obs(uu);
    x = qnorm(uu);
    norm_sim = simulate_normal(n, 2, true, { it, 5 }).array() * b;
    w = simulate_uniform(n, 1, true, { it, 55 });

    for (size_t i = 0; i < n; i++) {
      indices = covering.get_box_indices(lb.row(i), ub.row(i));
      double n_idx = static_cast<double>(indices.size());
      if (n_idx > 0) {
        size_t j = indices.at(static_cast<size_t>(w(i) * n_idx));
        x.row(i) = x.row(j) + norm_sim.row(i);
        uu.row(i) = pnorm(x.row(i));
        covering.swap_sample(i, uu.row(i));
      }
    }
  }

  return to_pseudo_obs(x);
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

//! window smoother
inline Eigen::VectorXd
win(const Eigen::VectorXd& x, size_t wl = 5)
{
  size_t n = x.size();
  // pad length to powers of 2 to force FFT to use its fastest algorithm
  size_t fftSize = next_power_of_two(n + 2 * wl);

  Eigen::VectorXd xx = Eigen::VectorXd::Zero(fftSize);
  Eigen::VectorXd yy = Eigen::VectorXd::Zero(fftSize);
  xx.segment(2 * wl, n) = x;
  yy.head(2 * wl + 1) = Eigen::VectorXd::Ones(2 * wl + 1);

  Eigen::FFT<double> fft;
  Eigen::VectorXcd tmp1 = fft.fwd(xx);
  Eigen::VectorXcd tmp2 = fft.fwd(yy);
  tmp2 = tmp2.conjugate();
  tmp1 = tmp1.cwiseProduct(tmp2);
  tmp2 = fft.inv(tmp1);

  Eigen::VectorXd result = tmp2.real().segment(wl, n);
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
    size_t wl = 5)
{
  Eigen::VectorXd cey = x(ind);
  cey = win(cey, wl);
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
        cef(phi.col(0).cwiseProduct(w), ind.col(1), ranks.col(1), wl);

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
    phi.col(0) = cef(phi.col(1).cwiseProduct(w), ind.col(0), ranks.col(0), wl);

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

//! calculates the pairwise maximum correlation coefficient.
inline double
pairwise_mcor(const Eigen::MatrixXd& x, const Eigen::VectorXd& weights)
{
  Eigen::MatrixXd phi = ace(x, weights);
  return wdm::wdm(phi, "cor", weights)(0, 1);
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
    V(i) = static_cast<size_t>(std::pow(2, 32 - (i + 1))); // all m's = 1
  }

  // Evalulate X scaled by pow(2,32)
  Eigen::Matrix<size_t, Eigen::Dynamic, 1> X(n);
  X(0) = static_cast<size_t>(scrambling(0) * std::pow(2.0, 32));
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

    // Evalulate X
    X(0) = static_cast<size_t>(scrambling(j + 1) * std::pow(2.0, 32));
    for (size_t i = 1; i < n; i++)
      X(i) = X(i - 1) ^ V(C(i - 1) - 1);
    output.block(0, j + 1, n, 1) = X.cast<double>();
  }

  // Scale output by pow(2,32)
  output /= std::pow(2.0, 32);

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

  auto f = [snu, nu, ors, rho](double h, double k) {
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
    hs = static_cast<int>(d1 >= 0 ? 1 : -1); // d_sign(&c_b91, &d1);
    d1 = k - rho * h;
    ks = static_cast<int>(d1 >= 0 ? 1 : -1);
    if (nu % 2 == 0) {
      bvt = atan2(sqrt(ors), -rho) / 6.2831853071795862;
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
//! with extensive modications for double precisions by
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

  // normal cdf
  boost::math::normal dist;
  auto phi = [&dist](double y) { return boost::math::cdf(dist, y); };

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

  auto f = [lg, rho, x, w, phi](double h, double k) {
    size_t i1;
    double d1, d2, hk, bvn;
    h = -h;
    k = -k;
    hk = h * k;
    bvn = 0.0;
    if (std::fabs(rho) < .925f) {
      double hs = (h * h + k * k) / 2;
      double asr = asin(rho);
      i1 = lg;
      for (size_t i = 0; i < i1; ++i) {
        double sn = std::sin(asr * (x(i) + 1) / 2);
        bvn += w(i) * std::exp((sn * hk - hs) / (1 - sn * sn));
        sn = std::sin(asr * (-x(i) + 1) / 2);
        bvn += w(i) * std::exp((sn * hk - hs) / (1 - sn * sn));
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
        double as = (1 - rho) * (rho + 1);
        double a = std::sqrt(as);
        /* Computing 2nd power */
        d1 = h - k;
        double bs = d1 * d1;
        double c = (4 - hk) / 8;
        double d = (12 - hk) / 16;
        bvn = a * std::exp(-(bs / as + hk) / 2) *
              (1 - c * (bs - as) * (1 - d * bs / 5) / 3 + c * d * as * as / 5);
        if (hk > -160.) {
          double b = std::sqrt(bs);
          d1 = -b / a;
          bvn -= std::exp(-hk / 2) * std::sqrt(6.283185307179586) * phi(d1) *
                 b * (1 - c * bs * (1 - d * bs / 5) / 3);
        }
        a /= 2;
        i1 = lg;
        for (size_t i = 0; i < i1; ++i) {
          /* Computing 2nd power */
          d1 = a * (x(i) + 1);
          double xs = d1 * d1;
          double rs = std::sqrt(1 - xs);
          bvn += a * w(i) *
                 (std::exp(-bs / (xs * 2) - hk / (rs + 1)) / rs -
                  std::exp(-(bs / xs + hk) / 2) * (c * xs * (d * xs + 1) + 1));
          /* Computing 2nd power */
          d1 = -x(i) + 1;
          xs = as * (d1 * d1) / 4;
          rs = std::sqrt(1 - xs);
          /* Computing 2nd power */
          d1 = rs + 1;
          bvn += a * w(i) * std::exp(-(bs / xs + hk) / 2) *
                 (std::exp(-hk * xs / (d1 * d1 * 2)) / rs -
                  (c * xs * (d * xs + 1) + 1));
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
