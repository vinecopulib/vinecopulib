// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <numeric>
#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {
inline TllBicop::TllBicop()
{
  family_ = BicopFamily::tll;
}

inline Eigen::VectorXd
TllBicop::gaussian_kernel_2d(const Eigen::MatrixXd& x)
{
  return tools_stats::dnorm(x).rowwise().prod();
}

//! selects the bandwidth matrix for local líkelihood estimator (covariance
//! times appropriate factor).
inline Eigen::Matrix2d
TllBicop::select_bandwidth(const Eigen::MatrixXd& x,
                           const std::string& method,
                           const Eigen::VectorXd& weights)
{
  size_t n = x.rows();
  double cor = wdm::wdm(x, "cor", weights)(0, 1);
  cor = std::min(std::max(cor, -0.95), 0.95);
  Eigen::Matrix2d cov = Eigen::MatrixXd::Identity(2, 2);
  cov(0, 1) = cor;
  cov(1, 0) = cor;

  double mult;
  if (method == "constant") {
    mult = std::pow(n, -1.0 / 3.0);
  } else {
    double degree;
    if (method == "linear") {
      degree = 1.0;
    } else {
      degree = 2.0;
    }
    mult = 1.5 * std::pow(n, -1.0 / (2.0 * degree + 1.0));
  }
  double mcor = tools_stats::pairwise_mcor(x, weights);
  double scale = std::pow(std::fabs(cor / mcor), 0.5 * mcor);

  return mult * cov * scale;
}

//! calculates the cholesky root of a 2x2 matrix.
inline Eigen::Matrix2d
chol22(const Eigen::Matrix2d& B)
{

  Eigen::Matrix2d rB;

  rB(0, 0) = std::sqrt(B(0, 0));
  rB(0, 1) = 0.0;
  rB(1, 0) = B(1, 0) / rB(0, 0);
  rB(1, 1) = std::sqrt(B(1, 1) - rB(1, 0) * rB(1, 0));

  return rB;
}

//! evaluates local likelihood density estimate.
//!
//! @details In the coordinates whitened by the bandwidth, let `f0` be the
//! kernel density estimate and `mu` and `S` the kernel-weighted mean and
//! covariance of the observations around an evaluation point. The estimate is
//! `f0` (constant), `f0 * exp(-mu' mu / 2)` (linear) or
//! `f0 * exp(-mu' S^{-1} mu / 2) / sqrt(det S)` (quadratic); see
//! `docs/tll/tll_and_interpolation.tex`.
//!
//! @param x Evaluation points.
//! @param x_data Observations.
//! @param B Bandwidth matrix.
//! @param method Order of local polynomial approximation; either `"constant"`,
//!   `"linear"`, or `"quadratic"`.
//! @param weights Vector of weights for the observations
//! @return a two-column matrix; first column is estimated density, second
//!    column is influence of evaluation point.
inline Eigen::MatrixXd
TllBicop::fit_local_likelihood(const Eigen::MatrixXd& x,
                               const Eigen::MatrixXd& x_data,
                               const Eigen::Matrix2d& B,
                               const std::string& method,
                               const Eigen::VectorXd& weights)
{
  size_t m = x.rows();      // number of evaluation points
  size_t n = x_data.rows(); // number of observations

  // pre-calculate inverse root of bandwidth matrix and determinant
  Eigen::Matrix2d irB = chol22(B).inverse();
  double det_irB = irB.determinant();

  // de-correlate data by applying B^{-1/2}; the bandwidth is then the identity
  Eigen::MatrixXd z = (irB * x.transpose()).transpose();
  Eigen::MatrixXd z_data = (irB * x_data.transpose()).transpose();

  Eigen::MatrixXd res(m, 2);
  Eigen::VectorXd kernels(n);
  Eigen::Vector2d mu = Eigen::Vector2d::Zero();
  Eigen::Matrix2d S = Eigen::Matrix2d::Identity();
  Eigen::MatrixXd zz(n, 2);
  for (size_t k = 0; k < m; ++k) {
    zz = z_data.rowwise() - z.row(k);
    kernels = gaussian_kernel_2d(zz) * det_irB;
    if (weights.size() > 0)
      kernels = kernels.cwiseProduct(weights);
    double f0 = kernels.mean();
    res(k, 0) = f0;
    if (method != "constant") {
      mu = zz.transpose() * kernels / kernels.sum();
      if (method == "quadratic") {
        zz.rowwise() -= mu.transpose();
        S = zz.transpose() * (zz.array().colwise() * kernels.array()).matrix() /
            kernels.sum();
      }
      double factor =
        std::exp(-0.5 * double(mu.transpose() * S.inverse() * mu)) /
        std::sqrt(S.determinant());
      if (!(std::isfinite)(factor)) {
        // the local covariance can be singular or not positive definite when
        // the true value is equal or close to zero
        factor = 0.0;
      }
      res(k, 0) *= factor;
    }
    if (weights.size() > 0) {
      // average weight in neighborhood of evaluation point (essentially a
      // kernel regression estimate);
      // kernels have already been multiplied with weights above
      double w = kernels.sum() / kernels.cwiseQuotient(weights).sum();
      res(k, 1) = calculate_infl(n, f0, mu, S, det_irB, method, w);
    } else {
      res(k, 1) = calculate_infl(n, f0, mu, S, det_irB, method, 1.0);
    }
  }

  if (weights.size() > 0) {
    // estimate can be negative if negative weights are used
    res.col(0) = res.col(0).array().max(0.0);
  }

  return res;
}

//! calculate influence for data point for density estimate based on
//! quantities pre-computed in `fit_local_likelihood()`.
//!
//! @details The influence is `W(0) / n * (1 + m' C^{-1} m) / f0`, with `W(0)`
//! the kernel at zero and `m` and `C` the mean and covariance of the
//! non-constant terms of the local polynomial under `N(mu, S)`, in the
//! whitened coordinates.
inline double
TllBicop::calculate_infl(const size_t& n,
                         const double& f0,
                         const Eigen::Vector2d& mu,
                         const Eigen::Matrix2d& S,
                         const double& det_irB,
                         const std::string& method,
                         const double& weight)
{
  // the Gaussian kernel at zero; static so it is computed only once
  static const double kernel0 =
    gaussian_kernel_2d(Eigen::MatrixXd::Zero(1, 2))(0);

  if (method == "constant") {
    // the 1x1 "matrix inverse" is just the reciprocal
    return kernel0 * det_irB / f0 * weight / static_cast<double>(n);
  }

  double q;
  if (method == "linear") {
    // the terms are v ~ N(mu, I)
    q = mu.squaredNorm();
  } else {
    // the terms are (v1, v2, v1^2 / 2, v1 v2, v2^2 / 2) with v ~ N(mu, S);
    // their moments follow from Isserlis' theorem
    const Eigen::Index ij[3][2] = { { 0, 0 }, { 0, 1 }, { 1, 1 } };
    const double coef[3] = { 0.5, 1.0, 0.5 };
    Eigen::Matrix<double, 5, 1> mean;
    Eigen::Matrix<double, 5, 5> cov;
    mean.head(2) = mu;
    cov.topLeftCorner(2, 2) = S;
    for (Eigen::Index a = 0; a < 3; ++a) {
      const Eigen::Index i = ij[a][0], j = ij[a][1];
      mean(2 + a) = coef[a] * (S(i, j) + mu(i) * mu(j));
      for (Eigen::Index r = 0; r < 2; ++r) {
        cov(r, 2 + a) = coef[a] * (S(r, i) * mu(j) + S(r, j) * mu(i));
        cov(2 + a, r) = cov(r, 2 + a);
      }
      for (Eigen::Index c = 0; c < 3; ++c) {
        const Eigen::Index k = ij[c][0], l = ij[c][1];
        cov(2 + a, 2 + c) = coef[a] * coef[c] *
                            (S(i, k) * S(j, l) + S(i, l) * S(j, k) +
                             mu(i) * mu(k) * S(j, l) + mu(i) * mu(l) * S(j, k) +
                             mu(j) * mu(k) * S(i, l) + mu(j) * mu(l) * S(i, k));
      }
    }
    q = mean.dot(cov.ldlt().solve(mean));
  }

  return kernel0 * det_irB / f0 * (1.0 + q) * weight / static_cast<double>(n);
}

//! The number of points close to each point, counted continuously:
//! \f$ m_i = \sum_j e^{-\|x_i - x_j\|^2 / (2 s^2)} \f$.
//!
//! @param x Points, one per row.
//! @param scale The distance \f$ s \f$.
//! @return The counts, each at least one.
inline Eigen::VectorXd
TllBicop::multiplicity(const Eigen::MatrixXd& x, double scale)
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
  // the distinct points in that order, and how often each occurs
  std::vector<Eigen::Index> point;
  std::vector<double> copies;
  std::vector<size_t> distinct(static_cast<size_t>(n));
  for (const Eigen::Index i : order) {
    if (point.empty() || (x(i, 0) != x(point.back(), 0)) ||
        (x(i, 1) != x(point.back(), 1))) {
      point.push_back(i);
      copies.push_back(0.0);
    }
    copies.back() += 1.0;
    distinct[static_cast<size_t>(i)] = point.size() - 1;
  }
  // where the run of distinct points sharing a first coordinate ends
  std::vector<size_t> run_end(point.size());
  for (size_t a = point.size(); a-- > 0;) {
    run_end[a] =
      ((a + 1 < point.size()) && (x(point[a + 1], 0) == x(point[a], 0)))
        ? run_end[a + 1]
        : a + 1;
  }
  const double reach = 9.0 * scale;
  std::vector<double> count(copies);
  for (size_t a = 0; a < point.size(); ++a) {
    const Eigen::Index i = point[a];
    for (size_t b = a + 1; b < point.size(); ++b) {
      const Eigen::Index j = point[b];
      if (x(j, 0) - x(i, 0) > reach) {
        break;
      }
      if (x(j, 1) - x(i, 1) > reach) {
        // the rest of this run is further still, in its second coordinate
        b = run_end[b] - 1;
        continue;
      }
      const double d2 = (x.row(i) - x.row(j)).squaredNorm() / (scale * scale);
      const double kernel = std::exp(-0.5 * d2);
      count[a] += copies[b] * kernel;
      count[b] += copies[a] * kernel;
    }
  }
  Eigen::VectorXd out(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    out(i) = count[distinct[static_cast<size_t>(i)]];
  }
  return out;
}

inline void
TllBicop::fit(const Eigen::MatrixXd& data,
              std::string method,
              double mult,
              size_t grid_size,
              const Eigen::VectorXd& weights)
{
  using namespace tools_interpolation;

  // fitted in its own order and transposed back: a pair and its flip are the
  // same fit, bit for bit
  if (tools_stats::swaps_pair(data)) {
    Eigen::MatrixXd swapped = data;
    swapped.col(0).swap(swapped.col(1));
    if (swapped.cols() == 4) {
      swapped.col(2).swap(swapped.col(3));
    }
    std::swap(var_types_[0], var_types_[1]);
    fit(swapped, method, mult, grid_size, weights);
    std::swap(var_types_[0], var_types_[1]);
    interp_grid_->flip();
    return;
  }

  // construct default grid (equally spaced on Gaussian scale)
  auto grid_points = this->make_normal_grid(grid_size);

  // expand the interpolation grid; a matrix with two columns where each row
  // contains one combination of the grid points
  auto grid_2d = tools_eigen::expand_grid(grid_points);

  // transform evaluation grid and data by inverse Gaussian cdf
  Eigen::MatrixXd z = tools_stats::qnorm(grid_2d);

  bool discrete = (var_types_[0] == "d") || (var_types_[1] == "d");
  // ties broken at random; on a discrete edge, soft ranks
  Eigen::MatrixXd psobs = tools_stats::pair_soft_pseudo_obs(
    data, weights, discrete ? tools_stats::default_soft_scale() : 0.0);
  Eigen::MatrixXd z_data = tools_stats::qnorm(psobs);

  // find bandwidth matrix
  Eigen::Matrix2d B = select_bandwidth(z_data, method, weights);
  B *= mult;

  // find latent sample in case observations are discrete
  if (discrete) {
    psobs =
      tools_stats::find_latent_sample(data, std::pow(B(0, 0) * B(1, 1), 0.25));
    z_data = tools_stats::qnorm(psobs);
  }

  // compute the density estimator (first column estimate, second influence)
  Eigen::MatrixXd ll_fit = fit_local_likelihood(z, z_data, B, method, weights);

  // transform density estimate to copula scale
  Eigen::VectorXd c =
    ll_fit.col(0).cwiseQuotient(tools_stats::dnorm(z).rowwise().prod());
  // store values in mxm grid
  Eigen::MatrixXd values(grid_size, grid_size);
  values =
    Eigen::Map<Eigen::MatrixXd>(c.data(), grid_size, grid_size).transpose();

  // create interpolation grid
  interp_grid_ = std::make_shared<InterpolationGrid>(grid_points, values);

  // compute effective degrees of freedom via interpolation ---------
  // stabilize interpolation by restricting to plausible range
  Eigen::VectorXd infl_vec = ll_fit.col(1).cwiseMin(1.3).cwiseMax(-0.2);
  Eigen::MatrixXd infl(grid_size, grid_size);
  infl = Eigen::Map<Eigen::MatrixXd>(infl_vec.data(), grid_size, grid_size)
           .transpose();
  // don't normalize margins of the EDF! (norm_times = 0)
  auto infl_grid = InterpolationGrid(grid_points, infl, 0);
  if (discrete) {
    // for discrete, use mid ranks to compute EDF and log-likelihood
    // (this is closer to "observations" than jittered or "upper" pseudo data);
    // an observation counts once however often it is repeated
    psobs = 0.5 * (data.leftCols(2) + data.rightCols(2)).array();
    npars_ =
      infl_grid.interpolate(psobs)
        .cwiseQuotient(multiplicity(psobs, tools_stats::default_soft_scale()))
        .sum();
    npars_ = std::max(npars_, 1.0);
  } else {
    npars_ = std::max(infl_grid.interpolate(data).sum(), 1.0);
  }
  set_loglik(pdf(data).array().log().sum());
}
}
