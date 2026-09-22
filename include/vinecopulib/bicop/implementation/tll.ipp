// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <array>
#include <boost/math/special_functions/bessel.hpp>
#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/tools_circular.hpp>
#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_var_types.hpp>
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

  // de-correlate data by applying B^{-1/2}
  Eigen::MatrixXd z = (irB * x.transpose()).transpose();
  Eigen::MatrixXd z_data = (irB * x_data.transpose()).transpose();

  Eigen::MatrixXd res(m, 2);
  res.col(0) = Eigen::VectorXd::Ones(m); // result will be a product
  Eigen::VectorXd kernels(n);
  Eigen::Vector2d f1;
  Eigen::Vector2d b;
  Eigen::Matrix2d S(B);
  Eigen::MatrixXd zz(n, 2), zz2(n, 2);
  for (size_t k = 0; k < m; ++k) {
    zz = z_data.rowwise() - z.row(k);
    kernels = gaussian_kernel_2d(zz) * det_irB;
    if (weights.size() > 0)
      kernels = kernels.cwiseProduct(weights);
    double f0 = kernels.mean();
    if (method != "constant") {
      zz = (irB * zz.transpose()).transpose();
      f1 = (zz.array().colwise() * kernels.array()).colwise().mean();
      b = f1 / f0;
      if (method == "quadratic") {
        zz2 = (zz.array().colwise() * kernels.array()).matrix() /
              (f0 * static_cast<double>(n));
        b = B * b;
        S = (B * (zz.transpose() * zz2) * B - b * b.transpose()).inverse();
        res(k) *= std::sqrt(S.determinant()) / det_irB;
      }
      res(k) *= std::exp(-0.5 * double(b.transpose() * S * b));
      if ((std::isnan)(res(k)) || (std::isinf)(res(k))) {
        // inverse operation might go wrong due to rounding when
        // true value is equal or close to zero
        res(k) = 0.0;
      }
    }
    res(k, 0) *= f0;
    if (weights.size() > 0) {
      // average weight in neighborhood of evaluation point (essentially a
      // kernel regression estimate);
      // kernels have already been multiplied with weights above
      double w = kernels.sum() / kernels.cwiseQuotient(weights).sum();
      res(k, 1) = calculate_infl(n, f0, b, B, det_irB, S, method, w);
    } else {
      res(k, 1) = calculate_infl(n, f0, b, B, det_irB, S, method, 1.0);
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
inline double
TllBicop::calculate_infl(const size_t& n,
                         const double& f0,
                         const Eigen::Vector2d& b,
                         const Eigen::Matrix2d& B,
                         const double& det_irB,
                         const Eigen::Matrix2d& S,
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

  double m_inv_00;
  if (method == "linear") {
    Eigen::Matrix3d M;
    M(0, 0) = f0;
    M.col(0).tail(2) = B * b * f0;
    M.row(0).tail(2) = M.col(0).tail(2);
    M.block(1, 1, 2, 2) = f0 * B + f0 * B * b * b.transpose() * B;
    m_inv_00 = M.inverse()(0, 0);
  } else {
    Eigen::Matrix<double, 6, 6> M = Eigen::Matrix<double, 6, 6>::Zero();
    M(0, 0) = f0;
    M.col(0).segment(1, 2) = f0 * b;
    M.row(0).segment(1, 2) = M.col(0).segment(1, 2);
    M.block(1, 1, 2, 2) = f0 * B + f0 * b * b.transpose();
    M(3, 0) = 0.5 * M(1, 1);
    M(4, 0) = 0.5 * M(2, 2);
    M(5, 0) = M(1, 2);
    M.row(0).tail(3) = M.col(0).tail(3);
    Eigen::MatrixXd Si = S.inverse();
    M(3, 1) = 0.5 * f0 * (3.0 * Si(0, 0) * b(0) + std::pow(b(0), 3));
    M(4, 2) = 0.5 * f0 * (3.0 * Si(1, 1) * b(1) + std::pow(b(1), 3));
    M(4, 1) = 0.5 * f0;
    M(4, 1) *= 2.0 * Si(0, 1) * b(1) + Si(1, 1) * b(0) + b(0) * b(1) * b(1);
    M(3, 2) = 0.5 * f0;
    M(3, 2) *= 2.0 * Si(0, 1) * b(0) + Si(0, 0) * b(1) + b(1) * b(0) * b(0);
    M(5, 1) = 2.0 * M(3, 2);
    M(5, 2) = 2.0 * M(4, 1);
    M.block(1, 3, 2, 3) = M.block(3, 1, 3, 2).transpose();
    M(3, 3) = 0.25 * f0;
    M(3, 3) *= 3.0 * Si(0, 0) * Si(0, 0) + 6.0 * Si(0, 0) * b(0) * b(0) +
               std::pow(b(0), 4);
    M(4, 4) = 0.25 * f0;
    M(4, 4) *= 3.0 * Si(1, 1) * Si(1, 1) + 6.0 * Si(1, 1) * b(1) * b(1) +
               std::pow(b(1), 4);
    M(5, 5) = Si(0, 0) * Si(1, 1) + 2.0 * S(0, 1) + b(0) * b(0) * b(1) * b(1);
    M(5, 5) += 4.0 * Si(0, 1) * b(0) * b(1);
    M(5, 5) += Si(0, 0) * b(1) * b(1) + Si(1, 1) * b(0) * b(0);
    M(5, 5) *= f0;
    M(4, 3) = M(5, 5) * 0.25;
    M(3, 4) = M(4, 3);
    M(5, 3) = 3.0 * Si(0, 0) * Si(0, 1) + 3.0 * Si(0, 1) * b(0) * b(0);
    M(5, 3) += 3.0 * Si(0, 0) * b(0) * b(1) + b(1) * std::pow(b(0), 3);
    M(5, 3) *= 0.5 * f0;
    M(3, 5) = M(5, 3);
    M(5, 4) = 3.0 * Si(1, 1) * Si(0, 1) + 3.0 * Si(0, 1) * b(1) * b(1);
    M(5, 4) += 3.0 * Si(1, 1) * b(0) * b(1) + b(0) * std::pow(b(1), 3);
    M(5, 4) *= 0.5 * f0;
    M(4, 5) = M(5, 4);
    m_inv_00 = M.inverse()(0, 0);
  }

  return kernel0 * det_irB * m_inv_00 * weight / static_cast<double>(n);
}

inline void
TllBicop::fit(const Eigen::MatrixXd& data,
              std::string method,
              double mult,
              size_t grid_size,
              const Eigen::VectorXd& weights)
{
  using namespace tools_interpolation;

  if (tools_var_types::any_circular(var_types_)) {
    fit_mixed(data, method, mult, grid_size, weights);
    return;
  }

  // construct default grid (equally spaced on Gaussian scale)
  auto grid_points = this->make_normal_grid(grid_size);

  // expand the interpolation grid; a matrix with two columns where each row
  // contains one combination of the grid points
  auto grid_2d = tools_eigen::expand_grid(grid_points);

  // transform evaluation grid and data by inverse Gaussian cdf
  Eigen::MatrixXd z = tools_stats::qnorm(grid_2d);

  // use jittering in case observations are discrete
  auto psobs =
    tools_stats::to_pseudo_obs(data.leftCols(2), "random", weights, { 5 });
  Eigen::MatrixXd z_data = tools_stats::qnorm(psobs);

  // find bandwidth matrix
  Eigen::Matrix2d B = select_bandwidth(z_data, method, weights);
  B *= mult;

  // find latent sample in case observations are discrete
  if (!tools_var_types::all_continuous(var_types_)) {
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
  if (!tools_var_types::all_continuous(var_types_)) {
    // for discrete, use mid ranks to compute EDF and log-likelihood
    // (this is closer to "observations" than jittered or "upper" pseudo data)
    psobs = 0.5 * (data.leftCols(2) + data.rightCols(2)).array();
    npars_ = tools_eigen::unique(infl_grid.interpolate(psobs)).sum();
    npars_ = std::max(npars_, 1.0);
  } else {
    npars_ = std::max(infl_grid.interpolate(data).sum(), 1.0);
  }
  set_loglik(pdf(data).array().log().sum());
}

// ---------------------------------------------------------------------------
// pairs with a circular variable

//! @brief Transforms one axis for the mixed-geometry estimator.
//!
//! A linear axis is mapped to the normal scale, where the default knots are
//! equally spaced; a circular axis is mapped to the angle \f$ 2\pi u \f$ on
//! equally spaced knots of the unit interval, so that the two ends of the
//! knot vector are the same point of the circle.
inline TllBicop::Axis
TllBicop::make_axis(const std::string& var_type,
                    const Eigen::VectorXd& u,
                    size_t grid_size)
{
  Axis axis;
  axis.circular = tools_var_types::is_circular(var_type);
  axis.knots = make_grid_points(var_type, grid_size);
  if (axis.circular) {
    const double two_pi = tools_circular::two_pi();
    axis.x = two_pi * u;
    axis.grid = two_pi * axis.knots;
    axis.jacobian = Eigen::VectorXd::Constant(grid_size, two_pi);
  } else {
    axis.x = tools_stats::qnorm(u);
    axis.grid = tools_stats::qnorm(axis.knots);
    axis.jacobian = tools_stats::dnorm(axis.grid).cwiseInverse();
  }
  return axis;
}

//! @brief The bandwidth of one axis.
//!
//! The variance of the kernel is the same fraction of the variance of the
//! transformed margin as in `select_bandwidth()`: \f$ n^{-1/3} \f$ for the
//! local constant fit and \f$ 1.5 n^{-1/3} \f$ for the local linear one, times
//! the multiplier and a dependence factor. The probit-transformed margin has
//! unit variance; the angle of a uniform circular variable has variance
//! \f$ (2\pi)^2 / 12 \f$. For a circular axis the result is the von Mises
//! concentration \f$ 1 / \sigma^2 \f$, for a linear axis the standard
//! deviation \f$ \sigma \f$.
inline double
TllBicop::select_bandwidth_mixed(const Axis& axis,
                                 size_t n,
                                 const std::string& method,
                                 double dependence)
{
  const double nn = static_cast<double>(n);
  double base = (method == "constant") ? std::pow(nn, -1.0 / 3.0)
                                       : 1.5 * std::pow(nn, -1.0 / 3.0);
  // stronger dependence concentrates the density, so it needs a finer kernel
  base *= std::max(1.0 - dependence, 0.1);
  double variance = base;
  if (axis.circular) {
    const double two_pi = tools_circular::two_pi();
    variance *= two_pi * two_pi / 12.0;
    return std::min(1.0 / variance, 500.0);
  }
  return std::sqrt(variance);
}

//! @brief The local likelihood fit at one knot of the mixed-geometry grid.
//!
//! The kernel is a product of a von Mises kernel on each circular axis and a
//! Gaussian kernel on each linear one. The local model is log-linear per
//! axis: linear in \f$ z \f$ on a linear axis, linear in \f$ (\cos\theta,
//! \sin\theta) \f$ on a circular one. With a product kernel and no
//! interaction term the local likelihood equations separate, and each axis
//! contributes a closed-form correction to the kernel density value: for a
//! linear axis \f$ \exp(-m^2 / 2\sigma^2) \f$ with the local mean \f$ m \f$,
//! for a circular axis \f$ \exp(r\, \hat m \cdot e_0 - \kappa)\, I_0(\kappa)
//! / I_0(r) \f$, where \f$ \hat m \f$ is the local mean direction, \f$ r =
//! A^{-1}(\bar R) \f$ the concentration matching the local mean resultant
//! length, and \f$ e_0 \f$ the direction of the knot.
//!
//! @return The density estimate at the knot and the influence of an
//!   observation at the knot on it.
inline std::pair<double, double>
TllBicop::local_fit_mixed(const std::vector<Axis>& axes,
                          const std::array<Eigen::Index, 2>& knot,
                          const std::string& method,
                          const Eigen::VectorXd& weights,
                          Eigen::VectorXd& kernels)
{
  const Eigen::Index n = axes[0].x.size();
  const double nn = static_cast<double>(n);
  const double two_pi = tools_circular::two_pi();
  const double sqrt_two_pi = std::sqrt(two_pi);

  // product kernel and its value at zero
  kernels.setOnes(n);
  double kernel0 = 1.0;
  std::array<Eigen::VectorXd, 2> diffs;
  for (size_t a = 0; a < 2; ++a) {
    const Axis& axis = axes[a];
    diffs[a] = axis.x.array() - axis.grid(knot[a]);
    if (axis.circular) {
      const double kappa = axis.scale;
      // scaled by exp(-kappa) so that large concentrations do not overflow
      const double norm =
        two_pi * boost::math::cyl_bessel_i(0, kappa) * std::exp(-kappa);
      kernels.array() *= (kappa * (diffs[a].array().cos() - 1.0)).exp() / norm;
      kernel0 /= norm;
    } else {
      const double sd = axis.scale;
      kernels.array() *=
        (-0.5 * (diffs[a].array() / sd).square()).exp() / (sd * sqrt_two_pi);
      kernel0 /= sd * sqrt_two_pi;
    }
  }
  if (weights.size() > 0) {
    kernels = kernels.cwiseProduct(weights);
  }
  const double f0 = kernels.sum() / nn;
  const double weight = (weights.size() > 0)
                          ? kernels.sum() / kernels.cwiseQuotient(weights).sum()
                          : 1.0;

  if (method == "constant") {
    return { f0, kernel0 / f0 * weight / nn };
  }

  // local log-linear fit: per-axis corrections and the moments of the local
  // model under the kernel, from which the influence follows
  double f = f0;
  std::vector<Eigen::VectorXd> first;  // E[psi_a]
  std::vector<Eigen::MatrixXd> second; // E[psi_a psi_a^T]
  for (size_t a = 0; a < 2; ++a) {
    const Axis& axis = axes[a];
    if (axis.circular) {
      const double kappa = axis.scale;
      Eigen::Vector2d m;
      m << (kernels.array() * diffs[a].array().cos()).sum() / kernels.sum(),
        (kernels.array() * diffs[a].array().sin()).sum() / kernels.sum();
      const double rbar = std::min(m.norm(), 0.9999);
      Eigen::Vector2d dir(1.0, 0.0);
      if (rbar > 1e-12) {
        dir = m / m.norm();
      }
      const double r = tools_circular::von_mises_a_inverse(rbar, 500.0);
      const double log_i0_kappa = std::log(boost::math::cyl_bessel_i(0, kappa));
      const double log_i0_r = std::log(boost::math::cyl_bessel_i(0, r));
      f *= std::exp(r * dir(0) - kappa + log_i0_kappa - log_i0_r);

      // moments of the local model, a von Mises law with direction `dir` and
      // concentration `r`, in the basis (cos - 1, sin) centered at the knot
      const double a1 = tools_circular::von_mises_a(r);
      const double a2 = (r > 0.0) ? boost::math::cyl_bessel_i(2, r) /
                                      boost::math::cyl_bessel_i(0, r)
                                  : 0.0;
      const double theta_m = std::atan2(dir(1), dir(0));
      Eigen::Vector2d e0(1.0, 0.0);
      Eigen::Vector2d mean = a1 * dir;
      Eigen::Matrix2d ee;
      ee << 0.5 * (1.0 + a2 * std::cos(2.0 * theta_m)),
        0.5 * a2 * std::sin(2.0 * theta_m), 0.5 * a2 * std::sin(2.0 * theta_m),
        0.5 * (1.0 - a2 * std::cos(2.0 * theta_m));
      first.push_back(mean - e0);
      second.push_back(ee - mean * e0.transpose() - e0 * mean.transpose() +
                       e0 * e0.transpose());
    } else {
      const double sd = axis.scale;
      const double m1 =
        (kernels.array() * diffs[a].array()).sum() / kernels.sum();
      f *= std::exp(-0.5 * m1 * m1 / (sd * sd));
      first.push_back(Eigen::VectorXd::Constant(1, m1));
      second.push_back(Eigen::MatrixXd::Constant(1, 1, sd * sd + m1 * m1));
    }
  }
  if (!std::isfinite(f)) {
    // the corrections can overflow where the true value is (close to) zero
    f = 0.0;
  }

  // the local information matrix M = f0 E[(1, psi)(1, psi)^T]; the axes are
  // independent under the product model, so the cross block is E[psi_a]
  // E[psi_b]^T
  const Eigen::Index p1 = first[0].size();
  const Eigen::Index p2 = first[1].size();
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(1 + p1 + p2, 1 + p1 + p2);
  M(0, 0) = 1.0;
  M.block(1, 0, p1, 1) = first[0];
  M.block(1 + p1, 0, p2, 1) = first[1];
  M.block(0, 1, 1, p1) = first[0].transpose();
  M.block(0, 1 + p1, 1, p2) = first[1].transpose();
  M.block(1, 1, p1, p1) = second[0];
  M.block(1 + p1, 1 + p1, p2, p2) = second[1];
  M.block(1, 1 + p1, p1, p2) = first[0] * first[1].transpose();
  M.block(1 + p1, 1, p2, p1) = first[1] * first[0].transpose();
  M *= f0;
  double infl = kernel0 * M.inverse()(0, 0) * weight / nn;
  if (!std::isfinite(infl)) {
    infl = 0.0;
  }
  return { f, infl };
}

//! @brief Fits the estimator for a pair with a circular variable.
inline void
TllBicop::fit_mixed(const Eigen::MatrixXd& data,
                    const std::string& method,
                    double mult,
                    size_t grid_size,
                    const Eigen::VectorXd& weights)
{
  using namespace tools_interpolation;
  if (method == "quadratic") {
    throw std::runtime_error(
      "nonparametric_method 'quadratic' is not available for pairs with a "
      "circular variable; use 'constant' or 'linear'");
  }
  const size_t n = data.rows();

  // a linear axis is rank-transformed as in the linear estimator; a circular
  // axis keeps its values, since ranks would move with the cut and the
  // estimate would no longer rotate with the data
  Eigen::MatrixXd psobs =
    tools_stats::to_pseudo_obs(data.leftCols(2), "random", weights, { 5 });
  for (Eigen::Index a = 0; a < 2; ++a) {
    if (tools_var_types::is_circular(var_types_[a])) {
      psobs.col(a) = data.col(a);
    }
  }
  std::vector<Axis> axes = { make_axis(var_types_[0], psobs.col(0), grid_size),
                             make_axis(
                               var_types_[1], psobs.col(1), grid_size) };
  const double dependence =
    tools_stats::pairwise_circular(psobs, var_types_, weights);
  for (auto& axis : axes) {
    axis.scale = select_bandwidth_mixed(axis, n, method, dependence);
    if (axis.circular) {
      axis.scale = std::min(axis.scale / mult, 500.0);
    } else {
      axis.scale *= std::sqrt(mult);
    }
  }

  // the estimate and the influence on every knot
  Eigen::MatrixXd values(grid_size, grid_size), infl(grid_size, grid_size);
  Eigen::VectorXd kernels(n);
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(grid_size); ++i) {
    for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(grid_size); ++j) {
      // the two ends of a circular axis are the same point
      const bool i_wraps =
        axes[0].circular && (i + 1 == static_cast<Eigen::Index>(grid_size));
      const bool j_wraps =
        axes[1].circular && (j + 1 == static_cast<Eigen::Index>(grid_size));
      if (i_wraps || j_wraps) {
        values(i, j) = values(i_wraps ? 0 : i, j_wraps ? 0 : j);
        infl(i, j) = infl(i_wraps ? 0 : i, j_wraps ? 0 : j);
        continue;
      }
      auto fit = local_fit_mixed(axes, { i, j }, method, weights, kernels);
      values(i, j) = fit.first * axes[0].jacobian(i) * axes[1].jacobian(j);
      infl(i, j) = fit.second;
    }
  }
  if (weights.size() > 0) {
    // estimate can be negative if negative weights are used
    values = values.array().max(0.0);
  }

  interp_grid_ =
    std::make_shared<InterpolationGrid>(axes[0].knots, axes[1].knots, values);

  // effective degrees of freedom via interpolation of the influence,
  // restricted to a plausible range; the margins of the influence are not
  // normalized
  infl = infl.array().min(1.3).max(-0.2);
  auto infl_grid = InterpolationGrid(axes[0].knots, axes[1].knots, infl, 0);
  npars_ = std::max(infl_grid.interpolate(data.leftCols(2)).sum(), 1.0);
  set_loglik(pdf(data).array().log().sum());
}
}
