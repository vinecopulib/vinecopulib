// Copyright © 2016-2023 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/math/distributions.hpp>
#include <memory>
#include <set>
#include <vinecopulib/misc/tools_eigen.hpp>
#include <unsupported/Eigen/SpecialFunctions>


namespace vinecopulib {

namespace tools_stats {

//! @brief Density function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Eigen::MatrixXd
dnorm(const Eigen::MatrixXd& x)
{
  static const double pi = 3.14159265358979323846;
  static const double sqrt_2pi = std::sqrt(2.0 * pi);
  return (1.0 / sqrt_2pi) * (-0.5 * x.array().square()).exp();
}

//! @brief Distribution function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Eigen::MatrixXd
pnorm(const Eigen::MatrixXd& x)
{
  static const double sqrt2 = std::sqrt(2.0);
  return 0.5 * (1.0 + (x.array() / sqrt2).erf());
}

//! @brief Quantile function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd
qnorm(const Eigen::MatrixXd& x)
{
  return x.array().ndtri();
}

//! @brief Quantile function of the Standard normal distribution with additional
//! bound checks for numerical stability.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd
safe_qnorm(const Eigen::MatrixXd& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) {
    if (y <= 0) {
      return -std::numeric_limits<double>::infinity();
    } else if (y >= 1) {
      return std::numeric_limits<double>::infinity();
    } else {
      return boost::math::quantile(dist, y);
    }
  };

  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Distribution function of the Standard normal distribution with
//! additional bound checks for numerical stability.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd
safe_pnorm(const Eigen::MatrixXd& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) {
    if (y >= std::numeric_limits<double>::max()) {
      return 1.0;
    } else if (y <= -std::numeric_limits<double>::max()) {
      return 0.0;
    } else {
      return boost::math::cdf(dist, y);
    }
  };

  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Density function of the Student t distribution.
//!
//! @param x Evaluation points.
//! @param nu Degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Eigen::MatrixXd
dt(const Eigen::MatrixXd& x, double nu)
{
  boost::math::students_t dist(nu);
  auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Distribution function of the Student t distribution.
//!
//! @param x Evaluation points.
//! @param nu Degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Eigen::MatrixXd
pt(const Eigen::MatrixXd& x, double nu)
{
  boost::math::students_t dist(nu);
  auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Quantile function of the Student t distribution.
//!
//! @param x Evaluation points.
//! @param nu Degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd
qt(const Eigen::MatrixXd& x, double nu)
{
  boost::math::students_t dist(nu);
  auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

Eigen::MatrixXd
simulate_uniform(const size_t& n,
                 const size_t& d,
                 bool qrng = false,
                 std::vector<int> seeds = std::vector<int>());

Eigen::MatrixXd
simulate_normal(const size_t& n,
                const size_t& d,
                std::vector<int> seeds = std::vector<int>());

Eigen::VectorXd
to_pseudo_obs_1d(Eigen::VectorXd x, const std::string& ties_method = "average");

Eigen::MatrixXd
to_pseudo_obs(Eigen::MatrixXd x, const std::string& ties_method = "average");

class BoxCovering
{
public:
  explicit BoxCovering(const Eigen::MatrixXd& u, uint16_t K = 40)
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

  std::vector<size_t> get_box_indices(const Eigen::VectorXd& lower,
                                      const Eigen::VectorXd& upper)
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

  void swap_sample(size_t i, const Eigen::VectorXd& new_sample)
  {
    auto k = static_cast<size_t>(std::floor(u_(i, 0) * K_));
    auto j = static_cast<size_t>(std::floor(u_(i, 1) * K_));
    boxes_[k][j]->indices_.erase(i);

    u_.row(i) = new_sample;
    k = static_cast<size_t>(std::floor(new_sample(0) * K_));
    j = static_cast<size_t>(std::floor(new_sample(1) * K_));
    boxes_[k][j]->indices_.insert(i);
  }

private:
  struct Box
  {
  public:
    Box(const std::vector<double>& lower, const std::vector<double>& upper)
      : lower_(lower)
      , upper_(upper)
    {
    }

    std::vector<double> lower_;
    std::vector<double> upper_;
    std::set<size_t> indices_;
  };

  Eigen::MatrixXd u_;
  size_t n_;
  uint16_t K_;
  std::vector<std::vector<std::unique_ptr<Box>>> boxes_;
};

Eigen::MatrixXd
find_latent_sample(const Eigen::MatrixXd& u, double b, size_t niter = 3);

double
pairwise_mcor(const Eigen::MatrixXd& x,
              const Eigen::VectorXd& weights = Eigen::VectorXd());

Eigen::MatrixXd
dependence_matrix(const Eigen::MatrixXd& x, const std::string& measure);

Eigen::MatrixXd
ghalton(const size_t& n,
        const size_t& d,
        const std::vector<int>& seeds = std::vector<int>());

Eigen::MatrixXd
sobol(const size_t& n,
      const size_t& d,
      const std::vector<int>& seeds = std::vector<int>());

Eigen::VectorXd
pbvt(const Eigen::MatrixXd& z, int nu, double rho);

Eigen::VectorXd
pbvnorm(const Eigen::MatrixXd& z, double rho);
}
}

#include <vinecopulib/misc/implementation/tools_stats.ipp>
