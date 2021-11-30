// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/math/distributions.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_stats {

//! @brief Density function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Matrix
dnorm(const Matrix& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Distribution function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Matrix
pnorm(const Matrix& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Quantile function of the Standard normal distribution.
//!
//! @param x Evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Matrix
qnorm(const Matrix& x)
{
  boost::math::normal dist;
  auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Density function of the Student t distribution.
//!
//! @param x Evaluation points.
//! @param nu Degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Matrix
dt(const Matrix& x, double nu)
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
inline Matrix
pt(const Matrix& x, double nu)
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
inline Matrix
qt(const Matrix& x, double nu)
{
  boost::math::students_t dist(nu);
  auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
  return tools_eigen::unaryExpr_or_nan(x, f);
}

Matrix
simulate_uniform(const size_t& n,
                 const size_t& d,
                 bool qrng = false,
                 std::vector<int> seeds = std::vector<int>());

Eigen::VectorXd
to_pseudo_obs_1d(Eigen::VectorXd x, const std::string& ties_method = "average");

Matrix
to_pseudo_obs(Matrix x, const std::string& ties_method = "average");

double
pairwise_mcor(const Matrix& x,
              const Eigen::VectorXd& weights = Eigen::VectorXd());

Matrix
dependence_matrix(const Matrix& x, const std::string& measure);

Matrix
ghalton(const size_t& n,
        const size_t& d,
        const std::vector<int>& seeds = std::vector<int>());

Matrix
sobol(const size_t& n,
      const size_t& d,
      const std::vector<int>& seeds = std::vector<int>());

Eigen::VectorXd
pbvt(const Matrix& z, int nu, double rho);

Eigen::VectorXd
pbvnorm(const Matrix& z, double rho);
}
}

#include <vinecopulib/misc/implementation/tools_stats.ipp>
