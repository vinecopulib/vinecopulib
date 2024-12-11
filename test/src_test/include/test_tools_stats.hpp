// Copyright © 2016-2023 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "test_vinecop_sanity_checks.hpp"
#include "gtest/gtest.h"
#include <vinecopulib.hpp>

namespace test_tools_stats {

using namespace vinecopulib;

TEST(test_tools_stats, to_pseudo_obs_is_correct)
{

  int n = 9;

  // X1 = (1,...,n) and X2 = (n, ..., 1)
  // X = (X1, X2)
  Eigen::MatrixXd X(n, 2);
  X.col(0) = Eigen::VectorXd::LinSpaced(n, 1, n);
  X.col(1) = Eigen::VectorXd::LinSpaced(n, n, 1);

  // U = pobs(X)
  Eigen::MatrixXd U = tools_stats::to_pseudo_obs(X);
  for (int i = 0; i < 9; i++) {
    EXPECT_NEAR(U(i, 0), (i + 1.0) * 0.1, 1e-2);
    EXPECT_NEAR(U(i, 1), 1.0 - (i + 1.0) * 0.1, 1e-2);
  }

  Eigen::MatrixXd X2 = tools_stats::simulate_uniform(100, 2);
  EXPECT_NO_THROW(tools_stats::to_pseudo_obs(X2, "random"));
  EXPECT_NO_THROW(tools_stats::to_pseudo_obs(X2, "first"));
  EXPECT_ANY_THROW(tools_stats::to_pseudo_obs(X2, "something"));

  X2.col(0).head(50) = Eigen::VectorXd::Constant(50, NAN);
  auto u = tools_stats::to_pseudo_obs(X2);
  EXPECT_TRUE(std::isnan(u(0, 0)));
  EXPECT_GE(u.col(0).tail(50).maxCoeff(), 0.98);
}

TEST(test_tools_stats, qrng_are_correct)
{

  size_t d = 2;
  size_t n = 10;
  size_t N = 1000;
  double Nd = static_cast<double>(N);

  auto cop = Bicop(BicopFamily::gaussian);
  auto u = cop.simulate(n);
  auto U = tools_stats::ghalton(N, d);
  auto U1 = tools_stats::sobol(N, d);
  auto U2 = tools_stats::simulate_uniform(N, d);

  Eigen::VectorXd x(N), p(n), p1(N), x2(N), p2(n);
  p2 = Eigen::VectorXd::Zero(n);
  for (size_t i = 0; i < n; i++) {
    auto f = [i, u](const double& u1, const double& u2) {
      return (u1 <= u(i, 0) && u2 <= u(i, 1)) ? 1.0 : 0.0;
    };
    x = U.col(0).binaryExpr(cop.hinv1(U), f);
    p(i) = x.sum() / Nd;
    x = U1.col(0).binaryExpr(cop.hinv1(U1), f);
    p1(i) = x.sum() / Nd;
    x2 = U2.col(0).binaryExpr(cop.hinv1(U2), f);
    p2(i) = x2.sum() / Nd;
  }

  x = cop.cdf(u);
  if (p2.isApprox(x, 1e-2)) {
    ASSERT_TRUE(p.isApprox(x, 1e-2));
    ASSERT_TRUE(p1.isApprox(x, 1e-2));
  }
}

TEST(test_tools_stats, mcor_works)
{
  std::vector<int> seeds = { 1, 2, 3, 4, 5 };
  Eigen::MatrixXd Z = tools_stats::simulate_uniform(10000, 2, true, seeds);
  Z = tools_stats::qnorm(Z);
  Z.block(0, 1, 5000, 1) =
    Z.block(0, 1, 5000, 1) + Z.block(0, 0, 5000, 1).cwiseAbs2();
  auto a1 = tools_stats::pairwise_mcor(Z);
  Eigen::VectorXd weights = Eigen::VectorXd::Ones(10000);
  auto a2 = tools_stats::pairwise_mcor(Z, weights);
  ASSERT_TRUE(std::fabs(a1 - a2) < 1e-4);

  a1 = tools_stats::pairwise_mcor(Z.block(0, 0, 5000, 2));
  weights.block(5000, 0, 5000, 1) = Eigen::VectorXd::Zero(5000);
  a2 = tools_stats::pairwise_mcor(Z, weights);
  ASSERT_TRUE(std::fabs(a1 - a2) < 0.05);
}

TEST(test_tools_stats, seed_works)
{
  size_t d = 2;
  size_t n = 10;
  std::vector<int> v = { 1, 2, 3 };

  auto U1 = tools_stats::simulate_uniform(n, d);
  auto U2 = tools_stats::simulate_uniform(n, d, false, v);
  auto U3 = tools_stats::simulate_uniform(n, d, false, v);

  ASSERT_TRUE(U1.cwiseNotEqual(U2).all());
  ASSERT_TRUE(U2.cwiseEqual(U3).all());
}

TEST(test_tools_stats, dpqnorm_are_nan_safe)
{
  Eigen::VectorXd X = Eigen::VectorXd::Random(10);
  X(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_NO_THROW(tools_stats::dnorm(X));
  EXPECT_NO_THROW(tools_stats::pnorm(X));
  EXPECT_NO_THROW(tools_stats::qnorm(tools_stats::pnorm(X)));
}

TEST(test_tools_stats, dpt_are_nan_safe)
{
  Eigen::VectorXd X = Eigen::VectorXd::Random(10);
  X(0) = std::numeric_limits<double>::quiet_NaN();
  double nu = 4.0;
  EXPECT_NO_THROW(tools_stats::dt(X, nu));
  EXPECT_NO_THROW(tools_stats::pt(X, nu));
  EXPECT_NO_THROW(tools_stats::qt(tools_stats::pt(X, nu), nu));
}

TEST(test_tools_stats, pbvt_and_pbvnorm_are_nan_safe)
{
  Eigen::MatrixXd X = Eigen::MatrixXd::Random(10, 2);
  X(0) = std::numeric_limits<double>::quiet_NaN();
  double rho = -0.95;
  int nu = 5;
  EXPECT_NO_THROW(tools_stats::pbvt(X, nu, rho));
  EXPECT_NO_THROW(tools_stats::pbvnorm(X, rho));
}
}
