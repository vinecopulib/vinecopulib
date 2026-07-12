// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "test_utils.hpp"
#include "test_vinecop_sanity_checks.hpp"
#include "gtest/gtest.h"
#include <vinecopulib.hpp>

namespace test_tools_stats {

using namespace vinecopulib;
using test_utils::all_close;

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

  Eigen::MatrixXd X2 = tools_stats::simulate_uniform(100, 2, false, { 1 });
  EXPECT_NO_THROW(tools_stats::to_pseudo_obs(X2, "random"));
  EXPECT_NO_THROW(tools_stats::to_pseudo_obs(X2, "first"));
  EXPECT_ANY_THROW(tools_stats::to_pseudo_obs(X2, "something"));

  auto weights = Eigen::VectorXd::Constant(100, 1.0);
  auto r1 = tools_stats::to_pseudo_obs(X2, "average");
  auto r2 = tools_stats::to_pseudo_obs(X2, "average", weights);
  EXPECT_TRUE(r1 == r2);

  r1 = tools_stats::to_pseudo_obs(X2, "first");
  r2 = tools_stats::to_pseudo_obs(X2, "first", weights);
  EXPECT_TRUE(r1 == r2);

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

  // seeded so the simulated evaluation points (and hence the test) are
  // deterministic across runs and platforms
  std::vector<int> seeds = { 1, 2, 3, 4, 5 };
  auto cop = Bicop(BicopFamily::gaussian);
  auto u = cop.simulate(n, false, seeds);
  auto U = tools_stats::ghalton(N, d);
  auto U1 = tools_stats::sobol(N, d);

  // Monte-Carlo CDF estimate at each simulated point using each low-discrepancy
  // sequence: p(i) = (1/N) sum_j 1{U_j1 <= u_i1, hinv1(U_j) <= u_i2}. ghalton
  // and sobol converge much faster than plain Monte Carlo, so at N = 1000 they
  // recover the analytical copula CDF to well within 1e-2.
  Eigen::VectorXd x(N), p(n), p1(n);
  for (size_t i = 0; i < n; i++) {
    auto f = [i, u](const double& u1, const double& u2) {
      return (u1 <= u(i, 0) && u2 <= u(i, 1)) ? 1.0 : 0.0;
    };
    x = U.col(0).binaryExpr(cop.hinv1(U), f);
    p(i) = x.sum() / Nd;
    x = U1.col(0).binaryExpr(cop.hinv1(U1), f);
    p1(i) = x.sum() / Nd;
  }

  x = cop.cdf(u);
  ASSERT_TRUE(all_close(p, x, 1e-2, 1e-2));  // ghalton
  ASSERT_TRUE(all_close(p1, x, 1e-2, 1e-2)); // sobol
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

TEST(test_tools_stats, dpqnorm_work)
{
  auto dnorm_boost = [](Eigen::MatrixXd x) {
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
  };

  auto pnorm_boost = [](Eigen::MatrixXd x) {
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
  };

  auto qnorm_boost = [](Eigen::MatrixXd x) {
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
  };

  // linspace from -5 to 5 (1000 points)
  Eigen::VectorXd X = Eigen::VectorXd::LinSpaced(1000, -5, 5);

  // tools_stats::dnorm is the same as dnorm_boost
  auto d1 = tools_stats::dnorm(X);
  auto d2 = dnorm_boost(X);
  ASSERT_TRUE(all_close(d1, d2, 1e-6, 1e-6));

  // tools_stats::pnorm is the same as pnorm_boost
  auto p1 = tools_stats::pnorm(X);
  auto p2 = pnorm_boost(X);
  ASSERT_TRUE(all_close(p1, p2, 1e-6, 1e-6));

  // tools_stats::qnorm is the same as qnorm_boost
  auto q1 = tools_stats::qnorm(p1);
  auto q2 = qnorm_boost(p1);
  ASSERT_TRUE(all_close(q1, q2, 1e-6, 1e-6));
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

TEST(test_tools_stats, find_latent_sample)
{
  Eigen::MatrixXd u(4, 4);
  u << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7;

  double bandwidth = 0.1;
  size_t niter = 10;

  Eigen::MatrixXd latent_sample =
    tools_stats::find_latent_sample(u, bandwidth, niter);

  EXPECT_EQ(latent_sample.rows(), u.rows());
  EXPECT_EQ(latent_sample.cols(), 2);

  u.resize(2, 8);
  EXPECT_THROW(tools_stats::find_latent_sample(u, bandwidth, niter),
               std::runtime_error);
}

// Golden-value regression guards for the performance work: these pin the
// exact output streams/values of deterministic numerical routines. The
// reference values were generated from the pre-optimization implementation;
// optimizations classified "bit-identical" must keep matching them.
TEST(test_tools_stats, golden_qrng_streams)
{
  Eigen::MatrixXd ghalton_exp(8, 5);
  // clang-format off
  ghalton_exp <<
    0.50134366401471198, 0.90017063378155804, 0.79033276262975316, 0.73629821272814733, 0.057817524712858233,
    0.0013436640147119761, 0.23350396711489133, 0.39033276262975314, 0.16486964129957585, 0.42145388834922187,
    0.75134366401471198, 0.56683730044822467, 0.99033276262975323, 0.59344106987100442, 0.78509025198558546,
    0.25134366401471198, 0.67794841155933572, 0.59033276262975309, 0.022012498442432991, 0.14872661562194914,
    0.62634366401471198, 0.011281744892669107, 0.19033276262975313, 0.45058392701386157, 0.51236297925831276,
    0.12634366401471198, 0.34461507822600246, 0.71033276262975309, 0.87915535558529012, 0.8759993428946764,
    0.87634366401471198, 0.78905952267044688, 0.31033276262975312, 0.30772678415671872, 0.23963570653104005,
    0.37634366401471198, 0.12239285600378021, 0.91033276262975316, 0.79752270252406565, 0.6032720701674037;
  // clang-format on
  EXPECT_TRUE(
    all_close(tools_stats::ghalton(8, 5, { 5 }), ghalton_exp, 1e-14, 1e-14));

  Eigen::MatrixXd sobol_exp(8, 5);
  // clang-format off
  sobol_exp <<
    0.8847676906734705, 0.88564104889519513, 0.78288133232854307, 0.82119966601021588, 0.067851732717826962,
    0.3847676906734705, 0.38564104889519513, 0.28288133232854307, 0.32119966601021588, 0.56785173271782696,
    0.1347676906734705, 0.63564104889519513, 0.53288133232854307, 0.57119966601021588, 0.81785173271782696,
    0.6347676906734705, 0.13564104889519513, 0.032881332328543067, 0.071199666010215878, 0.31785173271782696,
    0.5097676906734705, 0.51064104889519513, 0.40788133232854307, 0.19619966601021588, 0.44285173271782696,
    0.0097676906734704971, 0.010641048895195127, 0.90788133232854307, 0.69619966601021588, 0.94285173271782696,
    0.2597676906734705, 0.76064104889519513, 0.15788133232854307, 0.44619966601021588, 0.69285173271782696,
    0.7597676906734705, 0.26064104889519513, 0.65788133232854307, 0.94619966601021588, 0.19285173271782696;
  // clang-format on
  EXPECT_TRUE(
    all_close(tools_stats::sobol(8, 5, { 5 }), sobol_exp, 1e-14, 1e-14));

  Eigen::MatrixXd unif_exp(8, 3);
  // clang-format off
  unif_exp <<
    0.8847676906734705, 0.19036572566255927, 0.10514992498792708,
    0.88564104889519513, 0.60051596234552562, 0.79738970589824021,
    0.78288133232854307, 0.49468548316508532, 0.49710885854437947,
    0.82119966601021588, 0.29806191520765424, 0.97251652297563851,
    0.067851732717826962, 0.62834115885198116, 0.43421882297843695,
    0.10804224293678999, 0.0022314379457384348, 0.75123894470743835,
    0.97652392741292715, 0.94386291340924799, 0.82599193067289889,
    0.98571535851806402, 0.43902719020843506, 0.76561863790266216;
  // clang-format on
  EXPECT_TRUE(all_close(
    tools_stats::simulate_uniform(8, 3, false, { 5 }), unif_exp, 1e-14, 1e-14));

  Eigen::MatrixXd unif_qrng_exp(8, 3);
  // clang-format off
  unif_qrng_exp <<
    0.97411725879646838, 0.75434095744764029, 0.60310761896927689,
    0.47411725879646838, 0.087674290780973649, 0.20310761896927693,
    0.72411725879646838, 0.42100762411430698, 0.80310761896927685,
    0.22411725879646838, 0.86545206855875145, 0.40310761896927688,
    0.84911725879646838, 0.19878540189208474, 0.0031076189692769113,
    0.34911725879646838, 0.53211873522541808, 0.72310761896927689,
    0.59911725879646838, 0.9765631796698625, 0.32310761896927687,
    0.099117258796468377, 0.30989651300319582, 0.92310761896927696;
  // clang-format on
  EXPECT_TRUE(all_close(tools_stats::simulate_uniform(8, 3, true, { 5 }),
                        unif_qrng_exp,
                        1e-14,
                        1e-14));
}

TEST(test_tools_stats, golden_genz_kernels)
{
  Eigen::MatrixXd z(5, 2);
  // clang-format off
  z << -1.5, -0.5,
        0.0,  0.7,
        1.2, -2.1,
        2.5,  2.5,
       -0.3,  0.1;
  // clang-format on

  Eigen::VectorXd pbvnorm_05(5), pbvnorm_095(5), pbvt_4(5), pbvt_5(5);
  pbvnorm_05 << 0.046836527008394177, 0.44264120431058773, 0.01781391095522011,
    0.98825003409597401, 0.28443362319937548;
  pbvnorm_095 << 0.066790910076327439, 0.49943713401311252,
    0.017864420562816556, 0.99162723035221656, 0.37588082262768707;
  pbvt_4 << 0.071401026241030285, 0.43407568687068981, 0.048897595354819315,
    0.94393981082369438, 0.28708460935999386;
  pbvt_5 << 0.018620689838265002, 0.33309924773333083, 0.026605698240605276,
    0.94645806141370081, 0.16204495416733036;

  EXPECT_TRUE(
    all_close(tools_stats::pbvnorm(z, 0.5), pbvnorm_05, 1e-12, 1e-12));
  EXPECT_TRUE(
    all_close(tools_stats::pbvnorm(z, 0.95), pbvnorm_095, 1e-12, 1e-12));
  EXPECT_TRUE(all_close(tools_stats::pbvt(z, 4, 0.5), pbvt_4, 1e-12, 1e-12));
  EXPECT_TRUE(all_close(tools_stats::pbvt(z, 5, -0.3), pbvt_5, 1e-12, 1e-12));
}

TEST(test_tools_stats, golden_pseudo_obs)
{
  Eigen::MatrixXd x(10, 2);
  // clang-format off
  x << 0.1, 0.9,  0.3, 0.3,  0.3, 0.5,  0.7, 0.1,  0.2, 0.2,
       0.9, 0.6,  0.5, 0.6,  0.5, 0.8,  0.8, 0.4,  0.4, 0.7;
  // clang-format on
  Eigen::VectorXd w = Eigen::VectorXd::LinSpaced(10, 0.5, 1.5);

  Eigen::MatrixXd weighted_exp(10, 2);
  // clang-format off
  weighted_exp <<
    0.045454545454545456, 0.90909090909090906,
    0.2196969696969697, 0.21717171717171715,
    0.2196969696969697, 0.40909090909090912,
    0.68686868686868674, 0.075757575757575746,
    0.1313131313131313, 0.1616161616161616,
    0.90909090909090895, 0.55808080808080807,
    0.55303030303030309, 0.55808080808080807,
    0.55303030303030309, 0.86363636363636365,
    0.81313131313131304, 0.34343434343434343,
    0.3888888888888889, 0.7474747474747474;
  // clang-format on
  EXPECT_TRUE(all_close(
    tools_stats::to_pseudo_obs(x, "average", w), weighted_exp, 1e-14, 1e-14));

  Eigen::MatrixXd random_exp(10, 2);
  // clang-format off
  random_exp <<
    0.090909090909090912, 0.90909090909090906,
    0.27272727272727271, 0.27272727272727271,
    0.27272727272727271, 0.45454545454545453,
    0.72727272727272729, 0.090909090909090912,
    0.18181818181818182, 0.18181818181818182,
    0.90909090909090906, 0.54545454545454541,
    0.54545454545454541, 0.54545454545454541,
    0.54545454545454541, 0.81818181818181823,
    0.81818181818181823, 0.36363636363636365,
    0.45454545454545453, 0.72727272727272729;
  // clang-format on
  EXPECT_TRUE(all_close(
    tools_stats::to_pseudo_obs(x, "random", Eigen::VectorXd(), { 17 }),
    random_exp,
    1e-14,
    1e-14));
}
}
