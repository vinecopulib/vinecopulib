// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "include/test_utils.hpp"
#include "gtest/gtest.h"
#include <algorithm>
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

TEST(test_tools_stats, cxi_works)
{
  std::vector<int> seeds = { 1, 2, 3, 4, 5 };
  Eigen::MatrixXd Z = tools_stats::simulate_uniform(2000, 2, false, seeds);
  Z = tools_stats::qnorm(Z);

  // under independence, xi is centered at zero
  EXPECT_NEAR(tools_stats::pairwise_cxi(Z), 0.0, 0.1);

  // a deterministic but non-monotonic relationship: Kendall's tau is blind to
  // it, xi is not
  Eigen::MatrixXd V(Z.rows(), 2);
  V.col(0) = Z.col(0);
  V.col(1) = Z.col(0).cwiseAbs2();
  EXPECT_NEAR(std::fabs(wdm::wdm(V, "tau")(0, 1)), 0.0, 0.1);
  EXPECT_GT(tools_stats::pairwise_cxi(V), 0.9);

  // the two directions differ substantially here, so taking the larger of
  // them is what makes the criterion independent of the column order
  double xi12 = wdm::wdm(V.col(0), V.col(1), "cxi");
  double xi21 = wdm::wdm(V.col(1), V.col(0), "cxi");
  EXPECT_GT(xi12 - xi21, 0.1);
  Eigen::MatrixXd V_swapped(V.rows(), 2);
  V_swapped.col(0) = V.col(1);
  V_swapped.col(1) = V.col(0);
  EXPECT_DOUBLE_EQ(tools_stats::pairwise_cxi(V),
                   tools_stats::pairwise_cxi(V_swapped));

  // uniform weights must not change anything, zero weights must drop rows
  Eigen::VectorXd weights = Eigen::VectorXd::Ones(V.rows());
  EXPECT_NEAR(
    tools_stats::pairwise_cxi(V, weights), tools_stats::pairwise_cxi(V), 1e-10);
  weights.tail(V.rows() / 2) = Eigen::VectorXd::Zero(V.rows() / 2);
  EXPECT_NEAR(tools_stats::pairwise_cxi(V, weights),
              tools_stats::pairwise_cxi(V.topRows(V.rows() / 2)),
              0.05);

  // wdm breaks predictor ties at random but seeds that deterministically, so
  // the criterion must not vary between calls on the same data, nor between
  // the two column orders
  Eigen::MatrixXd tied(60, 2);
  for (Eigen::Index i = 0; i < tied.rows(); ++i) {
    tied(i, 0) = static_cast<double>(i / 12) / 5.0;
    tied(i, 1) = static_cast<double>(i % 12) / 11.0;
  }
  Eigen::MatrixXd tied_swapped(tied.rows(), 2);
  tied_swapped.col(0) = tied.col(1);
  tied_swapped.col(1) = tied.col(0);
  double tied_cxi = tools_stats::pairwise_cxi(tied);
  for (int rep = 0; rep < 5; ++rep) {
    EXPECT_DOUBLE_EQ(tools_stats::pairwise_cxi(tied), tied_cxi);
  }
  EXPECT_DOUBLE_EQ(tools_stats::pairwise_cxi(tied_swapped), tied_cxi);
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

TEST(test_tools_stats, degenerate_dimensions_throw)
{
  EXPECT_ANY_THROW(tools_stats::ghalton(0, 2));
  EXPECT_ANY_THROW(tools_stats::ghalton(2, 0));
  EXPECT_ANY_THROW(tools_stats::sobol(0, 2));
  EXPECT_ANY_THROW(tools_stats::sobol(2, 0));
  EXPECT_ANY_THROW(tools_stats::simulate_uniform(0, 2));
  EXPECT_ANY_THROW(tools_stats::simulate_uniform(0, 2, true));
}

TEST(test_tools_stats, dpqnorm_work)
{
  auto dnorm_boost = [](const Eigen::MatrixXd& x) {
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
  };

  auto pnorm_boost = [](const Eigen::MatrixXd& x) {
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
  };

  auto qnorm_boost = [](const Eigen::MatrixXd& x) {
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

// The normal CDF saturates at 0 and 1 far in the tails, including at
// infinite arguments, and passes NaN through.
TEST(test_tools_stats, pnorm_saturates_in_the_tails)
{
  const double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd X(6);
  X << inf, -inf, 1e300, -1e300, 50.0, -50.0;
  Eigen::VectorXd expected(6);
  expected << 1, 0, 1, 0, 1, 0;
  EXPECT_EQ(tools_stats::pnorm(X), expected);

  Eigen::VectorXd nan =
    Eigen::VectorXd::Constant(1, std::numeric_limits<double>::quiet_NaN());
  EXPECT_TRUE(std::isnan(tools_stats::pnorm(nan)(0)));
}

// A coordinate of exactly one and bounds beyond the unit square fall into the
// boundary cells instead of past the end of the covering.
TEST(test_tools_stats, box_covering_handles_the_boundary)
{
  Eigen::MatrixXd u(3, 2);
  u << 1.0, 1.0, 0.0, 0.0, 0.5, 1.0;
  tools_stats::BoxCovering covering(u, 4);

  Eigen::VectorXd lower(2), upper(2);
  lower << 0.9, 0.9;
  upper << 1.0, 1.0;
  EXPECT_EQ(covering.get_box_indices(lower, upper), std::vector<size_t>{ 0 });

  lower << 0.0, 0.0;
  upper << 1.5, 1.5;
  EXPECT_EQ(covering.get_box_indices(lower, upper).size(), 3u);

  covering.swap_sample(1, Eigen::Vector2d(1.0, 0.0));
  lower << 0.9, 0.0;
  upper << 1.0, 0.1;
  EXPECT_EQ(covering.get_box_indices(lower, upper), std::vector<size_t>{ 1 });
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

namespace {

//! Discretizes both variables to `levels` support points and returns the
//! four-column (u, u^-) layout `find_latent_sample` expects.
Eigen::MatrixXd
discretize_both(const Eigen::MatrixXd& u, double levels = 8.0)
{
  Eigen::MatrixXd out(u.rows(), 4);
  out.col(0) = (u.col(0).array() * levels).ceil() / levels;
  out.col(1) = (u.col(1).array() * levels).ceil() / levels;
  out.col(2) = (u.col(0).array() * levels).floor() / levels;
  out.col(3) = (u.col(1).array() * levels).floor() / levels;
  return out;
}

} // namespace

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

// The sweep counter has to be as wide as `niter`, or a count above its range
// wraps to zero and the loop never ends.
TEST(test_tools_stats, find_latent_sample_sweeps_more_than_a_short_can_count)
{
  Eigen::MatrixXd u(2, 4);
  u << 0.25, 0.75, 0.0, 0.5, 0.75, 0.25, 0.5, 0.0;

  const size_t niter = 70000;
  EXPECT_EQ(tools_stats::find_latent_sample(u, 0.1, niter).rows(), u.rows());
}

// Every random component is a fixed-seed quasi-random sequence, so the draw is
// reproducible.
TEST(test_tools_stats, find_latent_sample_is_deterministic)
{
  Bicop gauss(BicopFamily::gaussian, 0, Eigen::MatrixXd::Constant(1, 1, 0.7));
  const Eigen::MatrixXd u = discretize_both(gauss.simulate(300, false, { 5 }));
  EXPECT_EQ(tools_stats::find_latent_sample(u, 0.2),
            tools_stats::find_latent_sample(u, 0.2));
}

// The draws are indexed by column position, so without a canonical ordering
// the same pair of variables recovers two unrelated latent samples depending
// on which one is passed first -- and a vine reuses pair copulas in either
// orientation.
TEST(test_tools_stats, find_latent_sample_ignores_argument_order)
{
  Bicop gauss(BicopFamily::gaussian, 0, Eigen::MatrixXd::Constant(1, 1, 0.7));
  const Eigen::MatrixXd u = discretize_both(gauss.simulate(300, false, { 5 }));
  Eigen::MatrixXd u_swapped = u;
  u_swapped.col(0).swap(u_swapped.col(1));
  u_swapped.col(2).swap(u_swapped.col(3));

  Eigen::MatrixXd swapped_back =
    tools_stats::find_latent_sample(u_swapped, 0.2);
  swapped_back.col(0).swap(swapped_back.col(1));
  EXPECT_EQ(tools_stats::find_latent_sample(u, 0.2), swapped_back);
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
    0.50872809253633022, 0.88590955577006536, 0.17918644809253989, 0.87640347244944372, 0.23008731673957963,
    0.0087280925363302231, 0.21924288910339865, 0.7791864480925399, 0.30497490102087232, 0.59372368037594336,
    0.75872809253633022, 0.55257622243673199, 0.37918644809253987, 0.73354632959230082, 0.9573600440123069,
    0.25872809253633022, 0.99702066688117641, 0.97918644809253996, 0.16211775816372945, 0.32099640764867055,
    0.63372809253633022, 0.33035400021450972, 0.57918644809253994, 0.59068918673515802, 0.68463277128503419,
    0.13372809253633022, 0.66368733354784304, 0.099186448092539903, 0.019260615306586598, 0.048269134921397838,
    0.88372809253633022, 0.7747984446589542, 0.69918644809253983, 0.44783204387801517, 0.41190549855776148,
    0.38372809253633022, 0.10813177799228753, 0.29918644809253991, 0.93762796224536216, 0.77554186219412502;
  // clang-format on
  EXPECT_TRUE(
    all_close(tools_stats::ghalton(8, 5, { 5 }), ghalton_exp, 1e-14, 1e-14));

  Eigen::MatrixXd sobol_exp(8, 5);
  // clang-format off
  sobol_exp <<
    0.8847676906734705, 0.78288133232854307, 0.067851732717826962, 0.97652392741292715, 0.19036572566255927,
    0.3847676906734705, 0.28288133232854307, 0.56785173271782696, 0.47652392741292715, 0.69036572566255927,
    0.1347676906734705, 0.53288133232854307, 0.31785173271782696, 0.72652392741292715, 0.94036572566255927,
    0.6347676906734705, 0.032881332328543067, 0.81785173271782696, 0.22652392741292715, 0.44036572566255927,
    0.5097676906734705, 0.65788133232854307, 0.69285173271782696, 0.10152392741292715, 0.31536572566255927,
    0.0097676906734704971, 0.15788133232854307, 0.19285173271782696, 0.60152392741292715, 0.81536572566255927,
    0.2597676906734705, 0.90788133232854307, 0.94285173271782696, 0.35152392741292715, 0.56536572566255927,
    0.7597676906734705, 0.40788133232854307, 0.44285173271782696, 0.85152392741292715, 0.065365725662559271;
  // clang-format on
  EXPECT_TRUE(
    all_close(tools_stats::sobol(8, 5, { 5 }), sobol_exp, 1e-14, 1e-14));

  Eigen::MatrixXd unif_exp(8, 3);
  // clang-format off
  unif_exp <<
    0.88476769087967477, 0.10514992517358379, 0.97992718237261922,
    0.78288133251974346, 0.49710885877081101, 0.2111646659552131,
    0.067851732742982507, 0.43421882315334837, 0.69873678327599209,
    0.97652392764243179, 0.8259919308511583, 0.18496739747148594,
    0.19036572580237776, 0.48365020736264719, 0.13959733959910037,
    0.49468548323448325, 0.93642072859966208, 0.73217458711859484,
    0.62834115885250064, 0.86105880271030555, 0.84886235113604314,
    0.94386291351146689, 0.43112278247608815, 0.28077161281111263;
  // clang-format on
  EXPECT_TRUE(all_close(
    tools_stats::simulate_uniform(8, 3, false, { 5 }), unif_exp, 1e-14, 1e-14));

  Eigen::MatrixXd unif_qrng_exp(8, 3);
  // clang-format off
  unif_qrng_exp <<
    0.89059707918204367, 0.76444954341647342, 0.087747554984608137,
    0.39059707918204367, 0.097782876749806666, 0.68774755498460816,
    0.64059707918204367, 0.43111621008313999, 0.28774755498460813,
    0.14059707918204367, 0.87556065452758436, 0.88774755498460822,
    0.76559707918204367, 0.20889398786091776, 0.48774755498460809,
    0.26559707918204367, 0.5422273211942511, 0.0077475549846081366,
    0.51559707918204367, 0.98667176563869552, 0.60774755498460808,
    0.015597079182043672, 0.32000509897202889, 0.20774755498460812;
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
    0.36363636363636365, 0.45454545454545453,
    0.72727272727272729, 0.090909090909090912,
    0.18181818181818182, 0.18181818181818182,
    0.90909090909090906, 0.63636363636363635,
    0.63636363636363635, 0.54545454545454541,
    0.54545454545454541, 0.81818181818181823,
    0.81818181818181823, 0.36363636363636365,
    0.45454545454545453, 0.72727272727272729;
  // clang-format on
  auto random_obs =
    tools_stats::to_pseudo_obs(x, "random", Eigen::VectorXd(), { 17 });
  EXPECT_TRUE(all_close(random_obs, random_exp, 1e-14, 1e-14));

  // `"random"` breaks ties, so each column is a permutation of
  // 1 / (n + 1), ..., n / (n + 1); golden values alone would not catch a tie
  // group collapsing onto one rank.
  double n = static_cast<double>(x.rows());
  for (Eigen::Index j = 0; j < random_obs.cols(); ++j) {
    Eigen::VectorXd sorted = random_obs.col(j);
    std::sort(sorted.data(), sorted.data() + sorted.size());
    for (Eigen::Index i = 0; i < sorted.size(); ++i) {
      EXPECT_NEAR(sorted(i), (static_cast<double>(i) + 1.0) / (n + 1.0), 1e-14);
    }
  }
}

// the maximal correlation is symmetric in its two variables
TEST(tools_stats, pairwise_mcor_is_symmetric)
{
  for (int seed = 0; seed < 20; ++seed) {
    Eigen::MatrixXd u = tools_stats::simulate_uniform(1000, 2, false, { seed });
    Eigen::MatrixXd swapped(u.rows(), 2);
    swapped << u.col(1), u.col(0);
    EXPECT_EQ(tools_stats::pairwise_mcor(u),
              tools_stats::pairwise_mcor(swapped))
      << "seed = " << seed;
  }
}
}
