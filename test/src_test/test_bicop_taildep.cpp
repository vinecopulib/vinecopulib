// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include <cmath>
#include <vinecopulib.hpp>

namespace test_bicop_taildep {
using namespace vinecopulib;

// convenience: a one-parameter parameter vector
inline Eigen::VectorXd
par1(double x)
{
  return Eigen::VectorXd::Constant(1, x);
}

TEST(bicop_taildep, independence_has_no_dependence)
{
  auto cop = Bicop(BicopFamily::indep);
  Eigen::MatrixXd td = cop.get_taildep();
  EXPECT_EQ(td.rows(), 2);
  EXPECT_EQ(td.cols(), 2);
  EXPECT_TRUE(td.isZero());
  EXPECT_NEAR(cop.get_beta(), 0.0, 1e-10);
}

TEST(bicop_taildep, gaussian_and_frank_have_no_taildep)
{
  auto gauss = Bicop(BicopFamily::gaussian, 0, par1(0.5));
  EXPECT_TRUE(gauss.get_taildep().isZero());
  auto frank = Bicop(BicopFamily::frank, 0, par1(5.0));
  EXPECT_TRUE(frank.get_taildep().isZero());
}

// The lower tail dependence of the base Clayton copula moves to a different
// corner under each rotation. This exercises the 2x2 rotation logic.
TEST(bicop_taildep, clayton_rotations_move_the_corner)
{
  double theta = 2.0;
  double lambda = std::pow(2.0, -1.0 / theta); // ~ 0.7071

  // rotation 0: lower-left corner
  auto td0 = Bicop(BicopFamily::clayton, 0, par1(theta)).get_taildep();
  EXPECT_NEAR(td0(0, 0), lambda, 1e-12);
  EXPECT_NEAR(td0(0, 1) + td0(1, 0) + td0(1, 1), 0.0, 1e-12);

  // rotation 90: lower-right corner (U1 -> 1, U2 -> 0)
  auto td90 = Bicop(BicopFamily::clayton, 90, par1(theta)).get_taildep();
  EXPECT_NEAR(td90(1, 0), lambda, 1e-12);
  EXPECT_NEAR(td90(0, 0), 0.0, 1e-12);
  EXPECT_NEAR(td90(1, 1), 0.0, 1e-12);

  // rotation 180: upper-right corner
  auto td180 = Bicop(BicopFamily::clayton, 180, par1(theta)).get_taildep();
  EXPECT_NEAR(td180(1, 1), lambda, 1e-12);
  EXPECT_NEAR(td180(0, 0), 0.0, 1e-12);

  // rotation 270: upper-left corner (U1 -> 0, U2 -> 1)
  auto td270 = Bicop(BicopFamily::clayton, 270, par1(theta)).get_taildep();
  EXPECT_NEAR(td270(0, 1), lambda, 1e-12);
  EXPECT_NEAR(td270(0, 0), 0.0, 1e-12);
  EXPECT_NEAR(td270(1, 1), 0.0, 1e-12);
}

// The Gumbel copula has upper tail dependence; check the 180 rotation swaps it
// to the lower-left corner.
TEST(bicop_taildep, gumbel_upper_taildep_swaps_under_180)
{
  double theta = 2.0;
  double lambda = 2 - std::pow(2.0, 1.0 / theta);
  EXPECT_NEAR(Bicop(BicopFamily::gumbel, 0, par1(theta)).get_taildep()(1, 1),
              lambda,
              1e-12);
  EXPECT_NEAR(Bicop(BicopFamily::gumbel, 180, par1(theta)).get_taildep()(0, 0),
              lambda,
              1e-12);
}

// The Student t copula is the only family with dependence in all four corners;
// it is radially symmetric so the matrix is symmetric with equal diagonal and
// equal off-diagonal entries.
TEST(bicop_taildep, student_has_four_corners)
{
  Eigen::VectorXd par(2);
  par << 0.5, 4.0; // rho = 0.5, nu = 4

  auto td = Bicop(BicopFamily::student, 0, par).get_taildep();
  // symmetry
  EXPECT_NEAR(td(0, 0), td(1, 1), 1e-12); // concordant corners equal
  EXPECT_NEAR(td(0, 1), td(1, 0), 1e-12); // discordant corners equal
  // all four positive but < 1
  EXPECT_GT(td.minCoeff(), 0.0);
  EXPECT_LT(td.maxCoeff(), 1.0);
  // positive correlation => concordant dependence stronger than discordant
  EXPECT_GT(td(0, 0), td(0, 1));

  // negative correlation => discordant dependence stronger than concordant
  par(0) = -0.5;
  auto td_neg = Bicop(BicopFamily::student, 0, par).get_taildep();
  EXPECT_GT(td_neg(0, 1), td_neg(0, 0));
  // magnitudes mirror the positive-correlation case
  EXPECT_NEAR(td_neg(0, 1), td(0, 0), 1e-12);
  EXPECT_NEAR(td_neg(0, 0), td(0, 1), 1e-12);
}

// Blomqvist's beta flips sign under 90/270 rotations, like Kendall's tau.
TEST(bicop_taildep, beta_sign_flips_under_rotation)
{
  double theta = 2.0;
  double beta0 = Bicop(BicopFamily::clayton, 0, par1(theta)).get_beta();
  EXPECT_GT(beta0, 0.0);
  EXPECT_NEAR(
    Bicop(BicopFamily::clayton, 90, par1(theta)).get_beta(), -beta0, 1e-12);
  EXPECT_NEAR(
    Bicop(BicopFamily::clayton, 180, par1(theta)).get_beta(), beta0, 1e-12);
  EXPECT_NEAR(
    Bicop(BicopFamily::clayton, 270, par1(theta)).get_beta(), -beta0, 1e-12);
}

// The getters must agree with the parameter-based methods.
TEST(bicop_taildep, getters_match_parameter_methods)
{
  auto cop = Bicop(BicopFamily::clayton, 90, par1(2.0));
  auto pars = cop.get_parameters();
  EXPECT_TRUE(cop.get_taildep().isApprox(cop.parameters_to_taildep(pars)));
  EXPECT_NEAR(cop.get_beta(), cop.parameters_to_beta(pars), 1e-12);
}

// The nonparametric kernel estimator reports NaN tail dependence but a finite
// Blomqvist's beta (computable from the cdf).
TEST(bicop_taildep, tll_taildep_is_nan_beta_is_finite)
{
  auto u =
    Bicop(BicopFamily::clayton, 0, par1(2.0)).simulate(1000, false, { 1 });
  Bicop tll(BicopFamily::tll);
  tll.fit(u);
  auto td = tll.get_taildep();
  EXPECT_TRUE(td.array().isNaN().all());
  EXPECT_FALSE(std::isnan(tll.get_beta()));
}
}
