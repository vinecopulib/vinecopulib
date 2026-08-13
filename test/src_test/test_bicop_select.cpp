// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include <limits>
#include <vinecopulib.hpp>

namespace test_bicop_select {
using namespace vinecopulib;

TEST(bicop_select, works_in_parallel)
{
  Bicop cop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, -0.5));
  auto u = cop.simulate(15, false, { 1 });
  Bicop fit1, fit2;
  fit1.select(u);
  FitControlsBicop controls;
  controls.set_num_threads(2);
  fit2.select(u, controls);
  EXPECT_EQ(fit1.get_family(), fit2.get_family());
  EXPECT_EQ(fit1.get_parameters(), fit2.get_parameters());
}

TEST(bicop_select, allows_all_selcrits)
{
  Bicop cop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, -0.5));
  auto u = cop.simulate(15, false, { 1 });
  FitControlsBicop controls;
  controls.set_selection_criterion("loglik");
  cop.select(u, controls);
  controls.set_selection_criterion("aic");
  cop.select(u, controls);
  controls.set_selection_criterion("bic");
  cop.select(u, controls);
  controls.set_selection_criterion("mbic");
  cop.select(u, controls);
  controls.set_selection_criterion("mbicv");
  cop.select(u, controls);
}

TEST(bicop_select, allow_rotations_works)
{
  Bicop cop(BicopFamily::clayton, 180, Eigen::VectorXd::Constant(1, 5));
  auto u = cop.simulate(static_cast<size_t>(5e3), false, { 1 });
  FitControlsBicop controls;
  controls.set_family_set({ BicopFamily::clayton });
  cop.select(u, controls);
  EXPECT_EQ(cop.get_rotation(), 180);
  controls.set_allow_rotations(false);
  cop.select(u, controls);
  EXPECT_EQ(cop.get_rotation(), 0);
}

TEST(bicop_select, fit_stats_are_correct)
{
  Bicop cop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, -0.5));
  auto u = cop.simulate(15, false, { 1 });
  cop.select(u);
  EXPECT_EQ(cop.get_nobs(), 15);
  EXPECT_NEAR(cop.get_loglik(), cop.loglik(u), 1e-10);
  EXPECT_NEAR(cop.get_loglik(), cop.loglik(), 1e-10);
  EXPECT_NEAR(cop.get_aic(), cop.aic(u), 1e-10);
  EXPECT_NEAR(cop.get_aic(), cop.aic(), 1e-10);
  EXPECT_NEAR(cop.get_bic(), cop.bic(u), 1e-10);
  EXPECT_NEAR(cop.get_bic(), cop.bic(), 1e-10);
  EXPECT_NEAR(cop.get_mbic(), cop.mbic(u), 1e-10);
  EXPECT_NEAR(cop.get_mbic(), cop.mbic(), 1e-10);
}

TEST(bicop_select, mbic_and_bic_agree_on_sample_size_with_nans)
{
  Bicop cop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, -0.5));
  auto u = cop.simulate(30, false, { 1 });
  cop.select(u);

  // Append rows containing NaN; both criteria must count complete cases only,
  // so appending incomplete observations must not change either penalty.
  Eigen::MatrixXd u_nan(u.rows() + 3, 2);
  u_nan << u,
    Eigen::MatrixXd::Constant(3, 2, std::numeric_limits<double>::quiet_NaN());

  EXPECT_NEAR(cop.bic(u_nan), cop.bic(u), 1e-10);
  EXPECT_NEAR(cop.mbic(u_nan), cop.mbic(u), 1e-10);
}
}
