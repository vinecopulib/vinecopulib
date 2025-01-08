// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib.hpp>

namespace test_bicop_select {
using namespace vinecopulib;

TEST(bicop_select, works_in_parallel)
{
  Bicop cop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, -0.5));
  auto u = cop.simulate(15);
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
  auto u = cop.simulate(15);
  FitControlsBicop controls;
  controls.set_selection_criterion("loglik");
  cop.select(u, controls);
  controls.set_selection_criterion("aic");
  cop.select(u, controls);
  controls.set_selection_criterion("bic");
  cop.select(u, controls);
  controls.set_selection_criterion("mbic");
  cop.select(u, controls);
}

TEST(bicop_select, fit_stats_are_correct)
{
  Bicop cop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, -0.5));
  auto u = cop.simulate(15);
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
}
