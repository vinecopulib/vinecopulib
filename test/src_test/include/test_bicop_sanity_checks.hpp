// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib.hpp>

namespace test_bicop_sanity_checks {
using namespace vinecopulib;

TEST(bicop_sanity_checks, catches_wrong_parameter_size)
{
  for (auto family : bicop_families::nonparametric) {
    EXPECT_ANY_THROW(Bicop(family, 0, Eigen::VectorXd::Zero(1)));
  }
  for (auto family : bicop_families::one_par) {
    EXPECT_ANY_THROW(Bicop(family, 0, Eigen::VectorXd::Zero(2)));
  }
  for (auto family : bicop_families::two_par) {
    EXPECT_ANY_THROW(Bicop(family, 0, Eigen::VectorXd::Zero(1)));
  }
}

TEST(bicop_sanity_checks, catches_parameters_out_of_bounds)
{
  auto cop = Bicop(BicopFamily::gaussian);
  auto wrong_par = Eigen::VectorXd::Constant(1, 1.01);
  EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, 0, wrong_par));
  EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, 0, -wrong_par));
}

TEST(bicop_sanity_checks, catches_wrong_rotation)
{
  EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, -10));
  EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, 10));
}

TEST(bicop_sanity_checks, catches_var_types)
{
  auto rho = Eigen::VectorXd::Constant(1, 0.5);
  EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, 0, rho, { "c" }));
  EXPECT_ANY_THROW(Bicop(BicopFamily::gaussian, 0, rho, { "c", "u" }));
}

TEST(bicop_sanity_checks, catches_data_dim)
{
  Bicop bicop;
  auto u = tools_stats::simulate_uniform(10, 3);
  EXPECT_ANY_THROW(bicop.select(u));
  bicop.set_var_types({ "d", "d" });
  EXPECT_ANY_THROW(bicop.select(u));
}

TEST(bicop_sanity_checks, catches_not_fitted_to_data)
{
  auto bc = Bicop(BicopFamily::gaussian);
  EXPECT_ANY_THROW(bc.get_loglik());
  EXPECT_ANY_THROW(bc.get_nobs());
  EXPECT_ANY_THROW(bc.get_aic());
  EXPECT_ANY_THROW(bc.get_bic());
  EXPECT_ANY_THROW(bc.get_mbic(0.6));
}

TEST(bicop_sanity_checks, select_can_handle_zeros_and_ones)
{
  Bicop bicop;
  auto u = bicop.simulate(10);
  u(0, 0) = 0.0;
  u(1, 0) = 1.0;
  EXPECT_NO_THROW(bicop.select(u));
}

TEST(bicop_sanity_checks, controls_print)
{
  auto controls = FitControlsBicop();
  EXPECT_NO_THROW(controls.str());
}

TEST(bicop_sanity_checks, copy)
{
  auto rho = Eigen::VectorXd::Constant(1, 0.5);
  Bicop bc1(BicopFamily::gaussian, 0, rho);
  Bicop bc2 = bc1;
  bc2.set_parameters(rho.array() + 0.2);
  EXPECT_EQ(bc1.get_parameters(), rho);
  EXPECT_ANY_THROW(bc1.get_loglik());

  auto u = bc1.simulate(10);
  bc2.select(u);
  auto bc3 = bc2;
  EXPECT_EQ(bc2.get_loglik(), bc3.get_loglik());
  EXPECT_EQ(bc2.get_nobs(), bc3.get_nobs());
}
}