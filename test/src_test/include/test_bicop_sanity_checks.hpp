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

TEST(bicop_sanity_checks, controls_checks)
{
  auto controls = FitControlsBicop();
  EXPECT_ANY_THROW(controls.set_selection_criterion("foo"));
  EXPECT_ANY_THROW(controls.set_nonparametric_method("foo"));
  EXPECT_ANY_THROW(controls.set_parametric_method("foo"));
  EXPECT_ANY_THROW(controls.set_nonparametric_mult(0.0));
  EXPECT_ANY_THROW(controls.set_nonparametric_grid_size(2));
  EXPECT_ANY_THROW(controls.set_psi0(0.0));
  EXPECT_ANY_THROW(controls.set_psi0(1.0));
}

TEST(bicop_sanity_checks, fit_controls_config_works)
{
  // Some non-default controls for testing
  FitControlsBicop controls;
  controls.set_family_set(bicop_families::itau);
  controls.set_parametric_method("itau");
  controls.set_nonparametric_method("quadratic");
  controls.set_nonparametric_mult(2.0);
  controls.set_nonparametric_grid_size(100);
  controls.set_selection_criterion("bic");
  controls.set_weights(Eigen::VectorXd::Ones(10));
  controls.set_psi0(0.6);
  controls.set_preselect_families(false);
  controls.set_allow_rotations(false);
  // can't use non-default num_threads in CI

  // Create a config object from the controls
  FitControlsConfig config;
  config.family_set = controls.get_family_set();
  config.parametric_method = controls.get_parametric_method();
  config.nonparametric_method = controls.get_nonparametric_method();
  config.nonparametric_mult = controls.get_nonparametric_mult();
  config.nonparametric_grid_size = controls.get_nonparametric_grid_size();
  config.selection_criterion = controls.get_selection_criterion();
  config.weights = controls.get_weights();
  config.psi0 = controls.get_psi0();
  config.preselect_families = controls.get_preselect_families();
  config.allow_rotations = controls.get_allow_rotations();
  config.num_threads = controls.get_num_threads();

  // Create and test new controls from the config object
  FitControlsBicop controls2(config);
  EXPECT_EQ(controls.get_family_set(), controls2.get_family_set());
  EXPECT_EQ(controls.get_parametric_method(),
            controls2.get_parametric_method());
  EXPECT_EQ(controls.get_nonparametric_method(),
            controls2.get_nonparametric_method());
  EXPECT_EQ(controls.get_nonparametric_mult(),
            controls2.get_nonparametric_mult());
  EXPECT_EQ(controls.get_nonparametric_grid_size(),
            controls2.get_nonparametric_grid_size());
  EXPECT_EQ(controls.get_selection_criterion(),
            controls2.get_selection_criterion());
  EXPECT_EQ(controls.get_weights(), controls2.get_weights());
  EXPECT_EQ(controls.get_psi0(), controls2.get_psi0());
  EXPECT_EQ(controls.get_preselect_families(),
            controls2.get_preselect_families());
  EXPECT_EQ(controls.get_allow_rotations(), controls2.get_allow_rotations());
  EXPECT_EQ(controls.get_num_threads(), controls2.get_num_threads());
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