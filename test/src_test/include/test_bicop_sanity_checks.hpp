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

//! Pins the documented contract of the two tau maps: `parameters_to_tau` is
//! available for every family, `tau_to_parameters` only where tau determines
//! the parameters completely.
TEST(bicop_sanity_checks, tau_maps_are_available_where_documented)
{
  const std::vector<BicopFamily> invertible = {
    BicopFamily::indep,  BicopFamily::gaussian, BicopFamily::clayton,
    BicopFamily::gumbel, BicopFamily::frank,    BicopFamily::joe
  };
  for (auto family : bicop_families::all) {
    Bicop bicop(family);
    const bool expected = tools_stl::is_member(family, invertible);

    // parameters_to_tau: every family
    EXPECT_NO_THROW(bicop.parameters_to_tau(bicop.get_parameters()))
      << bicop.get_family_name();

    // tau_to_parameters: only the families listed above. Student is in
    // bicop_families::itau yet is not invertible here, because tau pins its
    // correlation and leaves the degrees of freedom free.
    if (expected) {
      EXPECT_NO_THROW(bicop.tau_to_parameters(0.4)) << bicop.get_family_name();
    } else {
      EXPECT_THROW(bicop.tau_to_parameters(0.4), std::runtime_error)
        << bicop.get_family_name();
    }
  }
}

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

TEST(bicop_sanity_checks, catches_transposed_parameter_shape)
{
  // Regression test for PR #700: a same-size but transposed shape (1 x 2
  // instead of the expected 2 x 1) must be rejected up front. It used to slip
  // past the size check (equal total size) and reach the coefficient-wise
  // bound comparisons with mismatched dimensions (an out-of-bounds read).
  auto par = (Eigen::MatrixXd(2, 1) << 0.5, 4.0).finished();
  auto bc = Bicop(BicopFamily::student, 0, par);
  Eigen::MatrixXd transposed = bc.get_parameters().transpose(); // 1 x 2
  ASSERT_EQ(transposed.size(), bc.get_parameters().size());
  EXPECT_ANY_THROW(bc.set_parameters(transposed));
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
  auto u = tools_stats::simulate_uniform(10, 3, false, { 1 });
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
  auto u = bicop.simulate(10, false, { 1 });
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

  auto u = bc1.simulate(10, false, { 1 });
  bc2.select(u);
  auto bc3 = bc2;
  EXPECT_EQ(bc2.get_loglik(), bc3.get_loglik());
  EXPECT_EQ(bc2.get_nobs(), bc3.get_nobs());
}

// Kendall's tau of the families whose tau is a numerical integral. The
// expectations are 50-digit references (Gauss-Kronrod on the same integrand),
// and each parameter set is one where the previous implementation lost
// precision or failed outright.
TEST(bicop_sanity_checks, tau_of_integral_families_is_accurate)
{
  struct Case
  {
    BicopFamily family;
    std::vector<double> parameters;
    double tau;
  };

  // References computed at 50 digits. The corners matter: each of these was
  // wrong by between 5e-4 and 100% before the tau integrands were reformulated.
  const std::vector<Case> cases = {
    // large theta: A'' is a narrow peak
    { BicopFamily::tawn, { 0.5, 0.5, 50.0 }, 0.33100293306647305 },
    { BicopFamily::tawn, { 0.5, 0.5, 60.0 }, 0.33140637575860532 },
    { BicopFamily::tawn, { 0.9, 0.9, 60.0 }, 0.80695186456976831 },
    { BicopFamily::tawn, { 0.999, 0.999, 60.0 }, 0.98140077647182433 },
    { BicopFamily::tawn, { 0.5, 0.5, 2.0 }, 0.21460183660255169 },
    // small psi1: the peak moves to psi1 / (psi1 + psi2), far from 1/2
    { BicopFamily::tawn, { 1e-4, 0.5, 60.0 }, 9.9989819367449941e-05 },
    { BicopFamily::tawn, { 1e-3, 0.5, 10.0 }, 0.00099884000269966056 },
    { BicopFamily::tawn, { 0.01, 0.5, 60.0 }, 0.0098992116497294558 },
    // BB7/BB8 near the upper limits, where the integrand used to cancel to
    // exactly zero
    { BicopFamily::bb7, { 6.0, 0.01 }, 0.72294689117590216 },
    { BicopFamily::bb7, { 1.5, 1.0 }, 0.42857142857142860 },
    { BicopFamily::bb8, { 8.0, 1.0 }, 0.78325404384180519 },
    { BicopFamily::bb8, { 3.0, 0.8 }, 0.34731888449523796 },
  };

  for (const auto& c : cases) {
    Eigen::VectorXd par(c.parameters.size());
    for (size_t i = 0; i < c.parameters.size(); ++i)
      par(i) = c.parameters[i];
    Bicop bicop(c.family, 0, par);
    EXPECT_NEAR(bicop.get_tau(), c.tau, 1e-9)
      << "family " << bicop.get_family_name() << ", parameters "
      << par.transpose();
  }
}

// tau is monotone in Tawn's theta, which the pre-fix implementation violated by
// collapsing to zero above theta = 50.
TEST(bicop_sanity_checks, tawn_tau_is_monotone_in_theta)
{
  double previous = -1.0;
  for (double theta : { 1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0 }) {
    Eigen::VectorXd par(3);
    par << 0.9, 0.9, theta;
    const double tau = Bicop(BicopFamily::tawn, 0, par).get_tau();
    EXPECT_GT(tau, previous) << "theta = " << theta;
    previous = tau;
  }
}
}
