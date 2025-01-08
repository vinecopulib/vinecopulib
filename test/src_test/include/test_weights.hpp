// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib.hpp>
#include <vinecopulib.hpp>

namespace test_weights {
using namespace vinecopulib;

TEST(test_weights, catches_incompatible_sizes)
{
  auto u = tools_stats::simulate_uniform(20, 2);
  auto w = tools_stats::simulate_uniform(9, 1);
  FitControlsBicop controls;
  controls.set_weights(w);
  EXPECT_ANY_THROW(Bicop(u, controls));
  u = tools_stats::simulate_uniform(20, 4);
  EXPECT_ANY_THROW(
    Vinecop(u, RVineStructure(), {}, FitControlsVinecop(controls)));
}

TEST(test_weights, allows_nans)
{
  auto u = tools_stats::simulate_uniform(20, 2);
  auto w = tools_stats::simulate_uniform(20, 1);
  w(1) = std::numeric_limits<double>::quiet_NaN();
  FitControlsBicop controls;
  controls.set_weights(w);
  EXPECT_NO_THROW(Bicop(u, controls));
}

TEST(test_weights, works_in_bicop_select)
{
  // if half of the weights are zero
  auto u = tools_stats::simulate_uniform(200, 2);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(200);
  w.head(100) = Eigen::VectorXd::Ones(100);

  FitControlsBicop controls(bicop_families::parametric);
  auto cop_uw = Bicop(u.block(0, 0, 100, 2), controls);
  controls.set_weights(w);
  auto cop_w = Bicop(u, controls);
  EXPECT_EQ(cop_uw.get_family(), cop_w.get_family());
  EXPECT_EQ(cop_uw.get_parameters(), cop_w.get_parameters());
}

TEST(test_weights, works_in_vinecop_select)
{
  // if half of the weights are zero
  auto u = tools_stats::simulate_uniform(200, 7);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(200);
  w.head(100) = Eigen::VectorXd::Ones(100);

  FitControlsBicop controls(bicop_families::parametric);
  auto cop_uw = Vinecop(
    u.block(0, 0, 100, 7), RVineStructure(), {}, FitControlsVinecop(controls));
  controls.set_weights(w);
  auto cop_w = Vinecop(u, RVineStructure(), {}, FitControlsVinecop(controls));
  EXPECT_EQ(cop_uw.get_all_families(), cop_w.get_all_families());
  EXPECT_EQ(cop_uw.get_all_parameters(), cop_w.get_all_parameters());
}
}
