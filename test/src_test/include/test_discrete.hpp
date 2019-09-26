// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/vinecop/class.hpp>

namespace test_discrete {

using namespace vinecopulib;

TEST(discrete, bicop)
{
  for (auto rot : { 0, 90, 180, 270 }) {
    auto bc = Bicop(BicopFamily::clayton, rot, Eigen::VectorXd::Constant(1, 3));
    auto u = bc.simulate(5000, true, { 1 });
    Eigen::MatrixXd u_new(u.rows(), 4);
    u_new.block(0, 0, u.rows(), 2) = u;
    u_new.block(0, 2, u.rows(), 2) = u;

    // c_c
    EXPECT_GE(bc.pdf(u).minCoeff(), 0);
    bc.fit(u);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);

    // d_c
    u_new.col(0) = (u.col(0).array() * 2).ceil() / 2;
    u_new.col(2) = (u.col(0).array() * 2).floor() / 2;
    bc.set_var_types({ "d", "c" });
    EXPECT_GE(bc.pdf(u_new).minCoeff(), 0);
    EXPECT_GE(bc.hfunc1(u_new).minCoeff(), 0);
    EXPECT_GE(bc.hfunc2(u_new).minCoeff(), 0);
    bc.fit(u_new);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);
    EXPECT_EQ(bc.cdf(u_new), bc.as_continuous().cdf(u_new.leftCols(2)));

    // c_d
    u_new.block(0, 0, u.rows(), 2) = u;
    u_new.block(0, 2, u.rows(), 2) = u;
    u_new.col(1) = (u.col(1).array() * 2).ceil() / 2;
    u_new.col(3) = (u.col(1).array() * 2).floor() / 2;
    bc.set_var_types({ "c", "d" });
    EXPECT_GE(bc.pdf(u_new).minCoeff(), 0);
    EXPECT_GE(bc.hfunc1(u_new).minCoeff(), 0);
    EXPECT_GE(bc.hfunc2(u_new).minCoeff(), 0);
    bc.fit(u_new);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);
    EXPECT_EQ(bc.cdf(u_new), bc.as_continuous().cdf(u_new.leftCols(2)));

    // d_d
    u_new.col(0) = (u_new.col(0).array() * 2).ceil() / 2.0;
    u_new.col(2) = (u_new.col(2).array() * 2).floor() / 2.0;
    bc.set_var_types({ "d", "d" });
    EXPECT_GE(bc.pdf(u_new).minCoeff(), 0);
    EXPECT_GE(bc.hfunc1(u_new).minCoeff(), 0);
    EXPECT_GE(bc.hfunc2(u_new).minCoeff(), 0);
    bc.fit(u_new);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);
    EXPECT_EQ(bc.cdf(u_new), bc.as_continuous().cdf(u_new.leftCols(2)));
    bc.select(u_new.topRows(20)); // all families
  }
}

TEST(discrete, vinecop)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(5);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 90, Eigen::VectorXd::Constant(1, 3));
    }
  }
  RVineStructure str(std::vector<size_t>{ 1, 2, 3, 4, 5 });
  Vinecop vc(pair_copulas, str);

  auto utmp = vc.simulate(5000, true, 1, {1});
  Eigen::MatrixXd u(utmp.rows(), 2 * utmp.cols());
  u.leftCols(5) = utmp;
  u.rightCols(5) = utmp;
  u.col(0) = (utmp.col(0).array() * 1000).ceil() / 1000;
  u.col(0 + utmp.cols()) = (utmp.col(0).array() * 1000).floor() / 1000;
  u.col(2) = (utmp.col(2).array() * 1000).ceil() / 1000;
  u.col(2 + utmp.cols()) = (utmp.col(2).array() * 1000).floor() / 1000;
  u.col(3) = (utmp.col(3).array() * 1000).ceil() / 1000;
  u.col(3 + utmp.cols()) = (utmp.col(3).array() * 1000).floor() / 1000;

  // fit vine
  vc.set_var_types({ "d", "c", "d", "d", "c" });
  auto controls = FitControlsVinecop({ BicopFamily::clayton });
  // controls.set_show_trace(true);
  vc.select_families(u, controls);

  // check output
  auto pcs = vc.get_all_pair_copulas();
  for (auto tree : pcs) {
    for (auto pc : tree) {
      EXPECT_EQ(pc.get_rotation(), 90);
      EXPECT_NEAR(pc.get_parameters()(0), 3, 1);
    }
  }
}
}
