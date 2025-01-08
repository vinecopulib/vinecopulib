// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib.hpp>
#include <wdm/eigen.hpp>

namespace test_discrete {

using namespace vinecopulib;

TEST(discrete, bicop)
{
  for (auto rot : { 0, 90, 180, 270 }) {
    auto bc = Bicop(BicopFamily::clayton, rot, Eigen::VectorXd::Constant(1, 3));
    auto u = bc.simulate(1000, true, { 1 });

    Eigen::MatrixXd u_disc(u.rows(), 4);
    u_disc.col(0) = (u.col(0).array() * 5).ceil() / 5;
    u_disc.col(2) = (u.col(0).array() * 5).floor() / 5;
    u_disc.col(1) = (u.col(1).array() * 5).ceil() / 5;
    u_disc.col(3) = (u.col(1).array() * 5).floor() / 5;

    // c_c
    EXPECT_GE(bc.pdf(u.topRows(20)).minCoeff(), 0);
    bc.fit(u);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);

    // d_c
    Eigen::MatrixXd uu(u.rows(), 3);
    uu.col(0) = u_disc.col(0);
    uu.col(1) = u.col(1);
    uu.col(2) = u_disc.col(2);
    bc.set_var_types({ "d", "c" });
    EXPECT_GE(bc.pdf(uu.topRows(20)).minCoeff(), 0);
    bc.fit(uu);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);
    EXPECT_EQ(bc.cdf(uu.topRows(20)),
              bc.as_continuous().cdf(uu.leftCols(2).topRows(20)));

    // c_d
    uu = Eigen::MatrixXd(u.rows(), 4);
    uu.col(0) = u.col(0);
    uu.col(2) = u.col(0);
    uu.col(1) = u_disc.col(1);
    uu.col(3) = u_disc.col(3);
    bc.set_var_types({ "c", "d" });
    EXPECT_GE(bc.pdf(uu.topRows(20)).minCoeff(), 0);
    bc.fit(uu);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);
    EXPECT_EQ(bc.cdf(uu.topRows(20)),
              bc.as_continuous().cdf(uu.leftCols(2).topRows(20)));

    // d_d
    uu = u_disc;
    bc.set_var_types({ "d", "d" });
    EXPECT_GE(bc.pdf(uu.topRows(20)).minCoeff(), 0);
    bc.fit(uu);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);
    EXPECT_EQ(bc.cdf(uu.topRows(20)),
              bc.as_continuous().cdf(uu.leftCols(2).topRows(20)));
    bc.select(uu.topRows(20)); // all families

    // tll
    bc.select(uu.topRows(20), FitControlsBicop({ BicopFamily::tll }));
    bc.parameters_to_tau(bc.get_parameters());
  }
}

TEST(zero_inflated, bicop)
{

  for (auto rot : { 0, 90, 180, 270 }) {
    auto bc = Bicop(BicopFamily::clayton, rot, Eigen::VectorXd::Constant(1, 3));
    auto tau = bc.parameters_to_tau(Eigen::VectorXd::Constant(1, 3));
    auto u = bc.simulate(1000, true, { 1 });

    auto thresh = Eigen::VectorXd::Constant(u.rows(), 0.1);
    auto zero = Eigen::VectorXd::Zero(u.rows());

    Eigen::MatrixXd u_disc(u.rows(), 4);
    u_disc.col(0) = u.col(0).cwiseMax(thresh);
    u_disc.col(2) = (u.col(0).array() < 0.1).select(zero, u.col(0));
    u_disc.col(1) = u.col(1).cwiseMax(thresh);
    u_disc.col(3) = (u.col(1).array() < 0.1).select(zero, u.col(1));

    // c_c
    EXPECT_GE(bc.pdf(u.topRows(20)).minCoeff(), 0);
    bc.fit(u);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);

    // d_c
    Eigen::MatrixXd uu(u.rows(), 3);
    uu.col(0) = u_disc.col(0);
    uu.col(1) = u.col(1);
    uu.col(2) = u_disc.col(2);
    bc = Bicop(BicopFamily::clayton, rot, Eigen::VectorXd::Constant(1, 3));
    bc.set_var_types({ "d", "c" });
    EXPECT_GE(bc.pdf(uu.topRows(20)).minCoeff(), 0);
    bc.fit(uu);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);
    EXPECT_EQ(bc.cdf(uu.topRows(20)),
              bc.as_continuous().cdf(uu.leftCols(2).topRows(20)));
    // tll
    bc.select(uu.topRows(20), FitControlsBicop({ BicopFamily::tll }));
    EXPECT_NEAR(bc.parameters_to_tau(bc.get_parameters()), tau, 0.15);

    // c_d
    uu = Eigen::MatrixXd(u.rows(), 4);
    uu.col(0) = u.col(0);
    uu.col(2) = u.col(0);
    uu.col(1) = u_disc.col(1);
    uu.col(3) = u_disc.col(3);
    bc = Bicop(BicopFamily::clayton, rot, Eigen::VectorXd::Constant(1, 3));
    bc.set_var_types({ "c", "d" });
    EXPECT_GE(bc.pdf(uu.topRows(20)).minCoeff(), 0);
    bc.fit(uu);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);
    EXPECT_EQ(bc.cdf(uu.topRows(20)),
              bc.as_continuous().cdf(uu.leftCols(2).topRows(20)));
    // tll
    bc.select(uu.topRows(20), FitControlsBicop({ BicopFamily::tll }));
    EXPECT_NEAR(bc.parameters_to_tau(bc.get_parameters()), tau, 0.15);

    // d_d
    uu = u_disc;
    bc = Bicop(BicopFamily::clayton, rot, Eigen::VectorXd::Constant(1, 3));
    bc.set_var_types({ "d", "d" });
    EXPECT_GE(bc.pdf(uu.topRows(20)).minCoeff(), 0);
    bc.fit(uu);
    EXPECT_NEAR(bc.get_parameters()(0), 3, 0.5);
    EXPECT_EQ(bc.cdf(uu.topRows(20)),
              bc.as_continuous().cdf(uu.leftCols(2).topRows(20)));
    bc.select(uu.topRows(20)); // all families

    // tll
    bc.select(uu.topRows(20), FitControlsBicop({ BicopFamily::tll }));
    EXPECT_NEAR(bc.parameters_to_tau(bc.get_parameters()), tau, 0.15);
  }
}

TEST(discrete, vinecop)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(5);
  for (size_t t = 0; t < 4; t++) {
    for (auto& pc : pair_copulas[t]) {
      auto par =
        Eigen::VectorXd::Constant(1, 2.0 / (static_cast<double>(t) + 1.0));
      pc = Bicop(BicopFamily::clayton, 90, par);
    }
  }
  RVineStructure str(std::vector<size_t>{ 1, 2, 3, 4, 5 });
  auto var_types = std::vector<std::string>{ "d", "c", "d", "d", "c" };
  Vinecop vc(str, pair_copulas, var_types);

  // simulate data with continuous and discrete variables
  size_t n = 500;
  auto utmp = vc.simulate(n, true, 1, { 1 });
  Eigen::MatrixXd u(n, 5 + 3); // 3 discrete vars
  u.leftCols(5) = utmp;

  u.col(0) = (utmp.col(0).array() * 10).ceil() / 10;
  u.col(5 + 0) = (utmp.col(0).array() * 10).floor() / 10;

  u.col(2) = (utmp.col(2).array() * 10).ceil() / 10;
  u.col(5 + 1) = (utmp.col(2).array() * 10).floor() / 10;

  u.col(3) = (utmp.col(3).array() * 10).ceil() / 10;
  u.col(5 + 2) = (utmp.col(3).array() * 10).floor() / 10;

  // fit vine
  auto controls = FitControlsVinecop({ BicopFamily::clayton });
  auto vc2 = vc;
  // controls.set_show_trace(true);
  vc2.select(u, controls);
  vc2.pdf(u);

  // check output
  auto pcs = vc2.get_all_pair_copulas();
  for (size_t t = 0; t < 4; t++) {
    for (auto pc : pcs[t]) {
      EXPECT_EQ(pc.get_rotation(), 90);
      EXPECT_NEAR(
        pc.get_parameters()(0), 2.0 / (static_cast<double>(t) + 1.0), 0.5);
    }
  }

  for (auto& pc : pcs[0])
    pc.set_parameters(Eigen::VectorXd::Constant(1, 1));
  Vinecop vc3(vc2.get_rvine_structure(), pcs, var_types);
  vc3.fit(u, controls);

  ASSERT_TRUE(vc2.str() == vc3.str());

  // test other input format
  u = Eigen::MatrixXd(n, 10);
  u.leftCols(5) = utmp;
  u.rightCols(5) = utmp;
  u.col(0) = (utmp.col(0).array() * 10).ceil() / 10;
  u.col(5) = (utmp.col(0).array() * 10).floor() / 10;
  u.col(2) = (utmp.col(2).array() * 10).ceil() / 10;
  u.col(7) = (utmp.col(2).array() * 10).floor() / 10;
  u.col(3) = (utmp.col(3).array() * 10).ceil() / 10;
  u.col(8) = (utmp.col(3).array() * 10).floor() / 10;
  vc2.select(u, controls);
  vc2.pdf(u);
  pcs = vc2.get_all_pair_copulas();
  for (size_t t = 0; t < 4; t++) {
    for (auto pc : pcs[t]) {
      EXPECT_EQ(pc.get_rotation(), 90);
      EXPECT_NEAR(
        pc.get_parameters()(0), 2.0 / (static_cast<double>(t) + 1.0), 0.5);
    }
  }

  // test for approximate uniformity of rosenblatt transformation
  auto u4 = vc.rosenblatt(u.topRows(100), 1, true, { 5 });
  for (int i = 0; i < 5; i++) {
    auto w = tools_stats::to_pseudo_obs(u4.col(i));
    // close to KS test with FWER ~ 0.001
    EXPECT_LE(std::sqrt(n) * (u4.col(i) - w).cwiseAbs().maxCoeff(), 3.5);
    // Kendall's tau test with FWER ~ 0.0001
    for (int j = i + 1; j < 5; j++) {
      EXPECT_GE(wdm::Indep_test(wdm::utils::convert_vec(u4.col(i)),
                                wdm::utils::convert_vec(u4.col(j)),
                                "kendall")
                  .p_value(),
                0.0001);
    }
  }

  // only check that it works
  vc.inverse_rosenblatt(u.topRows(10));
}

TEST(zero_inflated, vinecop)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(5);
  for (size_t t = 0; t < 4; t++) {
    for (auto& pc : pair_copulas[t]) {
      auto par =
        Eigen::VectorXd::Constant(1, 2.0 / (static_cast<double>(t) + 1.0));
      pc = Bicop(BicopFamily::clayton, 90, par);
    }
  }
  RVineStructure str(std::vector<size_t>{ 1, 2, 3, 4, 5 });
  Vinecop vc(str, pair_copulas, { "d", "c", "d", "d", "c" });

  // simulate data with continuous and zero-inflated variables
  auto utmp = vc.simulate(500, true, 1, { 1 });
  Eigen::MatrixXd u(utmp.rows(), 5 + 3); // 3 discrete vars
  u.leftCols(5) = utmp;

  auto thresh = Eigen::VectorXd::Constant(utmp.rows(), 0.1);
  auto zero = Eigen::VectorXd::Zero(utmp.rows());
  u.col(0) = utmp.col(0).cwiseMax(thresh);
  u.col(5 + 0) = (utmp.col(0).array() < 0.1).select(zero, u.col(0));

  u.col(2) = utmp.col(2).cwiseMax(thresh);
  u.col(5 + 1) = (utmp.col(2).array() < 0.1).select(zero, u.col(2));

  u.col(3) = utmp.col(3).cwiseMax(thresh);
  u.col(5 + 2) = (utmp.col(3).array() < 0.1).select(zero, u.col(3));

  // fit vine
  auto controls = FitControlsVinecop({ BicopFamily::clayton });
  // controls.set_show_trace(true);
  vc.select(u, controls);
  vc.pdf(u);

  // check output
  auto pcs = vc.get_all_pair_copulas();
  for (size_t t = 0; t < 4; t++) {
    for (auto pc : pcs[t]) {
      EXPECT_EQ(pc.get_rotation(), 90);
      EXPECT_NEAR(
        pc.get_parameters()(0), 2.0 / (static_cast<double>(t) + 1.0), 0.5);
    }
  }
}
}
