// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "parbicop_test.hpp"
#include "rscript.hpp"
#include <cmath>

namespace test_bicop_parametric {
using namespace vinecopulib;
using namespace tools_stl;
std::vector<int> rotations = { 0, 90, 180, 270 };

// Test that the serialization works
TEST_P(ParBicopTest, bicop_serialization_is_correct)
{

  bicop_.to_file(std::string("temp"));
  Bicop pc(std::string("temp"));

  // Remove temp file
  std::string cmd = rm + "temp";
  int sys_exit_code = system(cmd.c_str());
  if (sys_exit_code != 0) {
    throw std::runtime_error("error in system call");
  }

  EXPECT_EQ(bicop_.get_rotation(), pc.get_rotation());
  EXPECT_EQ(bicop_.get_family_name(), pc.get_family_name());
  EXPECT_EQ(bicop_.get_var_types(), pc.get_var_types());
  EXPECT_EQ(bicop_.get_npars(), pc.get_npars());
  ASSERT_TRUE(bicop_.get_parameters().isApprox(pc.get_parameters(), 1e-4));
}

// Test if the C++ implementation of the basic methods is correct
TEST_P(ParBicopTest, parametric_bicop_is_correct)
{
  if (needs_check_) {
    std::string cmd = std::string(RSCRIPT) + std::string(TEST_BICOP);
    cmd += " " + std::to_string(get_n());
    cmd += " " + std::to_string(get_family());
    cmd += " " + std::to_string(get_par());
    cmd += " " + std::to_string(get_par2());
    int sys_exit_code = system(cmd.c_str());
    if (sys_exit_code != 0) {
      throw std::runtime_error("error in system call");
    }

    int n = get_n();

    Eigen::MatrixXd results = tools_eigen::read_matxd("temp");

    // Remove temp file
    cmd = rm + "temp";
    sys_exit_code += system(cmd.c_str());
    if (sys_exit_code != 0) {
      throw std::runtime_error("error in system call");
    }

    auto absdiff = fabs(bicop_.get_tau() - results(0, 0));
    ASSERT_TRUE(absdiff < 1e-2) << bicop_.str();

    // Get u-data
    Eigen::MatrixXd u = results.block(0, 1, n, 2);

    // evaluate pdf in C++
    Eigen::VectorXd f = bicop_.pdf(u);
    // assert approximate equality
    ASSERT_TRUE(f.isApprox(results.block(0, 3, n, 1), 1e-4)) << bicop_.str();

    // evaluate cdf in C++
    f = bicop_.cdf(u);
    // assert approximate equality
    ASSERT_TRUE(f.isApprox(results.block(0, 4, n, 1), 1e-4)) << bicop_.str();

    // evaluate hfunc1 in C++
    f = bicop_.hfunc1(u);
    // assert approximate equality
    ASSERT_TRUE(f.isApprox(results.block(0, 5, n, 1), 1e-4)) << bicop_.str();

    // evaluate hfunc2 in C++
    f = bicop_.hfunc2(u);
    // assert approximate equality
    ASSERT_TRUE(f.isApprox(results.block(0, 6, n, 1), 1e-4)) << bicop_.str();

    // evaluate hinv1 in C++
    f = bicop_.hinv1(u);
    // assert approximate equality
    ASSERT_TRUE(f.isApprox(results.block(0, 7, n, 1), 1e-4)) << bicop_.str();

    // evaluate hinv2 in C++
    f = bicop_.hinv2(u);
    // assert approximate equality
    ASSERT_TRUE(f.isApprox(results.block(0, 8, n, 1), 1e-4)) << bicop_.str();

    u(0, 0) = std::numeric_limits<double>::quiet_NaN();
    u(1, 1) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_NO_THROW(bicop_.pdf(u.block(0, 0, 10, 2))) << bicop_.str();
    EXPECT_TRUE(bicop_.pdf(u.block(0, 0, 1, 2)).array().isNaN()(0))
      << bicop_.str();
    EXPECT_NO_THROW(bicop_.cdf(u.block(0, 0, 10, 2))) << bicop_.str();
    EXPECT_TRUE(bicop_.cdf(u.block(0, 0, 1, 2)).array().isNaN()(0))
      << bicop_.str();
    EXPECT_NO_THROW(bicop_.hfunc1(u.block(0, 0, 10, 2))) << bicop_.str();
    EXPECT_TRUE(bicop_.hfunc1(u.block(0, 0, 1, 2)).array().isNaN()(0))
      << bicop_.str();
    EXPECT_NO_THROW(bicop_.hinv1(u.block(0, 0, 10, 2))) << bicop_.str();
    EXPECT_TRUE(bicop_.hinv1(u.block(0, 0, 1, 2)).array().isNaN()(0))
      << bicop_.str();
    EXPECT_NO_THROW(bicop_.hfunc2(u.block(0, 0, 10, 2))) << bicop_.str();
    EXPECT_TRUE(bicop_.hfunc2(u.block(0, 0, 1, 2)).array().isNaN()(0))
      << bicop_.str();
    EXPECT_NO_THROW(bicop_.hinv2(u.block(0, 0, 10, 2))) << bicop_.str();
    EXPECT_TRUE(bicop_.hinv2(u.block(0, 0, 1, 2)).array().isNaN()(0))
      << bicop_.str();
    EXPECT_NO_THROW(bicop_.loglik(u.block(0, 0, 10, 2))) << bicop_.str();
    EXPECT_ANY_THROW(bicop_.loglik());
    EXPECT_ANY_THROW(bicop_.aic());
    EXPECT_ANY_THROW(bicop_.bic());
    EXPECT_ANY_THROW(bicop_.mbic());
    EXPECT_NO_THROW(bicop_.simulate(10, true));
    EXPECT_NO_THROW(bicop_.str());
    if ((bicop_.get_parameters().size() > 1) &&
        (bicop_.get_family() != BicopFamily::student)) {
      EXPECT_ANY_THROW(bicop_.tau_to_parameters(0.5));
    } else if (bicop_.get_parameters().size() == 1) {
      EXPECT_NO_THROW(bicop_.tau_to_parameters(1));
    }
  }
}

// Test if the C++ implementation of select method using the mle is correct
TEST_P(ParBicopTest, bicop_select_mle_bic_is_correct)
{
  std::vector<int> positive_rotations = { 0, 180 };
  auto true_family = bicop_.get_family_name();
  auto true_rotation = bicop_.get_rotation();
  FitControlsBicop controls(
    { BicopFamily::indep, BicopFamily::gaussian, bicop_.get_family() },
    "mle",
    "quadratic",
    1.0,
    30,
    "bic");

  if (needs_check_) {
    auto data = bicop_.simulate(get_n());
    auto bicop = Bicop(data, controls);
    EXPECT_EQ(bicop.loglik(data), bicop.get_loglik()) << bicop_.str();

    EXPECT_NO_THROW(bicop.loglik());
    EXPECT_NO_THROW(bicop.aic());
    EXPECT_NO_THROW(bicop.bic());
    EXPECT_NO_THROW(bicop.mbic());

    auto selected_family = bicop.get_family_name();
    EXPECT_EQ(selected_family, true_family)
      << bicop_.str() << std::endl
      << bicop.bic(data) << " " << bicop_.bic(data);

    if (is_member(bicop_.get_family(), bicop_families::bb)) {
      int rot_sel = bicop.get_rotation();
      if (is_member(true_rotation, positive_rotations)) {
        EXPECT_TRUE(is_member(rot_sel, positive_rotations));
      } else {
        EXPECT_FALSE(is_member(rot_sel, positive_rotations));
      }
    } else {
      EXPECT_EQ(bicop.get_rotation(), true_rotation)
        << bicop_.str() << std::endl
        << bicop.bic(data) << " " << bicop_.bic(data);
    }
  }
}

// Test if the C++ implementation of select method using itau is correct
TEST_P(ParBicopTest, bicop_select_itau_bic_is_correct)
{
  if (is_member(bicop_.get_family(), bicop_families::itau)) {
    auto true_family = bicop_.get_family_name();
    auto true_rotation = bicop_.get_rotation();
    FitControlsBicop controls(
      { BicopFamily::indep, BicopFamily::gaussian, bicop_.get_family() },
      "itau",
      "quadratic",
      1.0,
      30,
      "bic");

    if (needs_check_) {
      auto data = bicop_.simulate(get_n());
      auto bicop = Bicop(data, controls);
      auto selected_family = bicop.get_family_name();
      EXPECT_EQ(selected_family, true_family)
        << bicop.bic(data) << " " << bicop_.bic(data);
      EXPECT_EQ(bicop.get_rotation(), true_rotation)
        << bicop_.str() << std::endl
        << bicop.bic(data) << " " << bicop_.bic(data);
    }
  }
}

// Test that per-row-parameter evaluation matches a row-by-row loop over
// single-parameter Bicop objects (the unchanged, state-based path).
TEST_P(ParBicopTest, per_row_parameters_match_loop)
{
  if (!needs_check_)
    return;
  // per-row parameters are meaningless for the parameterless independence
  // copula
  if (bicop_.get_parameters().size() == 0)
    return;

  auto family = bicop_.get_family();
  auto rotation = bicop_.get_rotation();
  auto var_types = bicop_.get_var_types();
  Eigen::Index n = 100;
  Eigen::MatrixXd u = bicop_.simulate(n);

  // build n distinct, in-bounds parameter sets (one per row)
  Eigen::VectorXd lb = bicop_.get_parameters_lower_bounds();
  Eigen::VectorXd ub = bicop_.get_parameters_upper_bounds();
  Eigen::Index p = lb.size();
  Eigen::VectorXd a = lb + 0.2 * (ub - lb);
  Eigen::VectorXd c = lb + 0.6 * (ub - lb);
  Eigen::MatrixXd P(n, p);
  for (Eigen::Index i = 0; i < n; ++i) {
    double w = static_cast<double>(i % 5) / 4.0;
    P.row(i) = ((1 - w) * a + w * c).transpose();
  }

  // ground truth: loop over single-parameter Bicops (state-based path)
  Eigen::VectorXd ref_pdf(n), ref_cdf(n), ref_h1(n), ref_h2(n), ref_i1(n),
    ref_i2(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    Bicop bi(family, rotation, P.row(i).transpose().eval(), var_types);
    Eigen::MatrixXd ui = u.row(i);
    ref_pdf(i) = bi.pdf(ui)(0);
    ref_cdf(i) = bi.cdf(ui)(0);
    ref_h1(i) = bi.hfunc1(ui)(0);
    ref_h2(i) = bi.hfunc2(ui)(0);
    ref_i1(i) = bi.hinv1(ui)(0);
    ref_i2(i) = bi.hinv2(ui)(0);
  }

  ASSERT_TRUE(bicop_.pdf(u, P).isApprox(ref_pdf, 1e-8)) << bicop_.str();
  ASSERT_TRUE(bicop_.cdf(u, P).isApprox(ref_cdf, 1e-8)) << bicop_.str();
  ASSERT_TRUE(bicop_.hfunc1(u, P).isApprox(ref_h1, 1e-8)) << bicop_.str();
  ASSERT_TRUE(bicop_.hfunc2(u, P).isApprox(ref_h2, 1e-8)) << bicop_.str();
  ASSERT_TRUE(bicop_.hinv1(u, P).isApprox(ref_i1, 1e-8)) << bicop_.str();
  ASSERT_TRUE(bicop_.hinv2(u, P).isApprox(ref_i2, 1e-8)) << bicop_.str();

  // loglik matches the (NaN-ignoring) sum of the looped log-densities
  double ref_ll = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) {
    double lp = std::log(ref_pdf(i));
    if (!(std::isnan)(lp))
      ref_ll += lp;
  }
  ASSERT_NEAR(bicop_.loglik(u, P, 1), ref_ll, 1e-8) << bicop_.str();

  // threading parity (results must not depend on num_threads)
  ASSERT_TRUE(bicop_.pdf(u, P, 3).isApprox(bicop_.pdf(u, P, 1), 1e-12))
    << bicop_.str();
  ASSERT_TRUE(bicop_.hinv1(u, P, 3).isApprox(bicop_.hinv1(u, P, 1), 1e-12))
    << bicop_.str();

  // a single parameter set per row (broadcast) matches the single-arg path
  Eigen::MatrixXd Pb = bicop_.get_parameters().transpose().replicate(n, 1);
  ASSERT_TRUE(bicop_.pdf(u, Pb).isApprox(bicop_.pdf(u), 1e-8)) << bicop_.str();
  ASSERT_TRUE(bicop_.hfunc1(u, Pb).isApprox(bicop_.hfunc1(u), 1e-8))
    << bicop_.str();

  // validation errors
  EXPECT_ANY_THROW(bicop_.pdf(u, P.topRows(n - 1))); // wrong number of rows
  {
    Eigen::MatrixXd Pbad(n, p + 1);
    Pbad.leftCols(p) = P;
    Pbad.col(p).setConstant(0.5);
    EXPECT_ANY_THROW(bicop_.pdf(u, Pbad)); // wrong number of columns
  }

  // discrete parity (both variables discrete)
  {
    Bicop disc = bicop_;
    disc.set_var_types({ "d", "d" });
    Eigen::MatrixXd u4(n, 4);
    u4.leftCols(2) = u;
    u4.rightCols(2) = (u.array() * 0.9).matrix(); // left limits below the cdf
    Eigen::VectorXd ref_dpdf(n);
    for (Eigen::Index i = 0; i < n; ++i) {
      Bicop bi(family, rotation, P.row(i).transpose().eval(), { "d", "d" });
      ref_dpdf(i) = bi.pdf(u4.row(i))(0);
    }
    ASSERT_TRUE(disc.pdf(u4, P).isApprox(ref_dpdf, 1e-8)) << bicop_.str();
  }
}

// Test that nonparametric families reject per-row parameters
TEST(BicopPerRowParameters, tll_throws)
{
  Bicop tll(BicopFamily::tll);
  Eigen::MatrixXd u(2, 2);
  u << 0.3, 0.4, 0.5, 0.6;
  Eigen::MatrixXd parameters(2, 1);
  parameters << 1.0, 1.0;
  EXPECT_ANY_THROW(tll.pdf(u, parameters));
}

INSTANTIATE_TEST_SUITE_P(
  ParBicopTest,
  ParBicopTest,
  testing::Combine(testing::ValuesIn(bicop_families::parametric),
                   testing::ValuesIn(rotations)));
}
