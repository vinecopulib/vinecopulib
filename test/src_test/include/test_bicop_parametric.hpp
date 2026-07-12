// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "parbicop_test.hpp"
#include "rscript.hpp"
#include "test_utils.hpp"
#include <cmath>
#include <limits>

namespace test_bicop_parametric {
using namespace vinecopulib;
using namespace tools_stl;
std::vector<int> rotations = { 0, 90, 180, 270 };
using test_utils::all_close;

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
  ASSERT_TRUE(
    all_close(bicop_.get_parameters(), pc.get_parameters(), 1e-4, 1e-4));
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

    // check Blomqvist's beta against VineCopula
    absdiff = fabs(bicop_.get_beta() - results(0, 9));
    ASSERT_TRUE(absdiff < 1e-2) << bicop_.str();

    // check the two diagonal tail dependence coefficients against VineCopula
    // (VineCopula only reports the lower/upper, i.e. diagonal, corners)
    Eigen::MatrixXd taildep = bicop_.get_taildep();
    ASSERT_TRUE(fabs(taildep(0, 0) - results(0, 10)) < 1e-2) << bicop_.str();
    ASSERT_TRUE(fabs(taildep(1, 1) - results(0, 11)) < 1e-2) << bicop_.str();

    // Get u-data
    Eigen::MatrixXd u = results.block(0, 1, n, 2);

    // evaluate pdf in C++
    Eigen::VectorXd f = bicop_.pdf(u);
    // assert approximate equality
    ASSERT_TRUE(all_close(f, results.block(0, 3, n, 1), 1e-4, 1e-4))
      << bicop_.str();

    // evaluate cdf in C++
    f = bicop_.cdf(u);
    // assert approximate equality
    ASSERT_TRUE(all_close(f, results.block(0, 4, n, 1), 1e-4, 1e-4))
      << bicop_.str();

    // evaluate hfunc1 in C++
    f = bicop_.hfunc1(u);
    // assert approximate equality
    ASSERT_TRUE(all_close(f, results.block(0, 5, n, 1), 1e-4, 1e-4))
      << bicop_.str();

    // evaluate hfunc2 in C++
    f = bicop_.hfunc2(u);
    // assert approximate equality
    ASSERT_TRUE(all_close(f, results.block(0, 6, n, 1), 1e-4, 1e-4))
      << bicop_.str();

    // evaluate hinv1 in C++
    f = bicop_.hinv1(u);
    // assert approximate equality
    ASSERT_TRUE(all_close(f, results.block(0, 7, n, 1), 1e-4, 1e-4))
      << bicop_.str();

    // evaluate hinv2 in C++
    f = bicop_.hinv2(u);
    // assert approximate equality
    ASSERT_TRUE(all_close(f, results.block(0, 8, n, 1), 1e-4, 1e-4))
      << bicop_.str();

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
    EXPECT_NO_THROW(bicop_.simulate(10, true, { 1 }));
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
    auto data = bicop_.simulate(get_n(), false, { 1 });
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
      auto data = bicop_.simulate(get_n(), false, { 1 });
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
  // fixed seeds keep the test deterministic across runs and platforms
  Eigen::MatrixXd u =
    bicop_.simulate(static_cast<size_t>(n), false, { 1, 2, 3, 4, 5 });

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

  ASSERT_TRUE(all_close(bicop_.pdf(u, P), ref_pdf)) << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.cdf(u, P), ref_cdf)) << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hfunc1(u, P), ref_h1)) << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hfunc2(u, P), ref_h2)) << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hinv1(u, P), ref_i1)) << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hinv2(u, P), ref_i2)) << bicop_.str();

  // loglik matches the (NaN-ignoring) sum of the looped log-densities
  double ref_ll = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) {
    double lp = std::log(ref_pdf(i));
    if (!(std::isnan)(lp))
      ref_ll += lp;
  }
  ASSERT_NEAR(bicop_.loglik(u, P, 1), ref_ll, 1e-8 * (1.0 + std::abs(ref_ll)))
    << bicop_.str();

  // threading parity (results must not depend on num_threads)
  ASSERT_TRUE(all_close(bicop_.pdf(u, P, 3), bicop_.pdf(u, P, 1), 1e-12, 1e-12))
    << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.cdf(u, P, 3), bicop_.cdf(u, P, 1), 1e-12, 1e-12))
    << bicop_.str();
  ASSERT_TRUE(
    all_close(bicop_.hfunc2(u, P, 3), bicop_.hfunc2(u, P, 1), 1e-12, 1e-12))
    << bicop_.str();
  ASSERT_TRUE(
    all_close(bicop_.hinv1(u, P, 3), bicop_.hinv1(u, P, 1), 1e-12, 1e-12))
    << bicop_.str();

  // a single parameter set per row (broadcast) matches the single-arg path
  Eigen::MatrixXd Pb = bicop_.get_parameters().transpose().replicate(n, 1);
  ASSERT_TRUE(all_close(bicop_.pdf(u, Pb), bicop_.pdf(u))) << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.cdf(u, Pb), bicop_.cdf(u))) << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hfunc1(u, Pb), bicop_.hfunc1(u)))
    << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hfunc2(u, Pb), bicop_.hfunc2(u)))
    << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hinv1(u, Pb), bicop_.hinv1(u))) << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hinv2(u, Pb), bicop_.hinv2(u))) << bicop_.str();

  // validation errors
  EXPECT_ANY_THROW(bicop_.pdf(u, P.topRows(n - 1))); // wrong number of rows
  {
    Eigen::MatrixXd Pbad(n, p + 1);
    Pbad.leftCols(p) = P;
    Pbad.col(p).setConstant(0.5);
    EXPECT_ANY_THROW(bicop_.pdf(u, Pbad)); // wrong number of columns
  }
  {
    Eigen::MatrixXd Pnan = P;
    Pnan(0, 0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_ANY_THROW(bicop_.pdf(u, Pnan)); // NaN parameter
  }
  {
    Eigen::MatrixXd Poob = P;
    Poob(0, 0) = lb(0) - 1.0;              // below the lower bound
    EXPECT_ANY_THROW(bicop_.pdf(u, Poob)); // out-of-bounds parameter
  }

  // discrete parity: per-row parameters through the discrete density,
  // h-function and h-inverse code paths (pdf_c_d / pdf_d_d, the discrete
  // h-functions, and the numeric h-inverses) for all combinations of
  // discrete/continuous variables.
  {
    Eigen::MatrixXd u4(n, 4);
    u4.leftCols(2) = u;
    u4.rightCols(2) = (u.array() * 0.9).matrix(); // left limits below the cdf
    const std::vector<std::vector<std::string>> var_type_sets = {
      { "d", "d" }, { "c", "d" }, { "d", "c" }
    };
    for (const auto& vt : var_type_sets) {
      Bicop disc = bicop_;
      disc.set_var_types(vt);
      Eigen::VectorXd r_pdf(n), r_h1(n), r_h2(n), r_i1(n), r_i2(n);
      for (Eigen::Index i = 0; i < n; ++i) {
        Bicop bi(family, rotation, P.row(i).transpose().eval(), vt);
        Eigen::MatrixXd u4i = u4.row(i);
        r_pdf(i) = bi.pdf(u4i)(0);
        r_h1(i) = bi.hfunc1(u4i)(0);
        r_h2(i) = bi.hfunc2(u4i)(0);
        r_i1(i) = bi.hinv1(u4i)(0);
        r_i2(i) = bi.hinv2(u4i)(0);
      }
      ASSERT_TRUE(all_close(disc.pdf(u4, P), r_pdf)) << disc.str();
      ASSERT_TRUE(all_close(disc.hfunc1(u4, P), r_h1)) << disc.str();
      ASSERT_TRUE(all_close(disc.hfunc2(u4, P), r_h2)) << disc.str();
      ASSERT_TRUE(all_close(disc.hinv1(u4, P), r_i1)) << disc.str();
      ASSERT_TRUE(all_close(disc.hinv2(u4, P), r_i2)) << disc.str();
    }
  }

  // when the conditioning variable is continuous, the h-inverse strips the
  // discrete companion column and must reproduce the purely-continuous result
  // (the discrete left-limit provably cannot enter the conditional quantile).
  {
    Eigen::MatrixXd u4(n, 4);
    u4.leftCols(2) = u;
    u4.rightCols(2) = (u.array() * 0.9).matrix();
    Bicop dc = bicop_;
    dc.set_var_types({ "d", "c" }); // hinv2 conditions on the continuous var 2
    Bicop cd = bicop_;
    cd.set_var_types({ "c", "d" }); // hinv1 conditions on the continuous var 1
    ASSERT_TRUE(all_close(dc.hinv2(u4, P), bicop_.hinv2(u, P), 0.0, 1e-12))
      << bicop_.str();
    ASSERT_TRUE(all_close(cd.hinv1(u4, P), bicop_.hinv1(u, P), 0.0, 1e-12))
      << bicop_.str();
  }

  // the h-inverse inverts its own h-function (continuous round-trip):
  // hfunc1(u1, hinv1(u1, w)) == w and hfunc2(hinv2(w, u2), u2) == w
  {
    Eigen::MatrixXd u1inv = u;
    u1inv.col(1) = bicop_.hinv1(u, P);
    ASSERT_TRUE(all_close(bicop_.hfunc1(u1inv, P), u.col(1), 1e-5, 1e-5))
      << bicop_.str();
    Eigen::MatrixXd u2inv = u;
    u2inv.col(0) = bicop_.hinv2(u, P);
    ASSERT_TRUE(all_close(bicop_.hfunc2(u2inv, P), u.col(0), 1e-5, 1e-5))
      << bicop_.str();
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

// The independence copula has no parameters (p = 0); the per-row overloads
// accept an n x 0 matrix and fall back to the parameter-free leaves.
TEST(BicopPerRowParameters, independence_zero_parameters)
{
  Bicop ind(BicopFamily::indep);
  const Eigen::Index n = 20;
  Eigen::MatrixXd u = ind.simulate(static_cast<size_t>(n), false, { 1, 2, 3 });
  Eigen::MatrixXd P(n, 0); // no parameters
  EXPECT_TRUE(all_close(ind.pdf(u, P), ind.pdf(u)));
  EXPECT_TRUE(all_close(ind.cdf(u, P), ind.cdf(u)));
  EXPECT_TRUE(all_close(ind.hfunc1(u, P), ind.hfunc1(u)));
  EXPECT_TRUE(all_close(ind.hfunc2(u, P), ind.hfunc2(u)));
  EXPECT_TRUE(all_close(ind.hinv1(u, P), ind.hinv1(u)));
  EXPECT_TRUE(all_close(ind.hinv2(u, P), ind.hinv2(u)));
}

INSTANTIATE_TEST_SUITE_P(
  ParBicopTest,
  ParBicopTest,
  testing::Combine(testing::ValuesIn(bicop_families::parametric),
                   testing::ValuesIn(rotations)));
}
