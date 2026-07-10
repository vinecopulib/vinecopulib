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

// Analytic (or fallback) derivatives must match central finite differences
// of the public (rotation-aware) pdf/hfunc1/hfunc2. This validates both the
// derivative formulas and the facade's rotation chain rule; second-order
// selectors are checked against finite differences of the analytic first
// derivatives, in both component orders (a Schwarz-symmetry check).
TEST_P(ParBicopTest, derivatives_match_finite_differences)
{
  if (!needs_check_)
    return;

  auto family = bicop_.get_family();
  auto rotation = bicop_.get_rotation();

  // deterministic interior grid (derivatives explode at the boundary)
  const Eigen::Index m = 9;
  Eigen::MatrixXd u(m * m, 2);
  for (Eigen::Index i = 0; i < m; ++i) {
    for (Eigen::Index j = 0; j < m; ++j) {
      u(i * m + j, 0) = 0.05 + 0.1 * static_cast<double>(i);
      u(i * m + j, 1) = 0.05 + 0.1 * static_cast<double>(j);
    }
  }

  // keep parameters strictly inside the bounds so that the central
  // finite-difference reference never degrades to a one-sided difference
  // (0 x 0 matrices for the independence copula can't be mapped to vectors)
  const Eigen::Index p = bicop_.get_parameters().size();
  Eigen::VectorXd par(p), lb(p), ub(p);
  if (p > 0) {
    par = bicop_.get_parameters();
    lb = bicop_.get_parameters_lower_bounds();
    ub = bicop_.get_parameters_upper_bounds();
    for (Eigen::Index k = 0; k < p; ++k) {
      par(k) = std::min(std::max(par(k), lb(k) + 1e-2), ub(k) - 1e-2);
    }
  }
  Bicop cop(family, rotation, par);

  // finite-difference helpers on the public (rotated) methods
  auto fd_par = [&](const std::function<Eigen::VectorXd(const Bicop&)>& eval,
                    Eigen::Index k) {
    Eigen::VectorXd par_p = par, par_m = par;
    double h = 1e-4 * std::max(1.0, std::abs(par(k)));
    par_p(k) = par(k) + h;
    par_m(k) = par(k) - h;
    Bicop cop_p(family, rotation, par_p);
    Bicop cop_m(family, rotation, par_m);
    return ((eval(cop_p) - eval(cop_m)) / (2 * h)).eval();
  };
  auto fd_u = [&](const std::function<Eigen::VectorXd(
                    const Bicop&, const Eigen::MatrixXd&)>& eval,
                  Eigen::Index col) {
    Eigen::MatrixXd u_p = u, u_m = u;
    u_p.col(col).array() += 1e-5;
    u_m.col(col).array() -= 1e-5;
    return ((eval(cop, u_p) - eval(cop, u_m)) / 2e-5).eval();
  };

  // all first-order components: parameters, then u1, u2
  std::vector<std::string> comps;
  for (Eigen::Index k = 0; k < p; ++k) {
    comps.push_back("par" + std::to_string(k + 1));
  }
  comps.push_back("u1");
  comps.push_back("u2");

  auto fd_first = [&](const std::string& method, const std::string& comp) {
    auto eval = [&method](const Bicop& b,
                          const Eigen::MatrixXd& uu) -> Eigen::VectorXd {
      if (method == "pdf")
        return b.pdf(uu);
      if (method == "hfunc1")
        return b.hfunc1(uu);
      return b.hfunc2(uu);
    };
    if (comp == "u1")
      return fd_u(eval, 0);
    if (comp == "u2")
      return fd_u(eval, 1);
    Eigen::Index k = std::stoi(comp.substr(3)) - 1;
    return fd_par([&](const Bicop& b) { return eval(b, u); }, k);
  };
  auto eval_first = [&](const Bicop& b,
                        const Eigen::MatrixXd& uu,
                        const std::string& method,
                        const std::string& comp) -> Eigen::VectorXd {
    if (method == "pdf")
      return b.pdf_deriv(uu, comp);
    if (method == "hfunc1")
      return b.hfunc1_deriv(uu, comp);
    return b.hfunc2_deriv(uu, comp);
  };

  // first order: analytic vs finite differences of the public methods
  for (const std::string method : { "pdf", "hfunc1", "hfunc2" }) {
    for (const auto& comp : comps) {
      ASSERT_TRUE(all_close(
        eval_first(cop, u, method, comp), fd_first(method, comp), 1e-3, 1e-3))
        << cop.str() << "method: " << method << ", deriv: " << comp;
    }
  }

  // second order: analytic vs finite differences of the analytic first
  // derivatives, differencing each component of the pair in turn
  auto check_second = [&](const std::string& method) {
    for (const auto& c1 : comps) {
      for (const auto& c2 : comps) {
        Eigen::VectorXd ref;
        if (c2 == "u1" || c2 == "u2") {
          Eigen::Index col = (c2 == "u1") ? 0 : 1;
          ref = fd_u(
            [&](const Bicop& b, const Eigen::MatrixXd& uu) {
              return eval_first(b, uu, method, c1);
            },
            col);
        } else {
          Eigen::Index k = std::stoi(c2.substr(3)) - 1;
          ref = fd_par(
            [&](const Bicop& b) { return eval_first(b, u, method, c1); }, k);
        }
        Eigen::VectorXd val;
        if (method == "pdf") {
          val = cop.pdf_deriv2(u, c1 + c2);
        } else if (method == "hfunc1") {
          val = cop.hfunc1_deriv2(u, c1 + c2);
        } else {
          val = cop.hfunc2_deriv2(u, c1 + c2);
        }
        ASSERT_TRUE(all_close(val, ref, 5e-3, 5e-3))
          << cop.str() << "method: " << method << ", deriv: " << c1 + c2;
      }
    }
  };
  check_second("pdf");
  check_second("hfunc1");
  check_second("hfunc2");

  // log-density derivatives: quotient-rule identity against pdf derivatives
  Eigen::ArrayXd c = cop.pdf(u).array();
  for (const auto& comp : comps) {
    Eigen::VectorXd ref = (cop.pdf_deriv(u, comp).array() / c).matrix();
    ASSERT_TRUE(all_close(cop.logpdf_deriv(u, comp), ref, 1e-6, 1e-6))
      << cop.str() << "logpdf deriv: " << comp;
  }
  if (p > 0) {
    Eigen::ArrayXd c_1 = cop.pdf_deriv(u, "par1").array();
    Eigen::VectorXd ref =
      (cop.pdf_deriv2(u, "par1par1").array() / c - (c_1 / c).square()).matrix();
    ASSERT_TRUE(all_close(cop.logpdf_deriv2(u, "par1par1"), ref, 1e-6, 1e-6))
      << cop.str();
  }

  // NaN propagation without throwing
  Eigen::MatrixXd u_nan = u.topRows(5);
  u_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_NO_THROW(cop.pdf_deriv(u_nan, "u1")) << cop.str();
  EXPECT_TRUE(cop.pdf_deriv(u_nan.topRows(1), "u1").array().isNaN()(0))
    << cop.str();
}

// Per-row-parameter derivative overloads must match a loop over
// single-parameter-set copulas, be independent of the number of threads,
// and reduce to the stored-parameter path under broadcasting.
TEST_P(ParBicopTest, derivatives_per_row_parameters_match_loop)
{
  if (!needs_check_)
    return;
  if (bicop_.get_parameters().size() == 0)
    return;

  auto family = bicop_.get_family();
  auto rotation = bicop_.get_rotation();
  Eigen::Index n = 50;
  Eigen::MatrixXd u =
    bicop_.simulate(static_cast<size_t>(n), false, { 1, 2, 3, 4, 5 });

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

  const std::vector<std::string> derivs1 = { "par1", "u1", "u2" };
  const std::vector<std::string> derivs2 = { "par1par1", "par1u1", "u1u2" };
  for (const auto& d : derivs1) {
    Eigen::VectorXd r_pdf(n), r_h1(n), r_h2(n), r_lpdf(n);
    for (Eigen::Index i = 0; i < n; ++i) {
      Bicop bi(family, rotation, P.row(i).transpose().eval());
      Eigen::MatrixXd ui = u.row(i);
      r_pdf(i) = bi.pdf_deriv(ui, d)(0);
      r_h1(i) = bi.hfunc1_deriv(ui, d)(0);
      r_h2(i) = bi.hfunc2_deriv(ui, d)(0);
      r_lpdf(i) = bi.logpdf_deriv(ui, d)(0);
    }
    ASSERT_TRUE(all_close(bicop_.pdf_deriv(u, d, P), r_pdf))
      << bicop_.str() << "deriv: " << d;
    ASSERT_TRUE(all_close(bicop_.hfunc1_deriv(u, d, P), r_h1))
      << bicop_.str() << "deriv: " << d;
    ASSERT_TRUE(all_close(bicop_.hfunc2_deriv(u, d, P), r_h2))
      << bicop_.str() << "deriv: " << d;
    ASSERT_TRUE(all_close(bicop_.logpdf_deriv(u, d, P), r_lpdf))
      << bicop_.str() << "deriv: " << d;
  }
  for (const auto& d : derivs2) {
    Eigen::VectorXd r_pdf(n);
    for (Eigen::Index i = 0; i < n; ++i) {
      Bicop bi(family, rotation, P.row(i).transpose().eval());
      r_pdf(i) = bi.pdf_deriv2(u.row(i), d)(0);
    }
    // looser tolerance: the broadcast and per-row leaves may take different
    // (equivalent) evaluation paths, and nested finite differencing of the
    // fallback amplifies the resulting machine-precision noise
    ASSERT_TRUE(all_close(bicop_.pdf_deriv2(u, d, P), r_pdf, 1e-5, 1e-5))
      << bicop_.str() << "deriv: " << d;
  }

  // threading parity
  ASSERT_TRUE(all_close(bicop_.pdf_deriv(u, "par1", P, 3),
                        bicop_.pdf_deriv(u, "par1", P, 1),
                        1e-12,
                        1e-12))
    << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hfunc1_deriv2(u, "par1u1", P, 3),
                        bicop_.hfunc1_deriv2(u, "par1u1", P, 1),
                        1e-12,
                        1e-12))
    << bicop_.str();

  // broadcasting a single parameter set matches the stored-parameter path
  Eigen::MatrixXd Pb = bicop_.get_parameters().transpose().replicate(n, 1);
  ASSERT_TRUE(
    all_close(bicop_.pdf_deriv(u, "u1", Pb), bicop_.pdf_deriv(u, "u1")))
    << bicop_.str();
  ASSERT_TRUE(all_close(bicop_.hfunc2_deriv(u, "par1", Pb),
                        bicop_.hfunc2_deriv(u, "par1")))
    << bicop_.str();

  // validation errors propagate through format_parameters
  EXPECT_ANY_THROW(bicop_.pdf_deriv(u, "par1", P.topRows(n - 1)));
}

// Derivatives must match the VineCopula R implementation (BiCopDeriv,
// BiCopDeriv2, BiCopHfuncDeriv, BiCopHfuncDeriv2). VineCopula encodes
// rotations in the family code with negated parameters, so parameter
// derivatives flip sign once per parameter differentiation at 90/270
// (chain rule for theta -> -theta); argument derivatives need no fix.
// VineCopula's BiCopHfuncDeriv* differentiate h2(u1|u2); oracles for h1
// come from swapped arguments with the 90/270 code decades exchanged.
TEST_P(ParBicopTest, parametric_bicop_derivatives_match_R)
{
  if (!needs_check_)
    return;
  // the R oracle only implements derivatives for these families (indep is
  // trivial and covered by the finite-difference/identity tests)
  std::vector<BicopFamily> supported = {
    BicopFamily::gaussian, BicopFamily::student, BicopFamily::clayton,
    BicopFamily::gumbel,   BicopFamily::frank,   BicopFamily::joe
  };
  if (!is_member(bicop_.get_family(), supported))
    return;

  int n = 500;
  std::string cmd = std::string(RSCRIPT) + std::string(TEST_BICOP_DERIV);
  cmd += " " + std::to_string(n);
  cmd += " " + std::to_string(get_family());
  cmd += " " + std::to_string(get_par());
  cmd += " " + std::to_string(get_par2());
  int sys_exit_code = system(cmd.c_str());
  if (sys_exit_code != 0) {
    throw std::runtime_error("error in system call");
  }
  Eigen::MatrixXd results = tools_eigen::read_matxd("temp_deriv");
  cmd = rm + "temp_deriv";
  sys_exit_code += system(cmd.c_str());
  if (sys_exit_code != 0) {
    throw std::runtime_error("error in system call");
  }

  Eigen::MatrixXd u = results.block(0, 0, n, 2);
  bool st = (bicop_.get_family() == BicopFamily::student);
  bool r90 = is_member(bicop_.get_rotation(), { 90, 270 });
  double s = r90 ? -1.0 : 1.0;
  auto col = [&](int k) { return results.block(0, k - 1, n, 1).eval(); };
  auto check = [&](const Eigen::VectorXd& val,
                   int k,
                   double sign,
                   const std::string& what) {
    ASSERT_TRUE(all_close(val, (sign * col(k)).eval(), 1e-4, 1e-4))
      << bicop_.str() << what << " (oracle column " << k << ")";
  };

  // first derivatives of the density
  check(bicop_.pdf_deriv(u, "par1"), 3, s, "pdf_deriv par1");
  if (st)
    check(bicop_.pdf_deriv(u, "par2"), 4, s, "pdf_deriv par2");
  check(bicop_.pdf_deriv(u, "u1"), 5, 1.0, "pdf_deriv u1");
  check(bicop_.pdf_deriv(u, "u2"), 6, 1.0, "pdf_deriv u2");

  // first derivatives of the log-density
  check(bicop_.logpdf_deriv(u, "par1"), 7, s, "logpdf_deriv par1");
  if (st)
    check(bicop_.logpdf_deriv(u, "par2"), 8, s, "logpdf_deriv par2");

  // second derivatives of the density; VineCopula's pure second parameter
  // derivatives are inconsistent with its own pdf for 90/270 codes (the
  // 90-degree branch of diff2PDF_mod reflects (u, 1-v) where diffPDF_mod
  // uses (1-u, v), deriv2.c:53-66), so those columns are skipped there and
  // covered by the finite-difference self-consistency test instead.
  if (!r90)
    check(bicop_.pdf_deriv2(u, "par1par1"), 9, 1.0, "pdf_deriv2 par1par1");
  if (st)
    check(bicop_.pdf_deriv2(u, "par2par2"), 10, 1.0, "pdf_deriv2 par2par2");
  check(bicop_.pdf_deriv2(u, "u1u1"), 11, 1.0, "pdf_deriv2 u1u1");
  check(bicop_.pdf_deriv2(u, "u2u2"), 12, 1.0, "pdf_deriv2 u2u2");
  if (st)
    check(bicop_.pdf_deriv2(u, "par1par2"), 13, 1.0, "pdf_deriv2 par1par2");
  check(bicop_.pdf_deriv2(u, "par1u1"), 14, s, "pdf_deriv2 par1u1");
  if (st)
    check(bicop_.pdf_deriv2(u, "par2u1"), 15, s, "pdf_deriv2 par2u1");
  check(bicop_.pdf_deriv2(u, "par1u2"), 16, s, "pdf_deriv2 par1u2");
  if (st)
    check(bicop_.pdf_deriv2(u, "par2u2"), 17, s, "pdf_deriv2 par2u2");

  // derivatives of the second h-function h2(u1|u2)
  check(bicop_.hfunc2_deriv(u, "par1"), 18, s, "hfunc2_deriv par1");
  if (st)
    check(bicop_.hfunc2_deriv(u, "par2"), 19, s, "hfunc2_deriv par2");
  check(bicop_.hfunc2_deriv(u, "u2"), 20, 1.0, "hfunc2_deriv u2");
  if (!r90)
    check(
      bicop_.hfunc2_deriv2(u, "par1par1"), 21, 1.0, "hfunc2_deriv2 par1par1");
  if (st)
    check(
      bicop_.hfunc2_deriv2(u, "par2par2"), 22, 1.0, "hfunc2_deriv2 par2par2");
  check(bicop_.hfunc2_deriv2(u, "u2u2"), 23, 1.0, "hfunc2_deriv2 u2u2");
  if (st)
    check(
      bicop_.hfunc2_deriv2(u, "par1par2"), 24, 1.0, "hfunc2_deriv2 par1par2");
  check(bicop_.hfunc2_deriv2(u, "par1u2"), 25, s, "hfunc2_deriv2 par1u2");
  if (st)
    check(bicop_.hfunc2_deriv2(u, "par2u2"), 26, s, "hfunc2_deriv2 par2u2");

  // derivatives of the first h-function h1(u2|u1); the swapped oracle's
  // "u2" selector corresponds to our u1 (the conditioning argument)
  check(bicop_.hfunc1_deriv(u, "par1"), 27, s, "hfunc1_deriv par1");
  if (st)
    check(bicop_.hfunc1_deriv(u, "par2"), 28, s, "hfunc1_deriv par2");
  check(bicop_.hfunc1_deriv(u, "u1"), 29, 1.0, "hfunc1_deriv u1");
  if (!r90)
    check(
      bicop_.hfunc1_deriv2(u, "par1par1"), 30, 1.0, "hfunc1_deriv2 par1par1");
  if (st)
    check(
      bicop_.hfunc1_deriv2(u, "par2par2"), 31, 1.0, "hfunc1_deriv2 par2par2");
  check(bicop_.hfunc1_deriv2(u, "u1u1"), 32, 1.0, "hfunc1_deriv2 u1u1");
  if (st)
    check(
      bicop_.hfunc1_deriv2(u, "par1par2"), 33, 1.0, "hfunc1_deriv2 par1par2");
  check(bicop_.hfunc1_deriv2(u, "par1u1"), 34, s, "hfunc1_deriv2 par1u1");
  if (st)
    check(bicop_.hfunc1_deriv2(u, "par2u1"), 35, s, "hfunc1_deriv2 par2u1");
}

// Selector and family validation for the derivative methods
TEST(BicopDerivatives, selector_and_family_validation)
{
  Eigen::MatrixXd u(2, 2);
  u << 0.3, 0.4, 0.5, 0.6;

  // nonparametric families have no derivatives
  Bicop tll(BicopFamily::tll);
  EXPECT_ANY_THROW(tll.pdf_deriv(u, "u1"));
  EXPECT_ANY_THROW(tll.pdf_deriv2(u, "u1u1"));
  EXPECT_ANY_THROW(tll.hfunc1_deriv(u, "u1"));
  EXPECT_ANY_THROW(tll.hfunc2_deriv(u, "u2"));
  EXPECT_ANY_THROW(tll.logpdf_deriv(u, "u1"));

  Bicop cl(BicopFamily::clayton, 0, Eigen::VectorXd::Constant(1, 2.0));

  // discrete variable types are rejected
  Bicop cl_disc = cl;
  cl_disc.set_var_types({ "d", "c" });
  Eigen::MatrixXd u4(2, 4);
  u4 << 0.3, 0.4, 0.2, 0.4, 0.5, 0.6, 0.4, 0.6;
  EXPECT_ANY_THROW(cl_disc.pdf_deriv(u4, "u1"));
  EXPECT_ANY_THROW(cl_disc.hfunc1_deriv(u4, "par1"));

  // selectors referring to nonexistent parameters (or nonsense) are rejected
  EXPECT_ANY_THROW(cl.pdf_deriv(u, "par2"));
  EXPECT_ANY_THROW(cl.pdf_deriv(u, "bogus"));
  EXPECT_ANY_THROW(cl.pdf_deriv(u, ""));
  EXPECT_ANY_THROW(cl.pdf_deriv(u, "par1u1"));    // second order selector
  EXPECT_ANY_THROW(cl.pdf_deriv2(u, "par1u1u2")); // third order selector
  Bicop ind(BicopFamily::indep);
  EXPECT_ANY_THROW(ind.pdf_deriv(u, "par1")); // no parameters

  // aliases: "par" = "par1", components may come in any order
  EXPECT_TRUE(
    all_close(cl.pdf_deriv(u, "par"), cl.pdf_deriv(u, "par1"), 0.0, 0.0));
  EXPECT_TRUE(
    all_close(cl.pdf_deriv2(u, "par"), cl.pdf_deriv2(u, "par1par1"), 0.0, 0.0));
  EXPECT_TRUE(all_close(
    cl.pdf_deriv2(u, "u1par1"), cl.pdf_deriv2(u, "par1u1"), 0.0, 0.0));

  // independence copula: argument derivatives vanish, h-function derivative
  // w.r.t. the conditioned argument is the density (= 1)
  EXPECT_TRUE(
    all_close(ind.pdf_deriv(u, "u1"), Eigen::VectorXd::Zero(2), 1e-10, 1e-10));
  EXPECT_TRUE(all_close(
    ind.hfunc1_deriv(u, "u1"), Eigen::VectorXd::Zero(2), 1e-10, 1e-10));
  EXPECT_TRUE(all_close(
    ind.hfunc1_deriv(u, "u2"), Eigen::VectorXd::Ones(2), 1e-10, 1e-10));
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
