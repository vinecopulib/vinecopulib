// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include <boost/math/distributions/negative_binomial.hpp>
#include <cmath>
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
    // The whole sample, not topRows(20): now that the latent sample
    // reaches the fit, 20 observations no longer identify the density.
    bc.select(uu, FitControlsBicop({ BicopFamily::tll }));
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
    // The whole sample, not topRows(20): now that the latent sample
    // reaches the fit, 20 observations no longer identify the density.
    bc.select(uu, FitControlsBicop({ BicopFamily::tll }));
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
    // The whole sample, not topRows(20): now that the latent sample
    // reaches the fit, 20 observations no longer identify the density.
    bc.select(uu, FitControlsBicop({ BicopFamily::tll }));
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
  auto u_compact = u;

  // fit vine
  auto controls = FitControlsVinecop({ BicopFamily::clayton });
  auto vc2 = vc;
  // controls.set_show_trace(true);
  vc2.select(u, controls);
  vc2.pdf(u);

  // check output
  auto pcs = vc2.get_all_pair_copulas();
  for (size_t t = 0; t < 4; t++) {
    for (const auto& pc : pcs[t]) {
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
  EXPECT_TRUE(
    vc.pdf(u_compact.topRows(100)).isApprox(vc.pdf(u.topRows(100)), 1e-12));
  EXPECT_TRUE(
    vc.rosenblatt(u_compact.topRows(100), 1, true, { 5 })
      .isApprox(vc.rosenblatt(u.topRows(100), 1, true, { 5 }), 1e-12));
  vc2.select(u, controls);
  vc2.pdf(u);
  pcs = vc2.get_all_pair_copulas();
  for (size_t t = 0; t < 4; t++) {
    for (const auto& pc : pcs[t]) {
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
    for (const auto& pc : pcs[t]) {
      EXPECT_EQ(pc.get_rotation(), 90);
      EXPECT_NEAR(
        pc.get_parameters()(0), 2.0 / (static_cast<double>(t) + 1.0), 0.5);
    }
  }
}

TEST(discrete, check_d_d_stability)
{
  constexpr size_t n = 200;
  constexpr double rho_true = 0.6;

  auto cop =
    Bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, rho_true));
  auto u = cop.simulate(n, false, { 5 });

  auto nbinom1 =
    boost::math::negative_binomial_distribution<>(2, 2.0 / (302.0));
  auto nbinom2 = boost::math::negative_binomial_distribution<>(2, 1.0 / (51.0));

  Eigen::MatrixXd u_disc(n, 4);
  for (size_t i = 0; i < n; ++i) {
    const double x1 = boost::math::quantile(nbinom1, u(i, 0));
    const double x2 = boost::math::quantile(nbinom2, u(i, 1));

    u_disc(i, 0) = boost::math::cdf(nbinom1, x1);
    u_disc(i, 1) = boost::math::cdf(nbinom2, x2);
    u_disc(i, 2) = boost::math::cdf(nbinom1, std::max(x1 - 1, 0.0));
    u_disc(i, 3) = boost::math::cdf(nbinom2, std::max(x2 - 1, 0.0));
  }

  const std::vector<std::string> var_types = { "d", "d" };
  auto vine_controls = FitControlsVinecop({ BicopFamily::gaussian }, "mle");
  auto bicop_controls = FitControlsBicop({ BicopFamily::gaussian }, "mle");

  auto vinecop_fit =
    Vinecop(u_disc, RVineStructure(), var_types, vine_controls);
  auto vinecop_fit2 = Vinecop(
    u_disc, vinecop_fit.get_rvine_structure(), var_types, vine_controls);

  Bicop bicop_fit;
  bicop_fit.set_var_types(var_types);
  bicop_fit.select(u_disc, bicop_controls);

  const double vine1 =
    vinecop_fit.get_all_pair_copulas().at(0).at(0).get_parameters()(0);
  const double vine2 =
    vinecop_fit2.get_all_pair_copulas().at(0).at(0).get_parameters()(0);
  const double bicop = bicop_fit.get_parameters()(0);

  EXPECT_NEAR(vine1, vine2, 1e-2);
  EXPECT_NEAR(vine1, bicop, 1e-2);
  EXPECT_NEAR(vine2, bicop, 1e-2);
}

// When an atom is narrower than the 5e-5 cutoff, pdf_d_d replaces the
// difference quotient in that argument by an h-function evaluated at the atom's
// midpoint. Both h-function evaluations have to use the same midpoint, so the
// result is the difference quotient of the remaining argument.
TEST(discrete, d_d_degenerate_branch_uses_the_atom_midpoint)
{
  const auto par = Eigen::VectorXd::Constant(1, 0.6);
  const auto reference = Bicop(BicopFamily::gaussian, 0, par);
  auto bicop = Bicop(BicopFamily::gaussian, 0, par);
  bicop.set_var_types({ "d", "d" });

  const double narrow = 4e-5; // below the cutoff; the wide atom is 0.4
  Eigen::MatrixXd u(1, 4);
  Eigen::MatrixXd upper(1, 2), lower(1, 2);

  // narrow atom in the first argument: h1 collapsed in u1, differenced in u2
  u << 0.5, 0.7, 0.5 - narrow, 0.3;
  const double mid1 = (u(0, 0) + u(0, 2)) / 2;
  upper << mid1, u(0, 1);
  lower << mid1, u(0, 3);
  EXPECT_NEAR(bicop.pdf(u)(0),
              (reference.hfunc1(upper)(0) - reference.hfunc1(lower)(0)) /
                (u(0, 1) - u(0, 3)),
              1e-12);

  // and the mirror case: h2 collapsed in u2, differenced in u1
  u << 0.7, 0.5, 0.3, 0.5 - narrow;
  const double mid2 = (u(0, 1) + u(0, 3)) / 2;
  upper << u(0, 0), mid2;
  lower << u(0, 2), mid2;
  EXPECT_NEAR(bicop.pdf(u)(0),
              (reference.hfunc2(upper)(0) - reference.hfunc2(lower)(0)) /
                (u(0, 0) - u(0, 2)),
              1e-12);
}

// Regression test for PR #700: the per-row parameter evaluation in the
// discrete leaves (pdf_c_d, pdf_d_d and the discrete branches of hfunc1 /
// hfunc2) must not index the parameter matrix row-by-row for families whose
// parameter matrix is neither 1 x p nor n x p. Two such shapes exist:
//   - a parameterless family (independence, 0 x 0 parameters); and
//   - a nonparametric family (tll, a 30 x 30 interpolation grid).
// Both used to read out of bounds (Eigen abort with assertions on, silent UB
// under -DNDEBUG) whenever a discrete variable was involved.
TEST(discrete, edge_case_parameter_shapes)
{
  // Discretized sample from a Clayton copula.
  auto clayton =
    Bicop(BicopFamily::clayton, 0, Eigen::VectorXd::Constant(1, 3));
  auto u = clayton.simulate(200, true, { 1 });

  Eigen::MatrixXd u_disc(u.rows(), 4);
  u_disc.col(0) = (u.col(0).array() * 5).ceil() / 5;
  u_disc.col(2) = (u.col(0).array() * 5).floor() / 5;
  u_disc.col(1) = (u.col(1).array() * 5).ceil() / 5;
  u_disc.col(3) = (u.col(1).array() * 5).floor() / 5;

  // d_c / c_d layouts reuse the discrete corner columns.
  Eigen::MatrixXd u_dc(u.rows(), 3);
  u_dc.col(0) = u_disc.col(0);
  u_dc.col(1) = u.col(1);
  u_dc.col(2) = u_disc.col(2);

  Eigen::MatrixXd u_cd(u.rows(), 4);
  u_cd.col(0) = u.col(0);
  u_cd.col(2) = u.col(0);
  u_cd.col(1) = u_disc.col(1);
  u_cd.col(3) = u_disc.col(3);

  // Parameterless family: independence has a 0 x 0 parameter matrix. Its
  // density is 1 everywhere, discrete margins included, in every layout.
  auto check_indep = [](const std::vector<std::string>& var_types,
                        const Eigen::MatrixXd& uu) {
    auto indep = Bicop(BicopFamily::indep);
    indep.set_var_types(var_types);
    Eigen::VectorXd pdf;
    ASSERT_NO_THROW(pdf = indep.pdf(uu));
    EXPECT_LE((pdf.array() - 1.0).abs().maxCoeff(), 1e-10);
    EXPECT_NO_THROW(indep.hfunc1(uu));
    EXPECT_NO_THROW(indep.hfunc2(uu));
  };
  check_indep({ "d", "d" }, u_disc);
  check_indep({ "c", "d" }, u_cd);
  check_indep({ "d", "c" }, u_dc);

  // Nonparametric family: tll stores a 30 x 30 grid. Evaluating on more than
  // 30 rows must not index the grid row-by-row.
  ASSERT_GT(u_disc.rows(), 30);
  auto tll = Bicop();
  tll.set_var_types({ "d", "d" });
  tll.select(u_disc, FitControlsBicop({ BicopFamily::tll }));
  Eigen::VectorXd pdf_tll;
  ASSERT_NO_THROW(pdf_tll = tll.pdf(u_disc));
  EXPECT_GE(pdf_tll.minCoeff(), 0);
  EXPECT_TRUE(pdf_tll.allFinite());
  EXPECT_NO_THROW(tll.hfunc1(u_disc));
  EXPECT_NO_THROW(tll.hfunc2(u_disc));
}

// A discrete/continuous density must satisfy, for every u2,
//
//     sum_atoms c(u1, u2) * (u1 - u1^-)  =  h2(1, u2) - h2(0, u2)  =  1,
//
// because the difference quotients telescope over the atoms. KernelBicop used
// to override pdf / hfunc1 / hfunc2 and return the *continuous* density at the
// atom midpoint instead, which violates this by up to 10% and does not converge
// away as the atoms shrink. The parametric families satisfy it exactly, so the
// same check is run on gaussian as a control.
TEST(discrete, kernel_pdf_is_normalized)
{
  auto gauss =
    Bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.6));
  auto u = gauss.simulate(2000, true, { 7 });

  // Fit continuously and declare the types afterwards, so this exercises the
  // evaluation path alone and not the (separate) discrete fit.
  auto tll = Bicop();
  tll.select(u, FitControlsBicop({ BicopFamily::tll }));

  for (const auto& family : { BicopFamily::tll, BicopFamily::gaussian }) {
    auto bc = family == BicopFamily::tll ? tll : gauss;
    bc.set_var_types({ "d", "c" });
    for (size_t k : { 2, 8, 32 }) {
      for (double u2 : { 0.1, 0.5, 0.9 }) {
        Eigen::MatrixXd uu(k, 4);
        for (size_t i = 0; i < k; ++i) {
          uu(i, 0) = static_cast<double>(i + 1) / static_cast<double>(k);
          uu(i, 2) = static_cast<double>(i) / static_cast<double>(k);
          uu(i, 1) = u2;
          uu(i, 3) = u2;
        }
        const auto pdf = bc.pdf(uu);
        const double mass =
          (pdf.array() * (uu.col(0) - uu.col(2)).array()).sum();
        EXPECT_NEAR(mass, 1.0, 1e-6) << "family " << bc.get_family_name()
                                     << ", k = " << k << ", u2 = " << u2;
      }
    }
  }
}

// TllBicop::fit computes a latent sample for discrete data via
// find_latent_sample. Its result used to be assigned to `psobs` and never fed
// back into `z_data`, so the fitted grid was bit-identical to the one the same
// data produces when declared continuous -- i.e. declaring a variable discrete
// had no effect on the fit at all.
TEST(discrete, kernel_fit_uses_the_latent_sample)
{
  auto gauss =
    Bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.6));
  auto u = gauss.simulate(2000, true, { 3 });

  Eigen::MatrixXd u_disc(u.rows(), 4);
  u_disc.col(0) = (u.col(0).array() * 5).ceil() / 5;
  u_disc.col(2) = (u.col(0).array() * 5).floor() / 5;
  u_disc.col(1) = (u.col(1).array() * 5).ceil() / 5;
  u_disc.col(3) = (u.col(1).array() * 5).floor() / 5;

  auto cont = Bicop();
  cont.select(u_disc.leftCols(2), FitControlsBicop({ BicopFamily::tll }));

  auto disc = Bicop();
  disc.set_var_types({ "d", "d" });
  disc.select(u_disc, FitControlsBicop({ BicopFamily::tll }));

  const double diff =
    (disc.get_parameters() - cont.get_parameters()).array().abs().maxCoeff();
  EXPECT_GT(diff, 1e-6);
}

}
