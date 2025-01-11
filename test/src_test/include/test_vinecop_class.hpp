// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "vinecop_test.hpp"
#include <string>
#include <vinecopulib.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace test_vinecop_class {
using namespace vinecopulib;

TEST_F(VinecopTest, constructors_without_error)
{
  Vinecop vinecop(5);
  Vinecop vinecop_indep(model_matrix);

  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 90);
    }
  }

  Vinecop vinecop_parametrized(model_matrix, pair_copulas);
}

TEST_F(VinecopTest, copy)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 90);
    }
  }

  Vinecop vinecop1(model_matrix, pair_copulas);
  auto pc = vinecop1.get_pair_copula(0, 0);
  pc.set_parameters(pc.get_parameters().array() + 1);
  EXPECT_EQ(vinecop1.get_parameters(0, 0), vinecop1.get_parameters(0, 1));

  Vinecop vinecop2 = vinecop1;
  pair_copulas[0][0].set_parameters(pc.get_parameters().array() + 1);
  vinecop2.set_all_pair_copulas(pair_copulas);
  EXPECT_EQ(vinecop1.get_parameters(0, 0), vinecop1.get_parameters(0, 1));
}

TEST_F(VinecopTest, print)
{
  auto cvine = CVineStructure(std::vector<size_t>({ 5, 4, 3, 2, 1 }));
  auto vc1 = Vinecop(cvine);

  // check if first, second and last line are correct
  std::string expected_first_line = "Vinecop model with 5 variables";
  std::string expected_second_line =
    "tree edge conditioned variables conditioning variables var_types       "
    "family rotation parameters df tau ";
  std::string expected_last_line =
    "   4    1                  5, 4                3, 2, 1      c, c "
    "Independence                        0.0 ";

  std::istringstream input;
  input.str(vc1.str());

  std::string line;
  // get first, second and last line
  std::getline(input, line);
  EXPECT_EQ(line, expected_first_line);
  std::getline(input, line);
  EXPECT_EQ(line, expected_second_line);
  std::string last_line;
  while (std::getline(input, line)) {
    last_line = line;
  }
  EXPECT_EQ(last_line, expected_last_line);

  // create vine with 7 variables
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::tawn, 270);
    }
  }

  Vinecop vc2(model_matrix, pair_copulas);

  // check if first, second and last line are correct
  expected_first_line = "Vinecop model with 7 variables";
  expected_second_line =
    "tree edge conditioned variables conditioning variables var_types family "
    "rotation       parameters  df   tau ";
  expected_last_line = "   6    1                  4, 7          3, 1, 2, 6, 5 "
                       "     c, c   Tawn      270 1.00, 1.00, 1.00 3.0 -0.00 ";

  input.clear();
  input.str(vc2.str());

  // get first, second and last line
  std::getline(input, line);
  EXPECT_EQ(line, expected_first_line);
  std::getline(input, line);
  EXPECT_EQ(line, expected_second_line);
  while (std::getline(input, line)) {
    last_line = line;
  }
  EXPECT_EQ(last_line, expected_last_line);

  auto data = tools_stats::simulate_uniform(100, 5);
  auto controls = FitControlsVinecop({ BicopFamily::tll });
  vc1.select(data, controls);

  // check if first and second are correct
  // we don't check the last line as it is random
  // but at least we see it doesn't crash
  expected_first_line = "Vinecop model with 5 variables";
  expected_second_line =
    "tree edge conditioned variables conditioning variables var_types family rotation   parameters   df   tau ";

  input.clear();
  input.str(vc1.str());

  // get first, second and last line
  std::getline(input, line);
  EXPECT_EQ(line, expected_first_line) << vc1.str();
  std::getline(input, line);
  EXPECT_EQ(line, expected_second_line) << vc1.str();
}

TEST_F(VinecopTest, serialization)
{

  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
  mat << 5, 2, 6, 6, 6, 6, 6, 6, 6, 1, 2, 5, 5, 0, 2, 5, 2, 5, 2, 0, 0, 1, 1, 5,
    1, 0, 0, 0, 3, 7, 7, 0, 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0;

  // create vine with 7 variables, 2-truncated
  size_t d = 7;
  auto pc_store = Vinecop::make_pair_copula_store(d, 5);
  for (auto& tree : pc_store) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::bb1, 90);
    }
  }

  auto vc = Vinecop(mat, pc_store);
  vc.truncate(3);

  // serialize
  vc.to_file(std::string("temp"));

  // unserialize
  auto vc2 = Vinecop(std::string("temp"));

  // Remove temp file
  std::string cmd = rm + "temp";
  int sys_exit_code = system(cmd.c_str());
  if (sys_exit_code != 0) {
    throw std::runtime_error("error in system call");
  }

  EXPECT_EQ(vc.get_all_rotations(), vc2.get_all_rotations());
  EXPECT_EQ(vc.get_all_families(), vc2.get_all_families());
  EXPECT_EQ(vc.get_var_types(), vc2.get_var_types());
  EXPECT_EQ(vc.get_matrix(), vc2.get_matrix());
}

TEST_F(VinecopTest, getters_are_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 90);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  for (auto& tree : vinecop.get_all_families()) {
    for (auto& fam : tree) {
      EXPECT_EQ(fam, BicopFamily::clayton);
    }
  }

  for (auto& tree : vinecop.get_all_pair_copulas()) {
    for (auto& pc : tree) {
      EXPECT_EQ(pc.get_family(), BicopFamily::clayton);
      EXPECT_EQ(pc.get_rotation(), 90);
    }
  }

  for (auto& tree : vinecop.get_all_parameters()) {
    for (auto& par : tree) {
      EXPECT_EQ(par.size(), 1);
      EXPECT_EQ(par(0), 1e-10);
    }
  }

  for (auto& tree : vinecop.get_all_rotations()) {
    for (auto& rot : tree) {
      EXPECT_EQ(rot, 90);
    }
  }

  for (auto& tree : vinecop.get_all_taus()) {
    for (auto& tau : tree) {
      ASSERT_TRUE(fabs(tau) < 1e-4);
    }
  }

  EXPECT_NO_THROW(vinecop.get_dim());
  EXPECT_NO_THROW(vinecop.get_rvine_structure());
  EXPECT_ANY_THROW(vinecop.get_loglik());
  EXPECT_ANY_THROW(vinecop.get_nobs());
  EXPECT_ANY_THROW(vinecop.get_aic());
  EXPECT_ANY_THROW(vinecop.get_bic());
  EXPECT_ANY_THROW(vinecop.get_mbicv());
  EXPECT_ANY_THROW(vinecop.aic());
  EXPECT_ANY_THROW(vinecop.bic());
  EXPECT_ANY_THROW(vinecop.mbicv());
}

TEST_F(VinecopTest, 1dim)
{
  auto data = tools_stats::simulate_uniform(3, 1);
  auto vc = Vinecop(1);
  vc.select(data);
  EXPECT_TRUE((vc.pdf(data).array() == 1).all());
  EXPECT_TRUE((vc.rosenblatt(data).array() == data.array()).all());
  EXPECT_TRUE((vc.inverse_rosenblatt(data).array() == data.array()).all());
  vc.loglik();
  vc.aic();
  vc.simulate(3);
  Vinecop::make_pair_copula_store(1, 2);
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(1, 1);
  mat(0, 0) = 1;
  RVineStructure rvine_structure(mat);
}

TEST_F(VinecopTest, fit_statistics_getters_are_correct)
{
  auto data = tools_stats::simulate_uniform(100, 3);
  auto vc = Vinecop(
    data, RVineStructure(), {}, FitControlsVinecop({ BicopFamily::clayton }));
  EXPECT_NEAR(vc.get_loglik(), vc.loglik(data), 1e-10);
  EXPECT_NEAR(static_cast<double>(vc.get_nobs()), 100, 1e-10);
  EXPECT_NEAR(vc.get_aic(), vc.aic(data), 1e-10);
  EXPECT_NEAR(vc.get_bic(), vc.bic(data), 1e-10);
  EXPECT_NEAR(vc.get_mbicv(0.6), vc.mbicv(data, 0.6), 1e-10);
}

TEST_F(VinecopTest, truncate_methods_works)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);
  vinecop.truncate(2);
  EXPECT_EQ(vinecop.get_all_pair_copulas().size(), 2);
  EXPECT_EQ(vinecop.get_rvine_structure().get_trunc_lvl(), 2);
  vinecop.truncate(0);
  EXPECT_EQ(vinecop.get_all_pair_copulas().size(), 0);
  EXPECT_EQ(vinecop.get_rvine_structure().get_trunc_lvl(), 0);
}

TEST_F(VinecopTest, pdf_is_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  ASSERT_TRUE(vinecop.pdf(u).isApprox(f, 1e-4));
}

TEST_F(VinecopTest, cdf_is_correct)
{
  // Create a bivariate copula and a corresponding vine with two variables
  auto pair_copulas = Vinecop::make_pair_copula_store(2);
  auto par = Eigen::VectorXd::Constant(1, 0.5);
  auto bicop = Bicop(BicopFamily::gaussian, 0, par);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::gaussian, 0, par);
    }
  }
  Eigen::Matrix<size_t, 2, 2> matrix;
  matrix << 1, 1, 2, 0;
  Vinecop vinecop(matrix, pair_copulas);

  // Test whether the analytic and simulated versions are "close" enough
  auto u2 = vinecop.simulate(10);
  ASSERT_TRUE(vinecop.cdf(u2, 10000).isApprox(bicop.cdf(u2), 1e-2));

  // verify that qrng stuff works
  Vinecop vinecop2(301);
  vinecop.simulate(10, true);
}

TEST_F(VinecopTest, simulate_is_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  // only check if it works
  vinecop.simulate(10);
  // check the underlying transformation from independent samples
  ASSERT_TRUE(vinecop.inverse_rosenblatt(u).isApprox(sim, 1e-4));

  // verify that qrng stuff works
  vinecop.simulate(10, true);
  Vinecop vinecop2(301);
  vinecop.simulate(10, true);
}

TEST_F(VinecopTest, rosenblatt_is_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);
  auto u2 = vinecop.simulate(5);
  ASSERT_TRUE(
    vinecop.rosenblatt(vinecop.inverse_rosenblatt(u2)).isApprox(u2, 1e-6));

  // truncated multivariate
  pair_copulas = Vinecop::make_pair_copula_store(7, 2);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  vinecop = Vinecop(model_matrix, pair_copulas);
  ASSERT_TRUE(
    vinecop.rosenblatt(vinecop.inverse_rosenblatt(u)).isApprox(u, 1e-6));

  // bivariate case
  pair_copulas = Vinecop::make_pair_copula_store(2);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Eigen::Matrix<size_t, 2, 2> mat;
  mat << 1, 1, 2, 0;
  vinecop = Vinecop(mat, pair_copulas);
  u = vinecop.simulate(5);
  ASSERT_TRUE(
    vinecop.rosenblatt(vinecop.inverse_rosenblatt(u)).isApprox(u, 1e-6));
}

TEST_F(VinecopTest, aic_bic_are_correct)
{
  int d = 7;
  auto data = tools_stats::simulate_uniform(100, 7);
  Vinecop true_model(d);

  auto pair_copulas = Vinecop::make_pair_copula_store(d);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 0, Eigen::VectorXd::Constant(1, 3.0));
    }
  }
  Vinecop complex_model(model_matrix, pair_copulas);

  ASSERT_TRUE(true_model.aic(data) < complex_model.aic(data));
  ASSERT_TRUE(true_model.bic(data) < complex_model.bic(data));
}

TEST_F(VinecopTest, fit_parameters_is_correct)
{
  u.conservativeResize(50, 7);
  auto controls = FitControlsVinecop({ BicopFamily::clayton }, "itau");
  Vinecop vc(7);
  vc.select(u, controls);
  auto rvine_structure = vc.get_rvine_structure();

  auto pcs = vc.get_all_pair_copulas();
  for (auto& pc : pcs[0])
    pc.set_parameters(Eigen::VectorXd::Constant(1, 1));
  Vinecop vc2(rvine_structure, pcs);
  vc2.fit(u, controls);

  ASSERT_TRUE(vc.str() == vc2.str());

  Vinecop vc3(rvine_structure, pcs);
  controls.set_select_families(false);
  vc3.select(u, controls);
  ASSERT_TRUE(vc.str() == vc3.str());
}

TEST_F(VinecopTest, family_select_finds_true_rotations)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);
  auto data = vinecop.simulate(2000);

  auto controls = FitControlsVinecop({ BicopFamily::clayton }, "itau");
  // controls.set_show_trace(true);
  Vinecop fit(data, model_matrix, {}, controls);

  // don't check last two trees to avoid random failures because of
  // estimation uncertainty
  auto true_rots = vinecop.get_all_rotations();
  auto fitd_rots = fit.get_all_rotations();
  true_rots.erase(true_rots.end() - 2, true_rots.end());
  fitd_rots.erase(fitd_rots.end() - 2, fitd_rots.end());
  EXPECT_EQ(true_rots, fitd_rots);
}

TEST_F(VinecopTest, family_select_returns_pcs_in_right_order)
{
  u.conservativeResize(50, 7);
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  auto controls = FitControlsVinecop(bicop_families::itau, "itau");
  // controls.set_show_trace(true);
  Vinecop fit_struct(u, RVineStructure(), {}, controls);
  Vinecop fit_fam(u, fit_struct.get_matrix(), {}, controls);

  EXPECT_EQ(fit_struct.get_all_parameters(), fit_fam.get_all_parameters());
}

TEST_F(VinecopTest, trace_works)
{
  u.conservativeResize(10, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_show_trace(true);
  controls.set_select_threshold(true);
  controls.set_trunc_lvl(3);
  testing::internal::CaptureStdout();
  Vinecop fit(u, RVineStructure(), {}, controls);
  std::string output = testing::internal::GetCapturedStdout();
  EXPECT_TRUE(!output.empty());
}

TEST_F(VinecopTest, works_multi_threaded)
{
  u.conservativeResize(100, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_select_trunc_lvl(true);

  Vinecop fit1(u, RVineStructure(), {}, controls);
  controls.set_num_threads(2);
  Vinecop fit2(u, RVineStructure(), {}, controls);

  auto pcs = fit1.get_all_pair_copulas();
  for (auto& pc : pcs[0])
    pc.set_parameters(Eigen::VectorXd::Constant(1, 1));
  Vinecop fit3(fit1.get_rvine_structure(), pcs);
  fit3.fit(u, controls, 2);

  // check for equality in likelihood, since the pair copulas may be stored
  // in a different order when running in parallel
  EXPECT_NEAR(fit1.loglik(u), fit2.loglik(u), 1e-2);
  EXPECT_NEAR(fit1.loglik(u), fit3.loglik(u), 1e-2);

  // check if parallel evaluators have same output as single threaded ones
  EXPECT_TRUE(fit2.pdf(u, 2).isApprox(fit2.pdf(u), 1e-10));
  EXPECT_TRUE(
    fit2.inverse_rosenblatt(u, 2).isApprox(fit2.inverse_rosenblatt(u), 1e-10));
  EXPECT_TRUE(fit2.rosenblatt(u, 2).isApprox(fit2.rosenblatt(u), 1e-10));

  // just check that it works
  fit2.simulate(2, false, 3);
  fit2.cdf(u, 100, 2);
}

// check if the same conditioned sets appear for each tree
inline size_t
get_pairs_unequal(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix1,
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix2,
  size_t trunc_lvl)
{
  std::vector<std::vector<std::vector<size_t>>> vc_sets(trunc_lvl),
    vcl_sets(trunc_lvl);
  size_t pairs_unequal = 0;
  for (size_t tree = 0; tree < trunc_lvl; ++tree) {
    vc_sets[tree].resize(trunc_lvl - tree);
    vcl_sets[tree].resize(trunc_lvl - tree);
    for (size_t edge = 0; edge < trunc_lvl - tree; ++edge) {
      vc_sets[tree][edge].resize(2);
      vc_sets[tree][edge][0] = matrix1(tree, edge);
      vc_sets[tree][edge][1] = matrix1(trunc_lvl - edge, edge);
      vcl_sets[tree][edge].resize(2);
      vcl_sets[tree][edge][0] = matrix2(tree, edge);
      vcl_sets[tree][edge][1] = matrix2(trunc_lvl - edge, edge);
    }
    for (auto s1 : vc_sets[tree]) {
      bool is_in_both = false;
      for (auto s2 : vcl_sets[tree]) {
        if (tools_stl::is_same_set(s1, s2)) {
          is_in_both = true;
        }
      }
      if (!is_in_both) {
        ++pairs_unequal;
      }
    }
  }
  return pairs_unequal;
}

TEST_F(VinecopTest, select_finds_right_structure_prim)
{
  // check whether the same structure appears if we only allow for
  // independence (pair-copula estimates differ otherwise)

  // select structure and get matrix
  Vinecop fit(7);
  fit.select(u, FitControlsVinecop({ BicopFamily::indep }));
  auto vcl_matrix = fit.get_matrix();

  // check if the same conditioned sets appear for each tree
  size_t pairs_unequal = get_pairs_unequal(vc_matrix, vcl_matrix, 6);
  EXPECT_EQ(pairs_unequal, 0);
}

TEST_F(VinecopTest, select_finds_right_structure_kruskal)
{
  // check whether the same structure appears if we only allow for
  // independence (pair-copula estimates differ otherwise)
  FitControlsVinecop controls({ BicopFamily::indep });
  EXPECT_ANY_THROW(controls.set_mst_algorithm("foobar"));
  controls.set_mst_algorithm("kruskal");

  // select structure and get matrix
  Vinecop fit(7);
  fit.select(u, controls);
  auto vcl_matrix = fit.get_matrix();

  // check if the same conditioned sets appear for each tree
  size_t pairs_unequal = get_pairs_unequal(vc_matrix, vcl_matrix, 6);
  EXPECT_EQ(pairs_unequal, 0);
}

TEST_F(VinecopTest, fixed_truncation)
{
  u.conservativeResize(10, 7);
  FitControlsVinecop controls({ BicopFamily::indep });
  controls.set_trunc_lvl(2);
  // controls.set_show_trace(true);
  Vinecop fit(7);
  fit.select(u, controls);
  EXPECT_EQ(fit.get_all_pair_copulas().size(), 2);

  fit.select(u, controls);
  EXPECT_EQ(fit.get_all_pair_copulas().size(), 2);

  Vinecop fit2(u, fit.get_rvine_structure(), {}, controls);
  EXPECT_EQ(fit2.get_all_pair_copulas().size(), 2);

  Vinecop fit3(u, fit.get_rvine_structure());
  EXPECT_EQ(fit3.get_all_pair_copulas().size(), 6);

  fit3.truncate(2);
  size_t pairs_unequal =
    get_pairs_unequal(fit3.get_matrix(), fit.get_matrix(), 2);
  EXPECT_EQ(pairs_unequal, 0);
}

TEST_F(VinecopTest, sparse_threshold_selection)
{
  u.conservativeResize(20, 7);

  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_select_threshold(NAN);
  controls.set_threshold(true);
  // controls.set_show_trace(true);
  controls.set_selection_criterion("mbicv");

  Vinecop fit(7);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);

  u = tools_stats::simulate_uniform(100, 7);
  fit.select(u, controls);
}

TEST_F(VinecopTest, sparse_truncation_selection)
{
  u.conservativeResize(50, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_select_trunc_lvl(true);
  // controls.set_show_trace(true);
  u = tools_stats::simulate_uniform(100, 7);
  Vinecop fit(7);
  fit.select(u, controls);
  EXPECT_LE(fit.get_rvine_structure().get_trunc_lvl(), 6);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);
  fit.select(u, controls);
  EXPECT_LE(fit.get_rvine_structure().get_trunc_lvl(), 6);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);
}

TEST_F(VinecopTest, sparse_both_selection)
{
  u.conservativeResize(20, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_select_trunc_lvl(true);
  controls.set_select_threshold(true);
  controls.set_tree_criterion("joe");
  controls.set_selection_criterion("mbicv");
  // controls.set_show_trace(true);
  Vinecop fit(7);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);

  // 1d models
  Eigen::MatrixXd uu = u.col(0).matrix();
  fit = Vinecop(1);
  fit.select(uu, controls);
}

TEST_F(VinecopTest, partial_selection)
{
  u.conservativeResize(20, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  // controls.set_show_trace(true);
  auto fixed = CVineStructure(std::vector<size_t>{ 5, 4, 7, 1, 3, 6, 2 });
  fixed.truncate(1);
  Vinecop fit(fixed);
  fit.select(u, controls);

  // a C-vine with root node 2 contains 6 edges with vertex 2.
  auto rvm = fit.get_rvine_structure().get_matrix();
  int count2 = 0;
  for (int i = 0; i < 6; i++) {
    //  diagonal element         base element
    if ((rvm(6 - i, i) == 2) || (rvm(0, i) == 2))
      count2++;
  }
  EXPECT_EQ(count2, 6);
}
}
