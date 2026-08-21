// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "include/test_utils.hpp"
#include "include/vinecop_test.hpp"
#include <future>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <vinecopulib.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace test_vinecop_class {
using namespace vinecopulib;
using test_utils::all_close;

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

TEST(VinecopConstructor, omitted_pair_copulas_are_independence)
{
  const size_t d = 4;
  Vinecop model(DVineStructure({ 1, 2, 3, 4 }));

  ASSERT_TRUE(model.get_all_pair_copulas().empty());
  auto data = tools_stats::simulate_uniform(20, d, false, { 1 });
  EXPECT_TRUE(all_close(model.pdf(data), Eigen::VectorXd::Ones(data.rows())));

  std::vector<size_t> conditioning_set{ 1, 2 };
  Eigen::MatrixXd u_cond(20, conditioning_set.size());
  u_cond.col(0).setConstant(0.3);
  u_cond.col(1).setConstant(0.8);
  auto simulated =
    model.simulate_conditional(u_cond, conditioning_set, false, 1, { 2 });
  EXPECT_TRUE(all_close(simulated.col(0), u_cond.col(0)));
  EXPECT_TRUE(all_close(simulated.col(1), u_cond.col(1)));
}

// An empty pair-copula store means independence, whatever the structure's
// truncation level, so no structural operation may index it.
TEST(VinecopConstructor, omitted_pair_copulas_survive_structure_operations)
{
  const size_t d = 4;
  const auto data = tools_stats::simulate_uniform(20, d, false, { 1 });
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(data.rows());
  const Vinecop model(DVineStructure({ 1, 2, 3, 4 }));

  // {1} is not the order tail, so this does not take the identity fast path
  Vinecop reoriented = model;
  ASSERT_NO_THROW(reoriented.reorient({ 1 }));
  EXPECT_EQ(reoriented.get_order().back(), 1u);
  EXPECT_TRUE(reoriented.get_all_pair_copulas().empty());
  EXPECT_TRUE(all_close(reoriented.pdf(data), ones));
  EXPECT_NO_THROW(RVineStructure(reoriented.get_matrix()));

  // every edge of the decomposition is independence, and it still round-trips
  const auto trees = model.get_trees();
  ASSERT_EQ(trees.get_trunc_lvl(), d - 1);
  for (const auto& tree : trees.get_trees()) {
    for (const auto& edge : tree) {
      EXPECT_EQ(edge.pair_copula.get_family(), BicopFamily::indep);
    }
  }
  EXPECT_EQ(model.get_rvine_structure(), RVineStructure(trees));

  // truncation keeps the store empty rather than adding empty tree slots, and
  // the decomposition follows the truncated structure
  Vinecop truncated = model;
  truncated.truncate(2);
  EXPECT_EQ(truncated.get_trunc_lvl(), 2u);
  EXPECT_TRUE(truncated.get_all_pair_copulas().empty());
  EXPECT_TRUE(all_close(truncated.pdf(data), ones));
  EXPECT_EQ(truncated.get_trees().get_trunc_lvl(), 2u);

  // a truncated vine that does have pair copulas decomposes into exactly its
  // stored trees, each edge carrying its own copula
  auto pcs = Vinecop::make_pair_copula_store(d);
  for (auto& tree : pcs)
    for (auto& pc : tree)
      pc = Bicop(BicopFamily::clayton, 90, Eigen::VectorXd::Constant(1, 3.0));
  Vinecop fitted(DVineStructure({ 1, 2, 3, 4 }), pcs);
  fitted.truncate(2);
  const auto fitted_trees = fitted.get_trees();
  ASSERT_EQ(fitted_trees.get_trunc_lvl(), 2u);
  for (const auto& tree : fitted_trees.get_trees()) {
    for (const auto& edge : tree) {
      EXPECT_EQ(edge.pair_copula.get_family(), BicopFamily::clayton);
    }
  }

  // variable types are labelled in the original order and survive relabeling
  Vinecop discrete(DVineStructure({ 1, 2, 3, 4 }), {}, { "c", "d", "c", "d" });
  ASSERT_NO_THROW(discrete.reorient({ 1 }));
  EXPECT_EQ(discrete.get_var_types(),
            std::vector<std::string>({ "c", "d", "c", "d" }));
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

  auto data = tools_stats::simulate_uniform(100, 5, false, { 1 });
  auto controls = FitControlsVinecop({ BicopFamily::tll });
  vc1.select(data, controls);

  // check if first and second are correct
  // we don't check the last line as it is random
  // but at least we see it doesn't crash
  expected_first_line = "Vinecop model with 5 variables";
  expected_second_line =
    "tree edge conditioned variables conditioning variables var_types family "
    "rotation   parameters   df   tau ";

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

  test_r_parity::RTempDir dir;
  vc.to_file(dir.file("vinecop.json"));
  auto vc2 = Vinecop(dir.file("vinecop.json"));

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
  auto data = tools_stats::simulate_uniform(3, 1, false, { 1 });
  auto vc = Vinecop(1);
  vc.select(data);
  EXPECT_TRUE((vc.pdf(data).array() == 1).all());
  EXPECT_TRUE((vc.rosenblatt(data).array() == data.array()).all());
  EXPECT_TRUE((vc.inverse_rosenblatt(data).array() == data.array()).all());
  vc.loglik();
  vc.aic();
  vc.simulate(3, false, 1, { 1 });
  Vinecop::make_pair_copula_store(1, 2);
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(1, 1);
  mat(0, 0) = 1;
  RVineStructure rvine_structure(mat);
}

TEST_F(VinecopTest, fit_statistics_getters_are_correct)
{
  auto data = tools_stats::simulate_uniform(100, 3, false, { 1 });
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

  ASSERT_TRUE(all_close(vinecop.pdf(u), f, 1e-4, 1e-4));
}

TEST_F(VinecopTest, hfuncs_is_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  auto r = vinecop.pdf_full(u, 1, true);

  auto rvine_structure = vinecop.get_rvine_structure();
  size_t d = vinecop.get_dim();
  size_t trunc_lvl = rvine_structure.get_trunc_lvl();

  ASSERT_EQ(r.hfunc1.get_trunc_lvl(), trunc_lvl);
  ASSERT_EQ(r.hfunc2.get_trunc_lvl(), trunc_lvl);
  ASSERT_EQ(r.pdf_edges.get_trunc_lvl(), trunc_lvl);
  // the _sub arrays are only filled when at least one variable is discrete
  ASSERT_EQ(r.hfunc1_sub.get_trunc_lvl(), 0);
  ASSERT_EQ(r.hfunc2_sub.get_trunc_lvl(), 0);

  for (size_t tree = 0; tree < trunc_lvl; ++tree) {
    for (size_t edge = 0; edge < d - tree - 1; ++edge) {
      EXPECT_EQ(r.pdf_edges(tree, edge).size(), u.rows());
      if (rvine_structure.needed_hfunc1(tree, edge)) {
        EXPECT_EQ(r.hfunc1(tree, edge).size(), u.rows());
      } else {
        EXPECT_EQ(r.hfunc1(tree, edge).size(), 0);
      }
      if (rvine_structure.needed_hfunc2(tree, edge)) {
        EXPECT_EQ(r.hfunc2(tree, edge).size(), u.rows());
      } else {
        EXPECT_EQ(r.hfunc2(tree, edge).size(), 0);
      }
    }
  }

  // discrete model: _sub arrays must be allocated and filled
  std::vector<std::string> var_types(7, "c");
  var_types[0] = "d";
  Vinecop vinecop_disc(model_matrix, pair_copulas, var_types);
  Eigen::MatrixXd u_disc(u.rows(), 8);
  u_disc.leftCols(7) = u;
  u_disc.col(0) = (u.col(0).array() * 10).ceil() / 10;
  u_disc.col(7) = (u.col(0).array() * 10).floor() / 10;
  auto r_disc = vinecop_disc.pdf_full(u_disc, 1, true);
  ASSERT_EQ(r_disc.hfunc1_sub.get_trunc_lvl(), trunc_lvl);
  ASSERT_EQ(r_disc.hfunc2_sub.get_trunc_lvl(), trunc_lvl);
  for (size_t tree = 0; tree < trunc_lvl; ++tree) {
    for (size_t edge = 0; edge < d - tree - 1; ++edge) {
      if (rvine_structure.needed_hfunc1(tree, edge)) {
        EXPECT_EQ(r_disc.hfunc1_sub(tree, edge).size(), u.rows());
      } else {
        EXPECT_EQ(r_disc.hfunc1_sub(tree, edge).size(), 0);
      }
      if (rvine_structure.needed_hfunc2(tree, edge)) {
        EXPECT_EQ(r_disc.hfunc2_sub(tree, edge).size(), u.rows());
      } else {
        EXPECT_EQ(r_disc.hfunc2_sub(tree, edge).size(), 0);
      }
    }
  }
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
  auto u2 = vinecop.simulate(10, false, 1, { 1 });
  ASSERT_TRUE(
    all_close(vinecop.cdf(u2, 10000, 1, { 1 }), bicop.cdf(u2), 1e-2, 1e-2));

  // verify that qrng stuff works
  Vinecop vinecop2(301);
  vinecop.simulate(10, true, 1, { 1 });
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
  vinecop.simulate(10, false, 1, { 1 });
  // check the underlying transformation from independent samples
  ASSERT_TRUE(all_close(vinecop.inverse_rosenblatt(u), sim, 1e-4, 1e-4));

  // verify that qrng stuff works
  vinecop.simulate(10, true, 1, { 1 });
  Vinecop vinecop2(301);
  vinecop.simulate(10, true, 1, { 1 });
}

TEST_F(VinecopTest, inverse_rosenblatt_does_not_mutate_discrete_model)
{
  const size_t d = 5;
  const size_t n_discrete = 3;
  auto pair_copulas = Vinecop::make_pair_copula_store(d);
  auto par = Eigen::VectorXd::Constant(1, 2.0);
  for (size_t tree = 0; tree < pair_copulas.size(); ++tree) {
    for (auto& pc : pair_copulas[tree]) {
      pc = Bicop(BicopFamily::clayton, tree % 2 == 0 ? 90 : 270, par);
    }
  }
  auto structure = DVineStructure(std::vector<size_t>{ 1, 2, 3, 4, 5 });
  auto var_types = std::vector<std::string>{ "d", "c", "d", "d", "c" };
  Vinecop discrete(structure, pair_copulas, var_types);
  Vinecop continuous(structure, pair_copulas);
  auto w = tools_stats::simulate_uniform(100, d, false, { 17 });

  auto expected = discrete.inverse_rosenblatt(w);
  EXPECT_TRUE(
    all_close(expected, continuous.inverse_rosenblatt(w), 1e-12, 1e-12));

  Eigen::MatrixXd w_compact(w.rows(), d + n_discrete);
  w_compact.leftCols(d) = w;
  w_compact.rightCols(n_discrete) = w.leftCols(n_discrete);
  Eigen::MatrixXd w_expanded(w.rows(), 2 * d);
  w_expanded.leftCols(d) = w;
  w_expanded.rightCols(d) = w;
  EXPECT_TRUE(
    all_close(discrete.inverse_rosenblatt(w_compact), expected, 1e-12, 1e-12));
  EXPECT_TRUE(
    all_close(discrete.inverse_rosenblatt(w_expanded), expected, 1e-12, 1e-12));
  EXPECT_EQ(discrete.get_var_types(), var_types);

  auto pair_copulas_before = discrete.get_all_pair_copulas();
  EXPECT_ANY_THROW(discrete.inverse_rosenblatt(w.leftCols(d - 1)));
  EXPECT_EQ(discrete.get_var_types(), var_types);
  auto pair_copulas_after = discrete.get_all_pair_copulas();
  for (size_t tree = 0; tree < pair_copulas_before.size(); ++tree) {
    for (size_t edge = 0; edge < pair_copulas_before[tree].size(); ++edge) {
      EXPECT_EQ(pair_copulas_after[tree][edge].get_var_types(),
                pair_copulas_before[tree][edge].get_var_types());
    }
  }

  auto tll_data =
    Bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.5))
      .simulate(200, false, { 29 });
  Bicop tll(tll_data, FitControlsBicop({ BicopFamily::tll }));
  std::vector<std::vector<Bicop>> tll_pair_copulas{ { tll } };
  auto tll_structure = DVineStructure(std::vector<size_t>{ 1, 2 });
  Vinecop tll_continuous(tll_structure, tll_pair_copulas);
  Vinecop tll_discrete(tll_structure, tll_pair_copulas, { "d", "c" });
  auto w_tll = tools_stats::simulate_uniform(20, 2, false, { 31 });
  EXPECT_TRUE(all_close(tll_discrete.inverse_rosenblatt(w_tll),
                        tll_continuous.inverse_rosenblatt(w_tll),
                        1e-12,
                        1e-12));
}

TEST_F(VinecopTest, bicop_view_as_continuous_matches_materialized_copula)
{
  Eigen::VectorXd parameters(3);
  parameters << 0.3, 0.8, 2.0;
  auto eval_data = tools_stats::simulate_uniform(40, 2, false, { 37 });

  for (int rotation : { 0, 90, 180, 270 }) {
    Bicop bicop(BicopFamily::tawn, rotation, parameters, { "d", "c" });
    for (bool flipped : { false, true }) {
      Bicop materialized = bicop;
      if (flipped)
        materialized.flip();
      materialized = materialized.as_continuous();
      auto view = BicopView(bicop, flipped).as_continuous();

      EXPECT_EQ(view.get_var_types(), std::vector<std::string>({ "c", "c" }));
      EXPECT_TRUE(all_close(
        view.hfunc1(eval_data), materialized.hfunc1(eval_data), 1e-10, 1e-10));
      EXPECT_TRUE(all_close(
        view.hfunc2(eval_data), materialized.hfunc2(eval_data), 1e-10, 1e-10));
      EXPECT_TRUE(all_close(
        view.hinv2(eval_data), materialized.hinv2(eval_data), 1e-8, 1e-8));
    }
  }
}

TEST_F(VinecopTest, inverse_rosenblatt_is_safe_for_concurrent_calls)
{
  const size_t d = 5;
  auto pair_copulas = Vinecop::make_pair_copula_store(d);
  auto par = Eigen::VectorXd::Constant(1, 2.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 90, par);
    }
  }
  const Vinecop model(DVineStructure(std::vector<size_t>{ 1, 2, 3, 4, 5 }),
                      pair_copulas,
                      { "d", "c", "d", "d", "c" });
  auto w = tools_stats::simulate_uniform(100, d, false, { 23 });
  auto expected = model.inverse_rosenblatt(w);

  std::vector<std::future<Eigen::MatrixXd>> results;
  for (size_t i = 0; i < 4; ++i) {
    results.emplace_back(std::async(std::launch::async, [&model, &w]() {
      return model.inverse_rosenblatt(w);
    }));
  }
  for (auto& result : results) {
    EXPECT_TRUE(all_close(result.get(), expected, 1e-12, 1e-12));
  }
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
  auto u2 = vinecop.simulate(5, false, 1, { 1 });
  ASSERT_TRUE(all_close(
    vinecop.rosenblatt(vinecop.inverse_rosenblatt(u2)), u2, 1e-6, 1e-6));

  // truncated multivariate
  pair_copulas = Vinecop::make_pair_copula_store(7, 2);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  vinecop = Vinecop(model_matrix, pair_copulas);
  ASSERT_TRUE(all_close(
    vinecop.rosenblatt(vinecop.inverse_rosenblatt(u)), u, 1e-6, 1e-6));

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
  u = vinecop.simulate(5, false, 1, { 1 });
  ASSERT_TRUE(all_close(
    vinecop.rosenblatt(vinecop.inverse_rosenblatt(u)), u, 1e-6, 1e-6));
}

// helper: a Clayton D-vine on 1..d, whose order is 1..d so the sampling order
// (and hence the admissible conditioning tail) is well-defined and R-free.
namespace {
inline Vinecop
make_clayton_dvine(size_t d, double par_value = 4.0, int rotation = 0)
{
  auto pcs = Vinecop::make_pair_copula_store(d);
  auto par = Eigen::VectorXd::Constant(1, par_value);
  for (auto& tree : pcs)
    for (auto& pc : tree)
      pc = Bicop(BicopFamily::clayton, rotation, par);
  std::vector<size_t> order(d);
  for (size_t i = 0; i < d; ++i)
    order[i] = i + 1;
  return Vinecop(DVineStructure(order), pcs);
}
}

TEST_F(VinecopTest, reoriented_transforms_match_materialized_model)
{
  const size_t d = 5;
  auto pcs = Vinecop::make_pair_copula_store(d);
  auto clayton_par = Eigen::VectorXd::Constant(1, 2.0);
  for (size_t tree = 0; tree < pcs.size(); ++tree) {
    for (size_t edge = 0; edge < pcs[tree].size(); ++edge) {
      int rotation = static_cast<int>(90 * ((tree + edge) % 4));
      pcs[tree][edge] = Bicop(BicopFamily::clayton, rotation, clayton_par);
    }
  }
  Eigen::VectorXd tawn_par(3);
  tawn_par << 0.3, 0.8, 2.0;
  pcs[0][0] = Bicop(BicopFamily::tawn, 90, tawn_par);

  std::vector<size_t> order{ 1, 2, 3, 4, 5 };
  Vinecop model(DVineStructure(order), pcs);
  std::vector<size_t> conditioning_set{ 1, 2 };
  Vinecop materialized = model;
  materialized.reorient(conditioning_set);

  auto w = tools_stats::simulate_uniform(50, d, false, { 41 });
  auto expected_inverse = materialized.inverse_rosenblatt(w, 2);
  auto actual_inverse = model.inverse_rosenblatt(w, conditioning_set, 2);
  EXPECT_TRUE(all_close(actual_inverse, expected_inverse, 1e-10, 1e-10));

  auto expected = materialized.rosenblatt(expected_inverse, 2, false);
  auto actual = model.rosenblatt(expected_inverse, conditioning_set, 2, false);
  EXPECT_TRUE(all_close(actual, expected, 1e-10, 1e-10));
  EXPECT_TRUE(all_close(actual, w, 1e-6, 1e-6));

  // The view-based transform leaves the model representation untouched.
  EXPECT_EQ(model.get_order(), order);

  // Requesting the existing tail is exactly the original transform.
  std::vector<size_t> current_tail{ 5 };
  EXPECT_TRUE(all_close(model.inverse_rosenblatt(w, current_tail),
                        model.inverse_rosenblatt(w),
                        1e-12,
                        1e-12));
  EXPECT_TRUE(
    all_close(model.rosenblatt(expected_inverse, current_tail, 1, false),
              model.rosenblatt(expected_inverse, 1, false),
              1e-12,
              1e-12));
}

TEST_F(VinecopTest, reoriented_transforms_handle_discrete_variables)
{
  const size_t d = 4;
  auto pcs = Vinecop::make_pair_copula_store(d);
  auto par = Eigen::VectorXd::Constant(1, 1.5);
  for (size_t tree = 0; tree < pcs.size(); ++tree) {
    for (auto& pc : pcs[tree])
      pc = Bicop(BicopFamily::clayton, tree % 2 == 0 ? 90 : 270, par);
  }
  std::vector<std::string> var_types{ "d", "c", "d", "c" };
  Vinecop model(DVineStructure({ 1, 2, 3, 4 }), pcs, var_types);
  std::vector<size_t> conditioning_set{ 1 };
  Vinecop materialized = model;
  materialized.reorient(conditioning_set);

  auto values = tools_stats::simulate_uniform(60, d, false, { 43 });
  Eigen::MatrixXd data(values.rows(), d + 2);
  data.leftCols(d) = values;
  data.col(d) = 0.7 * values.col(0);
  data.col(d + 1) = 0.7 * values.col(2);

  EXPECT_TRUE(all_close(model.rosenblatt(data, conditioning_set, 2, false),
                        materialized.rosenblatt(data, 2, false),
                        1e-12,
                        1e-12));
  EXPECT_TRUE(
    all_close(model.rosenblatt(data, conditioning_set, 2, true, { 47 }),
              materialized.rosenblatt(data, 2, true, { 47 }),
              1e-12,
              1e-12));

  auto w = tools_stats::simulate_uniform(60, d, false, { 53 });
  EXPECT_TRUE(all_close(model.inverse_rosenblatt(w, conditioning_set, 2),
                        materialized.inverse_rosenblatt(w, 2),
                        1e-12,
                        1e-12));
  EXPECT_EQ(model.get_var_types(), var_types);
}

TEST_F(VinecopTest, reoriented_transforms_handle_tll_without_materializing_grid)
{
  auto training_data =
    Bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.6))
      .simulate(300, false, { 59 });
  Bicop tll(training_data, FitControlsBicop({ BicopFamily::tll }));
  auto pcs = Vinecop::make_pair_copula_store(3);
  for (auto& tree : pcs)
    for (auto& pc : tree)
      pc = tll;

  Vinecop model(DVineStructure({ 1, 2, 3 }), pcs);
  std::vector<size_t> conditioning_set{ 1 };
  Vinecop materialized = model;
  materialized.reorient(conditioning_set);

  auto w = tools_stats::simulate_uniform(20, 3, false, { 61 });
  auto expected_inverse = materialized.inverse_rosenblatt(w);
  auto actual_inverse = model.inverse_rosenblatt(w, conditioning_set);
  EXPECT_TRUE(all_close(actual_inverse, expected_inverse, 1e-8, 1e-8));
  EXPECT_TRUE(
    all_close(model.rosenblatt(actual_inverse, conditioning_set, 1, false),
              materialized.rosenblatt(expected_inverse, 1, false),
              1e-8,
              1e-8));

  auto u_cond = Eigen::MatrixXd::Constant(20, 1, 0.4);
  EXPECT_TRUE(all_close(
    model.simulate_conditional(u_cond, conditioning_set, false, 1, { 63 }),
    materialized.simulate_conditional(u_cond, false, 1, { 63 }),
    1e-8,
    1e-8));
}

TEST_F(VinecopTest, reoriented_transforms_validate_conditioning_set)
{
  auto model = make_clayton_dvine(4, 2.0);
  auto eval_data = tools_stats::simulate_uniform(10, 4, false, { 67 });
  EXPECT_ANY_THROW(
    model.rosenblatt(eval_data, std::vector<size_t>{}, 1, false));
  EXPECT_ANY_THROW(
    model.inverse_rosenblatt(eval_data, std::vector<size_t>{ 2, 2 }));
  EXPECT_ANY_THROW(
    model.rosenblatt(eval_data, std::vector<size_t>{ 5 }, 1, false));
  EXPECT_ANY_THROW(
    model.inverse_rosenblatt(eval_data, std::vector<size_t>{ 1, 3 }));

  // a truncated model is supported: the conditioning-set overload evaluates the
  // re-oriented model without materializing it, and agrees with reorient()
  model.truncate(2);
  Vinecop materialized = model;
  materialized.reorient(std::vector<size_t>{ 1 });
  EXPECT_TRUE(
    all_close(model.inverse_rosenblatt(eval_data, std::vector<size_t>{ 1 }),
              materialized.inverse_rosenblatt(eval_data),
              1e-10,
              1e-10));
}

TEST_F(VinecopTest, simulate_conditional_is_correct)
{
  size_t d = 4;
  auto vc = make_clayton_dvine(d);
  auto order = vc.get_order();
  // the conditioning variables are the last two of the order (drawn first);
  // u_cond column i corresponds to variable order[d - k + i].
  size_t c0 = order[d - 2] - 1; // <-> u_cond column 0
  size_t c1 = order[d - 1] - 1; // <-> u_cond column 1
  size_t n = 1000;

  // (a) fixed conditioning point repeated n times: the conditioning columns
  //     are reproduced exactly and n samples are drawn.
  Eigen::MatrixXd u0(1, 2);
  u0 << 0.3, 0.7;
  auto U = vc.simulate_conditional(u0.replicate(n, 1), false, 1, { 1 });
  ASSERT_EQ(static_cast<size_t>(U.rows()), n);
  ASSERT_EQ(static_cast<size_t>(U.cols()), d);
  EXPECT_TRUE(
    all_close(U.col(c0), Eigen::VectorXd::Constant(n, 0.3), 1e-9, 1e-9));
  EXPECT_TRUE(
    all_close(U.col(c1), Eigen::VectorXd::Constant(n, 0.7), 1e-9, 1e-9));

  // (b) per-sample form: one conditioning row per output sample
  Eigen::MatrixXd ucn(n, 2);
  ucn.col(0) = Eigen::VectorXd::LinSpaced(n, 0.1, 0.9);
  ucn.col(1).setConstant(0.5);
  auto Un = vc.simulate_conditional(ucn, false, 1, { 2 });
  EXPECT_TRUE(all_close(Un.col(c0), ucn.col(0), 1e-9, 1e-9));
  EXPECT_TRUE(all_close(Un.col(c1), ucn.col(1), 1e-9, 1e-9));

  // (c) multi-threaded == single-threaded (identical seeds)
  auto uc = u0.replicate(n, 1);
  auto U_mt = vc.simulate_conditional(uc, false, 2, { 1 });
  EXPECT_TRUE(all_close(U_mt, U, 1e-10, 1e-10));

  // (d) determinism: identical seeds -> identical output
  auto U_again = vc.simulate_conditional(uc, false, 1, { 1 });
  EXPECT_TRUE(all_close(U_again, U, 1e-12, 1e-12));

  // (e) statistical correctness: the conditional mean of a free variable
  //     matches a large unconditional sample filtered near the conditioning
  //     values, and clearly differs from the unconditional mean.
  size_t freevar = order[0] - 1; // drawn last, conditioned on all others
  Eigen::MatrixXd u02(1, 2);
  u02 << 0.8, 0.8;
  auto Uc = vc.simulate_conditional(u02.replicate(100000, 1), false, 1, { 7 });
  double cond_mean = Uc.col(freevar).mean();

  auto big = vc.simulate(1000000, false, 1, { 11 });
  double sum = 0.0, cnt = 0.0, eps = 0.03;
  for (int i = 0; i < big.rows(); ++i) {
    if ((std::abs(big(i, c0) - 0.8) < eps) &&
        (std::abs(big(i, c1) - 0.8) < eps)) {
      sum += big(i, freevar);
      cnt += 1.0;
    }
  }
  double ref_mean = sum / cnt;
  EXPECT_NEAR(cond_mean, ref_mean, 0.02);
  EXPECT_GT(std::abs(cond_mean - big.col(freevar).mean()), 0.1);
}

TEST(VinecopConditionalSimulation, accepts_custom_conditioning_set)
{
  const size_t d = 4;
  auto model = make_clayton_dvine(d, 2.0, 90);
  const auto original_order = model.get_order();
  const std::vector<size_t> conditioning_set{ 2, 1 };

  Eigen::MatrixXd u_cond(200, conditioning_set.size());
  u_cond.col(0).setConstant(0.3);
  u_cond.col(1).setConstant(0.7);
  auto actual =
    model.simulate_conditional(u_cond, conditioning_set, true, 2, { 71 });
  Eigen::MatrixXd u_cond_expanded(u_cond.rows(), 2 * u_cond.cols());
  u_cond_expanded.leftCols(u_cond.cols()) = u_cond;
  u_cond_expanded.rightCols(u_cond.cols()) = u_cond;
  auto actual_expanded = model.simulate_conditional(
    u_cond_expanded, conditioning_set, true, 2, { 71 });

  auto materialized = model;
  materialized.reorient(conditioning_set);
  auto order = materialized.get_order();
  Eigen::MatrixXd u_cond_tail(u_cond.rows(), u_cond.cols());
  for (size_t i = 0; i < conditioning_set.size(); ++i) {
    size_t var = order[d - conditioning_set.size() + i];
    auto it = std::find(conditioning_set.begin(), conditioning_set.end(), var);
    u_cond_tail.col(i) = u_cond.col(it - conditioning_set.begin());
  }
  auto expected =
    materialized.simulate_conditional(u_cond_tail, true, 2, { 71 });

  EXPECT_TRUE(all_close(actual, expected, 1e-12, 1e-12));
  EXPECT_TRUE(all_close(actual_expanded, actual, 1e-12, 1e-12));
  EXPECT_TRUE(all_close(actual.col(1), u_cond.col(0), 1e-12, 1e-12));
  EXPECT_TRUE(all_close(actual.col(0), u_cond.col(1), 1e-12, 1e-12));
  EXPECT_EQ(model.get_order(), original_order);
}

TEST(VinecopConditionalSimulation, custom_set_handles_discrete_variables)
{
  auto model = make_clayton_dvine(4, 2.0, 270);
  model.set_var_types({ "d", "c", "d", "c" });
  const std::vector<size_t> conditioning_set{ 2, 1 };

  Eigen::MatrixXd u_cond(100, 3);
  u_cond.col(0).setConstant(0.4);
  u_cond.col(1).setConstant(0.8);
  u_cond.col(2).setConstant(0.6);
  auto actual =
    model.simulate_conditional(u_cond, conditioning_set, false, 2, { 73 });

  Eigen::MatrixXd u_cond_expanded(100, 4);
  u_cond_expanded.col(0) = u_cond.col(0);
  u_cond_expanded.col(1) = u_cond.col(1);
  u_cond_expanded.col(2) = u_cond.col(0);
  u_cond_expanded.col(3) = u_cond.col(2);
  auto actual_expanded = model.simulate_conditional(
    u_cond_expanded, conditioning_set, false, 2, { 73 });

  EXPECT_TRUE(all_close(actual_expanded, actual, 1e-12, 1e-12));
  EXPECT_TRUE(all_close(actual.col(1), u_cond.col(0), 1e-12, 1e-12));
  EXPECT_TRUE((actual.col(0).array() >= 0.6 - 1e-12).all());
  EXPECT_TRUE((actual.col(0).array() <= 0.8 + 1e-12).all());
}

TEST(VinecopConditionalSimulation, custom_set_validates_inputs)
{
  auto model = make_clayton_dvine(4, 2.0);
  auto u_cond = Eigen::MatrixXd::Constant(3, 2, 0.5);

  EXPECT_ANY_THROW(
    model.simulate_conditional(u_cond, std::vector<size_t>{}, false, 1, { 1 }));
  EXPECT_ANY_THROW(model.simulate_conditional(
    u_cond, std::vector<size_t>{ 1, 1 }, false, 1, { 1 }));
  EXPECT_ANY_THROW(model.simulate_conditional(
    u_cond, std::vector<size_t>{ 1, 3 }, false, 1, { 1 }));
  try {
    model.simulate_conditional(Eigen::MatrixXd::Constant(3, 1, 0.5),
                               std::vector<size_t>{ 1, 2 },
                               false,
                               1,
                               { 1 });
    FAIL() << "expected invalid conditional column count to throw";
  } catch (const std::runtime_error& error) {
    EXPECT_EQ(error.what(),
              std::string("u_cond has wrong number of columns; expected: 2 "
                          "(n x k continuous layout), actual: 1."));
  }

  model.set_var_types({ "d", "c", "c", "c" });
  try {
    model.simulate_conditional(Eigen::MatrixXd::Constant(3, 2, 0.5),
                               std::vector<size_t>{ 1, 2 },
                               false,
                               1,
                               { 1 });
    FAIL() << "expected invalid discrete conditional column count to throw";
  } catch (const std::runtime_error& error) {
    EXPECT_EQ(error.what(),
              std::string("u_cond has wrong number of columns; expected: 4 "
                          "(n x 2k expanded layout) or 3 (n x (k + k_d) "
                          "compact layout), actual: 2."));
  }

  model.set_var_types({ "d", "d", "c", "c" });
  try {
    model.simulate_conditional(Eigen::MatrixXd::Constant(3, 3, 0.5),
                               std::vector<size_t>{ 1, 2 },
                               false,
                               1,
                               { 1 });
    FAIL() << "expected missing discrete conditional column to throw";
  } catch (const std::runtime_error& error) {
    EXPECT_EQ(error.what(),
              std::string("u_cond has wrong number of columns; expected: 4 "
                          "(n x 2k expanded or n x (k + k_d) compact layout), "
                          "actual: 3."));
  }
}

TEST_F(VinecopTest, simulate_conditional_throws)
{
  size_t d = 4;
  auto vc = make_clayton_dvine(d, 2.0);

  // empty conditioning set (0 columns)
  EXPECT_ANY_THROW(
    vc.simulate_conditional(Eigen::MatrixXd(3, 0), false, 1, { 1 }));
  // too many columns (matches no k)
  EXPECT_ANY_THROW(vc.simulate_conditional(
    Eigen::MatrixXd::Constant(3, d, 0.5), false, 1, { 1 }));
  // values outside the unit cube
  Eigen::MatrixXd uc_bad(1, 2);
  uc_bad << 0.3, 1.5;
  EXPECT_ANY_THROW(vc.simulate_conditional(uc_bad, false, 1, { 1 }));

  // discrete vine: a column count that matches no k must throw. With var_types
  // {c,c,d,c} (order 1..4), valid column counts are {1, 3, 4}; 2 is in a gap
  // (it would require variable 3, which is discrete and needs a left-limit).
  auto vc_gap = make_clayton_dvine(d, 2.0);
  vc_gap.set_var_types({ "c", "c", "d", "c" });
  EXPECT_ANY_THROW(vc_gap.simulate_conditional(
    Eigen::MatrixXd::Constant(3, 2, 0.5), false, 1, { 1 }));

  // discrete conditioning variable with F(x^-) > F(x) must throw
  auto vc_disc = make_clayton_dvine(d, 2.0);
  vc_disc.set_var_types({ "c", "c", "c", "d" }); // variable 4 = order tail
  Eigen::MatrixXd uc_bad_lim(1, 2);
  uc_bad_lim << 0.4, 0.7; // F(x) = 0.4 < F(x^-) = 0.7
  EXPECT_ANY_THROW(vc_disc.simulate_conditional(uc_bad_lim, false, 1, { 1 }));
}

TEST_F(VinecopTest, simulate_conditional_discrete_is_correct)
{
  size_t d = 4;
  auto vc = make_clayton_dvine(d, 3.0);
  vc.set_var_types({ "c", "c", "c", "d" }); // variable 4 is discrete
  auto order = vc.get_order();              // [1, 2, 3, 4]

  // k = 1: condition on order[d - 1] = 4 (the discrete variable). u_cond has
  // two columns: F(x) and the left-limit F(x^-), i.e. the atom interval [a, b].
  double a = 0.6, b = 0.8;
  size_t condvar = order[d - 1] - 1; // 3
  size_t freevar = order[0] - 1;     // 0, drawn last (conditioned on all)
  size_t n = 200000;

  Eigen::MatrixXd u_cond(1, 2);
  u_cond << b, a; // F(x) = b, F(x^-) = a
  auto U = vc.simulate_conditional(u_cond.replicate(n, 1), false, 1, { 7 });
  ASSERT_EQ(static_cast<size_t>(U.cols()), d);

  // (1) the conditioning column lands inside the atom interval [a, b] (it is
  //     reproduced up to the atom, not exactly)
  EXPECT_TRUE((U.col(condvar).array() >= a - 1e-9).all());
  EXPECT_TRUE((U.col(condvar).array() <= b + 1e-9).all());

  // (2) statistical correctness: E[free | U_condvar in [a, b]] matches a large
  //     unconditional sample filtered to the atom interval, and differs
  //     clearly from the unconditional mean.
  double cond_mean = U.col(freevar).mean();
  auto big = vc.simulate(2000000, false, 1, { 11 });
  double sum = 0.0, cnt = 0.0;
  for (int i = 0; i < big.rows(); ++i) {
    if ((big(i, condvar) >= a) && (big(i, condvar) <= b)) {
      sum += big(i, freevar);
      cnt += 1.0;
    }
  }
  double ref_mean = sum / cnt;
  EXPECT_NEAR(cond_mean, ref_mean, 0.02);
  EXPECT_GT(std::abs(cond_mean - big.col(freevar).mean()), 0.02);

  // (3) determinism: identical seeds -> identical output
  auto U2 = vc.simulate_conditional(u_cond.replicate(n, 1), false, 1, { 7 });
  EXPECT_TRUE(all_close(U2, U, 1e-12, 1e-12));
}

TEST_F(VinecopTest, simulate_conditional_multi_discrete)
{
  // Two discrete conditioning variables stacked in draw order (the conditioning
  // set is a self-contained sub-vine). The conditioned (free) draws must remain
  // correct; the first-drawn discrete conditioning column lands exactly in its
  // atom, while a later-drawn one may land slightly outside (hence not
  // asserted).
  size_t d = 4;
  auto vc = make_clayton_dvine(d, 3.0);
  vc.set_var_types({ "c", "c", "d", "d" }); // variables 3 and 4 discrete
  auto order = vc.get_order();              // [1, 2, 3, 4]
  size_t cv_first = order[d - 1] - 1;       // variable 4, drawn first
  size_t freevar = order[0] - 1;            // variable 1, drawn last
  size_t n = 50000;

  double a3 = 0.5, b3 = 0.7, a4 = 0.6, b4 = 0.8;
  // u_cond columns: [F(x3), F(x4), F(x3^-), F(x4^-)] (values then left-limits)
  Eigen::MatrixXd row(1, 4);
  row << b3, b4, a3, a4;
  Eigen::MatrixXd U;
  EXPECT_NO_THROW(
    U = vc.simulate_conditional(row.replicate(n, 1), false, 1, { 5 }));

  // the first-drawn discrete conditioning column is reproduced within its atom
  EXPECT_TRUE((U.col(cv_first).array() >= a4 - 1e-9).all());
  EXPECT_TRUE((U.col(cv_first).array() <= b4 + 1e-9).all());

  // conditioned draws are correct: the free mean matches a large unconditional
  // sample filtered to both atoms
  double cond_mean = U.col(freevar).mean();
  auto big = vc.simulate(500000, false, 1, { 6 });
  size_t cv3 = order[d - 2] - 1;
  double sum = 0.0, cnt = 0.0;
  for (int i = 0; i < big.rows(); ++i) {
    if ((big(i, cv3) >= a3) && (big(i, cv3) <= b3) &&
        (big(i, cv_first) >= a4) && (big(i, cv_first) <= b4)) {
      sum += big(i, freevar);
      cnt += 1.0;
    }
  }
  EXPECT_NEAR(cond_mean, sum / cnt, 0.03);
}

TEST_F(VinecopTest, conditional_select_places_conditioning_set_at_tail)
{
  size_t d = 6;
  auto data = make_clayton_dvine(d, 3.0).simulate(2000, false, 1, { 1 });
  auto fam = std::vector<BicopFamily>{ BicopFamily::gaussian,
                                       BicopFamily::clayton,
                                       BicopFamily::gumbel,
                                       BicopFamily::frank };
  std::vector<std::vector<size_t>> sets{ { 3 }, { 2, 5 }, { 1, 4, 6 } };
  for (const auto& B : sets) {
    FitControlsVinecop c;
    c.set_family_set(fam);
    c.set_conditioning_set(B);
    Vinecop vc(d);
    vc.select(data, c);

    // the conditioning set must be the tail of the order (drawn first)
    auto order = vc.get_order();
    std::vector<size_t> tail(order.end() - static_cast<ptrdiff_t>(B.size()),
                             order.end());
    std::sort(tail.begin(), tail.end());
    auto b_sorted = B;
    std::sort(b_sorted.begin(), b_sorted.end());
    EXPECT_EQ(tail, b_sorted);
    // and the produced structure is a valid R-vine (checked on construction)
    EXPECT_NO_THROW(RVineStructure(vc.get_matrix()));
  }
}

TEST_F(VinecopTest, conditional_select_empty_equals_plain)
{
  size_t d = 6;
  auto data = make_clayton_dvine(d, 3.0).simulate(2000, false, 1, { 1 });
  auto fam =
    std::vector<BicopFamily>{ BicopFamily::gaussian, BicopFamily::clayton };
  FitControlsVinecop c;
  c.set_family_set(fam);
  Vinecop a(d);
  a.select(data, c);

  FitControlsVinecop c2;
  c2.set_family_set(fam);
  c2.set_conditioning_set({}); // empty -> ordinary selection
  Vinecop b(d);
  b.select(data, c2);

  EXPECT_TRUE(a.get_matrix() == b.get_matrix());
  EXPECT_NEAR(a.get_loglik(), b.get_loglik(), 1e-10);
}

TEST_F(VinecopTest, conditional_select_then_simulate_conditional)
{
  size_t d = 5;
  auto data = make_clayton_dvine(d, 3.0).simulate(2000, false, 1, { 3 });
  FitControlsVinecop c;
  c.set_family_set(
    { BicopFamily::gaussian, BicopFamily::clayton, BicopFamily::frank });
  std::vector<size_t> B{ 2, 4 };
  c.set_conditioning_set(B);
  Vinecop vc(d);
  vc.select(data, c);

  // condition on the tail (== B); the tail columns must be reproduced exactly
  auto order = vc.get_order();
  Eigen::MatrixXd uc = Eigen::MatrixXd::Constant(400, B.size(), 0.4);
  auto U = vc.simulate_conditional(uc, false, 1, { 9 });
  for (size_t i = 0; i < B.size(); ++i) {
    size_t col = order[d - B.size() + i] - 1;
    EXPECT_TRUE(
      all_close(U.col(col), Eigen::VectorXd::Constant(400, 0.4), 1e-9, 1e-9));
  }
}

TEST_F(VinecopTest, reorient_preserves_model)
{
  size_t d = 6;
  // rotated Clayton D-vine (asymmetric + rotations) so the flip logic in
  // reorient() is genuinely exercised
  auto pcs = Vinecop::make_pair_copula_store(d);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (size_t t = 0; t < d - 1; ++t)
    for (auto& pc : pcs[t])
      pc = Bicop(BicopFamily::clayton, (t % 2 == 0) ? 90 : 270, par);
  Vinecop model(DVineStructure(std::vector<size_t>{ 1, 2, 3, 4, 5, 6 }), pcs);
  auto data = model.simulate(500, false, 1, { 7 });

  auto fam = std::vector<BicopFamily>{ BicopFamily::gaussian,
                                       BicopFamily::clayton,
                                       BicopFamily::gumbel };
  Vinecop vc(d);
  {
    FitControlsVinecop c;
    c.set_family_set(fam);
    vc.select(data, c);
  }
  auto pdf0 = vc.pdf(data);
  auto w = tools_stats::simulate_uniform(50, d, false, { 71 });

  size_t feasible = 0;
  std::vector<std::vector<size_t>> candidates{
    { 1 },    { 2 },    { 3 },    { 4 },       { 5 },      { 6 },
    { 1, 2 }, { 2, 5 }, { 3, 6 }, { 1, 4, 6 }, { 2, 3, 5 }
  };
  for (const auto& B : candidates) {
    Vinecop vr = vc;
    try {
      vr.reorient(B);
    } catch (const std::exception&) {
      continue; // infeasible for this structure -> allowed
    }
    ++feasible;
    // reorient is a pure relabeling: the density is unchanged (up to fp)
    EXPECT_TRUE(all_close(vr.pdf(data), pdf0, 1e-6, 1e-6));
    EXPECT_TRUE(all_close(
      vc.inverse_rosenblatt(w, B), vr.inverse_rosenblatt(w), 1e-10, 1e-10));
    EXPECT_TRUE(all_close(vc.rosenblatt(data, B, 1, false),
                          vr.rosenblatt(data, 1, false),
                          1e-10,
                          1e-10));
    // the conditioning set is now the tail of the order
    auto o = vr.get_order();
    std::vector<size_t> tail(o.end() - static_cast<ptrdiff_t>(B.size()),
                             o.end());
    std::sort(tail.begin(), tail.end());
    auto bs = B;
    std::sort(bs.begin(), bs.end());
    EXPECT_EQ(tail, bs);
    EXPECT_NO_THROW(RVineStructure(vr.get_matrix()));
  }
  EXPECT_GT(feasible, 0u); // at least some sets are admissible tails

  // validation throws
  EXPECT_ANY_THROW(vc.reorient({}));                   // empty
  EXPECT_ANY_THROW(vc.reorient({ 1, 2, 3, 4, 5, 6 })); // |B| == d
  EXPECT_ANY_THROW(vc.reorient({ 7 }));                // out of range
  EXPECT_ANY_THROW(vc.reorient({ 2, 2 }));             // duplicate
}

// A truncated model is relabeled like any other: the peeling works from the top
// stored tree down, so only the stored trees are permuted and the truncation
// level is preserved.
TEST(VinecopReorient, handles_truncated_vines)
{
  const size_t d = 6;
  auto model = make_clayton_dvine(d, 3.0, 90);
  model.truncate(2);
  ASSERT_EQ(model.get_trunc_lvl(), 2u);
  ASSERT_EQ(model.get_all_pair_copulas().size(), 2u);

  const auto data = tools_stats::simulate_uniform(200, d, false, { 7 });
  const auto pdf0 = model.pdf(data);
  const auto w = tools_stats::simulate_uniform(20, d, false, { 71 });

  size_t feasible = 0;
  const std::vector<std::vector<size_t>> candidates{ { 1 },    { 3 },
                                                     { 6 },    { 1, 2 },
                                                     { 2, 5 }, { 1, 3, 6 } };
  for (const auto& B : candidates) {
    Vinecop vr = model;
    try {
      vr.reorient(B);
    } catch (const std::exception&) {
      continue; // inadmissible tail -> allowed, as for a full vine
    }
    ++feasible;
    EXPECT_EQ(vr.get_trunc_lvl(), 2u);
    EXPECT_EQ(vr.get_all_pair_copulas().size(), 2u);
    EXPECT_TRUE(all_close(vr.pdf(data), pdf0, 1e-9, 1e-9));
    EXPECT_TRUE(all_close(
      model.inverse_rosenblatt(w, B), vr.inverse_rosenblatt(w), 1e-10, 1e-10));
    auto o = vr.get_order();
    std::vector<size_t> tail(o.end() - static_cast<ptrdiff_t>(B.size()),
                             o.end());
    std::sort(tail.begin(), tail.end());
    auto bs = B;
    std::sort(bs.begin(), bs.end());
    EXPECT_EQ(tail, bs);
    EXPECT_NO_THROW(RVineStructure(vr.get_matrix()));
  }
  EXPECT_GT(feasible, 0u);
}

// A vine truncated at zero has no trees to peel: every order represents the
// same independence model, so the tail is placed by permuting the diagonal.
TEST(VinecopReorient, handles_a_vine_truncated_at_zero)
{
  const size_t d = 5;
  const Vinecop model(d);
  ASSERT_EQ(model.get_trunc_lvl(), 0u);

  const std::vector<size_t> B{ 2, 4 };
  Vinecop vr = model;
  ASSERT_NO_THROW(vr.reorient(B));
  EXPECT_EQ(vr.get_order(), std::vector<size_t>({ 1, 3, 5, 2, 4 }));
  EXPECT_EQ(vr.get_trunc_lvl(), 0u);
  EXPECT_TRUE(vr.get_all_pair_copulas().empty());

  const auto data = tools_stats::simulate_uniform(20, d, false, { 1 });
  EXPECT_TRUE(all_close(vr.pdf(data), Eigen::VectorXd::Ones(data.rows())));

  Eigen::MatrixXd u_cond(20, B.size());
  u_cond.col(0).setConstant(0.3);
  u_cond.col(1).setConstant(0.7);
  const auto simulated = vr.simulate_conditional(u_cond, false, 1, { 3 });
  EXPECT_TRUE(all_close(simulated.col(1), u_cond.col(0))); // variable 2
  EXPECT_TRUE(all_close(simulated.col(3), u_cond.col(1))); // variable 4
}

// Conditioning-aware selection under truncation: the fit stops at the
// truncation level and the re-orientation that follows works from there.
TEST(VinecopConditionalSelection, supports_truncation)
{
  const size_t d = 6;
  const auto data = make_clayton_dvine(d, 2.0).simulate(400, false, 1, { 5 });
  const std::vector<size_t> B{ 1, 5 };

  for (const bool select_trunc : { false, true }) {
    FitControlsVinecop controls;
    controls.set_family_set({ BicopFamily::gaussian, BicopFamily::clayton });
    controls.set_conditioning_set(B);
    if (select_trunc) {
      controls.set_select_trunc_lvl(true);
    } else {
      controls.set_trunc_lvl(2);
    }

    Vinecop vc(d);
    ASSERT_NO_THROW(vc.select(data, controls));
    EXPECT_EQ(vc.get_all_pair_copulas().size(), vc.get_trunc_lvl());
    if (!select_trunc) {
      EXPECT_EQ(vc.get_trunc_lvl(), 2u);
    }

    auto o = vc.get_order();
    std::vector<size_t> tail(o.end() - static_cast<ptrdiff_t>(B.size()),
                             o.end());
    std::sort(tail.begin(), tail.end());
    EXPECT_EQ(tail, B);

    Eigen::MatrixXd u_cond(30, B.size());
    u_cond.col(0).setConstant(0.35);
    u_cond.col(1).setConstant(0.65);
    const auto simulated = vc.simulate_conditional(u_cond, false, 1, { 9 });
    for (size_t j = 0; j < B.size(); ++j) {
      const auto col = static_cast<long>(o[d - B.size() + j] - 1);
      EXPECT_TRUE(
        all_close(simulated.col(col), u_cond.col(static_cast<long>(j))));
    }
  }
}

TEST_F(VinecopTest, conditional_select_throws)
{
  size_t d = 5;
  auto data = make_clayton_dvine(d, 2.0).simulate(500, false, 1, { 1 });
  {
    // index out of range
    FitControlsVinecop c;
    c.set_conditioning_set({ 6 });
    Vinecop vc(d);
    EXPECT_ANY_THROW(vc.select(data, c));
  }
  {
    // |B| == d (no free variable)
    FitControlsVinecop c;
    c.set_conditioning_set({ 1, 2, 3, 4, 5 });
    Vinecop vc(d);
    EXPECT_ANY_THROW(vc.select(data, c));
  }
  {
    // random tree_algorithm is incompatible with conditioning-aware selection
    FitControlsVinecop c;
    c.set_conditioning_set({ 2 });
    c.set_tree_algorithm("random_weighted");
    c.set_seeds({ 1 });
    Vinecop vc(d);
    EXPECT_ANY_THROW(vc.select(data, c));
  }
}

TEST_F(VinecopTest, scores_stepwise)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  auto uu = vinecop.simulate(100, false, 1, { 1 });

  auto J = vinecop.hessian(uu, true);
  EXPECT_TRUE(J.isUpperTriangular());

  // hessian() is the observation-average of hessian_full()
  auto H = vinecop.hessian_full(uu, true);
  Eigen::MatrixXd J_manual(J.rows(), J.cols());
  size_t d = vinecop.get_dim();
  size_t trunc_lvl = vinecop.get_trunc_lvl();
  size_t ipar = 0;
  for (size_t t = 0; t < trunc_lvl; t++) {
    for (size_t e = 0; e < d - 1 - t; e++) {
      size_t np = static_cast<size_t>(
        vinecop.get_pair_copula(t, e).get_parameters().size());
      for (size_t p = 0; p < np; p++) {
        J_manual.row(ipar++) = H(t, e)[p].colwise().mean();
      }
    }
  }
  EXPECT_TRUE(all_close(J, J_manual, 1e-10, 1e-10));

  // scores_cov() is the mean-centered covariance of scores()
  auto s = vinecop.scores(uu, true);
  Eigen::MatrixXd sc = s.rowwise() - s.colwise().mean();
  Eigen::MatrixXd I_manual =
    (sc.adjoint() * sc) / static_cast<double>(s.rows());
  EXPECT_TRUE(all_close(vinecop.scores_cov(uu, true), I_manual, 1e-10, 1e-10));
}

TEST_F(VinecopTest, scores_joint)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  auto uu = vinecop.simulate(100, false, 1, { 1 });
  auto J = vinecop.hessian(uu, false);
  auto I = vinecop.scores_cov(uu, false);

  EXPECT_FALSE(J.isUpperTriangular());

  Eigen::MatrixXd Jinv = J.triangularView<Eigen::Upper>()
                           .solve(Eigen::MatrixXd::Identity(J.cols(), J.cols()))
                           .triangularView<Eigen::Upper>();
  // std::cout << J << std::endl << std::endl;
  // std::cout << Jinv * I << std::endl;
  I = vinecop.scores_cov(uu, true);
  // std::cout << (Jinv * I * Jinv.transpose() /
  // u.rows()).diagonal().cwiseSqrt()
  //           << std::endl
  //           << std::endl;
  // std::cout << vinecop.str() << std::endl;
}

// The analytic full gradient (step_wise = false) must match VineCopula's
// RVineGrad and the analytic joint Hessian must match RVineHessian, both
// element-wise (analytic vs analytic). Requires VineCopula >= 2.6.2, whose
// #101 fix makes RVineHessian depend on the data so it can be used directly
// (not finite-differenced). The R oracle (test_vinecop_derivatives.R)
// hardcodes the same model and reorders RVineGrad and RVineHessian into
// vinecopulib's parameter order.
//
// The model uses only 0/180 families: vinecopulib and VineCopula do not share
// the same structure / swapped-rotation conventions, so a 90/270 vine cannot
// be built identically in both (see the R oracle's header). 0/180 families
// are symmetric under the argument swap that distinguishes the conventions.
// Vine-level rotations are validated against brute force instead
// (hessian_matches_brute_force, below), and rotated bicop derivatives against
// BiCopDeriv* (ParBicopTest).
TEST(VinecopDerivatives, full_scores_match_RVineGrad_and_RVineHessian)
{
  SKIP_WITHOUT_RSCRIPT();
  test_r_parity::RTempDir dir;
  dir.run(TEST_VINECOP_DERIV, "300");
  auto mat = tools_eigen::read_matxs(dir.file("temp_vderiv_matrix").c_str())
               .colwise()
               .reverse()
               .eval();
  auto u = tools_eigen::read_matxd(dir.file("temp_vderiv_data").c_str());
  auto grads = tools_eigen::read_matxd(dir.file("temp_vderiv_grad").c_str());
  auto hess_r = tools_eigen::read_matxd(dir.file("temp_vderiv_hess").c_str());

  // must stay in sync with the R script (0/180 families only; pair copula
  // (t, e) <-> VineCopula matrix cell [d - t, e + 1])
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 0, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 0, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  pcs[3][0] = Bicop(BicopFamily::clayton, 0, P1(2.5));
  Vinecop vc(RVineStructure(mat), pcs);
  ASSERT_EQ(static_cast<size_t>(vc.get_npars()), grads.rows());

  // full gradient vs RVineGrad: sum over observations, and a single row
  Eigen::MatrixXd s = vc.scores(u, false, 1);
  EXPECT_TRUE(all_close(
    s.colwise().sum().transpose().eval(), grads.col(0).eval(), 1e-4, 1e-4));
  EXPECT_TRUE(
    all_close(s.row(0).transpose().eval(), grads.col(1).eval(), 1e-4, 1e-4));

  // gradient() is the observation-average of the scores
  EXPECT_TRUE(all_close(vc.gradient(u, false, 1),
                        s.colwise().mean().transpose().eval(),
                        1e-12,
                        1e-12));

  // scores_full(keep_all) returns the same scores plus non-empty caches
  auto full = vc.scores_full(u, false, 1, true);
  EXPECT_TRUE(all_close(full.scores, s, 1e-12, 1e-12));
  EXPECT_EQ(full.pdf_edges.get_dim(), vc.get_dim());
  EXPECT_EQ(full.logpdf_deriv_pars(0, 0).size(), 1u);
  EXPECT_EQ(full.logpdf_deriv_pars(0, 0)[0].size(), u.rows());

  // deterministic across the number of threads
  EXPECT_TRUE(all_close(vc.scores(u, false, 3), s, 1e-12, 1e-12));

  // joint Hessian vs RVineHessian, element-wise. n * hessian() is the summed
  // observed information, as RVineHessian returns.
  Eigen::MatrixXd H = vc.hessian(u, false, 1) * static_cast<double>(u.rows());
  EXPECT_TRUE(all_close(H, hess_r, 1e-4, 1e-4));
  EXPECT_TRUE(all_close(H, H.transpose().eval(), 1e-10, 1e-10)); // symmetric
}

// The analytic joint Hessian must match brute-force central finite
// differences of the whole-vine log-likelihood — the unforgiving ground
// truth, and (unlike the R oracle) valid for 90/270 rotations and
// truncation, which VineCopula's RVineHessian cannot provide.
TEST(VinecopDerivatives, hessian_matches_brute_force)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 90, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 270, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(-0.3));
  pcs[3][0] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  Vinecop vc(DVineStructure({ 1, 2, 3, 4, 5 }), pcs);
  auto u = vc.simulate(150, false, 1, { 42, 1, 2 });

  size_t npars = static_cast<size_t>(vc.get_npars());
  Eigen::MatrixXd H = vc.hessian(u, false, 1) * static_cast<double>(u.rows());

  // flat (t, e, p) list in hessian()'s column order
  std::vector<std::array<size_t, 3>> pl;
  for (size_t t = 0; t < 4; ++t) {
    for (size_t e = 0; e < 4 - t; ++e) {
      size_t np = static_cast<size_t>(pcs[t][e].get_parameters().size());
      for (size_t p = 0; p < np; ++p) {
        pl.push_back({ t, e, p });
      }
    }
  }
  auto structure = vc.get_rvine_structure();
  auto loglik = [&](const std::vector<std::vector<Bicop>>& pp) {
    return Vinecop(structure, pp).pdf(u).array().max(1e-300).log().sum();
  };
  auto bump = [&](std::vector<std::vector<Bicop>> pp, size_t a, double h) {
    auto pr = pp[pl[a][0]][pl[a][1]].get_parameters();
    pr(pl[a][2]) += h;
    pp[pl[a][0]][pl[a][1]].set_parameters(pr);
    return pp;
  };
  double base = loglik(pcs);
  double h = 1e-4;
  for (size_t a = 0; a < npars; ++a) {
    for (size_t bb = a; bb < npars; ++bb) {
      double fd;
      if (a == bb) {
        fd = (loglik(bump(pcs, a, h)) - 2 * base + loglik(bump(pcs, a, -h))) /
             (h * h);
      } else {
        fd = (loglik(bump(bump(pcs, a, h), bb, h)) -
              loglik(bump(bump(pcs, a, h), bb, -h)) -
              loglik(bump(bump(pcs, a, -h), bb, h)) +
              loglik(bump(bump(pcs, a, -h), bb, -h))) /
             (4 * h * h);
      }
      EXPECT_NEAR(H(a, bb), fd, 1e-3 * (1.0 + std::abs(fd)))
        << "H(" << a << ", " << bb << ")";
    }
  }
  EXPECT_TRUE(all_close(H, H.transpose().eval(), 1e-10, 1e-10));
  EXPECT_TRUE(all_close(vc.hessian(u, false, 3) * static_cast<double>(u.rows()),
                        H,
                        1e-12,
                        1e-12)); // threading determinism
}

// The analytic step-wise Hessian must match central finite differences of
// the step-wise scores (the quantity the finite-difference path produced),
// including a two-parameter (Student) edge's own-edge cross block.
TEST(VinecopDerivatives, stepwise_hessian_matches_fd_of_scores)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 90, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 270, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(-0.3));
  pcs[3][0] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  auto structure = DVineStructure({ 1, 2, 3, 4, 5 });
  Vinecop vc(structure, pcs);
  auto u = vc.simulate(120, false, 1, { 5, 6, 7 });

  auto Ha = vc.hessian_full(u, true, 1);

  std::vector<std::array<size_t, 3>> pl;
  for (size_t t = 0; t < 4; ++t) {
    for (size_t e = 0; e < 4 - t; ++e) {
      size_t np = static_cast<size_t>(pcs[t][e].get_parameters().size());
      for (size_t p = 0; p < np; ++p) {
        pl.push_back({ t, e, p });
      }
    }
  }
  auto bump = [&](std::vector<std::vector<Bicop>> pp, size_t a, double h) {
    auto pr = pp[pl[a][0]][pl[a][1]].get_parameters();
    pr(pl[a][2]) += h;
    pp[pl[a][0]][pl[a][1]].set_parameters(pr);
    return pp;
  };
  for (size_t a = 0; a < pl.size(); ++a) {
    double base = std::abs(pcs[pl[a][0]][pl[a][1]].get_parameters()(pl[a][2]));
    double h = 1e-4 * std::max(1.0, base);
    Eigen::MatrixXd sp = Vinecop(structure, bump(pcs, a, h)).scores(u, true, 1);
    Eigen::MatrixXd sm =
      Vinecop(structure, bump(pcs, a, -h)).scores(u, true, 1);
    Eigen::MatrixXd fd = (sp - sm) / (2 * h);
    EXPECT_TRUE(all_close(Ha(pl[a][0], pl[a][1])[pl[a][2]], fd, 1e-3, 1e-3))
      << "param " << a;
  }
  // threading determinism
  auto Ha3 = vc.hessian_full(u, true, 3);
  for (size_t t = 0; t < 4; ++t) {
    for (size_t e = 0; e < 4 - t; ++e) {
      for (size_t p = 0; p < Ha(t, e).size(); ++p) {
        EXPECT_TRUE(all_close(Ha(t, e)[p], Ha3(t, e)[p], 1e-12, 1e-12));
      }
    }
  }
}

// A vine with a family outside the analytic set (bb1) still gets an analytic
// joint Hessian: the second-order cascade calls the pair copula's derivative
// leaves, which fall back to finite differences inside the leaf for bb1.
TEST(VinecopDerivatives, hessian_bb1_edge_matches_brute_force)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(3);
  Eigen::VectorXd bb1_par(2);
  bb1_par << 1.0, 1.5;
  pcs[0][0] = Bicop(BicopFamily::bb1, 0, bb1_par);
  pcs[0][1] = Bicop(BicopFamily::gaussian, 0, P1(0.5));
  pcs[1][0] = Bicop(BicopFamily::gumbel, 0, P1(1.5));
  Vinecop vc(DVineStructure({ 1, 2, 3 }), pcs);
  auto u = vc.simulate(100, false, 1, { 3 });

  size_t npars = static_cast<size_t>(vc.get_npars());
  Eigen::MatrixXd H = vc.hessian(u, false, 1) * static_cast<double>(u.rows());
  EXPECT_TRUE(H.allFinite());

  std::vector<std::array<size_t, 3>> pl;
  for (size_t t = 0; t < 2; ++t) {
    for (size_t e = 0; e < 2 - t; ++e) {
      size_t np = static_cast<size_t>(pcs[t][e].get_parameters().size());
      for (size_t p = 0; p < np; ++p) {
        pl.push_back({ t, e, p });
      }
    }
  }
  auto structure = vc.get_rvine_structure();
  auto loglik = [&](const std::vector<std::vector<Bicop>>& pp) {
    return Vinecop(structure, pp).pdf(u).array().max(1e-300).log().sum();
  };
  auto bump = [&](std::vector<std::vector<Bicop>> pp, size_t a, double hh) {
    auto pr = pp[pl[a][0]][pl[a][1]].get_parameters();
    pr(pl[a][2]) += hh;
    pp[pl[a][0]][pl[a][1]].set_parameters(pr);
    return pp;
  };
  double base = loglik(pcs);
  double h = 1e-4;
  // bb1 now has analytic derivative leaves (like the other families), so the
  // cascade agrees with brute force to the finite-difference truncation error.
  for (size_t a = 0; a < npars; ++a) {
    for (size_t bb = a; bb < npars; ++bb) {
      double fd;
      if (a == bb) {
        fd = (loglik(bump(pcs, a, h)) - 2 * base + loglik(bump(pcs, a, -h))) /
             (h * h);
      } else {
        fd = (loglik(bump(bump(pcs, a, h), bb, h)) -
              loglik(bump(bump(pcs, a, h), bb, -h)) -
              loglik(bump(bump(pcs, a, -h), bb, h)) +
              loglik(bump(bump(pcs, a, -h), bb, -h))) /
             (4 * h * h);
      }
      EXPECT_NEAR(H(a, bb), fd, 1e-3 * (1.0 + std::abs(fd)))
        << "H(" << a << ", " << bb << ")";
    }
  }
}

// scores() and hessian() reject models with nonparametric pair copulas
// (differentiating w.r.t. the interpolation grid is meaningless; mirrors
// the guard added in the perf-overhaul PR).
TEST(VinecopDerivatives, reject_nonparametric)
{
  auto data = tools_stats::simulate_uniform(200, 3);
  Vinecop vc(DVineStructure({ 1, 2, 3 }));
  vc.select(data, FitControlsVinecop({ BicopFamily::tll }));
  EXPECT_THROW(vc.scores(data, true), std::runtime_error);
  EXPECT_THROW(vc.scores(data, false), std::runtime_error);
  EXPECT_THROW(vc.hessian(data, true), std::runtime_error);
  EXPECT_THROW(vc.hessian(data, false), std::runtime_error);
  EXPECT_THROW(vc.gradient(data, false), std::runtime_error);
}

// The step-wise scores must match per-edge finite differences of the edge
// log-densities, with the edge arguments reconstructed from pdf_full().
TEST(VinecopDerivatives, stepwise_scores_match_per_edge_reference)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 90, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 270, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  pcs[3][0] = Bicop(BicopFamily::gaussian, 0, P1(-0.3));
  Vinecop vc(DVineStructure({ 1, 2, 3, 4, 5 }), pcs);

  auto u = vc.simulate(200, false, 1, { 7, 8, 9 });
  Eigen::MatrixXd s = vc.scores(u, true, 1);

  // reconstruct each edge's arguments from the cached h-functions
  auto full = vc.pdf_full(u, 1, true);
  auto structure = vc.get_rvine_structure();
  auto order = structure.get_order();
  size_t d = vc.get_dim();
  size_t ipar = 0;
  for (size_t t = 0; t < d - 1; ++t) {
    for (size_t e = 0; e < d - 1 - t; ++e) {
      size_t m = structure.min_array(t, e);
      bool arg2_is_h2 = (m == structure.struct_array(t, e, true));
      Eigen::MatrixXd u_e(u.rows(), 2);
      if (t == 0) {
        u_e.col(0) = u.col(order[e] - 1);
        u_e.col(1) = u.col(order[m - 1] - 1);
      } else {
        u_e.col(0) = full.hfunc2(t - 1, e);
        u_e.col(1) =
          arg2_is_h2 ? full.hfunc2(t - 1, m - 1) : full.hfunc1(t - 1, m - 1);
      }

      // old-style central differences of the edge log-density
      Bicop ec = vc.get_pair_copula(t, e);
      auto pars = ec.get_parameters();
      auto ub = ec.get_parameters_upper_bounds();
      auto lb = ec.get_parameters_lower_bounds();
      for (size_t p = 0; p < static_cast<size_t>(pars.size()); ++p) {
        auto pars_tmp = pars;
        double eps = 0;
        pars_tmp(p) = std::min(pars(p) + 1e-3, ub(p));
        eps += pars_tmp(p) - pars(p);
        ec.set_parameters(pars_tmp);
        Eigen::VectorXd f1 = ec.pdf(u_e).array().max(1e-20).log();
        pars_tmp(p) = std::max(pars(p) - 1e-3, lb(p));
        eps -= pars_tmp(p) - pars(p);
        ec.set_parameters(pars_tmp);
        Eigen::VectorXd f2 = ec.pdf(u_e).array().max(1e-20).log();
        ec.set_parameters(pars);
        EXPECT_TRUE(
          all_close(s.col(ipar).eval(), ((f1 - f2) / eps).eval(), 1e-3, 1e-4))
          << "tree " << t << ", edge " << e << ", parameter " << p;
        ipar++;
      }
    }
  }
  EXPECT_EQ(ipar, static_cast<size_t>(vc.get_npars()));

  EXPECT_TRUE(all_close(vc.scores(u, true, 3), s, 1e-12, 1e-12));
}

// The analytic cascade must agree with brute-force finite differences of
// the whole-vine log-likelihood on a truncated vine with rotated pair
// copulas (covers 90/270 rotations and truncation, which the R oracle
// cannot: VineCopula's RVineGrad is broken for rotated families).
TEST_F(VinecopTest, full_scores_match_brute_force_on_truncated_vine)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);
  auto uu = vinecop.simulate(100, false, 1, { 17 });

  Eigen::VectorXd grad =
    vinecop.scores(uu, false, 1).colwise().sum().transpose();

  auto loglik_of = [&](const std::vector<std::vector<Bicop>>& pcs) {
    Vinecop v(model_matrix, pcs);
    return v.pdf(uu).array().max(1e-300).log().sum();
  };
  size_t ipar = 0;
  for (size_t t = 0; t < 3; ++t) {
    for (size_t e = 0; e < 7 - 1 - t; ++e) {
      double h = 1e-5 * 3.0;
      auto pp = pair_copulas;
      auto pm = pair_copulas;
      pp[t][e].set_parameters(Eigen::VectorXd::Constant(1, 3.0 + h));
      pm[t][e].set_parameters(Eigen::VectorXd::Constant(1, 3.0 - h));
      double fd = (loglik_of(pp) - loglik_of(pm)) / (2 * h);
      EXPECT_NEAR(grad(ipar), fd, 1e-4 * (1.0 + std::abs(fd)))
        << "tree " << t << ", edge " << e;
      ipar++;
    }
  }
  EXPECT_EQ(ipar, static_cast<size_t>(vinecop.get_npars()));
}

// Vines with families outside the analytic set or with discrete variables
// keep working through the fallback paths.
TEST(VinecopDerivatives, fallback_paths_run)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  // a bb1 edge is parametric, so scores/hessian still take the analytic
  // cascade — the bb1 pair copula's derivative leaves simply fall back to
  // finite differences internally
  auto pcs = Vinecop::make_pair_copula_store(3);
  Eigen::VectorXd bb1_par(2);
  bb1_par << 1.0, 1.5;
  pcs[0][0] = Bicop(BicopFamily::bb1, 0, bb1_par);
  pcs[0][1] = Bicop(BicopFamily::gaussian, 0, P1(0.5));
  pcs[1][0] = Bicop(BicopFamily::gumbel, 0, P1(1.5));
  Vinecop vc(DVineStructure({ 1, 2, 3 }), pcs);
  auto u = vc.simulate(50, false, 1, { 3 });

  Eigen::MatrixXd s_full = vc.scores(u, false, 1);
  Eigen::MatrixXd s_sw = vc.scores(u, true, 1);
  EXPECT_EQ(s_full.cols(), 4);
  EXPECT_EQ(s_sw.cols(), 4);
  EXPECT_TRUE(s_full.allFinite());
  EXPECT_TRUE(s_sw.allFinite());
  // full and step-wise gradients agree for the last tree's parameters (no
  // deeper edges depend on them), giving the fallback a value check
  EXPECT_TRUE(all_close(s_full.col(3).eval(), s_sw.col(3).eval(), 1e-3, 1e-4));

  // discrete variables keep the finite-difference paths
  Vinecop vc_disc(DVineStructure({ 1, 2, 3 }), pcs);
  vc_disc.set_var_types({ "c", "d", "c" });
  Eigen::MatrixXd u4(u.rows(), 4);
  u4.leftCols(3) = u;
  u4.col(3) = (u.col(1).array() * 0.9).matrix();
  Eigen::MatrixXd s_disc = vc_disc.scores(u4, true, 1);
  EXPECT_EQ(s_disc.cols(), 4);
  EXPECT_TRUE(s_disc.allFinite());
  Eigen::MatrixXd s_disc_full = vc_disc.scores(u4, false, 1);
  EXPECT_TRUE(s_disc_full.allFinite());
}

// The per-observation-parameter overloads must (a) reduce to the fixed-
// parameter path when a single parameter vector is broadcast to every row and
// (b) for heterogeneous per-observation parameters, reproduce row i from a
// vine built with observation i's parameter vector — the unforgiving ground
// truth that the whole per-observation cascade is correct, for both
// step_wise settings and for scores and the (per-observation) Hessian.
TEST(VinecopDerivatives, per_obs_scores_and_hessian)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 90, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 270, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(-0.3));
  pcs[3][0] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  Vinecop vc(DVineStructure({ 1, 2, 3, 4, 5 }), pcs);
  auto u = vc.simulate(50, false, 1, { 42, 1, 2 });
  const Eigen::Index n = u.rows();
  const size_t npars = static_cast<size_t>(vc.get_npars());

  // flatten the stored parameters (and their bounds) in the (t, e, p) column
  // order that scores() and the `parameters` argument share
  Eigen::VectorXd theta0(npars), lb_flat(npars), ub_flat(npars);
  {
    size_t idx = 0;
    for (size_t t = 0; t < 4; ++t) {
      for (size_t e = 0; e < 4 - t; ++e) {
        auto pr = pcs[t][e].get_parameters();
        auto lb = pcs[t][e].get_parameters_lower_bounds();
        auto ub = pcs[t][e].get_parameters_upper_bounds();
        for (Eigen::Index p = 0; p < pr.size(); ++p) {
          theta0(idx) = pr(p);
          lb_flat(idx) = lb(p);
          ub_flat(idx) = ub(p);
          idx++;
        }
      }
    }
  }

  // (a) broadcasting the stored parameters reduces to the fixed-parameter path
  Eigen::MatrixXd P0 = theta0.transpose().replicate(n, 1);
  for (bool sw : { true, false }) {
    EXPECT_TRUE(
      all_close(vc.scores(u, P0, sw, 1), vc.scores(u, sw, 1), 1e-10, 1e-10));
    EXPECT_TRUE(all_close(
      vc.gradient(u, P0, sw, 1), vc.gradient(u, sw, 1), 1e-10, 1e-10));
    EXPECT_TRUE(
      all_close(vc.hessian(u, P0, sw, 1), vc.hessian(u, sw, 1), 1e-10, 1e-10));
    EXPECT_TRUE(all_close(
      vc.scores_cov(u, P0, sw, 1), vc.scores_cov(u, sw, 1), 1e-10, 1e-10));
    auto hf_pr = vc.hessian_full(u, P0, sw, 1);
    auto hf_fx = vc.hessian_full(u, sw, 1);
    for (size_t t = 0; t < 4; ++t) {
      for (size_t e = 0; e < 4 - t; ++e) {
        for (size_t p = 0; p < hf_fx(t, e).size(); ++p) {
          EXPECT_TRUE(all_close(hf_pr(t, e)[p], hf_fx(t, e)[p], 1e-10, 1e-10));
        }
      }
    }
  }

  // threading determinism on the per-observation path
  EXPECT_TRUE(
    all_close(vc.scores(u, P0, false, 3), vc.scores(u, P0, false, 1), 1e-12));
  EXPECT_TRUE(
    all_close(vc.scores(u, P0, true, 3), vc.scores(u, P0, true, 1), 1e-12));

  // (b) heterogeneous per-observation parameters (interior-clamped), compared
  // row by row against a vine rebuilt with that observation's parameters
  Eigen::MatrixXd P(n, npars);
  for (Eigen::Index i = 0; i < n; ++i) {
    double w = 0.85 + 0.30 * static_cast<double>(i % 4) / 3.0; // [0.85, 1.15]
    for (size_t k = 0; k < npars; ++k) {
      double lo = lb_flat(k) + 1e-2 * (ub_flat(k) - lb_flat(k));
      double hi = ub_flat(k) - 1e-2 * (ub_flat(k) - lb_flat(k));
      P(i, k) = std::min(std::max(theta0(k) * w, lo), hi);
    }
  }

  for (bool sw : { true, false }) {
    Eigen::MatrixXd s = vc.scores(u, P, sw, 1);
    auto hf = vc.hessian_full(u, P, sw, 1);
    // threading determinism with heterogeneous parameters
    EXPECT_TRUE(all_close(vc.scores(u, P, sw, 3), s, 1e-12));
    for (Eigen::Index i = 0; i < n; ++i) {
      auto pcs_i = pcs;
      size_t idx = 0;
      for (size_t t = 0; t < 4; ++t) {
        for (size_t e = 0; e < 4 - t; ++e) {
          Eigen::VectorXd pr = pcs_i[t][e].get_parameters();
          for (Eigen::Index p = 0; p < pr.size(); ++p) {
            pr(p) = P(i, idx++);
          }
          pcs_i[t][e].set_parameters(pr);
        }
      }
      Vinecop vc_i(vc.get_rvine_structure(), pcs_i);
      Eigen::MatrixXd ui = u.row(i);
      EXPECT_TRUE(
        all_close(s.row(i), vc_i.scores(ui, sw, 1).row(0), 1e-9, 1e-9))
        << "row " << i;
      auto hi = vc_i.hessian_full(ui, sw, 1);
      for (size_t t = 0; t < 4; ++t) {
        for (size_t e = 0; e < 4 - t; ++e) {
          for (size_t p = 0; p < hi(t, e).size(); ++p) {
            EXPECT_TRUE(
              all_close(hf(t, e)[p].row(i), hi(t, e)[p].row(0), 1e-9, 1e-9))
              << "row " << i << " edge (" << t << ", " << e << ") par " << p;
          }
        }
      }
    }
  }
}

// The per-observation path rejects discrete variables (its finite-difference
// fallback would mutate shared parameters) and validates the n x npars shape.
TEST(VinecopDerivatives, per_obs_rejects_discrete_and_bad_shape)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(3);
  pcs[0][0] = Bicop(BicopFamily::clayton, 0, P1(2.0));
  pcs[0][1] = Bicop(BicopFamily::gaussian, 0, P1(0.5));
  pcs[1][0] = Bicop(BicopFamily::gumbel, 0, P1(1.5));
  Vinecop vc(DVineStructure({ 1, 2, 3 }), pcs);
  auto u = vc.simulate(30, false, 1, { 5 });
  const size_t npars = static_cast<size_t>(vc.get_npars());

  Eigen::MatrixXd P(u.rows(), npars);
  {
    size_t idx = 0;
    for (size_t t = 0; t < 2; ++t) {
      for (size_t e = 0; e < 2 - t; ++e) {
        auto pr = pcs[t][e].get_parameters();
        for (Eigen::Index p = 0; p < pr.size(); ++p) {
          P.col(idx++).setConstant(pr(p));
        }
      }
    }
  }
  EXPECT_NO_THROW(vc.scores(u, P, false, 1));
  EXPECT_NO_THROW(vc.pdf(u, P));
  EXPECT_NO_THROW(vc.loglik(u, P));

  // wrong number of columns / rows
  EXPECT_ANY_THROW(vc.scores(u, P.leftCols(npars - 1)));
  EXPECT_ANY_THROW(vc.scores(u, P.topRows(u.rows() - 1)));
  EXPECT_ANY_THROW(vc.pdf(u, P.leftCols(npars - 1)));
  EXPECT_ANY_THROW(vc.loglik(u, P.topRows(u.rows() - 1)));

  // discrete variables are rejected on the per-observation path
  Vinecop vc_disc(DVineStructure({ 1, 2, 3 }), pcs);
  vc_disc.set_var_types({ "c", "d", "c" });
  Eigen::MatrixXd u4(u.rows(), 4);
  u4.leftCols(3) = u;
  u4.col(3) = (u.col(1).array() * 0.9).matrix();
  EXPECT_ANY_THROW(vc_disc.scores(u4, P, true, 1));
  EXPECT_ANY_THROW(vc_disc.hessian(u4, P, false, 1));
  EXPECT_ANY_THROW(vc_disc.pdf(u4, P));
  EXPECT_ANY_THROW(vc_disc.loglik(u4, P));
}

// Per-observation-parameter pdf/pdf_full/loglik must (a) reduce to the fixed
// path when a single parameter vector is broadcast to every row and (b) for
// heterogeneous per-observation parameters reproduce row i from a vine built
// with observation i's parameter vector (and loglik = sum of those log
// densities) — the same ground truth used for the per-observation scores.
TEST(VinecopDerivatives, per_obs_pdf_and_loglik)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 90, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 270, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(-0.3));
  pcs[3][0] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  Vinecop vc(DVineStructure({ 1, 2, 3, 4, 5 }), pcs);
  auto u = vc.simulate(50, false, 1, { 42, 1, 2 });
  const Eigen::Index n = u.rows();
  const size_t npars = static_cast<size_t>(vc.get_npars());

  Eigen::VectorXd theta0(npars), lb_flat(npars), ub_flat(npars);
  {
    size_t idx = 0;
    for (size_t t = 0; t < 4; ++t) {
      for (size_t e = 0; e < 4 - t; ++e) {
        auto pr = pcs[t][e].get_parameters();
        auto lb = pcs[t][e].get_parameters_lower_bounds();
        auto ub = pcs[t][e].get_parameters_upper_bounds();
        for (Eigen::Index p = 0; p < pr.size(); ++p) {
          theta0(idx) = pr(p);
          lb_flat(idx) = lb(p);
          ub_flat(idx) = ub(p);
          idx++;
        }
      }
    }
  }

  // (a) broadcasting the stored parameters reduces to the fixed-parameter path
  Eigen::MatrixXd P0 = theta0.transpose().replicate(n, 1);
  EXPECT_TRUE(all_close(vc.pdf(u, P0), vc.pdf(u), 1e-10, 1e-10));
  EXPECT_TRUE(all_close(vc.pdf_full(u, P0, 1).pdf, vc.pdf(u), 1e-10, 1e-10));
  double ll_fixed = vc.loglik(u);
  EXPECT_NEAR(vc.loglik(u, P0), ll_fixed, 1e-8 * (1.0 + std::abs(ll_fixed)));
  // threading determinism on the per-observation path
  EXPECT_TRUE(all_close(vc.pdf(u, P0, 3), vc.pdf(u, P0, 1), 1e-12));

  // (b) heterogeneous per-observation parameters (interior-clamped), compared
  // row by row against a vine rebuilt with that observation's parameters
  Eigen::MatrixXd P(n, npars);
  for (Eigen::Index i = 0; i < n; ++i) {
    double w = 0.85 + 0.30 * static_cast<double>(i % 4) / 3.0; // [0.85, 1.15]
    for (size_t k = 0; k < npars; ++k) {
      double lo = lb_flat(k) + 1e-2 * (ub_flat(k) - lb_flat(k));
      double hi = ub_flat(k) - 1e-2 * (ub_flat(k) - lb_flat(k));
      P(i, k) = std::min(std::max(theta0(k) * w, lo), hi);
    }
  }

  Eigen::VectorXd pdf_po = vc.pdf(u, P, 1);
  double ll_po = vc.loglik(u, P, 1);
  EXPECT_TRUE(
    all_close(vc.pdf(u, P, 3), pdf_po, 1e-12)); // threading determinism
  double ll_ref = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) {
    auto pcs_i = pcs;
    size_t idx = 0;
    for (size_t t = 0; t < 4; ++t) {
      for (size_t e = 0; e < 4 - t; ++e) {
        Eigen::VectorXd pr = pcs_i[t][e].get_parameters();
        for (Eigen::Index p = 0; p < pr.size(); ++p) {
          pr(p) = P(i, idx++);
        }
        pcs_i[t][e].set_parameters(pr);
      }
    }
    Vinecop vc_i(vc.get_rvine_structure(), pcs_i);
    double pdf_i = vc_i.pdf(u.row(i))(0);
    EXPECT_NEAR(pdf_po(i), pdf_i, 1e-9 * (1.0 + std::abs(pdf_i)))
      << "row " << i;
    ll_ref += std::log(pdf_i);
  }
  EXPECT_NEAR(ll_po, ll_ref, 1e-8 * (1.0 + std::abs(ll_ref)));
}

TEST_F(VinecopTest, aic_bic_are_correct)
{
  int d = 7;
  auto data = tools_stats::simulate_uniform(100, 7, false, { 1 });
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
  true_model.select(data);

  FitControlsVinecop controls({ BicopFamily::gaussian, BicopFamily::tll });
  complex_model.select(data);
  ASSERT_NEAR(complex_model.get_aic(), complex_model.aic(data), 1e-2);
  ASSERT_NEAR(complex_model.get_bic(), complex_model.bic(data), 1e-2);
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

  ASSERT_NEAR(vc2.get_loglik(), vc2.loglik(u), 1e-2);
  ASSERT_NEAR(vc2.get_aic(), vc2.aic(u), 1e-2);
  ASSERT_NEAR(vc2.get_bic(), vc2.bic(u), 1e-2);
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
  auto data = vinecop.simulate(2000, false, 1, { 1 });

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
  EXPECT_TRUE(all_close(fit2.pdf(u, 2), fit2.pdf(u), 1e-10, 1e-10));
  EXPECT_TRUE(all_close(
    fit2.inverse_rosenblatt(u, 2), fit2.inverse_rosenblatt(u), 1e-10, 1e-10));
  EXPECT_TRUE(
    all_close(fit2.rosenblatt(u, 2), fit2.rosenblatt(u), 1e-10, 1e-10));

  // just check that it works
  fit2.simulate(2, false, 3, { 1 });
  fit2.cdf(u, 100, 2, { 1 });
}

// Records where a custom tree criterion is evaluated from. The bookkeeping is
// mutex-guarded so that a regression fails the assertions instead of racing on
// the recorder itself.
struct CriterionRecorder
{
  std::mutex m;
  std::set<std::thread::id> ids;
  size_t num_calls = 0;

  TreeCriterionFunction tau_function()
  {
    return [this](const Eigen::MatrixXd& data, const Eigen::VectorXd& weights) {
      {
        std::lock_guard<std::mutex> lk(m);
        ids.insert(std::this_thread::get_id());
        ++num_calls;
      }
      return wdm::wdm(data, "tau", weights)(0, 1);
    };
  }

  void reset()
  {
    ids.clear();
    num_calls = 0;
  }
};

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
    for (const auto& s1 : vc_sets[tree]) {
      bool is_in_both = false;
      for (const auto& s2 : vcl_sets[tree]) {
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
  EXPECT_EQ(controls.get_tree_algorithm(), "mst_prim");
  EXPECT_ANY_THROW(controls.set_tree_algorithm("foobar"));
  controls.set_tree_algorithm("mst_kruskal");

  // select structure and get matrix
  Vinecop fit(7);
  fit.select(u, controls);
  auto vcl_matrix = fit.get_matrix();

  // check if the same conditioned sets appear for each tree
  size_t pairs_unequal = get_pairs_unequal(vc_matrix, vcl_matrix, 6);
  EXPECT_EQ(pairs_unequal, 0);
}

TEST_F(VinecopTest, tree_criterion_custom_works)
{
  auto tau_fn = [](const Eigen::MatrixXd& data,
                   const Eigen::VectorXd& weights) {
    return wdm::wdm(data, "tau", weights)(0, 1);
  };

  // setter path: a custom criterion equal to Kendall's tau must recover the
  // same structure as the built-in "tau" (independence families so the
  // structure is determined solely by the criterion).
  Vinecop fit(7);
  FitControlsVinecop controls({ BicopFamily::indep });
  controls.set_tree_criterion("custom");
  controls.set_tree_criterion_function(tau_fn);
  fit.select(u, controls);
  EXPECT_EQ(get_pairs_unequal(vc_matrix, fit.get_matrix(), 6), 0);

  // FitControlsConfig construction path covers the new optional field.
  FitControlsConfig cfg;
  cfg.family_set = std::vector<BicopFamily>{ BicopFamily::indep };
  cfg.tree_criterion = "custom";
  cfg.tree_criterion_function = tau_fn;
  Vinecop fit_cfg(7);
  fit_cfg.select(u, FitControlsVinecop(cfg));
  EXPECT_EQ(get_pairs_unequal(vc_matrix, fit_cfg.get_matrix(), 6), 0);

  // "custom" without a callable must throw during selection.
  FitControlsVinecop bad({ BicopFamily::indep });
  bad.set_tree_criterion("custom");
  Vinecop fit_bad(7);
  EXPECT_ANY_THROW(fit_bad.select(u, bad));
  bad.set_num_threads(4);
  Vinecop fit_bad_mt(7);
  EXPECT_ANY_THROW(fit_bad_mt.select(u, bad));

  // a callable that no criterion would ever use must not be silently ignored
  FitControlsVinecop unused({ BicopFamily::indep });
  unused.set_tree_criterion_function(tau_fn);
  Vinecop fit_unused(7);
  EXPECT_ANY_THROW(fit_unused.select(u, unused));

  // the two fields may be set in either order
  FitControlsVinecop reversed({ BicopFamily::indep });
  reversed.set_tree_criterion_function(tau_fn);
  reversed.set_tree_criterion("custom");
  Vinecop fit_reversed(7);
  EXPECT_NO_THROW(fit_reversed.select(u, reversed));
}

// A custom criterion is arbitrary user code (an R closure or a Python callable
// in the wrappers), so it must run on the thread that starts the fit no matter
// how many threads the fit is allowed to use.
TEST_F(VinecopTest, tree_criterion_custom_is_serial)
{
  CriterionRecorder recorder;
  FitControlsVinecop controls({ BicopFamily::indep });
  controls.set_tree_criterion("custom");
  controls.set_tree_criterion_function(recorder.tau_function());

  Vinecop fit1(7);
  fit1.select(u, controls);
  size_t num_calls_serial = recorder.num_calls;
  EXPECT_EQ(recorder.ids.size(), static_cast<size_t>(1));
  EXPECT_EQ(*recorder.ids.begin(), std::this_thread::get_id());
  // tree 1 alone enumerates C(7, 2) = 21 candidate edges
  EXPECT_GT(num_calls_serial, static_cast<size_t>(20));

  recorder.reset();
  controls.set_num_threads(4);
  Vinecop fit4(7);
  fit4.select(u, controls);

  EXPECT_EQ(recorder.ids.size(), static_cast<size_t>(1));
  EXPECT_EQ(*recorder.ids.begin(), std::this_thread::get_id());
  // the candidate edges are enumerated serially, so the number of criterion
  // evaluations cannot depend on the thread count
  EXPECT_EQ(recorder.num_calls, num_calls_serial);
  EXPECT_EQ(get_pairs_unequal(fit1.get_matrix(), fit4.get_matrix(), 6), 0);
  EXPECT_EQ(get_pairs_unequal(vc_matrix, fit4.get_matrix(), 6), 0);
}

// Sparse selection re-enters add_allowed_edges once per threshold pass; the
// criterion must stay on the calling thread there too.
TEST_F(VinecopTest, tree_criterion_custom_is_serial_sparse)
{
  u.conservativeResize(200, 7); // above the 10-row floor of the criterion

  CriterionRecorder recorder;
  FitControlsVinecop controls({ BicopFamily::indep, BicopFamily::gaussian });
  controls.set_tree_criterion("custom");
  controls.set_tree_criterion_function(recorder.tau_function());
  controls.set_select_threshold(true);
  controls.set_select_trunc_lvl(true);
  controls.set_num_threads(4);

  Vinecop fit(7);
  fit.select(u, controls);

  EXPECT_GT(recorder.num_calls, static_cast<size_t>(20));
  EXPECT_EQ(recorder.ids.size(), static_cast<size_t>(1));
  EXPECT_EQ(*recorder.ids.begin(), std::this_thread::get_id());
}

TEST_F(VinecopTest, select_finds_different_structures_random)
{
  // Initialize the controls
  FitControlsVinecop controls_weighted({ BicopFamily::tll });
  controls_weighted.set_tree_algorithm("random_weighted");

  FitControlsVinecop controls_unweighted({ BicopFamily::tll });
  controls_unweighted.set_tree_algorithm("random_unweighted");

  // For reseeding the random number generator (deterministic, so the test
  // is reproducible in CI)
  std::vector<int> seeds(20);

  // To store the unique structures
  std::set<TriangularArray<size_t>> unique_structures_weighted;
  std::set<TriangularArray<size_t>> unique_structures_unweighted;

  // To store the first RNG value after each reseeding
  std::set<uint32_t> first_rng_outputs;

  const size_t num_trials = 10;

  for (size_t i = 0; i < num_trials; ++i) {
    // Deterministic but distinct seeds for each trial, so every run uses a
    // different RNG stream (and hence a different random structure) without
    // depending on std::random_device.
    for (size_t k = 0; k < seeds.size(); ++k)
      seeds[k] = static_cast<int>(i * seeds.size() + k + 1);
    controls_weighted.set_seeds(seeds);
    controls_unweighted.set_seeds(seeds);

    // Check RNG output changes
    auto rng_sample_weighted =
      controls_weighted.get_rng()(); // Get first sample
    first_rng_outputs.insert(rng_sample_weighted);

    auto rng_sample_unweighted =
      controls_unweighted.get_rng()(); // Get first sample
    first_rng_outputs.insert(rng_sample_unweighted);

    // Select a random structure for the weighted method
    Vinecop fit_weighted(u, RVineStructure(), {}, controls_weighted);
    auto struct_array_weighted = fit_weighted.get_struct_array();
    unique_structures_weighted.insert(struct_array_weighted);

    // Select a random structure for the unweighted method
    Vinecop fit_unweighted(u, RVineStructure(), {}, controls_unweighted);
    auto struct_array_unweighted = fit_unweighted.get_struct_array();
    unique_structures_unweighted.insert(struct_array_unweighted);
  }

  // With distinct seeds the RNG streams and selected structures should all
  // differ; being deterministic, this is stable across runs.
  EXPECT_EQ(first_rng_outputs.size(), num_trials);
  EXPECT_EQ(unique_structures_weighted.size(), num_trials);
  EXPECT_EQ(unique_structures_unweighted.size(), num_trials);
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
  controls.set_select_threshold(true);
  controls.set_threshold(NAN);
  // controls.set_show_trace(true);
  controls.set_selection_criterion("mbicv");

  Vinecop fit(7);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);

  u = tools_stats::simulate_uniform(100, 7, false, { 1 });
  fit.select(u, controls);
}

TEST_F(VinecopTest, sparse_truncation_selection)
{
  u.conservativeResize(50, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_select_trunc_lvl(true);
  // controls.set_show_trace(true);
  u = tools_stats::simulate_uniform(100, 7, false, { 1 });
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

TEST_F(VinecopTest, tawn_flipping)
{
  FitControlsVinecop controls({ BicopFamily::tawn });
  Vinecop fit1(7);
  fit1.select(u, controls);
  Vinecop fit2(fit1.get_rvine_structure());
  fit2.select(u, controls);

  for (size_t tree = 0; tree < fit1.get_trunc_lvl(); ++tree) {
    for (size_t edge = 0; edge < 6 - tree; ++edge) {
      auto pc1 = fit1.get_pair_copula(tree, edge);
      auto pc2 = fit2.get_pair_copula(tree, edge);
      ASSERT_EQ(pc1.get_family(), pc2.get_family());
      ASSERT_EQ(pc1.get_rotation(), pc2.get_rotation());
      ASSERT_TRUE(
        all_close(pc1.get_parameters(), pc2.get_parameters(), 0.01, 0.01));
    }
  }
}

// Selection reuses the pair copulas it fitted while building the trees, and
// places one flipped whenever the finalized slot's diagonal variable is not
// the edge's first argument. Refitting on the selected structure fits in the
// finalized order instead, so the two only agree if a `tll` fit is a function
// of the observations rather than of the order they arrive in.
TEST_F(VinecopTest, tll_selection_is_reproducible)
{
  FitControlsVinecop controls({ BicopFamily::tll });
  controls.set_num_threads(1);
  Vinecop fit1(7);
  fit1.select(u, controls);
  Vinecop fit2(fit1.get_rvine_structure());
  fit2.select(u, controls);

  for (size_t tree = 0; tree < fit1.get_trunc_lvl(); ++tree) {
    for (size_t edge = 0; edge < 6 - tree; ++edge) {
      auto pc1 = fit1.get_pair_copula(tree, edge);
      auto pc2 = fit2.get_pair_copula(tree, edge);
      ASSERT_EQ(pc1.get_family(), pc2.get_family());
      ASSERT_TRUE(
        all_close(pc1.get_parameters(), pc2.get_parameters(), 1e-10, 1e-10))
        << "tree " << tree << ", edge " << edge;
    }
  }
}

}
