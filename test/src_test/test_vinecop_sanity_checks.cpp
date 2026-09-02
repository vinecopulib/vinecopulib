// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include <vinecopulib.hpp>

namespace test_vinecop_sanity_checks {
using namespace vinecopulib;

// A model constructed from a structure alone has no families to fit; fitting a
// truncated-at-zero model is the independence fit, whose log-likelihood is 0.
TEST(vinecop_sanity_checks, fit_needs_pair_copulas)
{
  const auto data = tools_stats::simulate_uniform(50, 4, false, { 1 });

  Vinecop no_families(DVineStructure({ 1, 2, 3, 4 }));
  EXPECT_THROW(no_families.fit(data), std::runtime_error);

  Vinecop independence(4); // a D-vine truncated at 0
  ASSERT_EQ(independence.get_trunc_lvl(), static_cast<size_t>(0));
  EXPECT_NO_THROW(independence.fit(data));
  EXPECT_EQ(independence.get_loglik(), 0.0);
  EXPECT_EQ(independence.get_nobs(), static_cast<size_t>(50));
}

TEST(vinecop_sanity_checks, catches_wrong_tree)
{
  Vinecop vinecop(3);
  EXPECT_ANY_THROW(vinecop.get_pair_copula(3, 0));
}

TEST(vinecop_sanity_checks, catches_wrong_edge)
{
  Vinecop vinecop(3);
  EXPECT_ANY_THROW(vinecop.get_pair_copula(0, -1));
  EXPECT_ANY_THROW(vinecop.get_pair_copula(0, 2));
}

TEST(vinecop_sanity_checks, catches_wrong_size)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(3);
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(3, 3);
  mat << 2, 2, 2, 1, 1, 0, 3, 0, 0;
  Vinecop vinecop(mat, pair_copulas);

  Eigen::MatrixXd U = Eigen::MatrixXd::Constant(4, 4, 0.5);
  try {
    vinecop.pdf(U);
    FAIL() << "expected invalid column count to throw";
  } catch (const std::runtime_error& error) {
    EXPECT_EQ(error.what(),
              std::string("data has wrong number of columns; expected: 3 "
                          "(n x d continuous layout), actual: 4 (model "
                          "contains no discrete variables)."));
  }
  EXPECT_ANY_THROW(vinecop.cdf(U));
  try {
    vinecop.inverse_rosenblatt(U);
    FAIL() << "expected invalid inverse Rosenblatt column count to throw";
  } catch (const std::runtime_error& error) {
    EXPECT_EQ(error.what(),
              std::string("data has wrong number of columns; expected: 3 "
                          "(n x d input), actual: 4."));
  }
  EXPECT_ANY_THROW(vinecop.select(U));

  vinecop.set_var_types({ "d", "c", "c" });
  auto U_discrete = Eigen::MatrixXd::Constant(4, 3, 0.5);
  try {
    vinecop.pdf(U_discrete);
    FAIL() << "expected invalid discrete column count to throw";
  } catch (const std::runtime_error& error) {
    EXPECT_EQ(error.what(),
              std::string("data has wrong number of columns; expected: 6 "
                          "(n x 2d expanded layout) or 4 (n x (d + k) "
                          "compact layout), actual: 3 (model contains 1 "
                          "discrete variable)."));
  }

  vinecop.set_var_types({ "d", "d", "d" });
  try {
    vinecop.pdf(U_discrete);
    FAIL() << "expected missing discrete left-limit columns to throw";
  } catch (const std::runtime_error& error) {
    EXPECT_EQ(error.what(),
              std::string("data has wrong number of columns; expected: 6 "
                          "(n x 2d expanded or n x (d + k) compact layout), "
                          "actual: 3 (model contains 3 discrete variables)."));
  }

  vinecop = Vinecop(301);
  U = tools_stats::simulate_uniform(1, 301, false, { 1 });
  EXPECT_NO_THROW(vinecop.cdf(U, 10000, 1, { 1 }));
  vinecop = Vinecop(21202);
  U = tools_stats::simulate_uniform(1, 21202, false, { 1 });
  EXPECT_ANY_THROW(vinecop.cdf(U));

  pair_copulas.resize(4); // too many
  EXPECT_ANY_THROW(Vinecop(mat, pair_copulas));
  pair_copulas = Vinecop::make_pair_copula_store(3);
  pair_copulas[0].pop_back(); // to few
  EXPECT_ANY_THROW(Vinecop(mat, pair_copulas));
}

TEST(vinecop_sanity_checks, var_types_needs_one_entry_per_variable)
{
  // `check_var_types` rejected a vector longer than the dimension but accepted
  // a shorter one, which `set_var_types_internal` then stored and the pair-type
  // derivation indexed per variable -- an out-of-bounds read on every
  // evaluation. `Bicop::check_var_types` has always required exactly two.
  auto vinecop = Vinecop(RVineStructure(std::vector<size_t>{ 1, 2, 3, 4 }));
  EXPECT_ANY_THROW(vinecop.set_var_types({ "c", "c" }));
  EXPECT_ANY_THROW(vinecop.set_var_types({ "c" }));
  EXPECT_ANY_THROW(vinecop.set_var_types({ "c", "c", "c", "c", "c" }));
  EXPECT_NO_THROW(vinecop.set_var_types({ "c", "d", "c", "d" }));
}

TEST(vinecop_sanity_checks, controls_print)
{
  auto controls = FitControlsVinecop();
  EXPECT_NO_THROW(controls.str());

  // Each flag has to report itself. The two "Select ..." lines are adjacent
  // and otherwise identical, so one reading the other's getter still prints
  // something plausible.
  auto reports = [](const FitControlsVinecop& c, const std::string& label) {
    const auto str = c.str();
    const auto at = str.find(label + ": ");
    EXPECT_NE(at, std::string::npos) << label << " missing from str()";
    return str.substr(at + label.size() + 2, 3) == "yes";
  };

  controls.set_select_trunc_lvl(true);
  controls.set_select_threshold(false);
  EXPECT_TRUE(reports(controls, "Select truncation level"));
  EXPECT_FALSE(reports(controls, "Select threshold"));

  controls.set_select_trunc_lvl(false);
  controls.set_select_threshold(true);
  EXPECT_FALSE(reports(controls, "Select truncation level"));
  EXPECT_TRUE(reports(controls, "Select threshold"));
}

TEST(vinecop_sanity_checks, fit_controls_config_works)
{
  // Some controls for testing
  // Only the non FitControlsBicop fields are tested here
  FitControlsVinecop controls;
  controls.set_trunc_lvl(3);
  controls.set_tree_criterion("rho");
  controls.set_threshold(0.5);
  controls.set_select_threshold(true);
  controls.set_select_trunc_lvl(true);
  controls.set_select_families(false);
  controls.set_show_trace(true);
  controls.set_tree_algorithm("mst_kruskal");
  std::vector<int> seeds = { 1, 2, 3, 4, 5 };
  controls.set_seeds(seeds);
  controls.set_conditioning_set({ 2, 4 });

  // Create a config object from the controls
  FitControlsConfig config;
  config.trunc_lvl = controls.get_trunc_lvl();
  config.tree_criterion = controls.get_tree_criterion();
  config.threshold = controls.get_threshold();
  config.select_threshold = controls.get_select_threshold();
  config.select_trunc_lvl = controls.get_select_trunc_lvl();
  config.select_families = controls.get_select_families();
  config.show_trace = controls.get_show_trace();
  config.num_threads = controls.get_num_threads();
  config.tree_algorithm = controls.get_tree_algorithm();
  config.seeds = controls.get_seeds();
  config.conditioning_set = controls.get_conditioning_set();

  // Create and test new controls from the config object
  FitControlsVinecop controls2(config);
  EXPECT_EQ(controls.get_trunc_lvl(), controls2.get_trunc_lvl());
  EXPECT_EQ(controls.get_tree_criterion(), controls2.get_tree_criterion());
  EXPECT_EQ(controls.get_threshold(), controls2.get_threshold());
  EXPECT_EQ(controls.get_select_threshold(), controls2.get_select_threshold());
  EXPECT_EQ(controls.get_select_trunc_lvl(), controls2.get_select_trunc_lvl());
  EXPECT_EQ(controls.get_select_families(), controls2.get_select_families());
  EXPECT_EQ(controls.get_show_trace(), controls2.get_show_trace());
  EXPECT_EQ(controls.get_num_threads(), controls2.get_num_threads());
  EXPECT_EQ(controls.get_tree_algorithm(), controls2.get_tree_algorithm());
  EXPECT_EQ(controls.get_seeds(), controls2.get_seeds());
  EXPECT_EQ(controls.get_conditioning_set(), controls2.get_conditioning_set());
}

TEST(vinecop_sanity_checks, controls_check)
{
  auto controls = FitControlsVinecop();
  EXPECT_NO_THROW(controls.set_tree_criterion("cxi"));
  EXPECT_ANY_THROW(controls.set_tree_criterion("foo"));
  EXPECT_ANY_THROW(controls.set_threshold(-1.0));
  EXPECT_ANY_THROW(controls.set_threshold(2.0));
  EXPECT_ANY_THROW(controls.set_conditioning_set({ 0 }));
  EXPECT_ANY_THROW(controls.set_conditioning_set({ 2, 2 }));
}
}
