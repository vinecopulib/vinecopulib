// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "gtest/gtest.h"
#include <string>
#include <vector>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/bicop/tools_select.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/class.hpp>

namespace test_circular {

using namespace vinecopulib;
using tools_stl::is_member;

const std::vector<std::string> cc = { "c", "c" };
const std::vector<std::string> cd = { "c", "d" };
const std::vector<std::string> aa = { "a", "a" };
const std::vector<std::string> ac = { "a", "c" };
const std::vector<std::string> ca = { "c", "a" };
const std::vector<BicopFamily> every_family =
  tools_stl::cat(bicop_families::all, bicop_families::circular);

TEST(test_circular, family_names_round_trip)
{
  for (auto fam : bicop_families::circular) {
    EXPECT_EQ(get_family_enum(get_family_name(fam)), fam);
    // registered, but not implemented and therefore not in `all` yet
    EXPECT_FALSE(is_member(fam, bicop_families::all));
    EXPECT_FALSE(is_member(fam, bicop_families::parametric));
    EXPECT_FALSE(is_member(fam, bicop_families::itau));
    EXPECT_FALSE(is_member(fam, bicop_families::analytic_derivs));
  }
  EXPECT_EQ(get_family_name(BicopFamily::cardioid), "Cardioid");
  EXPECT_EQ(get_family_name(BicopFamily::wrapped_cauchy), "Wrapped Cauchy");
  EXPECT_EQ(get_family_name(BicopFamily::von_mises), "von Mises");
  EXPECT_EQ(get_family_name(BicopFamily::quad_sections), "Quadratic sections");
  EXPECT_EQ(get_family_name(BicopFamily::cubic_sections), "Cubic sections");
}

TEST(test_circular, rotation_groups_partition_the_circular_families)
{
  for (auto fam : bicop_families::two_rotations) {
    EXPECT_TRUE(is_member(fam, bicop_families::circular));
    EXPECT_FALSE(is_member(fam, bicop_families::rotationless));
    EXPECT_FALSE(is_member(fam, bicop_families::cylindrical));
  }
  for (auto fam : bicop_families::cylindrical) {
    EXPECT_TRUE(is_member(fam, bicop_families::circular));
    EXPECT_TRUE(is_member(fam, bicop_families::rotationless));
  }
}

TEST(test_circular, eligibility_table)
{
  // linear and discrete pairs: exactly the pre-existing families
  for (const auto& types : { cc, cd }) {
    for (auto fam : every_family) {
      EXPECT_EQ(family_accepts_var_types(fam, types),
                !is_member(fam, bicop_families::circular))
        << get_family_name(fam);
    }
  }
  // circular-circular: independence, the circulas, and tll
  for (auto fam : every_family) {
    bool expected = fam == BicopFamily::indep || fam == BicopFamily::tll ||
                    is_member(fam, bicop_families::two_rotations);
    EXPECT_EQ(family_accepts_var_types(fam, aa), expected)
      << get_family_name(fam);
  }
  // circular-linear, either order: additionally the cylindrical families
  for (const auto& types : { ac, ca }) {
    for (auto fam : every_family) {
      bool expected = fam == BicopFamily::indep || fam == BicopFamily::tll ||
                      is_member(fam, bicop_families::circular);
      EXPECT_EQ(family_accepts_var_types(fam, types), expected)
        << get_family_name(fam);
    }
  }
  // circular-discrete: nothing
  for (auto fam : every_family) {
    EXPECT_FALSE(family_accepts_var_types(fam, { "a", "d" }));
    EXPECT_FALSE(family_accepts_var_types(fam, { "d", "a" }));
  }
  // the default search for a linear pair is unchanged
  EXPECT_EQ(eligible_families(every_family, cc), bicop_families::all);
}

TEST(test_circular, bicop_accepts_the_token_and_rejects_bad_pairs)
{
  // independence models any geometry
  Bicop indep(BicopFamily::indep, 0, Eigen::MatrixXd(), aa);
  EXPECT_EQ(indep.get_var_types(), aa);
  indep.set_var_types(ac);
  EXPECT_EQ(indep.get_var_types(), ac);

  // a linear family cannot take a circular variable
  EXPECT_THROW(Bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Zero(1), aa),
               std::runtime_error);
  Bicop gauss(BicopFamily::gaussian, 0, Eigen::VectorXd::Zero(1));
  EXPECT_THROW(gauss.set_var_types(ca), std::runtime_error);

  // unknown tokens and circular-discrete pairs are rejected
  EXPECT_THROW(indep.set_var_types({ "a", "d" }), std::runtime_error);
  EXPECT_THROW(indep.set_var_types({ "x", "c" }), std::runtime_error);

  // the circular families are registered but not implemented yet
  for (auto fam : bicop_families::circular) {
    EXPECT_THROW(Bicop(fam, 0, Eigen::MatrixXd(), aa), std::runtime_error);
  }
}

TEST(test_circular, circular_types_round_trip_through_json)
{
  Bicop indep(BicopFamily::indep, 0, Eigen::MatrixXd(), ac);
  Bicop back(indep.to_json());
  EXPECT_EQ(back.get_var_types(), ac);
  auto u = tools_stats::simulate_uniform(20, 2, false, { 1 });
  EXPECT_TRUE(back.pdf(u).isApprox(indep.pdf(u)));
}

TEST(test_circular, candidate_generation_filters_by_geometry)
{
  auto u = tools_stats::simulate_uniform(200, 2, false, { 2 });
  FitControlsBicop controls;

  // linear pair: the candidate set is what it was before circular families
  auto linear = tools_select::create_candidate_bicops(u, controls, cc);
  for (const auto& bc : linear) {
    EXPECT_FALSE(is_member(bc.get_family(), bicop_families::circular));
  }
  auto linear_default = tools_select::create_candidate_bicops(u, controls);
  EXPECT_EQ(linear.size(), linear_default.size());

  // an explicit family set with nothing eligible for the pair throws and
  // names the eligible families
  controls.set_family_set({ BicopFamily::gaussian, BicopFamily::clayton });
  try {
    tools_select::get_candidate_families(controls, aa);
    FAIL() << "expected an exception";
  } catch (const std::runtime_error& e) {
    std::string msg = e.what();
    EXPECT_NE(msg.find("Independence"), std::string::npos);
    EXPECT_NE(msg.find("TLL"), std::string::npos);
    EXPECT_EQ(msg.find("Gaussian"), std::string::npos);
  }

  // an explicit set is intersected with the eligible families
  controls.set_family_set(
    { BicopFamily::gaussian, BicopFamily::indep, BicopFamily::quad_sections });
  auto fams = tools_select::get_candidate_families(controls, ac);
  EXPECT_EQ(fams,
            (std::vector<BicopFamily>{ BicopFamily::indep,
                                       BicopFamily::quad_sections }));
  fams = tools_select::get_candidate_families(controls, aa);
  EXPECT_EQ(fams, (std::vector<BicopFamily>{ BicopFamily::indep }));

  // selection on a circular pair with only independence eligible works
  controls.set_family_set({ BicopFamily::indep });
  Bicop bc;
  bc.set_var_types(aa);
  bc.select(u, controls);
  EXPECT_EQ(bc.get_family(), BicopFamily::indep);
  EXPECT_EQ(bc.get_var_types(), aa);
}

TEST(test_circular, vinecop_accepts_the_token_but_does_not_select_yet)
{
  auto u = tools_stats::simulate_uniform(50, 3, false, { 3 });
  Vinecop vc(3);
  vc.set_var_types({ "a", "c", "c" });
  EXPECT_EQ(vc.get_var_types(), (std::vector<std::string>{ "a", "c", "c" }));
  EXPECT_THROW(vc.set_var_types({ "a", "d", "c" }), std::runtime_error);
  EXPECT_THROW(vc.set_var_types({ "a", "x", "c" }), std::runtime_error);
  EXPECT_THROW(vc.select(u), std::runtime_error);

  // discrete-free linear vines are unaffected
  Vinecop linear(3);
  linear.set_var_types({ "c", "c", "c" });
  EXPECT_NO_THROW(linear.select(u));
}

} // namespace test_circular
