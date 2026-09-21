// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "include/circular_golden.hpp"
#include "include/test_utils.hpp"
#include "gtest/gtest.h"
#include <boost/math/constants/constants.hpp>
#include <string>
#include <vector>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/bicop/tools_select.hpp>
#include <vinecopulib/misc/nlohmann_json.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/class.hpp>

namespace test_circular {

using namespace vinecopulib;
using test_utils::all_close;
using test_utils::default_var_types;
using test_utils::make_bicop;
using tools_stl::is_member;

const std::vector<std::string> cc = { "c", "c" };
const std::vector<std::string> cd = { "c", "d" };
const std::vector<std::string> aa = { "a", "a" };
const std::vector<std::string> ac = { "a", "c" };
const std::vector<std::string> ca = { "c", "a" };
const double pi = boost::math::constants::pi<double>();

Eigen::VectorXd
to_vector(const nlohmann::json& j)
{
  Eigen::VectorXd v(j.size());
  for (size_t i = 0; i < j.size(); ++i) {
    v(i) = j[i].get<double>();
  }
  return v;
}

Eigen::MatrixXd
to_matrix(const Eigen::VectorXd& a, const Eigen::VectorXd& b)
{
  Eigen::MatrixXd m(a.size(), 2);
  m.col(0) = a;
  m.col(1) = b;
  return m;
}

BicopFamily
family_of(const std::string& key)
{
  static const std::vector<std::pair<std::string, BicopFamily>> table = {
    { "cardioid", BicopFamily::cardioid },
    { "wrapped_cauchy", BicopFamily::wrapped_cauchy },
    { "von_mises", BicopFamily::von_mises },
    { "quad_sections", BicopFamily::quad_sections },
    { "cubic_sections", BicopFamily::cubic_sections },
  };
  for (const auto& entry : table) {
    if (key.rfind(entry.first, 0) == 0) {
      return entry.second;
    }
  }
  throw std::runtime_error("unknown golden key " + key);
}

// -------------------------------------------------------------------------
// registration and eligibility

TEST(test_circular, family_names_round_trip)
{
  for (auto fam : bicop_families::circular) {
    EXPECT_EQ(get_family_enum(get_family_name(fam)), fam);
    EXPECT_TRUE(is_member(fam, bicop_families::all));
    EXPECT_TRUE(is_member(fam, bicop_families::parametric));
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
    EXPECT_TRUE(is_member(fam, bicop_families::two_par));
  }
  for (auto fam : bicop_families::cylindrical) {
    EXPECT_TRUE(is_member(fam, bicop_families::circular));
    EXPECT_TRUE(is_member(fam, bicop_families::rotationless));
  }
  EXPECT_TRUE(is_member(BicopFamily::quad_sections, bicop_families::two_par));
  EXPECT_TRUE(
    is_member(BicopFamily::cubic_sections, bicop_families::three_par));
}

TEST(test_circular, eligibility_table)
{
  // linear and discrete pairs: exactly the pre-existing families
  for (const auto& types : { cc, cd }) {
    for (auto fam : bicop_families::all) {
      EXPECT_EQ(family_accepts_var_types(fam, types),
                !is_member(fam, bicop_families::circular))
        << get_family_name(fam);
    }
  }
  // circular-circular: independence, the circulas, and tll
  for (auto fam : bicop_families::all) {
    bool expected = fam == BicopFamily::indep || fam == BicopFamily::tll ||
                    is_member(fam, bicop_families::two_rotations);
    EXPECT_EQ(family_accepts_var_types(fam, aa), expected)
      << get_family_name(fam);
  }
  // circular-linear, either order: additionally the cylindrical families
  for (const auto& types : { ac, ca }) {
    for (auto fam : bicop_families::all) {
      bool expected = fam == BicopFamily::indep || fam == BicopFamily::tll ||
                      is_member(fam, bicop_families::circular);
      EXPECT_EQ(family_accepts_var_types(fam, types), expected)
        << get_family_name(fam);
    }
  }
  // circular-discrete: nothing
  for (auto fam : bicop_families::all) {
    EXPECT_FALSE(family_accepts_var_types(fam, { "a", "d" }));
    EXPECT_FALSE(family_accepts_var_types(fam, { "d", "a" }));
  }
  // the default search for a linear pair is unchanged
  EXPECT_EQ(eligible_families(bicop_families::all, cc),
            tools_stl::set_diff(bicop_families::all, bicop_families::circular));
}

TEST(test_circular, bicop_accepts_the_token_and_rejects_bad_pairs)
{
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

  // circular families need a circular variable; cylindrical ones exactly one
  for (auto fam : bicop_families::circular) {
    EXPECT_THROW(Bicop(fam, 0, Eigen::MatrixXd(), cc), std::runtime_error);
    EXPECT_THROW((Bicop{ fam }), std::runtime_error);
    EXPECT_NO_THROW(make_bicop(fam));
  }
  for (auto fam : bicop_families::cylindrical) {
    EXPECT_THROW(Bicop(fam, 0, Eigen::MatrixXd(), aa), std::runtime_error);
    EXPECT_NO_THROW(Bicop(fam, 0, Eigen::MatrixXd(), ca));
  }
}

TEST(test_circular, rotations_follow_the_groups)
{
  for (auto fam : bicop_families::two_rotations) {
    EXPECT_NO_THROW(make_bicop(fam, 0));
    EXPECT_NO_THROW(make_bicop(fam, 90));
    EXPECT_THROW(make_bicop(fam, 180), std::runtime_error);
    EXPECT_THROW(make_bicop(fam, 270), std::runtime_error);
  }
  for (auto fam : bicop_families::cylindrical) {
    EXPECT_THROW(make_bicop(fam, 90), std::runtime_error);
    EXPECT_THROW(make_bicop(fam, 180), std::runtime_error);
  }
}

TEST(test_circular, candidate_generation_filters_by_geometry)
{
  auto u = tools_stats::simulate_uniform(200, 2, false, { 2 });
  FitControlsBicop controls;

  auto linear = tools_select::create_candidate_bicops(u, controls, cc);
  for (const auto& bc : linear) {
    EXPECT_FALSE(is_member(bc.get_family(), bicop_families::circular));
  }
  auto linear_default = tools_select::create_candidate_bicops(u, controls);
  EXPECT_EQ(linear.size(), linear_default.size());

  // both orientations of every circula, one of every cylindrical family
  auto circ = tools_select::create_candidate_bicops(u, controls, ac);
  size_t n_two_rot = 0, n_cyl = 0;
  for (const auto& bc : circ) {
    EXPECT_TRUE(family_accepts_var_types(bc.get_family(), ac));
    EXPECT_EQ(bc.get_var_types(), ac);
    n_two_rot += is_member(bc.get_family(), bicop_families::two_rotations);
    n_cyl += is_member(bc.get_family(), bicop_families::cylindrical);
  }
  EXPECT_EQ(n_two_rot, 2 * bicop_families::two_rotations.size());
  EXPECT_EQ(n_cyl, bicop_families::cylindrical.size());
  controls.set_allow_rotations(false);
  circ = tools_select::create_candidate_bicops(u, controls, aa);
  for (const auto& bc : circ) {
    EXPECT_EQ(bc.get_rotation(), 0);
  }
  controls.set_allow_rotations(true);

  controls.set_family_set({ BicopFamily::gaussian, BicopFamily::clayton });
  try {
    tools_select::get_candidate_families(controls, aa);
    FAIL() << "expected an exception";
  } catch (const std::runtime_error& e) {
    std::string msg = e.what();
    EXPECT_NE(msg.find("Independence"), std::string::npos);
    EXPECT_NE(msg.find("Cardioid"), std::string::npos);
    EXPECT_EQ(msg.find("Gaussian"), std::string::npos);
  }

  controls.set_family_set(
    { BicopFamily::gaussian, BicopFamily::indep, BicopFamily::quad_sections });
  auto fams = tools_select::get_candidate_families(controls, ac);
  EXPECT_EQ(fams,
            (std::vector<BicopFamily>{ BicopFamily::indep,
                                       BicopFamily::quad_sections }));
  fams = tools_select::get_candidate_families(controls, aa);
  EXPECT_EQ(fams, (std::vector<BicopFamily>{ BicopFamily::indep }));
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

  Vinecop linear(3);
  linear.set_var_types({ "c", "c", "c" });
  EXPECT_NO_THROW(linear.select(u));
}

// -------------------------------------------------------------------------
// the parametric families

TEST(test_circular, golden_values)
{
  const auto golden = nlohmann::json::parse(test_circular_golden::json);
  size_t n_cases = 0;
  for (auto it = golden.begin(); it != golden.end(); ++it) {
    if (!it.value().is_array()) {
      continue;
    }
    const auto family = family_of(it.key());
    const bool root_solved = family == BicopFamily::cardioid ||
                             family == BicopFamily::von_mises ||
                             is_member(family, bicop_families::cylindrical);
    const double tol = family == BicopFamily::von_mises ? 1e-8 : 1e-9;
    const double tol_inv = root_solved ? 1e-7 : 1e-9;
    for (const auto& c : it.value()) {
      const auto pars = to_vector(c["parameters"]);
      const int rotation = c["rotation"].get<int>();
      Bicop bc(family, rotation, pars, default_var_types(family));
      const auto uv = to_matrix(to_vector(c["u"]), to_vector(c["v"]));
      const std::string what = it.key() + " " + bc.str();

      EXPECT_TRUE(all_close(bc.pdf(uv), to_vector(c["pdf"]), tol, tol)) << what;
      EXPECT_TRUE(all_close(bc.cdf(uv), to_vector(c["cdf"]), 1e-7, 1e-7))
        << what;
      EXPECT_TRUE(all_close(bc.hfunc1(uv), to_vector(c["hfunc1"]), tol, tol))
        << what;
      EXPECT_TRUE(all_close(bc.hfunc2(uv), to_vector(c["hfunc2"]), tol, tol))
        << what;
      const auto w = to_vector(c["w"]);
      for (Eigen::Index k = 0; k < w.size(); ++k) {
        Eigen::MatrixXd uw = uv;
        uw.col(1).setConstant(w(k));
        EXPECT_TRUE(
          all_close(bc.hinv1(uw), to_vector(c["hinv1"][k]), tol_inv, tol_inv))
          << what << " w = " << w(k);
        Eigen::MatrixXd wv = uv;
        wv.col(0).setConstant(w(k));
        EXPECT_TRUE(
          all_close(bc.hinv2(wv), to_vector(c["hinv2"][k]), tol_inv, tol_inv))
          << what << " w = " << w(k);
      }
      EXPECT_NEAR(bc.parameters_to_tau(pars), c["tau"].get<double>(), 1e-6)
        << what;
      ++n_cases;
    }
  }
  EXPECT_GT(n_cases, 30u);
}

TEST(test_circular, identities_on_a_grid)
{
  Eigen::MatrixXd u(7 * 7, 2);
  Eigen::VectorXd g = Eigen::VectorXd::LinSpaced(7, 0.05, 0.95);
  for (long i = 0, k = 0; i < 7; ++i) {
    for (long j = 0; j < 7; ++j, ++k) {
      u(k, 0) = g(i);
      u(k, 1) = g(j);
    }
  }
  Eigen::MatrixXd u0 = u, u1 = u;
  u0.col(0).setConstant(1e-12);
  u1.col(0).setConstant(1 - 1e-12);

  struct Case
  {
    BicopFamily family;
    int rotation;
    std::vector<double> pars;
  };
  const std::vector<Case> cases = {
    { BicopFamily::cardioid, 0, { 0.4, 1.0 } },
    { BicopFamily::cardioid, 90, { 0.5, -2.0 } },
    { BicopFamily::wrapped_cauchy, 0, { 0.7, pi } },
    { BicopFamily::wrapped_cauchy, 90, { 0.95, 0.3 } },
    { BicopFamily::von_mises, 0, { 4.0, 2.5 } },
    { BicopFamily::von_mises, 90, { 30.0, -0.7 } },
    { BicopFamily::quad_sections, 0, { 0.8, 1.2 } },
    { BicopFamily::cubic_sections, 0, { 0.9, -0.6, -2.0 } },
  };
  for (const auto& c : cases) {
    Eigen::VectorXd pars =
      Eigen::Map<const Eigen::VectorXd>(c.pars.data(), c.pars.size());
    Bicop bc(c.family, c.rotation, pars, default_var_types(c.family));
    const std::string what = bc.str();

    // h(hinv) = id, hinv in [0, 1]
    Eigen::MatrixXd v = u;
    v.col(1) = bc.hinv1(u);
    EXPECT_TRUE((v.col(1).array() >= 0).all() && (v.col(1).array() <= 1).all())
      << what;
    EXPECT_TRUE(all_close(bc.hfunc1(v), u.col(1), 1e-7, 1e-7)) << what;
    v = u;
    v.col(0) = bc.hinv2(u);
    EXPECT_TRUE(all_close(bc.hfunc2(v), u.col(0), 1e-7, 1e-7)) << what;

    // periodic in the circular coordinate(s); the arguments are trimmed to
    // [1e-10, 1 - 1e-10] on the way in, so steep densities differ slightly
    EXPECT_TRUE(all_close(bc.pdf(u0), bc.pdf(u1), 1e-6, 1e-8)) << what;
    if (is_member(c.family, bicop_families::two_rotations)) {
      Eigen::MatrixXd v0 = u, v1 = u;
      v0.col(1).setConstant(1e-12);
      v1.col(1).setConstant(1 - 1e-12);
      EXPECT_TRUE(all_close(bc.pdf(v0), bc.pdf(v1), 1e-6, 1e-8)) << what;
    }

    // flip: the density of the flipped copula at swapped arguments
    Bicop flipped = bc;
    flipped.flip();
    EXPECT_TRUE(is_member(flipped.get_rotation(), { 0, 90 })) << what;
    EXPECT_TRUE(all_close(
      flipped.pdf(tools_eigen::swap_cols(u)), bc.pdf(u), 1e-10, 1e-10))
      << what;
    EXPECT_TRUE(all_close(
      flipped.hfunc1(tools_eigen::swap_cols(u)), bc.hfunc2(u), 1e-10, 1e-10))
      << what;

    // simulation returns uniform-ish margins in (0, 1)
    auto sim = bc.simulate(500, false, { 11 });
    EXPECT_TRUE((sim.array() > 0).all() && (sim.array() < 1).all()) << what;
    EXPECT_NEAR(sim.col(0).mean(), 0.5, 0.05) << what;
    EXPECT_NEAR(sim.col(1).mean(), 0.5, 0.05) << what;

    // per-row parameters broadcast to the single-row result
    Eigen::MatrixXd rows = pars.transpose().replicate(u.rows(), 1);
    EXPECT_TRUE(all_close(bc.pdf(u, rows), bc.pdf(u), 1e-12, 1e-12)) << what;
    EXPECT_TRUE(all_close(bc.hfunc1(u, rows), bc.hfunc1(u), 1e-12, 1e-12))
      << what;

    // no tau inversion, zero tail dependence
    EXPECT_THROW(bc.tau_to_parameters(0.3), std::runtime_error) << what;
    EXPECT_TRUE(bc.get_taildep().isZero()) << what;
  }
}

TEST(test_circular, phase_conventions)
{
  Eigen::VectorXd p1(2), p2(2);
  p1 << 0.4, 0.4;
  p2 << 0.4, 0.4 + 2 * pi;
  auto u = tools_stats::simulate_uniform(50, 2, false, { 4 });
  for (auto fam : bicop_families::two_rotations) {
    // phases differing by 2 pi define the same copula
    Bicop a(fam, 0, p1, aa), b(fam, 0, p2, aa);
    EXPECT_TRUE(all_close(a.pdf(u), b.pdf(u), 1e-10, 1e-10));
    // the rotation-90 copula is exchangeable: flipping keeps rotation and phase
    Bicop r(fam, 90, p1, aa);
    Bicop rf = r;
    rf.flip();
    EXPECT_EQ(rf.get_rotation(), 90);
    EXPECT_TRUE(all_close(rf.get_parameters(), r.get_parameters()));
    // flipping rotation 0 negates the phase
    Bicop af = a;
    af.flip();
    EXPECT_EQ(af.get_rotation(), 0);
    EXPECT_NEAR(af.get_parameters()(1), -p1(1), 1e-14);
  }
  // reflecting the linear axis maps (a, b, mu) to (-b, -a, mu)
  Eigen::VectorXd q1(3), q3(3);
  q1 << 0.5, -0.4, 1.0;
  q3 << 0.4, -0.5, 1.0;
  Bicop c1(BicopFamily::cubic_sections, 0, q1, ac);
  Bicop c3(BicopFamily::cubic_sections, 0, q3, ac);
  Eigen::MatrixXd reflected = u;
  reflected.col(1) = 1 - u.col(1).array();
  EXPECT_TRUE(all_close(c1.pdf(reflected), c3.pdf(u), 1e-10, 1e-10));
}

TEST(test_circular, half_turn_case_has_zero_tau_and_strong_dependence)
{
  Eigen::VectorXd pars(2);
  pars << 0.95, pi;
  Bicop bc(BicopFamily::wrapped_cauchy, 0, pars, aa);
  EXPECT_NEAR(bc.parameters_to_tau(pars), -0.054, 0.01);
  EXPECT_NEAR(bc.get_beta(), -0.90, 0.02);
  auto sim = bc.simulate(2000, false, { 8 });
  EXPECT_GT(bc.loglik(sim) / 2000, 1.0); // far above independence
}

TEST(test_circular, fit_recovers_the_parameters)
{
  struct Case
  {
    BicopFamily family;
    int rotation;
    std::vector<double> pars;
    std::vector<std::string> types;
  };
  const std::vector<Case> cases = {
    { BicopFamily::cardioid, 0, { 0.35, 1.0 }, aa },
    { BicopFamily::wrapped_cauchy, 0, { 0.7, -3.0 }, aa }, // phase near -pi
    { BicopFamily::wrapped_cauchy, 90, { 0.9, pi }, aa },  // half turn
    { BicopFamily::von_mises, 0, { 3.0, 2.0 }, ac },       // on the cylinder
    { BicopFamily::quad_sections, 0, { 0.8, 0.5 }, ac },
    { BicopFamily::quad_sections, 0, { 0.7, -1.0 }, ca },
    { BicopFamily::cubic_sections, 0, { 0.8, -0.7, 1.5 }, ac },
  };
  for (const auto& c : cases) {
    Eigen::VectorXd pars =
      Eigen::Map<const Eigen::VectorXd>(c.pars.data(), c.pars.size());
    Bicop truth(c.family, c.rotation, pars, c.types);
    auto data = truth.simulate(4000, false, { 21 });

    Bicop fitted(c.family, c.rotation, Eigen::MatrixXd(), c.types);
    fitted.fit(data);
    const auto est = fitted.get_parameters();
    const std::string what = truth.str() + fitted.str();
    for (Eigen::Index k = 0; k + 1 < pars.size(); ++k) {
      EXPECT_NEAR(est(k), pars(k), 0.15) << what;
    }
    // the phase is compared on the circle and reported in [-pi, pi)
    const double dphase = est(pars.size() - 1) - pars(pars.size() - 1);
    EXPECT_NEAR(std::sin(dphase), 0.0, 0.1) << what;
    EXPECT_NEAR(std::cos(dphase), 1.0, 0.1) << what;
    EXPECT_GE(est(pars.size() - 1), -pi) << what;
    EXPECT_LT(est(pars.size() - 1), pi) << what;
    EXPECT_GE(fitted.get_loglik(), truth.loglik(data) - 5.0) << what;
  }
}

TEST(test_circular, select_picks_a_circular_family)
{
  Eigen::VectorXd pars(2);
  pars << 0.8, 1.0;
  Bicop truth(BicopFamily::wrapped_cauchy, 90, pars, aa);
  auto data = truth.simulate(1500, false, { 5 });

  Bicop bc;
  bc.set_var_types(aa);
  bc.select(data);
  EXPECT_TRUE(is_member(bc.get_family(), bicop_families::two_rotations))
    << bc.str();
  EXPECT_EQ(bc.get_rotation(), 90) << bc.str();
  EXPECT_EQ(bc.get_var_types(), aa);

  // with an explicit set on the cylinder in the reversed order
  Bicop truth2(BicopFamily::quad_sections, 0, pars, ca);
  auto data2 = truth2.simulate(1500, false, { 6 });
  FitControlsBicop controls({ BicopFamily::indep, BicopFamily::quad_sections });
  Bicop bc2;
  bc2.set_var_types(ca);
  bc2.select(data2, controls);
  EXPECT_EQ(bc2.get_family(), BicopFamily::quad_sections) << bc2.str();
  EXPECT_EQ(bc2.get_var_types(), ca);
}

TEST(test_circular, json_round_trip_with_parameters)
{
  Eigen::VectorXd pars(3);
  pars << 0.6, -0.3, 2.0;
  Bicop bc(BicopFamily::cubic_sections, 0, pars, ca);
  Bicop back(bc.to_json());
  EXPECT_EQ(back.get_family(), bc.get_family());
  EXPECT_EQ(back.get_var_types(), ca);
  auto u = tools_stats::simulate_uniform(20, 2, false, { 1 });
  EXPECT_TRUE(all_close(back.pdf(u), bc.pdf(u)));

  Eigen::VectorXd p2(2);
  p2 << 5.0, -1.0;
  Bicop vm(BicopFamily::von_mises, 90, p2, aa);
  Bicop vm_back(vm.to_json());
  EXPECT_EQ(vm_back.get_rotation(), 90);
  EXPECT_TRUE(all_close(vm_back.hfunc1(u), vm.hfunc1(u)));
}

TEST(test_circular, derivatives_are_available)
{
  Eigen::VectorXd pars(2);
  pars << 0.5, 0.7;
  Bicop bc(BicopFamily::wrapped_cauchy, 0, pars, aa);
  auto u = tools_stats::simulate_uniform(10, 2, false, { 9 });
  // the finite-difference fallback against a manual difference quotient
  const double h = 1e-6;
  Eigen::VectorXd plus = pars, minus = pars;
  plus(1) += h;
  minus(1) -= h;
  Eigen::MatrixXd plus_rows = plus.transpose().replicate(u.rows(), 1);
  Eigen::MatrixXd minus_rows = minus.transpose().replicate(u.rows(), 1);
  Eigen::VectorXd fd = (bc.pdf(u, plus_rows) - bc.pdf(u, minus_rows)) / (2 * h);
  EXPECT_TRUE(all_close(bc.pdf_deriv(u, "par2"), fd, 1e-4, 1e-4));
}

} // namespace test_circular
