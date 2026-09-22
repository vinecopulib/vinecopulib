// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// The circular variable type and the parametric circular families: their
// registration and eligibility, golden values, identities, conventions,
// fitting, and selection at the pair level.

#include "include/circular_golden.hpp"
#include "include/test_utils.hpp"
#include "gtest/gtest.h"
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/bicop/tools_select.hpp>
#include <vinecopulib/misc/nlohmann_json.hpp>
#include <vinecopulib/misc/tools_circular.hpp>
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
    { "cubic_sections", BicopFamily::cubic_sections },
  };
  for (const auto& entry : table) {
    if (key.rfind(entry.first, 0) == 0) {
      return entry.second;
    }
  }
  throw std::runtime_error("unknown golden key " + key);
}

Bicop
make(BicopFamily family,
     int rotation,
     std::initializer_list<double> pars,
     const std::vector<std::string>& types)
{
  Eigen::VectorXd p(static_cast<Eigen::Index>(pars.size()));
  Eigen::Index i = 0;
  for (double v : pars) {
    p(i++) = v;
  }
  return Bicop(family, rotation, p, types);
}

// -------------------------------------------------------------------------
// registration and eligibility

TEST(test_circular, families_are_registered_with_their_groups)
{
  for (auto fam : bicop_families::circular) {
    EXPECT_EQ(get_family_enum(get_family_name(fam)), fam);
    EXPECT_TRUE(is_member(fam, bicop_families::all));
    EXPECT_TRUE(is_member(fam, bicop_families::parametric));
    EXPECT_FALSE(is_member(fam, bicop_families::itau));
    EXPECT_FALSE(is_member(fam, bicop_families::analytic_derivs));
    // the group decides the rotations and the parameter count
    const bool two_rot = is_member(fam, bicop_families::two_rotations);
    EXPECT_NE(two_rot, is_member(fam, bicop_families::cylindrical));
    EXPECT_NE(two_rot, is_member(fam, bicop_families::rotationless));
    EXPECT_NO_THROW(make_bicop(fam, 0));
    EXPECT_EQ(two_rot, [&] {
      try {
        make_bicop(fam, 90);
        return true;
      } catch (const std::runtime_error&) {
        return false;
      }
    }());
    EXPECT_THROW(make_bicop(fam, 180), std::runtime_error);
    EXPECT_THROW(make_bicop(fam, 270), std::runtime_error);
  }
  EXPECT_EQ(get_family_name(BicopFamily::cardioid), "Cardioid");
  EXPECT_EQ(get_family_name(BicopFamily::wrapped_cauchy), "Wrapped Cauchy");
  EXPECT_EQ(get_family_name(BicopFamily::von_mises), "von Mises");
  EXPECT_EQ(get_family_name(BicopFamily::cubic_sections), "Cubic sections");
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

TEST(test_circular, tokens_are_enforced_by_bicop_and_vinecop)
{
  Bicop indep(BicopFamily::indep, 0, Eigen::MatrixXd(), aa);
  EXPECT_EQ(indep.get_var_types(), aa);
  indep.set_var_types(ac);
  EXPECT_EQ(indep.get_var_types(), ac);
  EXPECT_THROW(indep.set_var_types({ "a", "d" }), std::runtime_error);
  EXPECT_THROW(indep.set_var_types({ "x", "c" }), std::runtime_error);

  // a linear family cannot take a circular variable; circular families need
  // one; cylindrical ones exactly one
  EXPECT_THROW(Bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Zero(1), aa),
               std::runtime_error);
  Bicop gauss(BicopFamily::gaussian, 0, Eigen::VectorXd::Zero(1));
  EXPECT_THROW(gauss.set_var_types(ca), std::runtime_error);
  for (auto fam : bicop_families::circular) {
    EXPECT_THROW(Bicop(fam, 0, Eigen::MatrixXd(), cc), std::runtime_error);
    EXPECT_THROW((Bicop{ fam }), std::runtime_error);
  }
  for (auto fam : bicop_families::cylindrical) {
    EXPECT_THROW(Bicop(fam, 0, Eigen::MatrixXd(), aa), std::runtime_error);
    EXPECT_NO_THROW(Bicop(fam, 0, Eigen::MatrixXd(), ca));
  }

  // the vine accepts the token and rejects mixes with discrete variables
  auto u = tools_stats::simulate_uniform(50, 3, false, { 3 });
  Vinecop vc(3);
  vc.set_var_types({ "a", "c", "c" });
  EXPECT_THROW(vc.set_var_types({ "a", "d", "c" }), std::runtime_error);
  EXPECT_THROW(vc.set_var_types({ "a", "x", "c" }), std::runtime_error);
  EXPECT_NO_THROW(vc.select(u));
  EXPECT_EQ(vc.get_var_types(), (std::vector<std::string>{ "a", "c", "c" }));
}

TEST(test_circular, candidate_generation_filters_by_geometry)
{
  auto u = tools_stats::simulate_uniform(200, 2, false, { 2 });
  FitControlsBicop controls;

  auto linear = tools_select::create_candidate_bicops(u, controls, cc);
  for (const auto& bc : linear) {
    EXPECT_FALSE(is_member(bc.get_family(), bicop_families::circular));
  }
  EXPECT_EQ(linear.size(),
            tools_select::create_candidate_bicops(u, controls).size());

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
  for (const auto& bc :
       tools_select::create_candidate_bicops(u, controls, aa)) {
    EXPECT_EQ(bc.get_rotation(), 0);
  }

  // an explicit set is intersected with the eligible families; an empty
  // intersection names the eligible ones
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
    { BicopFamily::gaussian, BicopFamily::indep, BicopFamily::cubic_sections });
  EXPECT_EQ(tools_select::get_candidate_families(controls, ac),
            (std::vector<BicopFamily>{ BicopFamily::indep,
                                       BicopFamily::cubic_sections }));
  EXPECT_EQ(tools_select::get_candidate_families(controls, aa),
            (std::vector<BicopFamily>{ BicopFamily::indep }));
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
      EXPECT_TRUE(std::isnan(bc.parameters_to_tau(pars))) << what;
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
  const std::vector<Bicop> cases = {
    make(BicopFamily::cardioid, 0, { 0.4, 1.0 }, aa),
    make(BicopFamily::cardioid, 90, { 0.5, -2.0 }, aa),
    make(BicopFamily::wrapped_cauchy, 0, { 0.7, pi }, aa),
    make(BicopFamily::wrapped_cauchy, 90, { 0.95, 0.3 }, aa),
    make(BicopFamily::von_mises, 0, { 4.0, 2.5 }, aa),
    make(BicopFamily::von_mises, 90, { 30.0, -0.7 }, aa),
    make(BicopFamily::cubic_sections, 0, { 0.8, 0.8, 1.2 }, ac),
    make(BicopFamily::cubic_sections, 0, { 0.9, -0.6, -2.0 }, ca),
  };
  for (const auto& bc : cases) {
    const std::string what = bc.str();
    const auto pars = bc.get_parameters();

    // h(hinv) = id, hinv in [0, 1]
    Eigen::MatrixXd v = u;
    v.col(1) = bc.hinv1(u);
    EXPECT_TRUE((v.col(1).array() >= 0).all() && (v.col(1).array() <= 1).all())
      << what;
    EXPECT_TRUE(all_close(bc.hfunc1(v), u.col(1), 1e-7, 1e-7)) << what;
    v = u;
    v.col(0) = bc.hinv2(u);
    EXPECT_TRUE(all_close(bc.hfunc2(v), u.col(0), 1e-7, 1e-7)) << what;

    // periodic in every circular coordinate; the arguments are trimmed to
    // [1e-10, 1 - 1e-10] on the way in, so steep densities differ slightly
    for (Eigen::Index j = 0; j < 2; ++j) {
      if (bc.get_var_types()[j] != "a") {
        continue;
      }
      Eigen::MatrixXd u0 = u, u1 = u;
      u0.col(j).setConstant(1e-12);
      u1.col(j).setConstant(1 - 1e-12);
      EXPECT_TRUE(all_close(bc.pdf(u0), bc.pdf(u1), 1e-6, 1e-8)) << what;
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

    // JSON keeps family, rotation, types, and values
    Bicop back(bc.to_json());
    EXPECT_EQ(back.get_rotation(), bc.get_rotation()) << what;
    EXPECT_EQ(back.get_var_types(), bc.get_var_types()) << what;
    EXPECT_TRUE(all_close(back.pdf(u), bc.pdf(u))) << what;

    // no Kendall's tau in either direction, zero tail dependence
    EXPECT_TRUE(std::isnan(bc.get_tau())) << what;
    EXPECT_THROW(bc.tau_to_parameters(0.3), std::runtime_error) << what;
    EXPECT_TRUE(bc.get_taildep().isZero()) << what;
  }
}

TEST(test_circular, circular_helpers_handle_degenerate_inputs)
{
  using namespace tools_circular;
  EXPECT_EQ(von_mises_a(0.0), 0.0);
  EXPECT_EQ(von_mises_a_inverse(0.0, 500.0), 0.0);
  EXPECT_EQ(von_mises_a_inverse(1.0, 500.0), 500.0); // clipped to the bound
  EXPECT_EQ(bessel_i_ratios(0.0, 3), (std::vector<double>{ 0.0, 0.0, 0.0 }));
  Eigen::VectorXd theta = Eigen::VectorXd::LinSpaced(4, 0.0, 3.0);
  EXPECT_EQ(weighted_resultant(theta, Eigen::VectorXd::Zero(4)),
            std::complex<double>(0.0, 0.0));
  EXPECT_NEAR(wrap_pi(3 * pi), -pi, 1e-12);
}

TEST(test_circular, phase_conventions_and_the_half_turn_case)
{
  auto u = tools_stats::simulate_uniform(50, 2, false, { 4 });
  for (auto fam : bicop_families::two_rotations) {
    // phases differing by 2 pi define the same copula
    Bicop a = make(fam, 0, { 0.4, 0.4 }, aa);
    Bicop b = make(fam, 0, { 0.4, 0.4 + 2 * pi }, aa);
    EXPECT_TRUE(all_close(a.pdf(u), b.pdf(u), 1e-10, 1e-10));
    // the rotation-90 copula is exchangeable: flipping keeps rotation and
    // phase; flipping rotation 0 negates the phase
    Bicop rf = make(fam, 90, { 0.4, 0.4 }, aa);
    rf.flip();
    EXPECT_EQ(rf.get_rotation(), 90);
    EXPECT_NEAR(rf.get_parameters()(1), 0.4, 1e-14);
    a.flip();
    EXPECT_EQ(a.get_rotation(), 0);
    EXPECT_NEAR(a.get_parameters()(1), -0.4, 1e-14);
  }
  // reflecting the linear axis maps (a, b, mu) to (-b, -a, mu)
  Bicop c1 = make(BicopFamily::cubic_sections, 0, { 0.5, -0.4, 1.0 }, ac);
  Bicop c3 = make(BicopFamily::cubic_sections, 0, { 0.4, -0.5, 1.0 }, ac);
  Eigen::MatrixXd reflected = u;
  reflected.col(1) = 1 - u.col(1).array();
  EXPECT_TRUE(all_close(c1.pdf(reflected), c3.pdf(u), 1e-10, 1e-10));

  // the half-turn wrapped Cauchy: strong dependence (Kendall's tau would be
  // near 0)
  Bicop half_turn = make(BicopFamily::wrapped_cauchy, 0, { 0.95, pi }, aa);
  EXPECT_NEAR(half_turn.get_beta(), -0.90, 0.02);
  auto sim = half_turn.simulate(2000, false, { 8 });
  EXPECT_GT(half_turn.loglik(sim) / 2000, 1.0); // far above independence
}

TEST(test_circular, fit_recovers_the_parameters)
{
  const std::vector<Bicop> truths = {
    make(BicopFamily::cardioid, 0, { 0.35, 1.0 }, aa),
    make(BicopFamily::wrapped_cauchy, 0, { 0.7, -3.0 }, aa), // phase near -pi
    make(BicopFamily::wrapped_cauchy, 90, { 0.9, pi }, aa),  // half turn
    make(BicopFamily::von_mises, 0, { 3.0, 2.0 }, ac),       // on the cylinder
    make(BicopFamily::cubic_sections, 0, { 0.8, 0.8, 0.5 }, ac),
    make(BicopFamily::cubic_sections, 0, { 0.7, 0.7, -1.0 }, ca),
    make(BicopFamily::cubic_sections, 0, { 0.8, -0.7, 1.5 }, ac),
    // a + b < 0: the moment start lands in the equivalent frame with a < 0
    // and is normalized to a >= 0 with the phase shifted by pi
    make(BicopFamily::cubic_sections, 0, { 0.3, -0.8, 1.0 }, ca),
  };
  for (const auto& truth : truths) {
    const Eigen::VectorXd pars = truth.get_parameters();
    auto data = truth.simulate(4000, false, { 21 });
    Bicop fitted(truth.get_family(),
                 truth.get_rotation(),
                 Eigen::MatrixXd(),
                 truth.get_var_types());
    fitted.fit(data);
    const auto est = fitted.get_parameters();
    const std::string what = truth.str() + fitted.str();
    for (Eigen::Index k = 0; k + 1 < pars.size(); ++k) {
      EXPECT_NEAR(est(k), pars(k), 0.15) << what;
    }
    // the phase is compared on the circle and reported in [-pi, pi)
    const double dphase = est(pars.size() - 1) - pars(pars.size() - 1);
    EXPECT_NEAR(std::sin(dphase), 0.0, 0.15) << what;
    EXPECT_NEAR(std::cos(dphase), 1.0, 0.15) << what;
    EXPECT_GE(est(pars.size() - 1), -pi) << what;
    EXPECT_LT(est(pars.size() - 1), pi) << what;
    EXPECT_GE(fitted.get_loglik(), truth.loglik(data) - 5.0) << what;
  }
  // only maximum likelihood is available
  FitControlsBicop itau;
  itau.set_parametric_method("itau");
  Bicop bc(BicopFamily::von_mises, 0, Eigen::MatrixXd(), aa);
  EXPECT_THROW(bc.fit(truths[0].simulate(100, false, { 1 }), itau),
               std::runtime_error);
}

TEST(test_circular, select_picks_a_circular_family)
{
  Bicop truth = make(BicopFamily::wrapped_cauchy, 90, { 0.8, 1.0 }, aa);
  Bicop bc;
  bc.set_var_types(aa);
  bc.select(truth.simulate(1500, false, { 5 }));
  EXPECT_TRUE(is_member(bc.get_family(), bicop_families::two_rotations))
    << bc.str();
  EXPECT_EQ(bc.get_rotation(), 90) << bc.str();
  EXPECT_EQ(bc.get_var_types(), aa);

  // with an explicit set on the cylinder in the reversed order
  Bicop truth2 = make(BicopFamily::cubic_sections, 0, { 0.8, 0.8, 1.0 }, ca);
  FitControlsBicop controls(
    { BicopFamily::indep, BicopFamily::cubic_sections });
  Bicop bc2;
  bc2.set_var_types(ca);
  bc2.select(truth2.simulate(1500, false, { 6 }), controls);
  EXPECT_EQ(bc2.get_family(), BicopFamily::cubic_sections) << bc2.str();
  EXPECT_EQ(bc2.get_var_types(), ca);
}

} // namespace test_circular
