// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// The nonparametric estimator on pairs with a circular variable: periodic
// kernels on circular axes, the probit-transformed estimator on linear ones,
// per-axis knots, and the grid's serialization.

#include "include/test_utils.hpp"
#include "gtest/gtest.h"
#include <boost/math/constants/constants.hpp>
#include <string>
#include <vector>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/nlohmann_json.hpp>
#include <vinecopulib/misc/tools_serialization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/vinecop/class.hpp>

namespace test_circular_tll {

using namespace vinecopulib;
using test_utils::all_close;

const std::vector<std::string> aa = { "a", "a" };
const std::vector<std::string> ac = { "a", "c" };
const std::vector<std::string> ca = { "c", "a" };
const std::vector<std::string> cc = { "c", "c" };
const double pi = boost::math::constants::pi<double>();

Eigen::VectorXd
par(std::initializer_list<double> values)
{
  Eigen::VectorXd p(static_cast<Eigen::Index>(values.size()));
  Eigen::Index i = 0;
  for (double v : values) {
    p(i++) = v;
  }
  return p;
}

FitControlsBicop
tll_controls(const std::string& method, double mult = 1.0, size_t m = 30)
{
  FitControlsBicop controls({ BicopFamily::tll });
  controls.set_nonparametric_method(method);
  controls.set_nonparametric_mult(mult);
  controls.set_nonparametric_grid_size(m);
  return controls;
}

Bicop
fit_tll(const Eigen::MatrixXd& u,
        const std::vector<std::string>& var_types,
        const std::string& method,
        double mult = 1.0,
        size_t m = 30)
{
  Bicop fit(BicopFamily::tll, 0, Eigen::MatrixXd(), var_types);
  fit.fit(u, tll_controls(method, mult, m));
  return fit;
}

//! the two circular test models and one circular-linear one
Bicop
von_mises()
{
  return Bicop(BicopFamily::von_mises, 0, par({ 2.0, 0.7 }), aa);
}
Bicop
cubic()
{
  return Bicop(BicopFamily::cubic_sections, 0, par({ 0.6, -0.4, 1.0 }), ac);
}

//! the knots recorded in a JSON representation
std::vector<Eigen::VectorXd>
read_knots(const nlohmann::json& json)
{
  std::vector<Eigen::VectorXd> knots;
  for (const auto& k : json["grid"]["knots"]) {
    knots.emplace_back(tools_serialization::json_to_matrix<double>(k));
  }
  return knots;
}

//! mean absolute error of a fit against the truth on a lattice
double
mean_abs_error(const Bicop& fit, const Bicop& truth, int m = 25)
{
  Eigen::MatrixXd u(m * m, 2);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      u(i * m + j, 0) = (i + 0.5) / m;
      u(i * m + j, 1) = (j + 0.5) / m;
    }
  }
  return (fit.pdf(u) - truth.pdf(u)).cwiseAbs().mean();
}

// -------------------------------------------------------------------------
// the grid follows the geometry

TEST(test_circular_tll, knots_follow_the_variable_types)
{
  Bicop linear(BicopFamily::tll);
  auto json = linear.to_json();
  EXPECT_FALSE(json.contains("grid"));

  // circular axes get uniform knots, and the independence values are kept
  Bicop circular(BicopFamily::tll, 0, Eigen::MatrixXd(), aa);
  json = circular.to_json();
  ASSERT_TRUE(json.contains("grid"));
  ASSERT_EQ(json["grid"]["knots"].size(), 2u);
  auto knots = read_knots(json);
  EXPECT_EQ(knots[0].size(), 30);
  EXPECT_NEAR(knots[0](0), 0.0, 1e-15);
  EXPECT_NEAR(knots[0](1), 1.0 / 29, 1e-15);
  EXPECT_NEAR(knots[0](15), 15.0 / 29, 1e-15);
  EXPECT_EQ(json["grid"]["types"], nlohmann::json(aa));
  auto u = tools_stats::simulate_uniform(20, 2, false, { 1 });
  EXPECT_TRUE(all_close(circular.pdf(u), Eigen::VectorXd::Ones(20), 1e-14));

  // mixed: only the circular axis is uniform
  Bicop mixed(BicopFamily::tll, 0, Eigen::MatrixXd(), ca);
  knots = read_knots(mixed.to_json());
  EXPECT_GT(knots[0](1), 1e-4); // normal grid, starts at 0
  EXPECT_LT(knots[0](1), 1e-2);
  EXPECT_NEAR(knots[1](1), 1.0 / 29, 1e-15);

  // changing the geometry rebuilds the knots
  Bicop retyped(BicopFamily::tll);
  retyped.set_var_types(aa);
  EXPECT_EQ(retyped.to_json()["grid"]["knots"],
            circular.to_json()["grid"]["knots"]);

  // values of a new shape land on the default knots; bad values, shapes,
  // and degrees of freedom are rejected
  Eigen::MatrixXd values = Eigen::MatrixXd::Constant(12, 20, 1.0);
  retyped.set_parameters(values);
  knots = read_knots(retyped.to_json());
  EXPECT_EQ(knots[0].size(), 12);
  EXPECT_EQ(knots[1].size(), 20);
  EXPECT_THROW(retyped.set_parameters(Eigen::MatrixXd::Constant(2, 5, 1.0)),
               std::runtime_error);
  EXPECT_THROW(retyped.set_parameters(Eigen::MatrixXd::Constant(5, 5, -1.0)),
               std::runtime_error);
}

TEST(test_circular_tll, json_round_trip_keeps_the_grid)
{
  auto u = von_mises().simulate(400, false, { 3 });
  auto fit = fit_tll(u, aa, "linear");
  auto json = fit.to_json();
  ASSERT_TRUE(json.contains("grid"));
  Bicop reloaded(json);
  EXPECT_EQ(reloaded.get_var_types(), aa);
  auto x = tools_stats::simulate_uniform(50, 2, false, { 5 });
  EXPECT_TRUE(all_close(reloaded.pdf(x), fit.pdf(x), 1e-14, 1e-14));
  EXPECT_TRUE(all_close(reloaded.hfunc1(x), fit.hfunc1(x), 1e-14, 1e-14));
  EXPECT_NEAR(reloaded.get_npars(), fit.get_npars(), 1e-12);

  // without the field the default knots of the types are rebuilt: the same
  // grid, since the estimator wrote the default knots
  json.erase("grid");
  Bicop rebuilt(json);
  EXPECT_TRUE(all_close(rebuilt.pdf(x), fit.pdf(x), 1e-14, 1e-14));

  // inconsistent fields are rejected with a message naming them
  auto bad = fit.to_json();
  bad["grid"]["types"] = cc;
  EXPECT_THROW(Bicop{ bad }, std::runtime_error);
  bad = fit.to_json();
  bad["grid"]["knots"][0]["shape"][0] = 29; // fewer knots than rows
  bad["grid"]["knots"][0]["data"].erase(0);
  EXPECT_THROW(Bicop{ bad }, std::runtime_error);
  bad = fit.to_json();
  bad["grid"]["knots"][0]["data"][3] = 0.0; // not increasing
  EXPECT_THROW(Bicop{ bad }, std::runtime_error);
}

// -------------------------------------------------------------------------
// the estimator

TEST(test_circular_tll, fit_is_a_periodic_copula_density)
{
  for (const char* method : { "constant", "linear", "quadratic" }) {
    for (const auto& truth : { von_mises(), cubic() }) {
      const auto var_types = truth.get_var_types();
      auto u = truth.simulate(800, false, { 11 });
      auto fit = fit_tll(u, var_types, method);
      EXPECT_EQ(fit.get_var_types(), var_types);
      EXPECT_EQ(fit.get_family(), BicopFamily::tll);

      // positive, periodic in every circular variable
      auto x = tools_stats::simulate_uniform(200, 2, false, { 13 });
      EXPECT_TRUE((fit.pdf(x).array() > 0).all());
      for (size_t j = 0; j < 2; ++j) {
        if (var_types[j] != "a") {
          continue;
        }
        Eigen::MatrixXd lo = x, hi = x;
        lo.col(j).setZero();
        hi.col(j).setOnes();
        // evaluation trims the arguments to [1e-10, 1 - 1e-10]
        EXPECT_TRUE(all_close(fit.pdf(lo), fit.pdf(hi), 1e-7, 1e-9))
          << method << " variable " << j + 1;
      }

      // uniform margins: the cdf at the boundary is the other coordinate
      Eigen::MatrixXd edge(9, 2);
      for (int k = 0; k < 9; ++k) {
        edge(k, 0) = 0.1 * (k + 1);
        edge(k, 1) = 1.0;
      }
      EXPECT_TRUE(all_close(fit.cdf(edge), edge.col(0), 1e-3, 1e-3)) << method;
      edge.col(1).swap(edge.col(0));
      EXPECT_TRUE(all_close(fit.cdf(edge), edge.col(1), 1e-3, 1e-3)) << method;

      // h-functions and inverses are anchored and invert each other
      Eigen::MatrixXd w = tools_stats::simulate_uniform(100, 2, false, { 17 });
      Eigen::MatrixXd v = w;
      v.col(1) = fit.hinv1(w);
      EXPECT_TRUE(all_close(fit.hfunc1(v), w.col(1), 1e-8, 1e-8)) << method;
      v = w;
      v.col(0) = fit.hinv2(w);
      EXPECT_TRUE(all_close(fit.hfunc2(v), w.col(0), 1e-8, 1e-8)) << method;

      // it is close to the truth and beats independence out of sample
      auto test = truth.simulate(2000, false, { 19 });
      EXPECT_GT(fit.pdf(test).array().log().mean(), 0.0) << method;
      EXPECT_LT(mean_abs_error(fit, truth), 0.25) << method;
      EXPECT_GT(fit.get_npars(), 1.0);
      EXPECT_LT(fit.get_npars(), 80.0);
      EXPECT_NEAR(fit.get_loglik(), fit.pdf(u).array().log().sum(), 1e-8);
    }
  }
}

TEST(test_circular_tll, fit_does_not_depend_on_the_cut)
{
  // shifting the cut of a circular variable shifts the fit with it, up to
  // the interpolation between knots
  auto u = von_mises().simulate(600, false, { 23 });
  auto fit = fit_tll(u, aa, "linear");
  Eigen::MatrixXd shifted = u;
  const double delta = 0.37;
  shifted.col(0) = (u.col(0).array() + delta).unaryExpr([](double x) {
    return x - std::floor(x);
  });
  auto fit_shifted = fit_tll(shifted, aa, "linear");

  auto x = tools_stats::simulate_uniform(300, 2, false, { 29 });
  Eigen::MatrixXd xs = x;
  xs.col(0) = (x.col(0).array() + delta).unaryExpr([](double v) {
    return v - std::floor(v);
  });
  EXPECT_TRUE(all_close(fit_shifted.pdf(xs), fit.pdf(x), 0.05, 0.02));
}

TEST(test_circular_tll, argument_order_and_flip_are_consistent)
{
  // circular-linear in either order gives the same density
  auto u = cubic().simulate(500, false, { 31 });
  auto fit_ac = fit_tll(u, ac, "linear");
  Eigen::MatrixXd u_swapped = u;
  u_swapped.col(0).swap(u_swapped.col(1));
  auto fit_ca = fit_tll(u_swapped, ca, "linear");
  EXPECT_EQ(fit_ca.get_var_types(), ca);
  auto x = tools_stats::simulate_uniform(100, 2, false, { 37 });
  Eigen::MatrixXd xs = x;
  xs.col(0).swap(xs.col(1));
  EXPECT_TRUE(all_close(fit_ca.pdf(xs), fit_ac.pdf(x), 1e-10, 1e-12));
  EXPECT_TRUE(all_close(fit_ca.hfunc2(xs), fit_ac.hfunc1(x), 1e-10, 1e-12));

  // flipping a fit is the fit of the swapped data, knots included
  Bicop flipped = fit_ac;
  flipped.flip();
  EXPECT_EQ(flipped.get_var_types(), ca);
  EXPECT_TRUE(all_close(flipped.pdf(xs), fit_ac.pdf(x), 1e-14, 1e-14));
  EXPECT_EQ(flipped.to_json()["grid"]["knots"],
            fit_ca.to_json()["grid"]["knots"]);
  flipped.flip();
  EXPECT_TRUE(all_close(flipped.pdf(x), fit_ac.pdf(x), 1e-14, 1e-14));
}

TEST(test_circular_tll, refinement_weights_and_bandwidth_extremes)
{
  auto truth = von_mises();
  auto u = truth.simulate(600, false, { 41 });
  auto x = tools_stats::simulate_uniform(100, 2, false, { 43 });

  // a finer grid changes the interpolant only a little
  auto coarse = fit_tll(u, aa, "constant", 1.0, 15);
  auto fine = fit_tll(u, aa, "constant", 1.0, 60);
  EXPECT_EQ(coarse.get_parameters().rows(), 15);
  EXPECT_EQ(fine.get_parameters().rows(), 60);
  EXPECT_TRUE(all_close(coarse.pdf(x), fine.pdf(x), 0.1, 0.05));

  // constant weights do not change the fit
  FitControlsBicop weighted = tll_controls("linear");
  weighted.set_weights(Eigen::VectorXd::Constant(u.rows(), 2.0));
  Bicop with_weights(BicopFamily::tll, 0, Eigen::MatrixXd(), aa);
  with_weights.fit(u, weighted);
  auto without = fit_tll(u, aa, "linear");
  EXPECT_TRUE(all_close(with_weights.pdf(x), without.pdf(x), 1e-10, 1e-10));

  // extreme bandwidths stay finite and normalized; the widest is nearly
  // independence, the narrowest very rough
  auto wide = fit_tll(u, aa, "constant", 50.0);
  auto narrow = fit_tll(u, aa, "linear", 0.05);
  for (const auto& fit : { wide, narrow }) {
    EXPECT_TRUE(fit.pdf(x).array().isFinite().all());
    EXPECT_TRUE((fit.pdf(x).array() >= 0).all());
    Eigen::MatrixXd edge(1, 2);
    edge << 0.4, 1.0;
    EXPECT_NEAR(fit.cdf(edge)(0), 0.4, 1e-2);
  }
  EXPECT_LT(
    mean_abs_error(wide, Bicop(BicopFamily::indep, 0, Eigen::MatrixXd(), aa)),
    0.1);
  EXPECT_GT(narrow.get_npars(), without.get_npars());
}

TEST(test_circular_tll, selection_and_vines_use_the_estimator)
{
  // selection with the nonparametric family alone
  auto truth = von_mises();
  auto u = truth.simulate(500, false, { 47 });
  Bicop selected(u, tll_controls("constant"), aa);
  EXPECT_EQ(selected.get_family(), BicopFamily::tll);
  EXPECT_EQ(selected.get_var_types(), aa);

  // ... and within the default family set it competes with the circulas
  Bicop any(u, FitControlsBicop(), aa);
  EXPECT_TRUE(family_accepts_var_types(any.get_family(), aa));

  // a mixed vine fitted with the nonparametric family on every edge
  std::vector<std::string> var_types{ "a", "a", "c" };
  Vinecop model(
    DVineStructure(std::vector<size_t>{ 1, 2, 3 }),
    { { von_mises(), cubic() },
      { Bicop(BicopFamily::cubic_sections, 0, par({ 0.8, 0.8, -0.5 }), ac) } },
    var_types);
  auto data = model.simulate(700, false, 1, { 53 });
  FitControlsVinecop controls;
  controls.set_family_set({ BicopFamily::tll });
  controls.set_nonparametric_method("linear");
  Vinecop fitted(data, RVineStructure(), var_types, controls);
  for (const auto& tree : fitted.get_all_pair_copulas()) {
    for (const auto& pc : tree) {
      EXPECT_EQ(pc.get_family(), BicopFamily::tll);
      EXPECT_TRUE(
        family_accepts_var_types(pc.get_family(), pc.get_var_types()));
    }
  }
  auto x = tools_stats::simulate_uniform(40, 3, false, { 59 });
  for (size_t j : { 0, 1 }) {
    Eigen::MatrixXd lo = x, hi = x;
    lo.col(j).setZero();
    hi.col(j).setOnes();
    EXPECT_TRUE(all_close(fitted.pdf(lo), fitted.pdf(hi), 1e-8, 1e-10));
  }
  auto w = fitted.rosenblatt(x);
  EXPECT_TRUE(all_close(fitted.inverse_rosenblatt(w), x, 1e-6, 1e-6));
  EXPECT_GT(fitted.loglik(data), 0.0);

  // the vine survives a JSON round trip with its grids
  Vinecop reloaded(fitted.to_json());
  EXPECT_TRUE(all_close(reloaded.pdf(x), fitted.pdf(x), 1e-12, 1e-12));
}
}
