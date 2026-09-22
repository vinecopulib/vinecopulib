// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Vine copulas with supplied structures whose variables are circular: the
// geometry of every edge, evaluation, transforms, simulation, refitting, and
// structure operations.

#include "include/test_utils.hpp"
#include "gtest/gtest.h"
#include <boost/math/constants/constants.hpp>
#include <map>
#include <string>
#include <vector>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/class.hpp>

namespace test_circular_vinecop {

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

//! A D-vine on (1, 2, 3) with variables 1 and 2 circular and variable 3
//! linear: a circula in tree 1, cylindrical families on the two edges that
//! pair a circular with a linear variable, including the tree-2 edge whose
//! circular conditioned variable arrives through an h-function.
struct MixedVine
{
  DVineStructure structure{ std::vector<size_t>{ 1, 2, 3 } };
  std::vector<std::string> var_types{ "a", "a", "c" };
  Bicop c12{ BicopFamily::von_mises, 90, par({ 2.0, 0.7 }), aa };
  Bicop c23{ BicopFamily::cubic_sections, 0, par({ 0.6, -0.4, 1.0 }), ac };
  Bicop c13_2{ BicopFamily::quad_sections, 0, par({ 0.8, -0.5 }), ac };

  std::vector<std::vector<Bicop>> pair_copulas() const
  {
    return { { c12, c23 }, { c13_2 } };
  }

  Vinecop model() const
  {
    return Vinecop(structure, pair_copulas(), var_types);
  }

  //! the joint density assembled by hand from the three pair copulas
  Eigen::VectorXd density(const Eigen::MatrixXd& u) const
  {
    Eigen::MatrixXd u12 = u.leftCols(2);
    Eigen::MatrixXd u23 = u.rightCols(2);
    Eigen::MatrixXd u13_2(u.rows(), 2);
    u13_2.col(0) = c12.hfunc2(u12); // F(u1 | u2)
    u13_2.col(1) = c23.hfunc1(u23); // F(u3 | u2)
    return c12.pdf(u12).array() * c23.pdf(u23).array() *
           c13_2.pdf(u13_2).array();
  }
};

//! A D-vine on (1, 2, 3) whose only circular variable, 2, is conditioned on in
//! tree 2, where both conditioned variables are linear.
struct ConditioningOnlyVine
{
  DVineStructure structure{ std::vector<size_t>{ 1, 2, 3 } };
  std::vector<std::string> var_types{ "c", "a", "c" };
  Bicop c12{ BicopFamily::cubic_sections, 0, par({ 0.5, 0.3, 0.4 }), ca };
  Bicop c23{ BicopFamily::quad_sections, 0, par({ 0.7, 2.0 }), ac };
  Bicop c13_2{ BicopFamily::gaussian, 0, par({ 0.5 }), cc };

  Vinecop model() const
  {
    return Vinecop(structure, { { c12, c23 }, { c13_2 } }, var_types);
  }
};

//! A D-vine on (1, 2, 3) with one circular variable, 1, whose transform is a
//! conditioned argument of the tree-2 edge; the other pair is linear.
struct CircularLinearLinearVine
{
  DVineStructure structure{ std::vector<size_t>{ 1, 2, 3 } };
  std::vector<std::string> var_types{ "a", "c", "c" };
  Bicop c12{ BicopFamily::quad_sections, 0, par({ 0.7, 0.2 }), ac };
  Bicop c23{ BicopFamily::gaussian, 0, par({ -0.6 }), cc };
  Bicop c13_2{ BicopFamily::cubic_sections, 0, par({ 0.3, 0.8, 2.5 }), ac };

  Vinecop model() const
  {
    return Vinecop(structure, { { c12, c23 }, { c13_2 } }, var_types);
  }
};

//! A five-dimensional vine on a simulated structure with three circular
//! variables, its pair copulas chosen by the geometry of each edge.
Vinecop
larger_mixed_vine()
{
  const size_t d = 5;
  std::vector<std::string> var_types{ "a", "c", "a", "c", "a" };
  auto structure = RVineStructure::simulate(d, false, { 7 });
  // independence is eligible on every edge, so this assigns every pair its
  // geometry, which then picks the family
  Vinecop vc(structure, Vinecop::make_pair_copula_store(d), var_types);
  auto pcs = vc.get_all_pair_copulas();
  for (size_t t = 0; t < pcs.size(); ++t) {
    for (size_t e = 0; e < pcs[t].size(); ++e) {
      const auto types = pcs[t][e].get_var_types();
      const double shift = 0.3 * static_cast<double>(t + e);
      if (types == aa) {
        pcs[t][e] =
          Bicop(BicopFamily::wrapped_cauchy, 0, par({ 0.5, 0.3 + shift }), aa);
      } else if (types == cc) {
        pcs[t][e] = Bicop(BicopFamily::gaussian, 0, par({ 0.4 }), cc);
      } else {
        pcs[t][e] = Bicop(
          BicopFamily::quad_sections, 0, par({ 0.5, 0.8 + shift }), types);
      }
    }
  }
  vc.set_all_pair_copulas(pcs);
  return vc;
}

std::vector<std::vector<std::string>>
edge_types(const Vinecop& vc)
{
  std::vector<std::vector<std::string>> types;
  for (const auto& tree : vc.get_all_pair_copulas()) {
    for (const auto& pc : tree) {
      types.push_back(pc.get_var_types());
    }
  }
  return types;
}

// -------------------------------------------------------------------------
// geometry of the edges

TEST(test_circular_vinecop, conditioned_variables_carry_their_geometry)
{
  MixedVine mixed;
  auto vc = mixed.model();
  EXPECT_EQ(vc.get_order(), (std::vector<size_t>{ 1, 2, 3 }));
  EXPECT_EQ(edge_types(vc),
            (std::vector<std::vector<std::string>>{ aa, ac, ac }));

  // a circular variable that is only conditioned on leaves the pair linear
  ConditioningOnlyVine conditioning;
  EXPECT_EQ(edge_types(conditioning.model()),
            (std::vector<std::vector<std::string>>{ ca, ac, cc }));
  EXPECT_EQ(edge_types(CircularLinearLinearVine().model()),
            (std::vector<std::vector<std::string>>{ ac, cc, ac }));

  // setting the types afterwards propagates the same way
  Vinecop relabeled(mixed.structure, Vinecop::make_pair_copula_store(3));
  relabeled.set_var_types(mixed.var_types);
  EXPECT_EQ(edge_types(relabeled),
            (std::vector<std::vector<std::string>>{ aa, ac, ac }));
}

TEST(test_circular_vinecop, ineligible_pair_copulas_are_rejected_by_edge)
{
  MixedVine mixed;
  auto pcs = mixed.pair_copulas();
  pcs[1][0] = Bicop(BicopFamily::gaussian, 0, par({ 0.5 }));
  try {
    Vinecop vc(mixed.structure, pcs, mixed.var_types);
    FAIL() << "a Gaussian on a circular-linear edge must be rejected";
  } catch (const std::runtime_error& err) {
    EXPECT_NE(std::string(err.what()).find("tree 2, edge 1"), std::string::npos)
      << err.what();
  }

  // the setter validates as well, and a valid store replaces the old one
  auto vc = mixed.model();
  EXPECT_THROW(vc.set_all_pair_copulas(pcs), std::runtime_error);
  pcs[1][0] = Bicop(BicopFamily::indep);
  vc.set_all_pair_copulas(pcs);
  EXPECT_EQ(vc.get_family(1, 0), BicopFamily::indep);
  EXPECT_EQ(vc.get_pair_copula(1, 0).get_var_types(), ac);

  // a cylindrical family on a circular-circular edge is equally rejected
  pcs[0][0] = Bicop(BicopFamily::quad_sections, 0, par({ 0.5, 0.0 }), ac);
  EXPECT_THROW(Vinecop(mixed.structure, pcs, mixed.var_types),
               std::runtime_error);
}

TEST(test_circular_vinecop, structure_selection_is_still_rejected)
{
  MixedVine mixed;
  auto u = mixed.model().simulate(100, false, 1, { 11 });
  EXPECT_THROW(Vinecop(u, mixed.structure, mixed.var_types),
               std::runtime_error);
  EXPECT_THROW(Vinecop(u, RVineStructure(), mixed.var_types),
               std::runtime_error);
}

// -------------------------------------------------------------------------
// evaluation

TEST(test_circular_vinecop, pdf_is_the_product_of_the_pair_densities)
{
  MixedVine mixed;
  auto vc = mixed.model();
  auto u = tools_stats::simulate_uniform(200, 3, false, { 13 });
  EXPECT_TRUE(all_close(vc.pdf(u), mixed.density(u), 1e-10, 1e-12));
  EXPECT_TRUE(all_close(vc.pdf(u, 2), vc.pdf(u), 1e-14, 1e-14));
  EXPECT_NEAR(vc.loglik(u), vc.pdf(u).array().log().sum(), 1e-10);

  auto full = vc.pdf_full(u);
  EXPECT_TRUE(all_close(full.pdf, vc.pdf(u), 1e-14, 1e-14));
}

TEST(test_circular_vinecop,
     integrating_out_the_linear_variable_recovers_the_circula)
{
  MixedVine mixed;
  auto vc = mixed.model();
  // Simpson's rule in u3 on a fine grid; the integrand is a smooth
  // polynomial-times-cosine in u3
  const int m = 400;
  Eigen::MatrixXd u(m + 1, 3);
  Eigen::MatrixXd u12(4, 2);
  u12 << 0.1, 0.2, 0.9, 0.35, 0.5, 0.5, 0.05, 0.95;
  Eigen::VectorXd c12 = mixed.c12.pdf(u12);
  for (int r = 0; r < 4; ++r) {
    for (int i = 0; i <= m; ++i) {
      u(i, 0) = u12(r, 0);
      u(i, 1) = u12(r, 1);
      u(i, 2) = static_cast<double>(i) / m;
    }
    Eigen::VectorXd f = vc.pdf(u);
    double integral = f(0) + f(m);
    for (int i = 1; i < m; ++i) {
      integral += (i % 2 == 1 ? 4.0 : 2.0) * f(i);
    }
    integral /= 3.0 * m;
    EXPECT_NEAR(integral, c12(r), 1e-5) << "row " << r;
  }
}

TEST(test_circular_vinecop, density_is_periodic_in_every_circular_variable)
{
  auto check = [](const Vinecop& vc, const std::vector<size_t>& circular) {
    const size_t d = vc.get_dim();
    auto u = tools_stats::simulate_uniform(50, d, false, { 17 });
    for (size_t j : circular) {
      Eigen::MatrixXd lo = u, hi = u;
      lo.col(j).setZero();
      hi.col(j).setOnes();
      EXPECT_TRUE(all_close(vc.pdf(lo), vc.pdf(hi), 1e-6, 1e-8))
        << "variable " << j + 1;
    }
  };
  check(MixedVine().model(), { 0, 1 });
  check(ConditioningOnlyVine().model(), { 1 });
  check(CircularLinearLinearVine().model(), { 0 });
  check(larger_mixed_vine(), { 0, 2, 4 });
}

TEST(test_circular_vinecop, transforms_round_trip_and_simulation_is_uniform)
{
  auto check = [](const Vinecop& vc) {
    const size_t n = 2000;
    const size_t d = vc.get_dim();
    auto u = vc.simulate(n, false, 1, { 19 });
    EXPECT_TRUE(all_close(vc.simulate(n, false, 2, { 19 }), u, 1e-12, 1e-12));
    EXPECT_TRUE((u.array() > 0).all() && (u.array() < 1).all());

    // every margin is uniform: the sorted sample stays close to the identity
    for (size_t j = 0; j < d; ++j) {
      std::vector<double> col(u.col(j).data(), u.col(j).data() + n);
      std::sort(col.begin(), col.end());
      double max_dev = 0;
      for (size_t i = 0; i < n; ++i) {
        const double target =
          (static_cast<double>(i) + 0.5) / static_cast<double>(n);
        max_dev = std::max(max_dev, std::abs(col[i] - target));
      }
      EXPECT_LT(max_dev, 0.05) << "variable " << j + 1;
    }

    auto w = vc.rosenblatt(u);
    EXPECT_TRUE(all_close(vc.inverse_rosenblatt(w), u, 1e-6, 1e-8));
    EXPECT_TRUE(all_close(vc.rosenblatt(u, 2), w, 1e-12, 1e-12));
    EXPECT_TRUE(all_close(
      vc.inverse_rosenblatt(w, 2), vc.inverse_rosenblatt(w), 1e-12, 1e-12));

    // the Rosenblatt transform of the model's own sample is uniform too
    for (size_t j = 0; j < d; ++j) {
      EXPECT_NEAR(w.col(j).mean(), 0.5, 0.03) << "variable " << j + 1;
    }
  };
  check(MixedVine().model());
  check(CircularLinearLinearVine().model());
  check(larger_mixed_vine());
}

TEST(test_circular_vinecop, cdf_is_a_distribution_function)
{
  auto vc = MixedVine().model();
  Eigen::MatrixXd u(4, 3);
  u << 0.3, 0.6, 0.2, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9;
  auto p = vc.cdf(u, 20000, 1, { 23 });
  EXPECT_GT(p(0), 0.0);
  EXPECT_LT(p(0), 0.2); // bounded by the smallest coordinate
  EXPECT_NEAR(p(1), 1.0, 1e-12);
  EXPECT_GT(p(2), 0.0);
  EXPECT_LT(p(2), 0.5);
  EXPECT_GT(p(3), p(2) + 0.2); // increasing in every coordinate
  EXPECT_LT(p(3), 0.9);
}

TEST(test_circular_vinecop, scores_match_a_finite_difference_of_the_loglik)
{
  MixedVine mixed;
  auto vc = mixed.model();
  auto u = vc.simulate(200, false, 1, { 29 });
  // the full (not step-wise) derivatives of the joint log-likelihood
  Eigen::MatrixXd scores = vc.scores(u, false);
  Eigen::MatrixXd hessian = vc.hessian(u, false);
  const auto pcs = mixed.pair_copulas();

  Eigen::Index k = 0;
  for (size_t t = 0; t < pcs.size(); ++t) {
    for (size_t e = 0; e < pcs[t].size(); ++e) {
      Eigen::VectorXd theta = pcs[t][e].get_parameters();
      for (Eigen::Index i = 0; i < theta.size(); ++i, ++k) {
        const double h = 1e-5;
        auto shifted = [&](double delta) {
          auto store = pcs;
          Eigen::VectorXd th = theta;
          th(i) += delta;
          store[t][e].set_parameters(th);
          return Vinecop(mixed.structure, store, mixed.var_types).loglik(u);
        };
        const double fd = (shifted(h) - shifted(-h)) / (2 * h);
        EXPECT_NEAR(scores.col(k).sum(), fd, 1e-3 * (1 + std::abs(fd)))
          << "tree " << t << ", edge " << e << ", parameter " << i;
      }
    }
  }
  EXPECT_EQ(scores.cols(), k);
  EXPECT_EQ(hessian.rows(), k);
  EXPECT_EQ(hessian.cols(), k);
}

// -------------------------------------------------------------------------
// fitting on a supplied structure

TEST(test_circular_vinecop, fit_on_a_fixed_structure_recovers_the_parameters)
{
  MixedVine mixed;
  auto truth = mixed.model();
  auto u = truth.simulate(3000, false, 1, { 31 });

  // start from the true families with default parameters
  auto pcs = mixed.pair_copulas();
  for (auto& tree : pcs) {
    for (auto& pc : tree) {
      pc = Bicop(pc.get_family(),
                 pc.get_rotation(),
                 Eigen::MatrixXd(),
                 pc.get_var_types());
    }
  }
  Vinecop fitted(mixed.structure, pcs, mixed.var_types);
  fitted.fit(u, FitControlsBicop(), 2);

  EXPECT_EQ(edge_types(fitted), edge_types(truth));
  auto on_circle = [](double a, double b) {
    double diff = std::fmod(std::abs(a - b), 2 * pi);
    return std::min(diff, 2 * pi - diff);
  };
  // von Mises: concentration, phase
  Eigen::VectorXd p12 = fitted.get_pair_copula(0, 0).get_parameters();
  EXPECT_NEAR(p12(0), 2.0, 0.3);
  EXPECT_LT(on_circle(p12(1), 0.7), 0.15);
  EXPECT_EQ(fitted.get_pair_copula(0, 0).get_rotation(), 90);
  // cubic sections: a, b, phase
  Eigen::VectorXd p23 = fitted.get_pair_copula(0, 1).get_parameters();
  EXPECT_NEAR(p23(0), 0.6, 0.2);
  EXPECT_NEAR(p23(1), -0.4, 0.2);
  EXPECT_LT(on_circle(p23(2), 1.0), 0.3);
  // quadratic sections in tree 2: a, phase
  Eigen::VectorXd p13 = fitted.get_pair_copula(1, 0).get_parameters();
  EXPECT_NEAR(p13(0), 0.8, 0.2);
  EXPECT_LT(on_circle(p13(1), -0.5), 0.3);

  EXPECT_GT(fitted.get_loglik(), 0.0);
  EXPECT_NEAR(fitted.get_loglik(), fitted.loglik(u), 1e-8);
  EXPECT_GT(truth.loglik(u) + 30, fitted.get_loglik()); // no overfit
}

// -------------------------------------------------------------------------
// structure operations, views, and serialization

TEST(test_circular_vinecop, truncation_keeps_the_geometry)
{
  MixedVine mixed;
  auto vc = mixed.model();
  vc.truncate(1);
  EXPECT_EQ(vc.get_all_pair_copulas().size(), 1u);
  EXPECT_EQ(vc.get_var_types(), mixed.var_types);
  EXPECT_EQ(edge_types(vc), (std::vector<std::vector<std::string>>{ aa, ac }));

  auto u = tools_stats::simulate_uniform(50, 3, false, { 37 });
  Eigen::VectorXd expected = mixed.c12.pdf(u.leftCols(2)).array() *
                             mixed.c23.pdf(u.rightCols(2)).array();
  EXPECT_TRUE(all_close(vc.pdf(u), expected, 1e-12, 1e-12));
}

TEST(test_circular_vinecop, reoriented_transforms_match_the_materialized_model)
{
  auto model = larger_mixed_vine();
  const size_t d = model.get_dim();
  std::vector<size_t> conditioning_set{ 1, 3 };
  Vinecop materialized = model;
  materialized.reorient(conditioning_set);
  EXPECT_EQ(materialized.get_var_types(), model.get_var_types());

  auto u = tools_stats::simulate_uniform(60, d, false, { 41 });
  // relabeling does not change the joint density
  EXPECT_TRUE(all_close(materialized.pdf(u), model.pdf(u), 1e-10, 1e-12));

  EXPECT_TRUE(all_close(model.rosenblatt(u, conditioning_set, 2, false),
                        materialized.rosenblatt(u, 2, false),
                        1e-12,
                        1e-12));
  auto w = tools_stats::simulate_uniform(60, d, false, { 43 });
  EXPECT_TRUE(all_close(model.inverse_rosenblatt(w, conditioning_set, 2),
                        materialized.inverse_rosenblatt(w, 2),
                        1e-12,
                        1e-12));
  EXPECT_EQ(model.get_var_types(),
            (std::vector<std::string>{ "a", "c", "a", "c", "a" }));
}

TEST(test_circular_vinecop,
     conditional_simulation_fixes_the_conditioning_values)
{
  auto model = larger_mixed_vine();
  const size_t d = model.get_dim();
  const size_t n = 200;
  std::vector<size_t> conditioning_set{ 1, 3 }; // both circular
  Vinecop materialized = model;
  materialized.reorient(conditioning_set);
  const auto order = materialized.get_order();

  // conditioning values per variable: variable 1 next to the cut
  std::map<size_t, double> value{ { 1, 0.98 }, { 3, 0.25 } };
  Eigen::MatrixXd u_cond(n, 2);
  for (size_t i = 0; i < 2; ++i) {
    u_cond.col(i).setConstant(value.at(conditioning_set[i]));
  }
  auto sim =
    model.simulate_conditional(u_cond, conditioning_set, false, 1, { 47 });
  EXPECT_EQ(sim.rows(), static_cast<Eigen::Index>(n));
  EXPECT_EQ(sim.cols(), static_cast<Eigen::Index>(d));
  EXPECT_TRUE(all_close(sim.col(0), u_cond.col(0), 0, 1e-14));
  EXPECT_TRUE(all_close(sim.col(2), u_cond.col(1), 0, 1e-14));
  EXPECT_TRUE((sim.array() >= 0).all() && (sim.array() <= 1).all());

  // the materialized model takes the conditioning values in the order of
  // its tail and gives the same draw
  Eigen::MatrixXd u_cond_tail(n, 2);
  for (size_t i = 0; i < 2; ++i) {
    u_cond_tail.col(i).setConstant(value.at(order[d - 2 + i]));
  }
  EXPECT_TRUE(
    all_close(materialized.simulate_conditional(u_cond_tail, false, 1, { 47 }),
              sim,
              1e-12,
              1e-12));
}

TEST(test_circular_vinecop, json_round_trip_preserves_the_model)
{
  auto vc = larger_mixed_vine();
  Vinecop reloaded(vc.to_json());
  EXPECT_EQ(reloaded.get_var_types(), vc.get_var_types());
  EXPECT_EQ(reloaded.get_all_families(), vc.get_all_families());
  EXPECT_EQ(edge_types(reloaded), edge_types(vc));
  auto u = tools_stats::simulate_uniform(40, vc.get_dim(), false, { 53 });
  EXPECT_TRUE(all_close(reloaded.pdf(u), vc.pdf(u), 1e-12, 1e-12));
  EXPECT_TRUE(
    all_close(reloaded.rosenblatt(u), vc.rosenblatt(u), 1e-12, 1e-12));
}

TEST(test_circular_vinecop, omitted_pair_copulas_are_independence)
{
  MixedVine mixed;
  Vinecop vc(mixed.structure, {}, mixed.var_types);
  auto u = tools_stats::simulate_uniform(20, 3, false, { 59 });
  EXPECT_TRUE(all_close(vc.pdf(u), Eigen::VectorXd::Ones(20), 1e-14, 1e-14));
  EXPECT_EQ(vc.get_pair_copula(1, 0).get_var_types(), ac);
  EXPECT_EQ(vc.get_pair_copula(0, 0).get_var_types(), aa);
}
}
