// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Vine copulas with circular variables: the geometry of every edge,
// evaluation, transforms, simulation, refitting, structure operations, and
// structure and family selection.

#include "include/test_utils.hpp"
#include "gtest/gtest.h"
#include <boost/math/constants/constants.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>
#include <wdm/eigen.hpp>

namespace test_circular_vinecop {

using namespace vinecopulib;
using test_utils::all_close;
using Types = std::vector<std::string>;
using EdgeTypes = std::vector<Types>;

const Types aa = { "a", "a" };
const Types ac = { "a", "c" };
const Types ca = { "c", "a" };
const Types cc = { "c", "c" };
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

//! shifts one column around the circle
Eigen::MatrixXd
shift(Eigen::MatrixXd x, Eigen::Index j, double by)
{
  x.col(j) = (x.col(j).array() + by).unaryExpr([](double v) {
    return v - std::floor(v);
  });
  return x;
}

//! a D-vine on (1, 2, 3) with the given types and pair copulas
struct DVine3
{
  Types var_types;
  Bicop c12, c23, c13_2;

  DVineStructure structure() const
  {
    return DVineStructure(std::vector<size_t>{ 1, 2, 3 });
  }
  std::vector<std::vector<Bicop>> pair_copulas() const
  {
    return { { c12, c23 }, { c13_2 } };
  }
  Vinecop model() const
  {
    return Vinecop(structure(), pair_copulas(), var_types);
  }
  //! the joint density assembled by hand from the three pair copulas
  Eigen::VectorXd density(const Eigen::MatrixXd& u) const
  {
    Eigen::MatrixXd u13_2(u.rows(), 2);
    u13_2.col(0) = c12.hfunc2(u.leftCols(2));  // F(u1 | u2)
    u13_2.col(1) = c23.hfunc1(u.rightCols(2)); // F(u3 | u2)
    return c12.pdf(u.leftCols(2)).array() * c23.pdf(u.rightCols(2)).array() *
           c13_2.pdf(u13_2).array();
  }
};

//! variables 1 and 2 circular, 3 linear: a circula in tree 1, cylindrical
//! families on the mixed edges, including the tree-2 edge whose circular
//! argument arrives through an h-function
DVine3
mixed()
{
  return { { "a", "a", "c" },
           Bicop(BicopFamily::von_mises, 90, par({ 2.0, 0.7 }), aa),
           Bicop(BicopFamily::cubic_sections, 0, par({ 0.6, -0.4, 1.0 }), ac),
           Bicop(BicopFamily::cubic_sections, 0, par({ 0.8, 0.8, -0.5 }), ac) };
}
//! the only circular variable, 2, is conditioned on in tree 2
DVine3
conditioning_only()
{
  return { { "c", "a", "c" },
           Bicop(BicopFamily::cubic_sections, 0, par({ 0.5, 0.3, 0.4 }), ca),
           Bicop(BicopFamily::cubic_sections, 0, par({ 0.7, 0.7, 2.0 }), ac),
           Bicop(BicopFamily::gaussian, 0, par({ 0.5 }), cc) };
}
//! the circular variable 1 is a conditioned argument of the tree-2 edge
DVine3
circular_linear_linear()
{
  return { { "a", "c", "c" },
           Bicop(BicopFamily::cubic_sections, 0, par({ 0.7, 0.7, 0.2 }), ac),
           Bicop(BicopFamily::gaussian, 0, par({ -0.6 }), cc),
           Bicop(BicopFamily::cubic_sections, 0, par({ 0.3, 0.8, 2.5 }), ac) };
}

//! A five-dimensional vine on a simulated structure with three circular
//! variables, its pair copulas chosen by the geometry of each edge.
Vinecop
larger_mixed_vine()
{
  const size_t d = 5;
  Types var_types{ "a", "c", "a", "c", "a" };
  auto structure = RVineStructure::simulate(d, false, { 7 });
  // independence is eligible on every edge, so this assigns every pair its
  // geometry, which then picks the family
  Vinecop vc(structure, Vinecop::make_pair_copula_store(d), var_types);
  auto pcs = vc.get_all_pair_copulas();
  for (size_t t = 0; t < pcs.size(); ++t) {
    for (size_t e = 0; e < pcs[t].size(); ++e) {
      const auto types = pcs[t][e].get_var_types();
      const double phase = 0.3 * static_cast<double>(t + e);
      if (types == aa) {
        pcs[t][e] =
          Bicop(BicopFamily::wrapped_cauchy, 0, par({ 0.5, 0.3 + phase }), aa);
      } else if (types == cc) {
        pcs[t][e] = Bicop(BicopFamily::gaussian, 0, par({ 0.4 }), cc);
      } else {
        pcs[t][e] = Bicop(BicopFamily::cubic_sections,
                          0,
                          par({ 0.5, 0.5, 0.8 + phase }),
                          types);
      }
    }
  }
  vc.set_all_pair_copulas(pcs);
  return vc;
}

EdgeTypes
edge_types(const Vinecop& vc)
{
  EdgeTypes types;
  for (const auto& tree : vc.get_all_pair_copulas()) {
    for (const auto& pc : tree) {
      types.push_back(pc.get_var_types());
    }
  }
  return types;
}

//! every pair copula supports the geometry of its edge
void
expect_eligible_everywhere(const Vinecop& vc)
{
  for (const auto& tree : vc.get_all_pair_copulas()) {
    for (const auto& pc : tree) {
      EXPECT_TRUE(family_accepts_var_types(pc.get_family(), pc.get_var_types()))
        << get_family_name(pc.get_family());
    }
  }
}

//! the joint density agrees at both ends of every circular variable (the
//! arguments are trimmed to [1e-10, 1 - 1e-10] on the way in)
void
expect_periodic(const Vinecop& vc, const Eigen::MatrixXd& u)
{
  const auto types = vc.get_var_types();
  for (size_t j = 0; j < types.size(); ++j) {
    if (types[j] != "a") {
      continue;
    }
    Eigen::MatrixXd lo = u, hi = u;
    lo.col(j).setZero();
    hi.col(j).setOnes();
    EXPECT_TRUE(all_close(vc.pdf(lo), vc.pdf(hi), 1e-6, 1e-8))
      << "variable " << j + 1;
  }
}

//! the pair copulas of the first tree (a copy, so it can be iterated)
std::vector<Bicop>
first_tree(const Vinecop& vc)
{
  return vc.get_all_pair_copulas()[0];
}

//! the conditioned pairs of the first tree, each sorted
std::set<std::pair<size_t, size_t>>
first_tree_pairs(const Vinecop& vc)
{
  std::set<std::pair<size_t, size_t>> pairs;
  const auto& structure = vc.get_rvine_structure();
  const auto order = structure.get_order();
  for (size_t e = 0; e + 1 < order.size(); ++e) {
    size_t a = order[e];
    size_t b = order[structure.struct_array(0, e, true) - 1];
    pairs.emplace(std::min(a, b), std::max(a, b));
  }
  return pairs;
}

// -------------------------------------------------------------------------
// supplied structures

TEST(test_circular_vinecop, edges_carry_the_geometry_of_conditioned_variables)
{
  // a conditioned circular variable stays circular through the trees; a
  // variable that is only conditioned on leaves the pair linear
  const DVine3 vine = mixed();
  auto vc = vine.model();
  EXPECT_EQ(edge_types(vc), (EdgeTypes{ aa, ac, ac }));
  EXPECT_EQ(edge_types(conditioning_only().model()), (EdgeTypes{ ca, ac, cc }));
  EXPECT_EQ(edge_types(circular_linear_linear().model()),
            (EdgeTypes{ ac, cc, ac }));

  // setting the types afterwards, and omitted pair copulas, follow the rule
  Vinecop relabeled(vine.structure(), Vinecop::make_pair_copula_store(3));
  relabeled.set_var_types(vine.var_types);
  EXPECT_EQ(edge_types(relabeled), (EdgeTypes{ aa, ac, ac }));
  Vinecop omitted(vine.structure(), {}, vine.var_types);
  EXPECT_EQ(omitted.get_pair_copula(1, 0).get_var_types(), ac);
  auto u = tools_stats::simulate_uniform(20, 3, false, { 59 });
  EXPECT_TRUE(all_close(omitted.pdf(u), Eigen::VectorXd::Ones(20), 1e-14));

  // truncation keeps the geometry and drops the tree-2 factor
  vc.truncate(1);
  EXPECT_EQ(edge_types(vc), (EdgeTypes{ aa, ac }));
  Eigen::VectorXd expected =
    vine.c12.pdf(u.leftCols(2)).array() * vine.c23.pdf(u.rightCols(2)).array();
  EXPECT_TRUE(all_close(vc.pdf(u), expected, 1e-12, 1e-12));

  // a family that does not support its edge is rejected by tree and edge,
  // by the constructor and by the setter
  auto pcs = vine.pair_copulas();
  pcs[1][0] = Bicop(BicopFamily::gaussian, 0, par({ 0.5 }));
  try {
    Vinecop bad(vine.structure(), pcs, vine.var_types);
    FAIL() << "a Gaussian on a circular-linear edge must be rejected";
  } catch (const std::runtime_error& err) {
    EXPECT_NE(std::string(err.what()).find("tree 2, edge 1"), std::string::npos)
      << err.what();
  }
  auto model = vine.model();
  EXPECT_THROW(model.set_all_pair_copulas(pcs), std::runtime_error);
  pcs[1][0] = Bicop(BicopFamily::indep);
  model.set_all_pair_copulas(pcs);
  EXPECT_EQ(model.get_pair_copula(1, 0).get_var_types(), ac);
  pcs[0][0] = Bicop(BicopFamily::cubic_sections, 0, par({ 0.5, 0.5, 0.0 }), ac);
  EXPECT_THROW(Vinecop(vine.structure(), pcs, vine.var_types),
               std::runtime_error);
}

TEST(test_circular_vinecop, density_is_the_product_of_the_pair_densities)
{
  const DVine3 vine = mixed();
  auto vc = vine.model();
  auto u = tools_stats::simulate_uniform(200, 3, false, { 13 });
  EXPECT_TRUE(all_close(vc.pdf(u), vine.density(u), 1e-10, 1e-12));
  EXPECT_TRUE(all_close(vc.pdf(u, 2), vc.pdf(u), 1e-14, 1e-14));
  EXPECT_TRUE(all_close(vc.pdf_full(u).pdf, vc.pdf(u), 1e-14, 1e-14));
  EXPECT_NEAR(vc.loglik(u), vc.pdf(u).array().log().sum(), 1e-10);

  // integrating out the linear variable (Simpson's rule) leaves the circula
  const int m = 400;
  Eigen::MatrixXd grid(m + 1, 3);
  Eigen::MatrixXd u12(3, 2);
  u12 << 0.1, 0.2, 0.9, 0.35, 0.05, 0.95;
  Eigen::VectorXd c12 = vine.c12.pdf(u12);
  for (int r = 0; r < 3; ++r) {
    for (int i = 0; i <= m; ++i) {
      grid.row(i) << u12(r, 0), u12(r, 1), static_cast<double>(i) / m;
    }
    Eigen::VectorXd f = vc.pdf(grid);
    double integral = f(0) + f(m);
    for (int i = 1; i < m; ++i) {
      integral += (i % 2 == 1 ? 4.0 : 2.0) * f(i);
    }
    EXPECT_NEAR(integral / (3.0 * m), c12(r), 1e-5) << "row " << r;
  }

  // the Monte Carlo cdf is anchored and increasing
  Eigen::MatrixXd corners(3, 3);
  corners << 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9;
  auto p = vc.cdf(corners, 20000, 1, { 23 });
  EXPECT_NEAR(p(0), 1.0, 1e-12);
  EXPECT_GT(p(1), 0.0);
  EXPECT_GT(p(2), p(1) + 0.2);
  EXPECT_LT(p(2), 0.9);
}

TEST(test_circular_vinecop, evaluation_is_periodic_and_transforms_round_trip)
{
  for (const auto& vc : { mixed().model(),
                          conditioning_only().model(),
                          circular_linear_linear().model(),
                          larger_mixed_vine() }) {
    const size_t n = 1000;
    const size_t d = vc.get_dim();
    auto u = vc.simulate(n, false, 1, { 19 });
    expect_periodic(vc, u.topRows(50));
    EXPECT_TRUE(all_close(vc.simulate(n, false, 2, { 19 }), u, 1e-12, 1e-12));
    EXPECT_TRUE((u.array() > 0).all() && (u.array() < 1).all());

    // uniform margins: the sorted sample stays close to the identity
    for (size_t j = 0; j < d; ++j) {
      std::vector<double> col(u.col(j).data(), u.col(j).data() + n);
      std::sort(col.begin(), col.end());
      double max_dev = 0;
      for (size_t i = 0; i < n; ++i) {
        const double target =
          (static_cast<double>(i) + 0.5) / static_cast<double>(n);
        max_dev = std::max(max_dev, std::abs(col[i] - target));
      }
      EXPECT_LT(max_dev, 0.06) << "variable " << j + 1;
    }

    // Rosenblatt and its inverse, serial and threaded
    auto w = vc.rosenblatt(u);
    EXPECT_TRUE(all_close(vc.inverse_rosenblatt(w), u, 1e-6, 1e-8));
    EXPECT_TRUE(all_close(vc.rosenblatt(u, 2), w, 1e-12, 1e-12));
    EXPECT_TRUE(all_close(
      vc.inverse_rosenblatt(w, 2), vc.inverse_rosenblatt(w), 1e-12, 1e-12));
    for (size_t j = 0; j < d; ++j) {
      EXPECT_NEAR(w.col(j).mean(), 0.5, 0.04) << "variable " << j + 1;
    }
  }
}

TEST(test_circular_vinecop, scores_match_a_finite_difference_of_the_loglik)
{
  const DVine3 vine = mixed();
  auto vc = vine.model();
  auto u = vc.simulate(200, false, 1, { 29 });
  // the full (not step-wise) derivatives of the joint log-likelihood
  Eigen::MatrixXd scores = vc.scores(u, false);
  Eigen::MatrixXd hessian = vc.hessian(u, false);
  const auto pcs = vine.pair_copulas();

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
          return Vinecop(vine.structure(), store, vine.var_types).loglik(u);
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

TEST(test_circular_vinecop, fit_on_a_fixed_structure_recovers_the_parameters)
{
  const DVine3 vine = mixed();
  auto truth = vine.model();
  auto u = truth.simulate(3000, false, 1, { 31 });

  // start from the true families with default parameters
  auto pcs = vine.pair_copulas();
  for (auto& tree : pcs) {
    for (auto& pc : tree) {
      pc = Bicop(pc.get_family(),
                 pc.get_rotation(),
                 Eigen::MatrixXd(),
                 pc.get_var_types());
    }
  }
  Vinecop fitted(vine.structure(), pcs, vine.var_types);
  fitted.fit(u, FitControlsBicop(), 2);
  EXPECT_EQ(edge_types(fitted), edge_types(truth));

  // amplitudes and concentrations within 0.3, phases on the circle
  auto on_circle = [](double a, double b) {
    double diff = std::fmod(std::abs(a - b), 2 * pi);
    return std::min(diff, 2 * pi - diff);
  };
  for (size_t t = 0; t < 2; ++t) {
    for (size_t e = 0; e + t < 2; ++e) {
      Eigen::VectorXd est = fitted.get_pair_copula(t, e).get_parameters();
      Eigen::VectorXd truth_par = truth.get_pair_copula(t, e).get_parameters();
      for (Eigen::Index k = 0; k + 1 < est.size(); ++k) {
        EXPECT_NEAR(est(k), truth_par(k), 0.3)
          << "tree " << t << ", edge " << e;
      }
      EXPECT_LT(on_circle(est(est.size() - 1), truth_par(truth_par.size() - 1)),
                0.3)
        << "tree " << t << ", edge " << e;
    }
  }
  EXPECT_EQ(fitted.get_pair_copula(0, 0).get_rotation(), 90);
  EXPECT_NEAR(fitted.get_loglik(), fitted.loglik(u), 1e-8);
  EXPECT_GT(fitted.get_loglik(), truth.loglik(u) - 30);
}

TEST(test_circular_vinecop, views_conditional_simulation_and_json)
{
  auto model = larger_mixed_vine();
  const size_t d = model.get_dim();
  std::vector<size_t> conditioning_set{ 1, 3 }; // both circular
  Vinecop materialized = model;
  materialized.reorient(conditioning_set);
  EXPECT_EQ(materialized.get_var_types(), model.get_var_types());
  const auto order = materialized.get_order();

  // relabeling keeps the density; the views match the materialized model
  auto u = tools_stats::simulate_uniform(60, d, false, { 41 });
  EXPECT_TRUE(all_close(materialized.pdf(u), model.pdf(u), 1e-10, 1e-12));
  EXPECT_TRUE(all_close(model.rosenblatt(u, conditioning_set, 2, false),
                        materialized.rosenblatt(u, 2, false),
                        1e-12,
                        1e-12));
  EXPECT_TRUE(all_close(model.inverse_rosenblatt(u, conditioning_set, 2),
                        materialized.inverse_rosenblatt(u, 2),
                        1e-12,
                        1e-12));

  // conditional simulation fixes the conditioning values (variable 1 next to
  // the cut); the materialized model takes them in the order of its tail
  const size_t n = 200;
  std::map<size_t, double> value{ { 1, 0.98 }, { 3, 0.25 } };
  Eigen::MatrixXd u_cond(n, 2), u_cond_tail(n, 2);
  for (size_t i = 0; i < 2; ++i) {
    u_cond.col(i).setConstant(value.at(conditioning_set[i]));
    u_cond_tail.col(i).setConstant(value.at(order[d - 2 + i]));
  }
  auto sim =
    model.simulate_conditional(u_cond, conditioning_set, false, 1, { 47 });
  EXPECT_EQ(sim.cols(), static_cast<Eigen::Index>(d));
  EXPECT_TRUE(all_close(sim.col(0), u_cond.col(0), 0, 1e-14));
  EXPECT_TRUE(all_close(sim.col(2), u_cond.col(1), 0, 1e-14));
  EXPECT_TRUE((sim.array() >= 0).all() && (sim.array() <= 1).all());
  EXPECT_TRUE(
    all_close(materialized.simulate_conditional(u_cond_tail, false, 1, { 47 }),
              sim,
              1e-12,
              1e-12));

  // JSON keeps types, families, and the evaluation
  Vinecop reloaded(model.to_json());
  EXPECT_EQ(reloaded.get_var_types(), model.get_var_types());
  EXPECT_EQ(reloaded.get_all_families(), model.get_all_families());
  EXPECT_EQ(edge_types(reloaded), edge_types(model));
  EXPECT_TRUE(all_close(reloaded.pdf(u), model.pdf(u), 1e-12, 1e-12));
  EXPECT_TRUE(
    all_close(reloaded.rosenblatt(u), model.rosenblatt(u), 1e-12, 1e-12));
}

// -------------------------------------------------------------------------
// structure and family selection

TEST(test_circular_vinecop, circular_criterion_is_cut_invariant)
{
  auto u = mixed().model().simulate(400, false, 1, { 61 });
  Eigen::MatrixXd u12 = u.leftCols(2), u23 = u.rightCols(2);
  double aa_crit = tools_stats::pairwise_circular(u12, aa);
  double ac_crit = tools_stats::pairwise_circular(u23, ac);
  EXPECT_GT(aa_crit, 0.3);
  EXPECT_GT(ac_crit, 0.1);
  EXPECT_LE(std::max(aa_crit, ac_crit), 1.0);

  // invariant to the cut of either circular variable, to reflecting the
  // linear one, to the column order, and to constant weights
  EXPECT_NEAR(
    tools_stats::pairwise_circular(shift(u12, 0, 0.37), aa), aa_crit, 1e-12);
  EXPECT_NEAR(
    tools_stats::pairwise_circular(shift(u12, 1, 0.81), aa), aa_crit, 1e-12);
  EXPECT_NEAR(
    tools_stats::pairwise_circular(shift(u23, 0, 0.37), ac), ac_crit, 1e-12);
  Eigen::MatrixXd u23_reflected = u23;
  u23_reflected.col(1) = 1.0 - u23.col(1).array();
  EXPECT_NEAR(
    tools_stats::pairwise_circular(u23_reflected, ac), ac_crit, 1e-12);
  Eigen::MatrixXd u32 = u23;
  u32.col(0).swap(u32.col(1));
  EXPECT_NEAR(tools_stats::pairwise_circular(u32, ca), ac_crit, 1e-12);
  Eigen::VectorXd w = Eigen::VectorXd::Constant(u.rows(), 2.0);
  EXPECT_NEAR(tools_stats::pairwise_circular(u12, aa, w), aa_crit, 1e-12);
  EXPECT_EQ(tools_stats::pairwise_circular(u12, aa, 0.0 * w), 0.0);

  // one for a rotation, small under independence, rejected for linear pairs
  Eigen::MatrixXd rot(200, 2);
  rot.col(0) = tools_stats::simulate_uniform(200, 1, false, { 67 });
  rot.col(1) = rot.col(0);
  EXPECT_NEAR(
    tools_stats::pairwise_circular(shift(rot, 1, 0.3), aa), 1.0, 1e-12);
  auto indep = tools_stats::simulate_uniform(2000, 2, false, { 71 });
  EXPECT_LT(tools_stats::pairwise_circular(indep, aa), 0.1);
  EXPECT_LT(tools_stats::pairwise_circular(indep, ac), 0.1);
  EXPECT_THROW(tools_stats::pairwise_circular(indep, cc), std::runtime_error);

  // the selector routes pairs with a circular variable to it under every
  // built-in criterion, and leaves linear pairs alone
  for (const char* crit : { "tau", "rho", "hoeffd", "cxi", "mcor" }) {
    EXPECT_NEAR(
      tools_select::calculate_criterion(u12, crit, Eigen::VectorXd(), {}, aa),
      aa_crit,
      1e-12);
  }
  EXPECT_NEAR(
    tools_select::calculate_criterion(u23, "tau", Eigen::VectorXd(), {}, cc),
    std::fabs(wdm::wdm(u23, "tau")(0, 1)),
    1e-12);
}

TEST(test_circular_vinecop, select_recovers_the_mixed_vine)
{
  // strong tree-1 edges, a weak tree-2 edge, so that the dependence induced
  // between variables 1 and 3 is clearly weaker than either tree-1 edge
  const DVine3 vine{
    { "a", "a", "c" },
    Bicop(BicopFamily::von_mises, 0, par({ 1.0, 0.7 }), aa),
    Bicop(BicopFamily::cubic_sections, 0, par({ 1.0, 1.0, 1.0 }), ac),
    Bicop(BicopFamily::cubic_sections, 0, par({ 0.5, 0.5, -0.5 }), ac)
  };
  auto truth = vine.model();
  auto u = truth.simulate(1500, false, 1, { 73 });

  Vinecop selected(u, RVineStructure(), vine.var_types);
  EXPECT_EQ(selected.get_var_types(), vine.var_types);
  EXPECT_EQ(first_tree_pairs(selected),
            (std::set<std::pair<size_t, size_t>>{ { 1, 2 }, { 2, 3 } }));
  expect_eligible_everywhere(selected);
  size_t n_aa = 0;
  for (const auto& pc : first_tree(selected)) {
    if (pc.get_var_types() == aa) {
      ++n_aa;
      EXPECT_TRUE(
        tools_stl::is_member(pc.get_family(), bicop_families::two_rotations));
    }
  }
  EXPECT_EQ(n_aa, 1u);
  EXPECT_GT(selected.loglik(u), truth.loglik(u) - 30);

  // the same data and a thread pool give the same model
  FitControlsVinecop controls;
  controls.set_num_threads(2);
  Vinecop threaded(u, RVineStructure(), vine.var_types, controls);
  EXPECT_EQ(threaded.get_all_families(), selected.get_all_families());
  EXPECT_EQ(threaded.get_order(), selected.get_order());
  EXPECT_TRUE(all_close(threaded.pdf(u), selected.pdf(u), 1e-10, 1e-12));
}

TEST(test_circular_vinecop, half_turn_pair_survives_selection_and_thresholding)
{
  // a perfectly dependent half-turn pair (Kendall's tau zero) and a weakly
  // dependent linear variable
  Types var_types{ "a", "a", "c" };
  Vinecop truth(
    DVineStructure(std::vector<size_t>{ 1, 2, 3 }),
    { { Bicop(BicopFamily::wrapped_cauchy, 0, par({ 0.95, pi }), aa),
        Bicop(BicopFamily::cubic_sections, 0, par({ 0.3, 0.3, 0.5 }), ac) } },
    var_types);
  auto u = truth.simulate(1000, false, 1, { 79 });
  EXPECT_LT(std::fabs(wdm::wdm(u.leftCols(2), "tau")(0, 1)), 0.15);

  Vinecop selected(u, RVineStructure(), var_types);
  EXPECT_TRUE(first_tree_pairs(selected).count({ 1, 2 }));

  // with a threshold the half-turn edge stays and the weak one goes
  FitControlsVinecop controls;
  controls.set_threshold(0.5);
  Vinecop thresholded(u, RVineStructure(), var_types, controls);
  EXPECT_TRUE(first_tree_pairs(thresholded).count({ 1, 2 }));
  for (const auto& pc : first_tree(thresholded)) {
    EXPECT_EQ(pc.get_family() == BicopFamily::indep, pc.get_var_types() != aa);
  }

  // a custom criterion is honored as is: Kendall's tau misses the pair
  FitControlsVinecop tau_controls;
  tau_controls.set_tree_criterion("custom");
  tau_controls.set_tree_criterion_function(
    [](const Eigen::MatrixXd& x, const Eigen::VectorXd& w) {
      return std::fabs(wdm::wdm(x, "tau", w)(0, 1));
    });
  tau_controls.set_threshold(0.5);
  Vinecop by_tau(u, RVineStructure(), var_types, tau_controls);
  for (const auto& pc : first_tree(by_tau)) {
    EXPECT_EQ(pc.get_family(), BicopFamily::indep);
  }
}

TEST(test_circular_vinecop, selection_controls_apply_per_edge)
{
  const DVine3 vine = mixed();
  auto u = vine.model().simulate(300, false, 1, { 83 });

  // no eligible family for a circular pair; independence rescues every edge
  FitControlsVinecop only_linear;
  only_linear.set_family_set({ BicopFamily::gaussian });
  EXPECT_THROW(Vinecop(u, RVineStructure(), vine.var_types, only_linear),
               std::runtime_error);
  only_linear.set_family_set({ BicopFamily::gaussian, BicopFamily::indep });
  Vinecop indep(u, RVineStructure(), vine.var_types, only_linear);
  for (const auto& tree : indep.get_all_families()) {
    for (auto fam : tree) {
      EXPECT_EQ(fam, BicopFamily::indep);
    }
  }

  // an explicit mixed set restricts each edge to its eligible members
  FitControlsVinecop mixed_set;
  mixed_set.set_family_set({ BicopFamily::von_mises,
                             BicopFamily::cubic_sections,
                             BicopFamily::gaussian });
  Vinecop restricted(u, RVineStructure(), vine.var_types, mixed_set);
  expect_eligible_everywhere(restricted);
  for (const auto& tree : restricted.get_all_families()) {
    for (auto fam : tree) {
      EXPECT_NE(fam, BicopFamily::gaussian);
    }
  }

  // the threshold and truncation searches and a random spanning tree run on
  // a larger mixed vine and keep every edge eligible and periodic
  auto truth = larger_mixed_vine();
  auto v = truth.simulate(600, false, 1, { 89 });
  FitControlsVinecop sparse_controls;
  sparse_controls.set_select_trunc_lvl(true);
  sparse_controls.set_select_threshold(true);
  Vinecop sparse(v, RVineStructure(), truth.get_var_types(), sparse_controls);
  EXPECT_LE(sparse.get_trunc_lvl(), 4u);
  EXPECT_TRUE(std::isfinite(sparse.loglik(v)));
  expect_eligible_everywhere(sparse);
  expect_periodic(sparse, v.topRows(20));
  FitControlsVinecop random;
  random.set_tree_algorithm("random_weighted");
  random.set_seeds({ 97 });
  expect_eligible_everywhere(
    Vinecop(v, RVineStructure(), truth.get_var_types(), random));
}
}
