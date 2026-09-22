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
  Bicop c13_2{ BicopFamily::cubic_sections, 0, par({ 0.8, 0.8, -0.5 }), ac };

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
  Bicop c23{ BicopFamily::cubic_sections, 0, par({ 0.7, 0.7, 2.0 }), ac };
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
  Bicop c12{ BicopFamily::cubic_sections, 0, par({ 0.7, 0.7, 0.2 }), ac };
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
        pcs[t][e] = Bicop(BicopFamily::cubic_sections,
                          0,
                          par({ 0.5, 0.5, 0.8 + shift }),
                          types);
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
  pcs[0][0] = Bicop(BicopFamily::cubic_sections, 0, par({ 0.5, 0.5, 0.0 }), ac);
  EXPECT_THROW(Vinecop(mixed.structure, pcs, mixed.var_types),
               std::runtime_error);
}

// -------------------------------------------------------------------------
// structure and family selection

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

TEST(test_circular_vinecop, circular_criterion_is_cut_invariant)
{
  auto u = MixedVine().model().simulate(400, false, 1, { 61 });
  Eigen::MatrixXd u12 = u.leftCols(2), u23 = u.rightCols(2);
  double aa_crit = tools_stats::pairwise_circular(u12, aa);
  double ac_crit = tools_stats::pairwise_circular(u23, ac);
  EXPECT_GT(aa_crit, 0.3);
  EXPECT_GT(ac_crit, 0.1);
  EXPECT_LE(aa_crit, 1.0);
  EXPECT_LE(ac_crit, 1.0);

  // shifting the cut of a circular variable leaves the measure unchanged
  auto shift = [](Eigen::MatrixXd x, Eigen::Index j, double by) {
    x.col(j) = (x.col(j).array() + by).unaryExpr([](double v) {
      return v - std::floor(v);
    });
    return x;
  };
  EXPECT_NEAR(
    tools_stats::pairwise_circular(shift(u12, 0, 0.37), aa), aa_crit, 1e-12);
  EXPECT_NEAR(
    tools_stats::pairwise_circular(shift(u12, 1, 0.81), aa), aa_crit, 1e-12);
  EXPECT_NEAR(
    tools_stats::pairwise_circular(shift(u23, 0, 0.37), ac), ac_crit, 1e-12);
  // ... and so does reflecting the linear variable or swapping the columns
  Eigen::MatrixXd u23_reflected = u23;
  u23_reflected.col(1) = 1.0 - u23.col(1).array();
  EXPECT_NEAR(
    tools_stats::pairwise_circular(u23_reflected, ac), ac_crit, 1e-12);
  Eigen::MatrixXd u32 = u23;
  u32.col(0).swap(u32.col(1));
  EXPECT_NEAR(tools_stats::pairwise_circular(u32, ca), ac_crit, 1e-12);

  // a rotation is a perfectly dependent pair
  Eigen::MatrixXd rot(200, 2);
  rot.col(0) = tools_stats::simulate_uniform(200, 1, false, { 67 });
  rot = shift(rot, 1, 0.0);
  rot.col(1) = rot.col(0);
  rot = shift(rot, 1, 0.3);
  EXPECT_NEAR(tools_stats::pairwise_circular(rot, aa), 1.0, 1e-12);

  // weights are honored; independence gives a small value
  Eigen::VectorXd w = Eigen::VectorXd::Constant(u.rows(), 2.0);
  EXPECT_NEAR(tools_stats::pairwise_circular(u12, aa, w), aa_crit, 1e-12);
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
  std::vector<std::string> var_types{ "a", "a", "c" };
  Vinecop truth(
    DVineStructure(std::vector<size_t>{ 1, 2, 3 }),
    { { Bicop(BicopFamily::von_mises, 0, par({ 1.0, 0.7 }), aa),
        Bicop(BicopFamily::cubic_sections, 0, par({ 1.0, 1.0, 1.0 }), ac) },
      { Bicop(BicopFamily::cubic_sections, 0, par({ 0.5, 0.5, -0.5 }), ac) } },
    var_types);
  auto u = truth.simulate(1500, false, 1, { 73 });

  Vinecop selected(u, RVineStructure(), var_types);
  EXPECT_EQ(selected.get_var_types(), var_types);
  EXPECT_EQ(first_tree_pairs(selected),
            (std::set<std::pair<size_t, size_t>>{ { 1, 2 }, { 2, 3 } }));
  // every edge carries its geometry and a family that supports it; the
  // circular-circular edge gets a circula
  size_t n_aa = 0;
  for (size_t t = 0; t < 2; ++t) {
    for (size_t e = 0; e + t < 2; ++e) {
      const auto pc = selected.get_pair_copula(t, e);
      EXPECT_TRUE(family_accepts_var_types(pc.get_family(), pc.get_var_types()))
        << get_family_name(pc.get_family());
      if (pc.get_var_types() == aa) {
        ++n_aa;
        EXPECT_TRUE(
          tools_stl::is_member(pc.get_family(), bicop_families::two_rotations));
      }
    }
  }
  EXPECT_EQ(n_aa, 1u);
  // the selected model fits about as well as the truth
  EXPECT_GT(selected.loglik(u), truth.loglik(u) - 30);

  // the same data and a thread pool give the same model
  FitControlsVinecop controls;
  controls.set_num_threads(2);
  Vinecop threaded(u, RVineStructure(), var_types, controls);
  EXPECT_EQ(threaded.get_all_families(), selected.get_all_families());
  EXPECT_EQ(threaded.get_rvine_structure().get_order(),
            selected.get_rvine_structure().get_order());
  EXPECT_TRUE(all_close(threaded.pdf(u), selected.pdf(u), 1e-10, 1e-12));
}

TEST(test_circular_vinecop, half_turn_pair_survives_selection_and_thresholding)
{
  // a perfectly dependent half-turn pair (Kendall's tau zero) and a weakly
  // dependent linear variable
  std::vector<std::string> var_types{ "a", "a", "c" };
  Bicop half_turn(BicopFamily::wrapped_cauchy, 0, par({ 0.95, pi }), aa);
  Bicop weak(BicopFamily::cubic_sections, 0, par({ 0.3, 0.3, 0.5 }), ac);
  Vinecop truth(DVineStructure(std::vector<size_t>{ 1, 2, 3 }),
                { { half_turn, weak } },
                var_types);
  truth.truncate(1);
  auto u = truth.simulate(1000, false, 1, { 79 });
  EXPECT_LT(std::fabs(wdm::wdm(u.leftCols(2), "tau")(0, 1)), 0.15);

  Vinecop selected(u, RVineStructure(), var_types);
  EXPECT_TRUE(first_tree_pairs(selected).count({ 1, 2 }));
  EXPECT_NE(selected.get_family(0, 0), BicopFamily::indep);

  // with a threshold the half-turn edge stays, and only the weak one may go
  FitControlsVinecop controls;
  controls.set_threshold(0.5);
  Vinecop thresholded(u, RVineStructure(), var_types, controls);
  auto pairs = first_tree_pairs(thresholded);
  EXPECT_TRUE(pairs.count({ 1, 2 }));
  for (size_t e = 0; e < 2; ++e) {
    const auto pc = thresholded.get_pair_copula(0, e);
    if (pc.get_var_types() == aa) {
      EXPECT_NE(pc.get_family(), BicopFamily::indep);
    } else {
      EXPECT_EQ(pc.get_family(), BicopFamily::indep);
    }
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
  for (size_t e = 0; e < 2; ++e) {
    EXPECT_EQ(by_tau.get_family(0, e), BicopFamily::indep);
  }
}

TEST(test_circular_vinecop, family_restrictions_apply_per_edge)
{
  MixedVine mixed;
  auto u = mixed.model().simulate(300, false, 1, { 83 });

  // no eligible family for a circular pair
  FitControlsVinecop only_linear;
  only_linear.set_family_set({ BicopFamily::gaussian });
  EXPECT_THROW(Vinecop(u, RVineStructure(), mixed.var_types, only_linear),
               std::runtime_error);

  // independence rescues every edge
  FitControlsVinecop with_indep;
  with_indep.set_family_set({ BicopFamily::gaussian, BicopFamily::indep });
  Vinecop indep(u, RVineStructure(), mixed.var_types, with_indep);
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
  Vinecop restricted(u, RVineStructure(), mixed.var_types, mixed_set);
  for (size_t t = 0; t < 2; ++t) {
    for (size_t e = 0; e + t < 2; ++e) {
      const auto pc = restricted.get_pair_copula(t, e);
      EXPECT_TRUE(tools_stl::is_member(
        pc.get_family(),
        { BicopFamily::von_mises, BicopFamily::cubic_sections }));
      EXPECT_TRUE(
        family_accepts_var_types(pc.get_family(), pc.get_var_types()));
    }
  }
}

TEST(test_circular_vinecop, truncation_and_threshold_searches_run)
{
  auto truth = larger_mixed_vine();
  auto u = truth.simulate(600, false, 1, { 89 });
  FitControlsVinecop controls;
  controls.set_select_trunc_lvl(true);
  controls.set_select_threshold(true);
  Vinecop sparse(u, RVineStructure(), truth.get_var_types(), controls);
  EXPECT_LE(sparse.get_trunc_lvl(), 4u);
  EXPECT_TRUE(std::isfinite(sparse.loglik(u)));
  EXPECT_EQ(sparse.get_var_types(), truth.get_var_types());
  for (const auto& tree : sparse.get_all_pair_copulas()) {
    for (const auto& pc : tree) {
      EXPECT_TRUE(
        family_accepts_var_types(pc.get_family(), pc.get_var_types()));
    }
  }
  // the density is still periodic in the circular variables
  Eigen::MatrixXd lo = u.topRows(20), hi = lo;
  lo.col(0).setZero();
  hi.col(0).setOnes();
  EXPECT_TRUE(all_close(sparse.pdf(lo), sparse.pdf(hi), 1e-6, 1e-8));

  // a random spanning tree also respects geometry
  FitControlsVinecop random;
  random.set_tree_algorithm("random_weighted");
  random.set_seeds({ 97 });
  Vinecop rnd(u, RVineStructure(), truth.get_var_types(), random);
  for (const auto& tree : rnd.get_all_pair_copulas()) {
    for (const auto& pc : tree) {
      EXPECT_TRUE(
        family_accepts_var_types(pc.get_family(), pc.get_var_types()));
    }
  }
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
  // quadratic sections (a = b) in tree 2: a, a, phase
  Eigen::VectorXd p13 = fitted.get_pair_copula(1, 0).get_parameters();
  EXPECT_NEAR(p13(0), 0.8, 0.2);
  EXPECT_NEAR(p13(1), 0.8, 0.2);
  EXPECT_LT(on_circle(p13(2), -0.5), 0.3);

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
