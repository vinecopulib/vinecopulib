// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Snippets for docs/conditional.dox. Compiled as part of the documentation
// build so the examples cannot drift from the API.

#include <iostream>
#include <vinecopulib.hpp>

using namespace vinecopulib;

void
snippet_select_conditioning()
{
  //! [select]
  size_t d = 5;
  auto data = tools_stats::simulate_uniform(500, d, false, { 1 });

  // Ask the selection to end up with variables 3 and 4 last in the order, so
  // that they can be conditioned on afterwards. Indices are 1-based, matching
  // the vine order.
  FitControlsVinecop controls;
  controls.set_conditioning_set({ 3, 4 });
  Vinecop model(data, RVineStructure(), {}, controls);

  const auto order = model.get_order();
  std::cout << "order: ";
  for (auto v : order)
    std::cout << v << " ";
  std::cout << "(the conditioning set is the tail)" << std::endl;
  //! [select]
}

void
snippet_simulate_conditional()
{
  //! [simulate]
  size_t d = 4;
  auto data = tools_stats::simulate_uniform(500, d, false, { 2 });

  FitControlsVinecop controls;
  controls.set_conditioning_set({ 3, 4 });
  Vinecop model(data, RVineStructure(), {}, controls);

  // One conditioning point per output row; the columns correspond, left to
  // right, to the last k entries of the vine order -- here variables 3 and 4.
  size_t n = 100;
  Eigen::MatrixXd u_cond(n, 2);
  u_cond.col(0).setConstant(0.3);
  u_cond.col(1).setConstant(0.8);

  // n x d: the conditioning columns reproduce u_cond, the rest are draws from
  // the conditional distribution.
  Eigen::MatrixXd sim = model.simulate_conditional(u_cond, false, 1, { 7 });
  std::cout << "conditional draws: " << sim.rows() << " x " << sim.cols()
            << std::endl;
  //! [simulate]
}

void
snippet_reorient()
{
  //! [reorient]
  auto data = tools_stats::simulate_uniform(500, 5, false, { 3 });
  Vinecop model(data);

  const double loglik_before = model.loglik(data);

  // Relabel an already-fitted vine so that a given pair comes last and can be
  // conditioned on. Only sets admissible as a sampling-order tail of this vine
  // can be reached; swapping the two variables already at the tail always is.
  // To guarantee a particular set, fit with
  // FitControlsVinecop::set_conditioning_set() instead.
  const auto order = model.get_order();
  const size_t d = order.size();
  model.reorient({ order[d - 1], order[d - 2] });

  // The model itself is unchanged: only its sampling-order representation
  // moves, so the density and the log-likelihood are invariant.
  std::cout << "loglik before: " << loglik_before
            << ", after: " << model.loglik(data) << std::endl;

  // An arbitrary set generally is not admissible.
  try {
    model.reorient({ 1, 2 });
  } catch (const std::runtime_error& e) {
    std::cout << "not admissible: " << e.what() << std::endl;
  }
  //! [reorient]
}

int
main()
{
  // The snippets above are the documentation; this wrapper is not shown.
  // Catching here keeps `main` non-throwing, which clang-tidy requires.
  try {
    snippet_select_conditioning();
    snippet_simulate_conditional();
    snippet_reorient();
  } catch (const std::exception& e) {
    std::cerr << "snippet failed: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
