// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Snippets for docs/discrete.dox. Compiled as part of the documentation build
// so the examples cannot drift from the API.

#include <iostream>
#include <vinecopulib.hpp>

using namespace vinecopulib;

//! Discretizes a column to a grid of `k` atoms, returning F(x) and F(x-).
std::pair<Eigen::VectorXd, Eigen::VectorXd>
discretize(const Eigen::VectorXd& u, size_t k)
{
  Eigen::VectorXd upper(u.size()), lower(u.size());
  for (Eigen::Index i = 0; i < u.size(); ++i) {
    const double atom = std::ceil(u(i) * static_cast<double>(k));
    upper(i) = atom / static_cast<double>(k);
    lower(i) = (atom - 1.0) / static_cast<double>(k);
  }
  return { upper, lower };
}

void
snippet_bicop_discrete()
{
  //! [bicop]
  // A discrete variable needs two columns: F(x) and the left limit F(x-). The
  // data matrix therefore carries one extra column per discrete variable, all
  // the left limits following the d ordinary columns.
  auto continuous = tools_stats::simulate_uniform(500, 2, false, { 1 });
  auto disc = discretize(continuous.col(0), 5);

  Eigen::MatrixXd data(500, 3);
  data.col(0) = disc.first; // F(x) of the discrete variable
  data.col(1) = continuous.col(1);
  data.col(2) = disc.second; // F(x-) of the discrete variable

  // Declare the types, then fit as usual.
  FitControlsBicop controls({ BicopFamily::gaussian, BicopFamily::clayton });
  Bicop model(data, controls, { "d", "c" });

  std::cout << "family: " << model.get_family_name()
            << ", var_types: " << model.get_var_types()[0]
            << model.get_var_types()[1] << std::endl;

  // Evaluation takes data in the same layout.
  auto pdf = model.pdf(data);
  std::cout << "loglik: " << model.loglik(data) << std::endl;
  //! [bicop]
}

void
snippet_vinecop_discrete()
{
  //! [vinecop]
  // Same convention for vines: n x (d + k) columns, where k is the number of
  // discrete variables.
  size_t d = 4;
  auto continuous = tools_stats::simulate_uniform(500, d, false, { 2 });
  auto disc0 = discretize(continuous.col(0), 4);
  auto disc2 = discretize(continuous.col(2), 6);

  Eigen::MatrixXd data(500, d + 2);
  data.col(0) = disc0.first;
  data.col(1) = continuous.col(1);
  data.col(2) = disc2.first;
  data.col(3) = continuous.col(3);
  data.col(4) = disc0.second; // left limits, in the order the discrete
  data.col(5) = disc2.second; // variables appear among the first d columns

  Vinecop model(data, RVineStructure(), { "d", "c", "d", "c" });

  std::cout << "loglik: " << model.get_loglik() << std::endl;

  // simulate() returns the d "F(x)" columns only
  auto sim = model.simulate(10, false, 1, { 3 });
  std::cout << "simulated: " << sim.rows() << " x " << sim.cols() << std::endl;

  // set_var_types can also change the types of an existing model; the pair
  // copulas are re-typed to match their position in the vine.
  Vinecop continuous_model(continuous);
  continuous_model.set_var_types({ "d", "c", "d", "c" });
  //! [vinecop]
}

int
main()
{
  // The snippets above are the documentation; this wrapper is not shown.
  // Catching here keeps `main` non-throwing, which clang-tidy requires.
  try {
    snippet_bicop_discrete();
    snippet_vinecop_discrete();
  } catch (const std::exception& e) {
    std::cerr << "snippet failed: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
