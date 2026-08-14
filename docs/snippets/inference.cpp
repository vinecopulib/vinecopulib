// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Snippets for docs/inference.dox. Compiled as part of the documentation build
// so the examples cannot drift from the API.

#include <iostream>
#include <vinecopulib.hpp>

using namespace vinecopulib;

//! Data from a known 5-dimensional D-vine, so that the fits below have
//! something to find.
Eigen::MatrixXd
dependent_data(size_t n, size_t d, int seed)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(d);
  const std::vector<BicopFamily> families = { BicopFamily::clayton,
                                              BicopFamily::gaussian,
                                              BicopFamily::gumbel,
                                              BicopFamily::frank };
  size_t k = 0;
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      const auto family = families[k++ % families.size()];
      const double par = (family == BicopFamily::gaussian) ? 0.6 : 2.0;
      pc = Bicop(family, 0, Eigen::VectorXd::Constant(1, par));
    }
  }
  Vinecop truth(DVineStructure(tools_stl::seq_int(1, d)), pair_copulas);
  return truth.simulate(n, false, 1, { seed });
}

void
snippet_bicop_scores()
{
  //! [bicop]
  Bicop model(
    BicopFamily::student, 0, (Eigen::VectorXd(2) << 0.5, 4.0).finished());
  auto u = model.simulate(500, false, { 1 });

  // n x npars matrix of per-observation scores: d log c / d theta
  Eigen::MatrixXd s = model.scores(u);

  // the summed score, i.e. the gradient of the log-likelihood
  Eigen::VectorXd g = model.gradient(u);

  // npars x npars Hessian of the log-likelihood, and the per-observation
  // outer-product covariance of the scores
  Eigen::MatrixXd h = model.hessian(u);
  Eigen::MatrixXd cov = model.scores_cov(u);

  std::cout << "npars: " << s.cols() << ", gradient norm: " << g.norm()
            << std::endl;

  // At the MLE the gradient is (numerically) zero, which is the cheapest
  // available check that a fit converged.
  Bicop fitted(BicopFamily::student);
  fitted.fit(u);
  std::cout << "gradient norm at the MLE: " << fitted.gradient(u).norm()
            << std::endl;

  // A sandwich standard error: (-H)^-1 cov (-H)^-1
  Eigen::MatrixXd hinv = (-h).inverse();
  Eigen::MatrixXd sandwich = hinv * cov * hinv;
  std::cout << "std. errors: " << sandwich.diagonal().cwiseSqrt().transpose()
            << std::endl;
  //! [bicop]
}

void
snippet_vinecop_scores()
{
  //! [vinecop]
  auto data = dependent_data(500, 4, 2);
  FitControlsVinecop controls(bicop_families::parametric);
  Vinecop model(data, RVineStructure(), {}, controls);

  // Columns are ordered (tree, edge, parameter): tree 0's edges first, and
  // within an edge the pair copula's own parameters in order.
  Eigen::MatrixXd s = model.scores(data);
  Eigen::VectorXd g = model.gradient(data);
  Eigen::MatrixXd h = model.hessian(data);

  std::cout << "total parameters: " << s.cols()
            << ", gradient norm: " << g.norm() << std::endl;

  // step_wise = true (the default) is the score of the step-wise MLE: each
  // pair copula's gradient treats its pseudo-observations as fixed. Pass false
  // for the gradient of the full log-likelihood, which propagates through the
  // h-function cascade.
  Eigen::VectorXd g_full = model.gradient(data, false);
  std::cout << "step-wise norm: " << g.norm()
            << ", full norm: " << g_full.norm() << std::endl;
  //! [vinecop]
}

void
snippet_reuse()
{
  //! [reuse]
  auto data = dependent_data(300, 4, 3);
  FitControlsVinecop controls(bicop_families::parametric);
  Vinecop model(data, RVineStructure(), {}, controls);

  // scores_full() returns the intermediate per-edge quantities alongside the
  // scores, so a caller needing several of them pays the cascade once.
  auto full = model.scores_full(data);
  std::cout << "scores: " << full.scores.rows() << " x " << full.scores.cols()
            << std::endl;

  // hessian_full() keeps the per-edge blocks rather than assembling them
  auto blocks = model.hessian_full(data);
  //! [reuse]
}

void
snippet_per_row()
{
  //! [per-row]
  auto data = dependent_data(200, 3, 4);
  FitControlsVinecop controls(bicop_families::parametric);
  Vinecop model(data, RVineStructure(), {}, controls);

  // Flatten the stored parameters in the same (tree, edge, parameter) order
  // that scores() reports and the `parameters` argument expects.
  std::vector<double> flat;
  for (const auto& tree : model.get_all_parameters()) {
    for (const auto& edge : tree) {
      for (Eigen::Index p = 0; p < edge.size(); ++p) {
        flat.push_back(edge(p));
      }
    }
  }
  Eigen::VectorXd theta = Eigen::Map<Eigen::VectorXd>(
    flat.data(), static_cast<Eigen::Index>(flat.size()));

  // One parameter vector per observation. Broadcasting the stored parameters
  // reproduces the fixed-parameter results; varying the rows is what makes
  // covariate-dependent parameters possible. Continuous, all-parametric models
  // only.
  Eigen::MatrixXd parameters = theta.transpose().replicate(data.rows(), 1);

  auto pdf = model.pdf(data, parameters);
  auto s = model.scores(data, parameters);
  std::cout << "parameters per row: " << parameters.cols()
            << ", pdf matches the fixed-parameter path: "
            << pdf.isApprox(model.pdf(data)) << std::endl;
  //! [per-row]
}

int
main()
{
  // The snippets above are the documentation; this wrapper is not shown.
  // Catching here keeps `main` non-throwing, which clang-tidy requires.
  try {
    snippet_bicop_scores();
    snippet_vinecop_scores();
    snippet_reuse();
    snippet_per_row();
  } catch (const std::exception& e) {
    std::cerr << "snippet failed: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
