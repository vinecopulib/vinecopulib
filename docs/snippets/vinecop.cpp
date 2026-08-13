// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Snippets for docs/overview-vinecop.dox. Compiled as part of the documentation
// build so the examples cannot drift from the API.

#include <fstream>
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
snippet_structure()
{
  //! [structure]
  // an R-vine structure from a matrix (columns are trees, the diagonal holds
  // the conditioned variable of each column)
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(4, 4);
  mat << 4, 4, 4, 4, 3, 3, 3, 0, 2, 2, 0, 0, 1, 0, 0, 0;
  RVineStructure structure(mat);

  // C- and D-vines are determined by an order alone
  CVineStructure cvine({ 1, 2, 3, 4 });
  DVineStructure dvine({ 1, 2, 3, 4 });

  // a random structure on 5 variables, truncated after 2 trees
  auto random = RVineStructure::simulate(5, false, { 7 });

  std::cout << "dimension: " << structure.get_dim()
            << ", truncation level: " << structure.get_trunc_lvl() << std::endl;
  std::cout << structure.str() << std::endl;
  //! [structure]
}

void
snippet_custom()
{
  //! [custom]
  // a 4-dimensional D-vine
  DVineStructure structure({ 1, 2, 3, 4 });

  // all pair copulas independence by default; fill in the ones you want
  auto pair_copulas = Vinecop::make_pair_copula_store(4);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 0, Eigen::VectorXd::Constant(1, 2.0));
    }
  }

  Vinecop model(structure, pair_copulas);
  //! [custom]
}

void
snippet_fit()
{
  //! [fit]
  size_t d = 5;
  auto data = dependent_data(500, d, 1);

  // select the structure along with the families
  Vinecop model(data);

  // keep a known structure and select only the families and parameters
  Vinecop model2(data, DVineStructure({ 1, 2, 3, 4, 5 }));

  // refit an existing model's families and parameters, structure unchanged
  model2.fit(data);

  std::cout << "loglik: " << model.get_loglik() << ", bic: " << model.get_bic()
            << std::endl;
  //! [fit]
}

void
snippet_controls()
{
  //! [controls]
  auto data = dependent_data(500, 5, 2);

  FitControlsVinecop controls;
  controls.set_family_set(bicop_families::itau);
  controls.set_parametric_method("itau");
  controls.set_trunc_lvl(2);           // force independence above tree 2
  controls.set_tree_criterion("mcor"); // edge weights for the spanning tree
  controls.set_threshold(0.05);        // weaker dependencies -> independence
  controls.set_tree_algorithm("mst_kruskal");
  controls.set_num_threads(2);

  Vinecop model(data, RVineStructure(), {}, controls);

  // or let the fit choose the truncation level and threshold by mBICv
  FitControlsVinecop selected;
  selected.set_select_trunc_lvl(true);
  selected.set_select_threshold(true);
  Vinecop sparse(data, RVineStructure(), {}, selected);
  std::cout << "selected truncation level: "
            << sparse.get_rvine_structure().get_trunc_lvl() << std::endl;
  //! [controls]
}

void
snippet_model()
{
  //! [model]
  auto data = dependent_data(200, 4, 3);
  Vinecop model(data);

  // density, and the Monte-Carlo distribution function
  auto pdf = model.pdf(data);
  auto cdf = model.cdf(data, 1000);

  // pdf_full also returns the pair-copula densities and h-functions computed
  // along the way, which the derivative routines reuse
  auto full = model.pdf_full(data);

  // simulate, and the Rosenblatt transform with its inverse
  auto sim = model.simulate(100, false, 1, { 5 });
  auto z = model.rosenblatt(data);
  auto back = model.inverse_rosenblatt(z);

  // fit statistics
  std::cout << "loglik: " << model.get_loglik() << ", aic: " << model.get_aic()
            << ", bic: " << model.get_bic()
            << ", mbicv: " << model.get_mbicv(0.9) << std::endl;

  // a truncated copy
  auto truncated = model;
  truncated.truncate(1);
  //! [model]
}

void
snippet_trees()
{
  //! [trees]
  auto data = dependent_data(300, 5, 4);
  Vinecop model(data);

  // the structure as an edge list, tree by tree: each edge carries the
  // conditioned pair (a, b), the conditioning set C, and its pair copula
  RVineTrees trees = model.get_trees();
  for (const auto& tree : trees.get_trees()) {
    for (const auto& edge : tree) {
      std::cout << edge.a << "," << edge.b;
      if (!edge.C.empty()) {
        std::cout << " | ";
        for (auto c : edge.C)
          std::cout << c << " ";
      }
      std::cout << "(" << edge.pair_copula.get_family_name() << ")  ";
    }
    std::cout << std::endl;
  }

  // and back again
  RVineStructure round_trip(trees);
  //! [trees]
}

void
snippet_serialization()
{
  //! [serialization]
  auto data = dependent_data(200, 3, 6);
  Vinecop model(data);

  // JSON, via a file or an nlohmann::json object
  model.to_file(std::string("myvine.json"));
  Vinecop model2(std::string("myvine.json"));

  nlohmann::json j = model.to_json();
  Vinecop model3(j);

  // the .cbor extension selects binary CBOR automatically
  model.to_file(std::string("myvine.cbor"));
  Vinecop model4(std::string("myvine.cbor"));
  //! [serialization]
}

int
main()
{
  // The snippets above are the documentation; this wrapper is not shown.
  // Catching here keeps `main` non-throwing, which clang-tidy requires.
  try {
    snippet_structure();
    snippet_custom();
    snippet_fit();
    snippet_controls();
    snippet_model();
    snippet_trees();
    snippet_serialization();
  } catch (const std::exception& e) {
    std::cerr << "snippet failed: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
