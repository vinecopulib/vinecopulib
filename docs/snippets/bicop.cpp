// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Snippets for docs/overview-bicop.dox. Compiled as part of the documentation
// build so the examples cannot drift from the API.

#include <fstream>
#include <iostream>
#include <vinecopulib.hpp>

using namespace vinecopulib;

void
snippet_families()
{
  //! [families]
  // print all available families
  std::cout << "Available families: ";
  for (auto family : bicop_families::all) {
    std::cout << get_family_name(family) << " ";
  }
  std::cout << std::endl;
  //! [families]
}

void
snippet_custom()
{
  //! [custom]
  // 90 degree rotated Clayton with default parameter (independence)
  Bicop clayton(BicopFamily::clayton, 90);

  // Gauss copula with parameter 0.5
  Bicop gauss(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.5));
  //! [custom]
}

void
snippet_fit()
{
  //! [fit]
  // create a Gauss copula with parameter 0.5 and simulate 1000 observations
  Bicop model(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.5));
  auto data = model.simulate(1000);

  // instantiate a Gaussian copula with default parameters and fit to data
  Bicop fitted(BicopFamily::gaussian);
  fitted.fit(data);
  std::cout << "estimated parameter: " << fitted.get_parameters() << std::endl;

  // assign another family to the same variable and fit to data
  fitted = Bicop(BicopFamily::student);
  fitted.fit(data);
  std::cout << "estimated parameters: " << fitted.get_parameters() << std::endl;

  // alternatively, select the family automatically
  fitted.select(data);
  std::cout << "family: " << fitted.get_family_name()
            << ", rotation: " << fitted.get_rotation() << std::endl;
  //! [fit]
}

void
snippet_select_controls(const Eigen::MatrixXd& data)
{
  //! [select-controls]
  // select the "best" archimedean family by the default criterion (AIC), with
  // parameters from maximum likelihood
  Bicop best_archimedean(data, FitControlsBicop(bicop_families::archimedean));
  std::cout << "family: " << best_archimedean.get_family_name()
            << ", rotation: " << best_archimedean.get_rotation() << ", "
            << best_archimedean.get_parameters() << std::endl;

  // select by BIC, with parameters from Kendall's tau inversion
  FitControlsBicop controls(bicop_families::itau, "itau");
  controls.set_selection_criterion("bic");
  Bicop best_itau(data, controls);
  std::cout << "family: " << best_itau.get_family_name()
            << ", rotation: " << best_itau.get_rotation() << ", "
            << best_itau.get_parameters() << std::endl;
  //! [select-controls]
}

void
snippet_model()
{
  //! [model]
  // Gauss copula with parameter 0.5
  Bicop bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.5));

  // simulate 100 observations
  auto sim_data = bicop.simulate(100);

  // evaluate the pdf and cdf
  auto pdf = bicop.pdf(sim_data);
  auto cdf = bicop.cdf(sim_data);

  // evaluate the two h-functions and their inverses
  auto h1 = bicop.hfunc1(sim_data);
  auto h2 = bicop.hfunc2(sim_data);
  auto hi1 = bicop.hinv1(sim_data);
  auto hi2 = bicop.hinv2(sim_data);

  // evaluate the log-likelihood and the information criteria
  auto ll = bicop.loglik(sim_data);
  auto aic = bicop.aic(sim_data);
  auto bic = bicop.bic(sim_data);
  auto mbic = bicop.mbic(sim_data);
  std::cout << "loglik: " << ll << ", aic: " << aic << ", bic: " << bic
            << ", mbic: " << mbic << std::endl;
  //! [model]
}

void
snippet_dependence()
{
  //! [dependence]
  Bicop bicop(BicopFamily::bb1, 0, (Eigen::VectorXd(2) << 1.0, 2.0).finished());

  // Kendall's tau of the current parameters
  double tau = bicop.get_tau();

  // Blomqvist's beta, and the 2 x 2 matrix of tail-dependence coefficients
  // (lower in (0, 0), upper in (1, 1))
  double beta = bicop.get_beta();
  Eigen::MatrixXd taildep = bicop.get_taildep();

  // The parameter <-> tau maps are only defined for the families in
  // bicop_families::itau, which BB1 is not.
  Bicop clayton(BicopFamily::clayton);
  Eigen::MatrixXd par = clayton.tau_to_parameters(0.5);
  double tau_again = clayton.parameters_to_tau(par);
  std::cout << "BB1 tau: " << tau << ", beta: " << beta
            << ", lower/upper tail dependence: " << taildep(0, 0) << "/"
            << taildep(1, 1) << "; Clayton tau round trip: " << tau_again
            << std::endl;
  //! [dependence]
}

void
snippet_per_row()
{
  //! [per-row]
  Bicop bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.5));

  // One parameter set per observation: simulate 100 draws, each from its own
  // Gauss copula, then evaluate at the same per-observation parameters.
  Eigen::MatrixXd parameters = Eigen::VectorXd::LinSpaced(100, 0.1, 0.8);
  auto sim = bicop.simulate(parameters);
  auto pdf = bicop.pdf(sim, parameters);
  auto loglik = bicop.loglik(sim, parameters);

  std::cout << "per-row draws: " << sim.rows() << " x " << sim.cols()
            << ", loglik: " << loglik << std::endl;
  //! [per-row]
}

void
snippet_serialization()
{
  //! [serialization]
  // Gauss copula with parameter 0.5
  Bicop bicop(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.5));

  // save as a nlohmann::json object
  nlohmann::json bicop_json = bicop.to_json();

  // write into a JSON file
  std::string filename = "myfile.json";
  std::ofstream file(filename);
  file << bicop_json << std::endl;
  file.close();

  // equivalently
  bicop.to_file(filename);

  // a new Bicop can be constructed from the nlohmann::json object
  Bicop bicop2(bicop_json);

  // or from the JSON file
  Bicop bicop3(std::string("myfile.json"));

  // the .cbor extension selects binary CBOR automatically
  bicop.to_file(std::string("myfile.cbor"));
  Bicop bicop4(std::string("myfile.cbor"));
  //! [serialization]
}

int
main()
{
  // The snippets above are the documentation; this wrapper is not shown.
  // Catching here keeps `main` non-throwing, which clang-tidy requires.
  try {
    snippet_families();
    snippet_custom();
    snippet_fit();
    Bicop model(BicopFamily::gaussian, 0, Eigen::VectorXd::Constant(1, 0.5));
    snippet_select_controls(model.simulate(500));
    snippet_model();
    snippet_dependence();
    snippet_per_row();
    snippet_serialization();
  } catch (const std::exception& e) {
    std::cerr << "snippet failed: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
