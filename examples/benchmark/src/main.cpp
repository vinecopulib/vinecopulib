#include "benchmark.hpp"
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vinecopulib.hpp>

#include "benchmark.hpp"

using namespace std;
using namespace vinecopulib;
using Eigen::MatrixXd;

// Utility function to generate random data using Eigen
Eigen::MatrixXd
generate_data(int n, int d, unsigned seed)
{

  std::mt19937 gen(seed);
  std::normal_distribution<> dist(0.0, 1.0);

  MatrixXd data(n, d);
  for (int i = 0; i < n; ++i) {
    double x = dist(gen);
    for (int j = 0; j < d; ++j) {
      data(i, j) = x * 1.0 + 0.5 * dist(gen);
    }
  }
  return tools_stats::to_pseudo_obs(data);
}

void
benchmark_vinecop_fitting(int n = 1000, int d = 5, unsigned int repeats = 10)
{

  // Seeds for benchmarking
  Eigen::VectorXi seeds = Eigen::VectorXi::LinSpaced(repeats, 1, repeats);

  // Benchmark data generation
  Eigen::VectorXd times_generate_data = benchmark_func(
    [&](unsigned seed) { auto u = generate_data(n, d, seed); }, seeds);

  // Define different FitControls configurations
  std::map<std::string, FitControlsVinecop> controls_configs = {
    { "itau", FitControlsVinecop(bicop_families::itau) },
    { "itau_par_method", FitControlsVinecop(bicop_families::itau, "itau") },
    { "tll", FitControlsVinecop({ BicopFamily::tll }) }
  };

  // Wrapper function to generate data and fit a Vinecop
  auto generate_and_fit = [&](unsigned seed,
                              const FitControlsVinecop& controls) {
    auto u = generate_data(n, d, seed);
    Vinecop vc(u, RVineStructure(), {}, controls);
  };

  // Benchmark different configurations
  cout << "Benchmark Results for Vinecop Fitting (ms):" << endl;
  for (const auto& config : controls_configs) {
    Eigen::VectorXd time_fit = benchmark_func(
      [&](unsigned seed) { generate_and_fit(seed, config.second); }, seeds);
    cout << config.first << ": "
         << benchmark_stats(time_fit - times_generate_data).transpose() << endl;
  }
}

void
benchmark_bicop_tll(int n = 1000, unsigned int repeats = 10)
{
  FitControlsBicop controls_tll({ BicopFamily::tll }); // Define FitControls

  // Seeds for benchmarking
  Eigen::VectorXi seeds = Eigen::VectorXi::LinSpaced(repeats, 1, repeats);

  // Warmup
  auto warmup = benchmark_func(
    [&](unsigned seed) {
      auto u = generate_data(n, 2, seed);
      Bicop bc(u, controls_tll);
    },
    seeds);


  // Benchmark data generation
  Eigen::VectorXd times_generate_data = benchmark_func(
    [&](unsigned seed) { auto u = generate_data(n, 2, seed); }, seeds);

  // Benchmark fitting
  Eigen::VectorXd time_fit_tll = benchmark_func(
    [&](unsigned seed) {
      auto u = generate_data(n, 2, seed);
      Bicop bc(u, controls_tll);
    },
    seeds);

  cout << "Benchmark Results for Bicop tll (ms):" << endl;
  cout << "fit: "
       << benchmark_stats(time_fit_tll - times_generate_data).transpose()
       << endl;

  // Benchmark different methods
  std::vector<std::string> methods = { "pdf",    "cdf",   "hfunc1",
                                       "hfunc2", "hinv1", "hinv2" };
  for (const auto& method : methods) {
    Eigen::VectorXd time_method = benchmark_func(
      [&](unsigned seed) {
        auto u = generate_data(n, 2, seed);
        Bicop bc(u, controls_tll);
        if (method == "pdf") {
          auto d = bc.pdf(u);
        } else if (method == "cdf") {
          auto p = bc.cdf(u);
        } else if (method == "hfunc1") {
          auto h1 = bc.hfunc1(u);
        } else if (method == "hfunc2") {
          auto h2 = bc.hfunc2(u);
        } else if (method == "hinv1") {
          auto u1 = bc.hinv1(u);
        } else if (method == "hinv2") {
          auto u2 = bc.hinv2(u);
        }
      },
      seeds);
    cout << method << ": "
         << benchmark_stats(time_method - time_fit_tll).transpose() << endl;
  }
}

int
main()
{

  benchmark_vinecop_fitting();
  benchmark_bicop_tll();

  return 0;
}
