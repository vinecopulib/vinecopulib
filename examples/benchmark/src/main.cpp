#include <Eigen/Dense>
#include <chrono>
#include <iostream>
#include <random>
#include <vinecopulib.hpp>

using namespace std;
using namespace vinecopulib;
using Eigen::MatrixXd;

// Utility function to generate random data using Eigen
MatrixXd
generate_data(int n, int d)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> dist(0.0, 1.0);

  MatrixXd data(n, d);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < d; ++j) {
      data(i, j) = dist(gen) * 1.0 + 0.5 * dist(gen);
    }
  }
  return data;
}

// Benchmarking function
template<typename Function>
double
benchmark(Function func, int repeats = 10, int n = 1000, int d = 5)
{
  double total_time = 0.0;
  for (int i = 0; i < repeats; ++i) {
    // Generate and preprocess data
    MatrixXd x = generate_data(n, d);

    // Convert to pseudo-observations
    auto u = tools_stats::to_pseudo_obs(x);

    auto start = chrono::high_resolution_clock::now();
    func(u);
    auto end = chrono::high_resolution_clock::now();
    total_time += chrono::duration<double, std::milli>(end - start).count();
  }
  return total_time / repeats;
}

int
main()
{
  // Define different FitControls configurations
  FitControlsVinecop controls_itau(bicop_families::itau);
  FitControlsVinecop controls_itau_par_method(bicop_families::itau, "itau");
  FitControlsVinecop controls_tll({ BicopFamily::tll });



  // Benchmark different configurations
  double time_itau = benchmark([&](const Eigen::MatrixXd& u) {
    Vinecop vc(u, RVineStructure(), {}, controls_itau);
  });
  double time_itau_par_method = benchmark([&](const Eigen::MatrixXd& u) {
    Vinecop vc(u, RVineStructure(), {}, controls_itau_par_method);
  });
  double time_tll = benchmark([&](const Eigen::MatrixXd& u) {
    Vinecop vc(u, RVineStructure(), {}, controls_tll);
  });

  // Output results
  cout << "Benchmark Results (ms):" << endl;
  cout << "itau: " << time_itau << " ms" << endl;
  cout << "itau_par_method: " << time_itau_par_method << " ms" << endl;
  cout << "tll: " << time_tll << " ms" << endl;

  return 0;
}