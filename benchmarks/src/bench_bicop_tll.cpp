// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "helpers.hpp"
#include <benchmark/benchmark.h>
#include <memory>

using namespace vinecopulib;

namespace {

using DataPtr = std::shared_ptr<const Eigen::MatrixXd>;

void
register_tll_fit(const std::string& method, size_t n)
{
  auto data = std::make_shared<const Eigen::MatrixXd>(
    bench::sim_data(BicopFamily::gaussian, 0, n));
  FitControlsBicop controls({ BicopFamily::tll });
  controls.set_nonparametric_method(method);
  benchmark::RegisterBenchmark(
    ("tll/fit/" + method + "/n=" + std::to_string(n)).c_str(),
    [data, controls](benchmark::State& st) {
      for (auto _ : st) {
        Bicop bc(BicopFamily::tll);
        bc.fit(*data, controls);
        benchmark::DoNotOptimize(bc);
      }
    });
}

//! Grid size is the one axis on which the fit and the margin normalization
//! scale together, so the ratio between them is only measurable if both are
//! benchmarked at the same `m`.
void
register_tll_fit_grid_size(const std::string& method, size_t m, size_t n)
{
  auto data = std::make_shared<const Eigen::MatrixXd>(
    bench::sim_data(BicopFamily::gaussian, 0, n));
  FitControlsBicop controls({ BicopFamily::tll });
  controls.set_nonparametric_method(method);
  controls.set_nonparametric_grid_size(m);
  benchmark::RegisterBenchmark(("tll/fit/" + method + "/m=" +
                                std::to_string(m) + "/n=" + std::to_string(n))
                                 .c_str(),
                               [data, controls](benchmark::State& st) {
                                 for (auto _ : st) {
                                   Bicop bc(BicopFamily::tll);
                                   bc.fit(*data, controls);
                                   benchmark::DoNotOptimize(bc);
                                 }
                               });
}

void
register_tll_eval(const std::string& method)
{
  const size_t n = 1000;
  Eigen::MatrixXd fit_data = bench::sim_data(BicopFamily::gaussian, 0, n);
  FitControlsBicop controls({ BicopFamily::tll });
  controls.set_nonparametric_method(method);
  Bicop bc(BicopFamily::tll);
  bc.fit(fit_data, controls);
  auto data = std::make_shared<const Eigen::MatrixXd>(std::move(fit_data));

  benchmark::RegisterBenchmark(("tll/pdf/" + method + "/n=1000").c_str(),
                               [bc, data](benchmark::State& st) {
                                 for (auto _ : st)
                                   benchmark::DoNotOptimize(bc.pdf(*data));
                               });
  benchmark::RegisterBenchmark(("tll/hfunc1/" + method + "/n=1000").c_str(),
                               [bc, data](benchmark::State& st) {
                                 for (auto _ : st)
                                   benchmark::DoNotOptimize(bc.hfunc1(*data));
                               });
  benchmark::RegisterBenchmark(("tll/hinv1/" + method + "/n=1000").c_str(),
                               [bc, data](benchmark::State& st) {
                                 for (auto _ : st)
                                   benchmark::DoNotOptimize(bc.hinv1(*data));
                               });
  benchmark::RegisterBenchmark(("tll/cdf/" + method + "/n=100").c_str(),
                               [bc, data](benchmark::State& st) {
                                 const Eigen::MatrixXd q = data->topRows(100);
                                 for (auto _ : st)
                                   benchmark::DoNotOptimize(bc.cdf(q));
                               });
}

void
register_tll_discrete_fit()
{
  const size_t n = 1000;
  auto data = std::make_shared<const Eigen::MatrixXd>(
    bench::discretize_first(bench::sim_data(BicopFamily::gaussian, 0, n)));
  FitControlsBicop controls({ BicopFamily::tll });
  benchmark::RegisterBenchmark("tll/fit_disc/constant/n=1000",
                               [data, controls](benchmark::State& st) {
                                 for (auto _ : st) {
                                   Bicop bc(BicopFamily::tll);
                                   bc.set_var_types({ "d", "c" });
                                   bc.fit(*data, controls);
                                   benchmark::DoNotOptimize(bc);
                                 }
                               });
}

struct Registrar
{
  Registrar()
  {
    for (const auto& method : { std::string("constant"),
                                std::string("linear"),
                                std::string("quadratic") }) {
      register_tll_fit(method, 1000);
      register_tll_eval(method);
    }
    register_tll_fit("constant", 5000);
    register_tll_fit("quadratic", 5000);
    // normalization is O(m^2) and n-independent while the fit is O(m^2 n),
    // so the smallest sample is where its share of the fit is largest
    register_tll_fit("constant", 200);
    register_tll_fit("quadratic", 200);
    register_tll_fit_grid_size("quadratic", 50, 1000);
    register_tll_discrete_fit();
  }
};
const Registrar registrar;

} // namespace
