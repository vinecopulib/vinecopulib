// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "helpers.hpp"
#include <benchmark/benchmark.h>
#include <memory>
#include <vinecopulib/misc/tools_stats.hpp>

using namespace vinecopulib;

namespace {

void
register_pseudo_obs()
{
  for (size_t n : { size_t(1000), size_t(10000) }) {
    for (size_t d : { size_t(2), size_t(10) }) {
      auto x = std::make_shared<const Eigen::MatrixXd>(
        tools_stats::simulate_uniform(n, d, false, { 5 }));
      const std::string suffix =
        "n=" + std::to_string(n) + "/d=" + std::to_string(d);
      benchmark::RegisterBenchmark(
        ("stats/to_pseudo_obs/" + suffix).c_str(), [x](benchmark::State& st) {
          for (auto _ : st)
            benchmark::DoNotOptimize(tools_stats::to_pseudo_obs(*x));
        });
    }
  }
  // weighted variant
  auto x = std::make_shared<const Eigen::MatrixXd>(
    tools_stats::simulate_uniform(1000, 10, false, { 5 }));
  auto w = std::make_shared<const Eigen::VectorXd>(
    Eigen::VectorXd::LinSpaced(1000, 0.5, 1.5));
  benchmark::RegisterBenchmark(
    "stats/to_pseudo_obs_weighted/n=1000/d=10", [x, w](benchmark::State& st) {
      for (auto _ : st)
        benchmark::DoNotOptimize(tools_stats::to_pseudo_obs(*x, "average", *w));
    });
}

void
register_genz()
{
  auto u = std::make_shared<const Eigen::MatrixXd>(
    tools_stats::simulate_uniform(10000, 2, false, { 5 }));
  auto z = std::make_shared<const Eigen::MatrixXd>(tools_stats::qnorm(*u));
  for (double rho : { 0.1, 0.5, 0.9 }) {
    benchmark::RegisterBenchmark(
      ("stats/pbvnorm/rho=" + std::to_string(rho).substr(0, 3) + "/n=10000")
        .c_str(),
      [z, rho](benchmark::State& st) {
        for (auto _ : st)
          benchmark::DoNotOptimize(tools_stats::pbvnorm(*z, rho));
      });
  }
  auto zt = std::make_shared<const Eigen::MatrixXd>(tools_stats::qt(*u, 4.0));
  benchmark::RegisterBenchmark(
    "stats/pbvt/nu=4/rho=0.5/n=10000", [zt](benchmark::State& st) {
      for (auto _ : st)
        benchmark::DoNotOptimize(tools_stats::pbvt(*zt, 4, 0.5));
    });
}

void
register_qrng()
{
  for (size_t d : { size_t(5), size_t(50) }) {
    const std::string suffix = "n=10000/d=" + std::to_string(d);
    benchmark::RegisterBenchmark(
      ("stats/ghalton/" + suffix).c_str(), [d](benchmark::State& st) {
        for (auto _ : st)
          benchmark::DoNotOptimize(tools_stats::ghalton(10000, d, { 5 }));
      });
    benchmark::RegisterBenchmark(
      ("stats/sobol/" + suffix).c_str(), [d](benchmark::State& st) {
        for (auto _ : st)
          benchmark::DoNotOptimize(tools_stats::sobol(10000, d, { 5 }));
      });
  }
  benchmark::RegisterBenchmark(
    "stats/simulate_uniform/n=10000/d=5", [](benchmark::State& st) {
      for (auto _ : st)
        benchmark::DoNotOptimize(
          tools_stats::simulate_uniform(10000, 5, false, { 5 }));
    });
}

void
register_mcor()
{
  auto x = std::make_shared<const Eigen::MatrixXd>(
    bench::sim_data(BicopFamily::gaussian, 0, 1000));
  benchmark::RegisterBenchmark(
    "stats/pairwise_mcor/n=1000", [x](benchmark::State& st) {
      for (auto _ : st)
        benchmark::DoNotOptimize(tools_stats::pairwise_mcor(*x));
    });
}

struct Registrar
{
  Registrar()
  {
    register_pseudo_obs();
    register_genz();
    register_qrng();
    register_mcor();
  }
};
const Registrar registrar;

} // namespace
