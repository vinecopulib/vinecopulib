// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "helpers.hpp"
#include <benchmark/benchmark.h>
#include <memory>
#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

using namespace vinecopulib;

namespace {

//! Builds a grid mimicking a fitted TLL surface: pnorm of an equidistant
//! z-grid with a smooth positive density-like surface.
tools_interpolation::InterpolationGrid
make_grid(size_t m)
{
  Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(m, -3.25, 3.25);
  Eigen::VectorXd grid_points = tools_stats::pnorm(z);
  Eigen::MatrixXd values(m, m);
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      values(i, j) = std::exp(-0.5 * (z(i) * z(i) + z(j) * z(j)) / 2.0);
    }
  }
  return tools_interpolation::InterpolationGrid(grid_points, values, 3);
}

void
register_interpolation(size_t m)
{
  auto grid =
    std::make_shared<tools_interpolation::InterpolationGrid>(make_grid(m));
  auto queries = std::make_shared<const Eigen::MatrixXd>(
    tools_stats::simulate_uniform(10000, 2, false, { 5 }));
  const std::string suffix = "m=" + std::to_string(m) + "/n=10000";

  benchmark::RegisterBenchmark(("interp/interpolate/" + suffix).c_str(),
                               [grid, queries](benchmark::State& st) {
                                 for (auto _ : st)
                                   benchmark::DoNotOptimize(
                                     grid->interpolate(*queries));
                               });
  benchmark::RegisterBenchmark(("interp/integrate_1d/" + suffix).c_str(),
                               [grid, queries](benchmark::State& st) {
                                 for (auto _ : st)
                                   benchmark::DoNotOptimize(
                                     grid->integrate_1d(*queries, 1));
                               });
  benchmark::RegisterBenchmark(
    ("interp/integrate_2d/m=" + std::to_string(m) + "/n=1000").c_str(),
    [grid, queries](benchmark::State& st) {
      const Eigen::MatrixXd q = queries->topRows(1000);
      for (auto _ : st)
        benchmark::DoNotOptimize(grid->integrate_2d(q));
    });
}

struct Registrar
{
  Registrar()
  {
    register_interpolation(30);
    register_interpolation(50);
  }
};
const Registrar registrar;

} // namespace
