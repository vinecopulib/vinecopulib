// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "helpers.hpp"
#include <benchmark/benchmark.h>
#include <cmath>
#include <memory>
#include <vector>
#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

using namespace vinecopulib;

namespace {

//! The grid points a fitted TLL model uses: pnorm of an equidistant z-grid.
Eigen::VectorXd
make_grid_points(size_t m)
{
  Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(m, -3.25, 3.25);
  return tools_stats::pnorm(z);
}

//! A real TLL fit at strong dependence, pushed back off the uniform-margin
//! manifold by a row/column rescaling. Normalizing an already-normalized grid
//! would measure the convergence test rather than the work, and strong
//! dependence is what keeps the iteration from converging early.
Eigen::MatrixXd
unnormalized_values(size_t m)
{
  FitControlsBicop controls({ BicopFamily::tll });
  controls.set_nonparametric_grid_size(m);
  Bicop gauss(BicopFamily::gaussian, 0, Eigen::MatrixXd::Constant(1, 1, 0.9));
  Bicop bc(BicopFamily::tll);
  bc.fit(gauss.simulate(1000, false, { 5 }), controls);

  Eigen::MatrixXd values = bc.get_parameters();
  Eigen::VectorXd s(m);
  for (size_t i = 0; i < m; ++i)
    s(i) = std::exp(std::sin(3.0 * static_cast<double>(i)));
  values.array().colwise() *= s.array();
  values.array().rowwise() *= s.transpose().array();
  return values;
}

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

//! Isolates `normalize_margins`: `set_values` differs from the constructor
//! by not rebuilding the cell-lookup table, whose 1024 binary searches cost
//! about as much as a sweep. `norm=0` is the baseline to subtract.
void
register_normalize(size_t m, const std::vector<int>& norm_times)
{
  auto grid_points = make_grid_points(m);
  auto values = std::make_shared<const Eigen::MatrixXd>(unnormalized_values(m));
  auto grid = std::make_shared<tools_interpolation::InterpolationGrid>(
    grid_points, *values, 0);

  for (int times : norm_times) {
    benchmark::RegisterBenchmark(("interp/set_values/m=" + std::to_string(m) +
                                  "/norm=" + std::to_string(times))
                                   .c_str(),
                                 [grid, values, times](benchmark::State& st) {
                                   for (auto _ : st)
                                     grid->set_values(*values, times);
                                 });
  }

  if (m == 30) {
    benchmark::RegisterBenchmark(
      "interp/ctor/m=30/norm=3", [grid_points, values](benchmark::State& st) {
        for (auto _ : st)
          benchmark::DoNotOptimize(
            tools_interpolation::InterpolationGrid(grid_points, *values, 3));
      });
  }
}

struct Registrar
{
  Registrar()
  {
    register_interpolation(30);
    register_interpolation(50);
    register_normalize(30, { 0, 1, 3, 10, 25 });
    register_normalize(50, { 0, 3, 25 });
  }
};
const Registrar registrar;

} // namespace
