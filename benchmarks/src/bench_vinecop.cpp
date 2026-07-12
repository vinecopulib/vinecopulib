// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "helpers.hpp"
#include <benchmark/benchmark.h>
#include <memory>

using namespace vinecopulib;

namespace {

void
register_eval(size_t d)
{
  const size_t n = 1000;
  auto vc = std::make_shared<const Vinecop>(bench::make_gaussian_vine(d));
  auto u =
    std::make_shared<const Eigen::MatrixXd>(vc->simulate(n, false, 1, { 5 }));
  const std::string suffix = "d=" + std::to_string(d) + "/n=1000";

  for (size_t threads : { size_t(1), size_t(4) }) {
    benchmark::RegisterBenchmark(
      ("vinecop/pdf/" + suffix + "/threads=" + std::to_string(threads)).c_str(),
      [vc, u, threads](benchmark::State& st) {
        for (auto _ : st)
          benchmark::DoNotOptimize(vc->pdf(*u, threads));
      });
  }
  benchmark::RegisterBenchmark(("vinecop/rosenblatt/" + suffix).c_str(),
                               [vc, u](benchmark::State& st) {
                                 for (auto _ : st)
                                   benchmark::DoNotOptimize(vc->rosenblatt(*u));
                               });
  benchmark::RegisterBenchmark(("vinecop/inverse_rosenblatt/" + suffix).c_str(),
                               [vc, u](benchmark::State& st) {
                                 for (auto _ : st)
                                   benchmark::DoNotOptimize(
                                     vc->inverse_rosenblatt(*u));
                               });
  benchmark::RegisterBenchmark(
    ("vinecop/simulate/" + suffix).c_str(), [vc](benchmark::State& st) {
      for (auto _ : st)
        benchmark::DoNotOptimize(vc->simulate(1000, false, 1, { 5 }));
    });
}

void
register_scores(size_t d,
                BicopFamily family = BicopFamily::gaussian,
                double par = 0.5,
                const std::string& fam_label = "gaussian")
{
  const size_t n = 500;
  auto vc = std::make_shared<const Vinecop>(bench::make_vine(family, par, d));
  auto u =
    std::make_shared<const Eigen::MatrixXd>(vc->simulate(n, false, 1, { 5 }));
  const std::string suffix = fam_label + "/d=" + std::to_string(d) + "/n=500";

  for (bool step_wise : { true, false }) {
    const std::string mode = step_wise ? "stepwise" : "full";
    benchmark::RegisterBenchmark(
      ("vinecop/scores_" + mode + "/" + suffix).c_str(),
      [vc, u, step_wise](benchmark::State& st) {
        // scores() restores var_types internally, so it is safe to call
        // repeatedly on the same object; copying per iteration would time a
        // Vinecop deep-copy instead of scores()
        auto vc_copy = *vc;
        for (auto _ : st) {
          benchmark::DoNotOptimize(vc_copy.scores(*u, step_wise, 1));
        }
      });
  }
  benchmark::RegisterBenchmark(("vinecop/hessian_stepwise/" + suffix).c_str(),
                               [vc, u](benchmark::State& st) {
                                 auto vc_copy = *vc;
                                 for (auto _ : st) {
                                   benchmark::DoNotOptimize(
                                     vc_copy.hessian(*u, true, 1));
                                 }
                               });
  benchmark::RegisterBenchmark(
    ("vinecop/hessian_joint/" + suffix).c_str(), [vc, u](benchmark::State& st) {
      auto vc_copy = *vc;
      for (auto _ : st) {
        benchmark::DoNotOptimize(vc_copy.hessian(*u, false, 1));
      }
    });
}

void
register_select(size_t d, const std::string& label, FitControlsVinecop controls)
{
  const size_t n = 1000;
  auto vc = bench::make_gaussian_vine(d);
  auto u =
    std::make_shared<const Eigen::MatrixXd>(vc.simulate(n, false, 1, { 5 }));
  const std::string suffix = label + "/d=" + std::to_string(d) + "/n=1000";

  for (size_t threads : { size_t(1), size_t(4) }) {
    FitControlsVinecop ctrl = controls;
    ctrl.set_num_threads(threads);
    benchmark::RegisterBenchmark(
      ("vinecop/select/" + suffix + "/threads=" + std::to_string(threads))
        .c_str(),
      [u, ctrl](benchmark::State& st) {
        for (auto _ : st) {
          Vinecop fitted(*u, RVineStructure(), {}, ctrl);
          benchmark::DoNotOptimize(fitted);
        }
      });
  }
}

void
register_cdf()
{
  auto vc = std::make_shared<const Vinecop>(bench::make_gaussian_vine(5));
  auto u =
    std::make_shared<const Eigen::MatrixXd>(vc->simulate(100, false, 1, { 5 }));
  benchmark::RegisterBenchmark(
    "vinecop/cdf/d=5/n=100/N=10000", [vc, u](benchmark::State& st) {
      for (auto _ : st)
        benchmark::DoNotOptimize(vc->cdf(*u, 10000, 1, { 5 }));
    });
}

struct Registrar
{
  Registrar()
  {
    for (size_t d : { size_t(5), size_t(10), size_t(25) }) {
      register_eval(d);
    }
    register_scores(5);
    register_scores(5, BicopFamily::clayton, 3.0, "clayton");
    register_scores(5, BicopFamily::frank, 5.0, "frank");
    register_cdf();

    register_select(
      5, "itau", FitControlsVinecop(bicop_families::itau, "itau"));
    register_select(
      10, "itau", FitControlsVinecop(bicop_families::itau, "itau"));
    register_select(5, "tll", FitControlsVinecop({ BicopFamily::tll }));
  }
};
const Registrar registrar;

} // namespace
