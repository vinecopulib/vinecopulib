// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "helpers.hpp"
#include <algorithm>
#include <benchmark/benchmark.h>
#include <memory>

using namespace vinecopulib;

namespace {

using DataPtr = std::shared_ptr<const Eigen::MatrixXd>;

void
register_eval_benchmarks(BicopFamily family)
{
  const std::string fam_name = get_family_name(family);
  for (size_t n : { size_t(1000), size_t(10000) }) {
    auto data =
      std::make_shared<const Eigen::MatrixXd>(bench::sim_data(family, 0, n));
    const Bicop bc(family, 0, bench::family_parameters(family));
    const std::string suffix = fam_name + "/n=" + std::to_string(n);

    benchmark::RegisterBenchmark(("bicop/pdf/" + suffix).c_str(),
                                 [bc, data](benchmark::State& st) {
                                   for (auto _ : st)
                                     benchmark::DoNotOptimize(bc.pdf(*data));
                                 });
    benchmark::RegisterBenchmark(("bicop/cdf/" + suffix).c_str(),
                                 [bc, data](benchmark::State& st) {
                                   for (auto _ : st)
                                     benchmark::DoNotOptimize(bc.cdf(*data));
                                 });
    benchmark::RegisterBenchmark(("bicop/hfunc1/" + suffix).c_str(),
                                 [bc, data](benchmark::State& st) {
                                   for (auto _ : st)
                                     benchmark::DoNotOptimize(bc.hfunc1(*data));
                                 });
    benchmark::RegisterBenchmark(("bicop/hfunc2/" + suffix).c_str(),
                                 [bc, data](benchmark::State& st) {
                                   for (auto _ : st)
                                     benchmark::DoNotOptimize(bc.hfunc2(*data));
                                 });
    benchmark::RegisterBenchmark(("bicop/hinv1/" + suffix).c_str(),
                                 [bc, data](benchmark::State& st) {
                                   for (auto _ : st)
                                     benchmark::DoNotOptimize(bc.hinv1(*data));
                                 });
    benchmark::RegisterBenchmark(("bicop/hinv2/" + suffix).c_str(),
                                 [bc, data](benchmark::State& st) {
                                   for (auto _ : st)
                                     benchmark::DoNotOptimize(bc.hinv2(*data));
                                 });
    if (n == 1000) {
      benchmark::RegisterBenchmark(
        ("bicop/loglik/" + suffix).c_str(), [bc, data](benchmark::State& st) {
          for (auto _ : st)
            benchmark::DoNotOptimize(bc.loglik(*data));
        });
    }
  }
}

void
register_fit_benchmarks(BicopFamily family)
{
  const std::string fam_name = get_family_name(family);
  const size_t n = 1000;
  auto data =
    std::make_shared<const Eigen::MatrixXd>(bench::sim_data(family, 0, n));

  FitControlsBicop controls_mle({ family });
  controls_mle.set_parametric_method("mle");
  benchmark::RegisterBenchmark(
    ("bicop/fit_mle/" + fam_name + "/n=1000").c_str(),
    [family, data, controls_mle](benchmark::State& st) {
      for (auto _ : st) {
        Bicop bc(family);
        bc.fit(*data, controls_mle);
        benchmark::DoNotOptimize(bc);
      }
    });

  const auto& itau_families = bicop_families::itau;
  if (std::find(itau_families.begin(), itau_families.end(), family) !=
      itau_families.end()) {
    FitControlsBicop controls_itau({ family });
    controls_itau.set_parametric_method("itau");
    benchmark::RegisterBenchmark(
      ("bicop/fit_itau/" + fam_name + "/n=1000").c_str(),
      [family, data, controls_itau](benchmark::State& st) {
        for (auto _ : st) {
          Bicop bc(family);
          bc.fit(*data, controls_itau);
          benchmark::DoNotOptimize(bc);
        }
      });
  }
}

void
register_discrete_benchmarks(BicopFamily family)
{
  const std::string fam_name = get_family_name(family);
  const size_t n = 1000;
  auto data = std::make_shared<const Eigen::MatrixXd>(
    bench::discretize_first(bench::sim_data(family, 0, n)));
  Bicop bc(family, 0, bench::family_parameters(family));
  bc.set_var_types({ "d", "c" });

  benchmark::RegisterBenchmark(
    ("bicop/pdf_disc/" + fam_name + "/n=1000").c_str(),
    [bc, data](benchmark::State& st) {
      for (auto _ : st)
        benchmark::DoNotOptimize(bc.pdf(*data));
    });
  benchmark::RegisterBenchmark(
    ("bicop/hfunc1_disc/" + fam_name + "/n=1000").c_str(),
    [bc, data](benchmark::State& st) {
      for (auto _ : st)
        benchmark::DoNotOptimize(bc.hfunc1(*data));
    });
}

// Times the analytic parameter scores d log c / d theta_k. Independent of the
// optimizer; isolates the cost of the `logpdf_deriv_raw` leaves (the Student
// df score in particular).
void
register_score_benchmarks(BicopFamily family)
{
  const std::string fam_name = get_family_name(family);
  const auto par = bench::family_parameters(family);
  const int npars = static_cast<int>(par.rows());
  for (size_t n : { size_t(1000), size_t(10000) }) {
    auto data =
      std::make_shared<const Eigen::MatrixXd>(bench::sim_data(family, 0, n));
    const Bicop bc(family, 0, par);
    for (int k = 1; k <= npars; ++k) {
      const std::string sel = "par" + std::to_string(k);
      benchmark::RegisterBenchmark(("bicop/logpdf_deriv/" + fam_name + "/" +
                                    sel + "/n=" + std::to_string(n))
                                     .c_str(),
                                   [bc, data, sel](benchmark::State& st) {
                                     for (auto _ : st)
                                       benchmark::DoNotOptimize(
                                         bc.logpdf_deriv(*data, sel));
                                   });
    }
  }
}

void
register_smoke_benchmark()
{
  auto data = std::make_shared<const Eigen::MatrixXd>(
    bench::sim_data(BicopFamily::gaussian, 0, 100));
  const Bicop bc(
    BicopFamily::gaussian, 0, bench::family_parameters(BicopFamily::gaussian));
  benchmark::RegisterBenchmark("smoke/bicop_pdf",
                               [bc, data](benchmark::State& st) {
                                 for (auto _ : st)
                                   benchmark::DoNotOptimize(bc.pdf(*data));
                               });
}

struct Registrar
{
  Registrar()
  {
    const std::vector<BicopFamily> families = {
      BicopFamily::gaussian, BicopFamily::student, BicopFamily::clayton,
      BicopFamily::gumbel,   BicopFamily::frank,   BicopFamily::joe,
      BicopFamily::bb1,      BicopFamily::bb6,     BicopFamily::bb7,
      BicopFamily::bb8,      BicopFamily::tawn
    };
    register_smoke_benchmark();
    for (auto family : families) {
      register_eval_benchmarks(family);
      register_fit_benchmarks(family);
      register_score_benchmarks(family);
    }
    for (auto family : { BicopFamily::gaussian, BicopFamily::clayton }) {
      register_discrete_benchmarks(family);
    }
  }
};
const Registrar registrar;

} // namespace
