# Benchmarks

A [google/benchmark](https://github.com/google/benchmark) suite (`bench_all`)
and a numerical parity harness (`parity_dump`), plus the quadrature study that
chose the Kendall's-tau integration rule (`quadrature_study`).

> **Not** `examples/benchmark/`. That is an unrelated, pre-existing `chrono`
> harness and has nothing to do with this suite.

## Building and running

```bash
cmake --preset bench          # Release + VINECOPULIB_BUILD_BENCHMARKS + native arch
cmake --build build-bench -j
build-bench/bin/bench_all
```

Or spelled out:

```bash
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release \
      -DVINECOPULIB_BUILD_BENCHMARKS=ON -DBUILD_TESTING=OFF
cmake --build build-bench -j
```

Useful flags:

| flag | effect |
| --- | --- |
| `--benchmark_filter=<regex>` | run a subset, e.g. `--benchmark_filter='bicop/pdf/Clayton'` |
| `--benchmark_repetitions=10` | repeat, for the aggregate statistics |
| `--benchmark_report_aggregates_only=true` | print only mean/median/stddev |
| `--benchmark_out=x.json --benchmark_out_format=json` | machine-readable output |

## Layout

| file | covers |
| --- | --- |
| `src/bench_bicop_families.cpp` | per-family evaluation, fitting, discrete dispatch |
| `src/bench_bicop_tll.cpp` | the nonparametric `tll` family |
| `src/bench_interpolation.cpp` | `InterpolationGrid` |
| `src/bench_tools_stats.cpp` | pseudo-observations, QRNG, distributions, dependence measures |
| `src/bench_vinecop.cpp` | vine selection, evaluation, simulation, Rosenblatt |
| `src/bench_main.cpp` | hand-rolled `main` |
| `src/helpers.hpp` | seeded fixtures shared by all of the above |

`bench_main.cpp` is hand-written rather than `BENCHMARK_MAIN()` only so that the
file stays a normal translation unit; it does exactly what the macro does.

## Adding a case

Cases are registered at static-initialization time from a file-local
`Registrar`, not with the `BENCHMARK()` macro, because the names are built at
runtime from the family list:

```cpp
struct Registrar
{
  Registrar()
  {
    for (auto family : families) {
      register_eval_benchmarks(family);
    }
  }
};
const Registrar registrar;   // registers when the TU is loaded
```

Follow the slash-namespaced naming convention — `bicop/pdf/Clayton/n=1000` — so
that `--benchmark_filter` can select by surface, family or size. Use the seeded
fixtures in `helpers.hpp`: every input must be deterministic, or run-to-run
comparisons mean nothing.

## Comparing two runs

`scripts/bench_baseline.sh` is the one-command driver; see
[`../scripts/README.md`](../scripts/README.md) for the whole workflow, including
the numerical parity harness that guards against a "speedup" that quietly
changes results.

## The quadrature study

`src/quadrature_study.cpp` measures accuracy *and* time for several quadrature
rules on the four `parameters_to_tau` integrands (BB6, BB7, BB8, Tawn), against
50-digit references computed in the same harness with
`boost::multiprecision`. It cross-checks its restated integrands against
`Bicop::get_tau()` first, so a mismatch means the study is wrong rather than the
library.

It is a study, not a regression test: run it when changing an integrand or the
integration rule.

```bash
cmake --build build-bench -j --target quadrature_study
build-bench/bin/quadrature_study
```
