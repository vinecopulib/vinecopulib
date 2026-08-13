# Scripts

Tooling for the two things that are easy to break silently: numerical results
and the release version.

## Numerical parity

A performance change is only acceptable if it does not move the numbers. The
parity harness makes that checkable.

`benchmarks/src/parity_dump.cpp` writes one JSON file of outputs across the
surfaces the optimization work touches — bivariate evaluation and fitting, TLL,
`tools_stats`, vine evaluation and selection, and the Kendall's-tau sweeps — on a
fixed 21x21 lattice with fixed seeds. It is fully deterministic, so two dumps
differ only where the code does.

```bash
# on the baseline
git checkout main
cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release \
      -DVINECOPULIB_BUILD_BENCHMARKS=ON -DBUILD_TESTING=OFF
cmake --build build-bench -j --target parity_dump
build-bench/bin/parity_dump bench_results/parity_main.json

# on the branch
git checkout my-branch
cmake --build build-bench -j --target parity_dump
build-bench/bin/parity_dump bench_results/parity_branch.json

python3 scripts/compare_parity.py \
    bench_results/parity_main.json bench_results/parity_branch.json
```

`compare_parity.py` flattens both files to slash-joined paths
(`bicop_eval/Clayton/rot0/pdf`) and compares elementwise, choosing a tolerance by
**longest matching prefix** from `parity_gates.json`. Anything unmatched must be
bit-identical. `--strict` demands that everywhere, which is the right setting for
a change that claims to be a pure refactor.

### Re-baselining a gate

When a change is *meant* to move numbers, relax the gate and say why. The
`_comment*` keys exist for that and are skipped when the gates are read:

```json
"_comment_tau": "tau integrands reformulated and integrated with tanh-sinh ...",
"tau/": 1e-07,
```

A gate without a recorded reason is indistinguishable from a regression someone
gave up on. Keep the relaxation as narrow as the affected keys.

Note that the lattice covers interiors. Near-boundary corners live in separate
`*_edge` keys, added after a round of τ bugs that the interior sweeps missed
entirely — if a change affects extreme parameters, check that those moved as
intended.

The golden-value tests in the suite
(`test/src_test/test_tools_stats.cpp`, `test_bicop_kernel.cpp`,
`test_bicop_sanity_checks.cpp`) are the part of this that CI enforces on every
run; the dump comparison is a manual, before/after tool.

## Benchmark baselines

`bench_baseline.sh` builds the suite and captures one run:

```bash
scripts/bench_baseline.sh bench_results/baseline.json --benchmark_repetitions=10
# ... make the change ...
scripts/bench_baseline.sh bench_results/contender.json --benchmark_repetitions=10
python3 build-bench/_deps/googlebenchmark-src/tools/compare.py \
    benchmarks bench_results/baseline.json bench_results/contender.json
```

Output goes to `bench_results/` (gitignored). Pin the CPU governor and use
repetitions if you intend to trust small differences. See
[`../benchmarks/README.md`](../benchmarks/README.md) for the suite itself.

## Version consistency

`check_version.py` cross-checks the project version everywhere it appears —
`CMakeLists.txt`, both macros in `version.hpp`, `CITATION.cff`, `.zenodo.json`,
and the top heading of `NEWS.md`:

```bash
python3 scripts/check_version.py            # every PR (CI runs this)
python3 scripts/check_version.py --released # also requires a dated NEWS heading
python3 scripts/check_version.py --tag v1.0.0
```

`--released` is what `release.yml` runs against a tag. The weekly
`release-drift.yml` catches the opposite failure: a version dated in `NEWS.md`
with no tag pushed, which is exactly how 0.8.0 ended up written up but never
released.
