#!/usr/bin/env bash
# Build the benchmark suite and capture a baseline (or contender) run.
#
# Usage:
#   scripts/bench_baseline.sh [output.json] [extra --benchmark_* flags]
#
# Output defaults to bench_results/ (gitignored). Compare two runs with
# google/benchmark's compare.py:
#   python3 <build>/_deps/googlebenchmark-src/tools/compare.py \
#     benchmarks bench_results/baseline.json bench_results/contender.json
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${BENCH_BUILD_DIR:-$repo_root/build-bench}"
out="${1:-$repo_root/bench_results/benchmark_results.json}"
shift || true
mkdir -p "$(dirname "$out")"

cmake -S "$repo_root" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release \
      -DVINECOPULIB_BUILD_BENCHMARKS=ON -DBUILD_TESTING=OFF
cmake --build "$build_dir" -j"$(nproc)" --target bench_all parity_dump

"$build_dir/bin/bench_all" \
  --benchmark_repetitions="${BENCH_REPETITIONS:-5}" \
  --benchmark_report_aggregates_only=true \
  --benchmark_out="$out" --benchmark_out_format=json "$@"
echo "wrote $out"
