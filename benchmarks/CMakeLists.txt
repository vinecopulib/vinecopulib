# Benchmark suite (google/benchmark) + numerical parity harness.
#
# Not built by default and never installed. Enable with:
#   cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release \
#         -DVINECOPULIB_BUILD_BENCHMARKS=ON -DBUILD_TESTING=OFF
#   cmake --build build-bench -j
#
# Baseline capture (reference machine, performance governor pinned):
#   build-bench/bin/bench_all --benchmark_repetitions=10 \
#     --benchmark_report_aggregates_only=true \
#     --benchmark_out=baseline.json --benchmark_out_format=json
# Then after changes, write contender.json the same way and compare with
#   _deps/googlebenchmark-src/tools/compare.py benchmarks baseline.json contender.json
#
# Numerical parity:
#   build-bench/bin/parity_dump parity_master.json     (on master)
#   build-bench/bin/parity_dump parity_branch.json     (on the branch)
#   python3 scripts/compare_parity.py parity_master.json parity_branch.json

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_executable(bench_all
        src/bench_bicop_families.cpp
        src/bench_bicop_tll.cpp
        src/bench_interpolation.cpp
        src/bench_tools_stats.cpp
        src/bench_vinecop.cpp
        src/bench_main.cpp)
target_link_libraries(bench_all PRIVATE vinecopulib benchmark::benchmark)
target_include_directories(bench_all PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(parity_dump src/parity_dump.cpp)
target_link_libraries(parity_dump PRIVATE vinecopulib)

add_executable(quadrature_study src/quadrature_study.cpp)
target_link_libraries(quadrature_study PRIVATE vinecopulib)
