# Benchmark suite (google/benchmark), numerical parity harness, and the
# quadrature study. Not built by default and never installed:
#   cmake --preset bench && cmake --build build-bench -j
# See benchmarks/README.md and scripts/README.md.

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
