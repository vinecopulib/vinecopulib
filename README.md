# vinecopulib

[![Build Status](https://github.com/vinecopulib/vinecopulib/actions/workflows/continuous_integration.yml/badge.svg?branch=main)](https://github.com/vinecopulib/vinecopulib/actions/workflows/continuous_integration.yml)
[![Coverage Status](https://img.shields.io/codecov/c/github/vinecopulib/vinecopulib/main.svg)](https://codecov.io/github/vinecopulib/vinecopulib?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-website-blue.svg)](https://vinecopulib.github.io/vinecopulib/)
[![DOI](https://zenodo.org/badge/76354683.svg)](https://zenodo.org/badge/latestdoi/76354683)

A header-only C++ library for **vine copula models**, built on
[Eigen](https://eigen.tuxfamily.org/). It provides high-performance
implementations of the core features of the
[VineCopula R package](https://github.com/tnagler/VineCopula) — inference
algorithms for both vine copula and bivariate copula models — with

- nonparametric and multi-parameter families,
- shorter runtimes and lower memory consumption, especially in high dimensions,
- a sleeker, more modern API,
- interfaces for [R](https://github.com/vinecopulib/rvinecopulib) and
  [Python](https://github.com/vinecopulib/pyvinecopulib).

## Requirements

| | |
| --- | --- |
| Compiler | anything with **C++17** |
| CMake | **3.14** or later |
| [Eigen](https://eigen.tuxfamily.org/) | 3.3 or later |
| [Boost](https://www.boost.org/) | **1.75** or later (headers only) |
| [wdm](https://github.com/tnagler/wdm) | **0.3.0** or later — found via `find_package`, otherwise fetched at configure time |
| R (optional) | only for the parity tests, with VineCopula >= 2.6.2 |

## Quickstart

```bash
git clone https://github.com/vinecopulib/vinecopulib.git
cd vinecopulib
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build -j
sudo cmake --install build
```

Then, from your own project:

```cmake
find_package(vinecopulib REQUIRED)
target_link_libraries(my_target PRIVATE vinecopulib::vinecopulib)
```

```cpp
#include <vinecopulib.hpp>

using namespace vinecopulib;

int main() {
  auto data = tools_stats::simulate_uniform(500, 5);

  // select the structure, the families, and the parameters
  Vinecop model(data);
  std::cout << model.str() << std::endl;

  auto pdf = model.pdf(data);
  auto sim = model.simulate(100);
}
```

The library is header-only by default; `-DVINECOPULIB_PRECOMPILED=ON` compiles it
into a real library instead.

## Documentation

The [website](https://vinecopulib.github.io/vinecopulib/) is the reference:

- [Setup](https://vinecopulib.github.io/vinecopulib/setup.html) — building, installing, and every CMake option
- [Bivariate copulas](https://vinecopulib.github.io/vinecopulib/overview-bicop.html) and [vine copulas](https://vinecopulib.github.io/vinecopulib/overview-vinecop.html)
- [Discrete variables](https://vinecopulib.github.io/vinecopulib/discrete.html), [conditional simulation](https://vinecopulib.github.io/vinecopulib/conditional.html), [scores and Hessians](https://vinecopulib.github.io/vinecopulib/inference.html)

Every code example on those pages is compiled and run in CI, from the sources
under [`docs/snippets/`](docs/snippets).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the workflow and
[AGENTS.md](AGENTS.md) for the engineering invariants. The benchmark and
numerical-parity tooling is documented in
[benchmarks/README.md](benchmarks/README.md) and
[scripts/README.md](scripts/README.md); releases in
[docs/releasing.md](docs/releasing.md).

## Citation

If you use vinecopulib in your work, please cite it. [CITATION.cff](CITATION.cff)
carries the machine-readable metadata, which GitHub renders as a citation in the
sidebar.

## Contact

[Open an issue](https://github.com/vinecopulib/vinecopulib/issues/new) or write to
[info@vinecopulib.org](mailto:info@vinecopulib.org).
