# AGENTS.md

Normative engineering spec for contributors and coding agents working on
this repository: scope, API-stability policy, module boundaries, build
modes, conventions, and where to look for what.

For the **user-facing pitch** and install pointers see
[README.md](README.md); for the **release-by-release history** see
[NEWS.md](NEWS.md); for **prose documentation** of the maths and the API
see the `docs/` folder (`*.dox`) and the rendered
[website](https://vinecopulib.github.io/vinecopulib/). This file does not
duplicate those — it concentrates on engineering invariants that survive
across PRs.

## Project overview

`vinecopulib` is a **header-only C++14 library** for vine-copula and
bivariate-copula inference, built on
[Eigen](https://eigen.tuxfamily.org/). It provides high-performance
implementations of the core features of the
[VineCopula R package](https://github.com/tnagler/VineCopula): inference
algorithms for both vine-copula and bivariate-copula models, with
nonparametric and multi-parameter families, lower runtimes, and lower
memory consumption in high dimensions.

It depends on one sister library, [`wdm`](https://github.com/tnagler/wdm)
(weighted Kendall's τ / Spearman's ρ / Pearson etc.), pulled in at
configure time (see [Build & tooling](#build--tooling)). It is **not** a
git submodule, and there is no vendored `lib/` tree.

This repository is the **upstream**. Two interface libraries wrap it and
track its tags: [`rvinecopulib`](https://github.com/vinecopulib/rvinecopulib)
(R) and [`pyvinecopulib`](https://github.com/vinecopulib/pyvinecopulib)
(Python). Behaviour and API changes here ripple downstream.

Three design principles inform the rest of this file:

- **This repo is upstream.** A signature change, a renamed CMake option,
  or a shifted numerical default propagates to R and Python. Weigh every
  public-API change against that blast radius and record it in
  [NEWS.md](NEWS.md).
- **Header-only, with an optional precompiled mode.** The `.hpp`
  (declaration) / `.ipp` (inline implementation) split is the single
  source of truth; there are no generated, gitignored artefacts to keep
  in sync. The `VINECOPULIB_PRECOMPILED` build merely *derives* `.cpp`
  translation units from the same `.ipp` files — see
  [Build modes](#build-modes).
- **Code is quantitatively sensitive.** Pseudo-observation transforms,
  h-functions, the Rosenblatt cascade, family parameterisations,
  Kendall's-τ ↔ parameter maps, TLL interpolation grids, and JSON
  round-trips all encode mathematically precise behaviour. Small
  "obvious-looking" changes can silently break copula identities. Treat
  numerical paths as correctness-critical and prefer round-trip / parity
  tests over structural ones.

## Scope

### Included

- **Bivariate copula modelling** — every family in `BicopFamily`
  ([bicop/family.hpp](include/vinecopulib/bicop/family.hpp)): `indep`,
  `gaussian`, `student`, `clayton`, `gumbel`, `frank`, `joe`,
  `bb1/6/7/8`, `tawn`, and the nonparametric `tll` (transformation local
  likelihood); with rotations (`0/90/180/270`), mixed-discrete handling
  (`var_types` of `"c"`/`"d"`), and family-set constraints via
  `FitControlsBicop`.
- **Vine copula modelling** — `Vinecop` with Dissmann selection over
  maximum spanning trees (`mst_prim` / `mst_kruskal`) and random
  spanning trees (`random_weighted` / `random_unweighted`), and
  user-supplied `RVineStructure` / `CVineStructure` / `DVineStructure`.
  Truncation, dependence thresholding, sparse (mBICv) selection,
  threading, and family sets are exposed through `FitControlsVinecop`.
- **Vine structures** — `RVineStructure` and its `CVineStructure` /
  `DVineStructure` specialisations, the underlying `TriangularArray<T>`
  container, and random-structure simulation.
- **Dependence measures** — via `wdm` (Kendall's τ, Spearman's ρ, etc.),
  used as edge-weighting criteria during structure selection.
- **Quasi-random sampling** — `sobol` (up to 21201 dimensions) and
  `ghalton`, behind `simulate_uniform`.
- **Pseudo-observations** — `to_pseudo_obs` (rank-normalise to the unit
  hypercube), the canonical input transform for copula fitting.
- **Numerics infrastructure** — bound-constrained MLE (Powell's BOBYQA),
  bilinear interpolation grids for TLL, 1-d numerical integration, a
  thread pool, and JSON (de)serialisation.

### Excluded (explicit)

- **Copula families outside the bound set.** New parametric families are
  added *here*, in `bicop/` (see [Extension points](#extension-points)) —
  but the set is otherwise fixed; don't reach for ad-hoc families.
- **1-d marginal estimation (`kde1d`).** Marginal KDE lives only in the
  downstream `pyvinecopulib` wheel. `vinecopulib` operates on copula
  data (uniform margins); it does not ship `kde1d`.
- **Binding / interface code.** The R (`rvinecopulib`) and Python
  (`pyvinecopulib`) wrappers are separate repositories. The only
  binding-aware hook here is the `INTERFACED_FROM_R` conditional in
  [misc/tools_interface.hpp](include/vinecopulib/misc/tools_interface.hpp)
  (RcppThread + R printing/interrupts when compiled inside R).
- **Higher-level estimators.** Scikit-learn-style estimators, PyTorch
  evaluators, and forest ensembles are entirely a `pyvinecopulib`
  concern. Nothing of the sort belongs here.

## API stability & releases

`vinecopulib` follows **semantic versioning** (`MAJOR.MINOR.PATCH`;
currently 0.7.3). Because the R and Python interfaces pin a tag of this
repo, public-API changes are real breaks for downstream users.

- **Record every change in [NEWS.md](NEWS.md).** Each release has
  `### BREAKING API CHANGES`, `### NEW FEATURES`, and `### BUG FIXES`
  sections, with every item tagged by its PR number (e.g. `(#637)`).
  There is no separate `CHANGELOG.md`.
- **Bump both version macros in
  [version.hpp](include/vinecopulib/version.hpp)** on release:
  `VINECOPULIB_VERSION` (encoded integer, `000703` for 0.7.3 —
  `major*100000 + minor*100 + patch`) and `VINECOPULIB_LIB_VERSION`
  (string, `"0_7_3"`).
- **Prefer deprecation over a hard break.** The `DEPRECATED` macro
  (defined in
  [vinecop/fit_controls.hpp](include/vinecopulib/vinecop/fit_controls.hpp))
  marks superseded API while keeping it callable — e.g.
  `Vinecop::select_all` / `select_families` (use `select`),
  `FitControlsVinecop::get_truncation_level` (use `get_trunc_lvl`). When
  a rename is unavoidable, ship both for a cycle and document the
  migration. Recent examples: `mst_algorithm` → `tree_algorithm` and the
  `"prim"`/`"kruskal"` → `"mst_prim"`/`"mst_kruskal"` value renames
  (#637); the CMake option `VINECOPULIB_BUILD_SHARED_LIBS` →
  `VINECOPULIB_PRECOMPILED` (#641).

## Project layout

```text
vinecopulib/
  AGENTS.md, CLAUDE.md           # this file + thin pointer (`@AGENTS.md`)
  README.md, NEWS.md, LICENSE, .zenodo.json
  CMakeLists.txt                 # 19-line entry point; includes the cmake/ modules
  .clang-format, .clang-tidy     # Mozilla style; advisory tidy checks
  codecov.yml, .codacy.yml       # coverage / static-analysis service config

  include/
    vinecopulib.hpp              # umbrella header — the canonical user include
    vinecopulib/
      version.hpp                # VINECOPULIB_VERSION / VINECOPULIB_LIB_VERSION
      bicop/
        class.hpp                # Bicop (the public facade)
        family.hpp               # BicopFamily enum + bicop_families:: group lists
        fit_controls.hpp         # FitControlsBicop
        abstract.hpp             # AbstractBicop (internal base)
        parametric.hpp, elliptical.hpp, archimedean.hpp,
        extreme_value.hpp, kernel.hpp           # family base classes
        gaussian.hpp ... tll.hpp                # the 13 families
        tools_select.hpp         # bicop candidate generation / pre-selection
        implementation/*.ipp     # inline definitions (one per .hpp)
      vinecop/
        class.hpp                # Vinecop (the public facade)
        rvine_structure.hpp      # RVineStructure + C/DVineStructure
        fit_controls.hpp         # FitControlsVinecop (: public FitControlsBicop)
        tools_select.hpp         # VinecopSelector (Dissmann + random trees)
        implementation/*.ipp
      misc/
        tools_stats.hpp          # pseudo-obs, QRNG, distributions, dep. measures
        tools_stats_sobol.hpp, tools_stats_ghalton.hpp
        tools_eigen.hpp, tools_stl.hpp
        tools_optimization.hpp, tools_bobyqa.hpp   # BOBYQA MLE
        tools_interpolation.hpp  # InterpolationGrid (TLL)
        tools_serialization.hpp, nlohmann_json.hpp # JSON
        tools_thread.hpp, tools_interface.hpp, tools_batch.hpp  # parallelism
        tools_integration.hpp, tools_constants.hpp, tools_optional.hpp
        triangular_array.hpp     # TriangularArray<T> (vine container)
        fit_controls.hpp         # FitControlsConfig (optional-field builder)
        implementation/*.ipp

  cmake/                         # options, dependency/header discovery, targets
    options.cmake, findDependencies.cmake, findHeaders.cmake,
    compilerDefOpt.cmake, buildTargets.cmake, codeCoverage.cmake,
    findR.cmake, printInfo.cmake
    templates/                   # Config.cmake.in, rscript.hpp.in, *.R parity tests
  test/                          # GoogleTest; thin test_*.cpp + src_test/ logic
  examples/                      # bicop/, vinecop/, benchmark/ — standalone projects
  docs/                          # Doxyfile(.in), Doxyfile-mcss, *.dox, create-docs.md
  .github/workflows/             # continuous_integration.yml
```

There is **no** `src/` of library `.cpp` files, **no** `lib/` /
`contrib/` vendored tree, **no** git submodules, and **no**
`CONTRIBUTING.md` — the workflow lives in this file and in CI.

## Build & tooling

CMake **≥ 3.10**, C++ **14**. The 19-line root
[CMakeLists.txt](CMakeLists.txt) sets the standard, the project version,
and includes the `cmake/` modules in order.

### Build modes

The `VINECOPULIB_PRECOMPILED` option (default **OFF**) selects between two
ways of building the `vinecopulib` target — see
[cmake/buildTargets.cmake](cmake/buildTargets.cmake) and
[cmake/findHeaders.cmake](cmake/findHeaders.cmake):

- **OFF (header-only, the default):**
  `add_library(vinecopulib INTERFACE)`. The install ships the `.hpp` *and*
  `.ipp` trees; users just `#include <vinecopulib.hpp>`.
- **ON (precompiled):** `findHeaders.cmake` derives `.cpp` translation
  units from the `.ipp` files (it strips the `inline` keyword and the
  license header and prepends the matching `#include`), compiles them into
  a real library (PIC, `WINDOWS_EXPORT_ALL_SYMBOLS`), and installs the
  binary instead of the `.ipp` files. This is what CI builds and installs.

The `.ipp` files are the source of truth in both modes — never hand-author
a `.cpp` under the generated tree.

### Options

| Option | Default | Effect |
| --- | --- | --- |
| `VINECOPULIB_PRECOMPILED` | `OFF` | Compile a real library vs. header-only `INTERFACE`. |
| `BUILD_TESTING` | `ON` | Build the GoogleTest suite (fetches GoogleTest + finds Rscript). |
| `OPT_ASAN` | `ON` | Add `-fsanitize=address,undefined` to Debug builds. |
| `CODE_COVERAGE` | `OFF` | Add `--coverage` flags (Debug + `BUILD_TESTING`, non-Windows). |
| `STRICT_COMPILER` | `OFF` | GCC: `-Werror` plus a long list of `-W…` flags. |
| `BUILD_DOC` | `OFF` | Add the Doxygen `doc` target. |

Default `CMAKE_BUILD_TYPE` is `Release` when unset.

### Dependencies

Discovered in
[cmake/findDependencies.cmake](cmake/findDependencies.cmake):

- **Eigen3** — `find_package(Eigen3 REQUIRED CONFIG)`; or set
  `EIGEN3_INCLUDE_DIR`. Target: `Eigen3::Eigen`.
- **Boost ≥ 1.56** — CONFIG mode then MODULE fallback; or set
  `Boost_INCLUDE_DIRS`. Used for math distributions (Student-t etc.), RNG
  (`boost::random`), `boost::optional` (C++14), and `odeint`
  (integration). Target: `Boost::boost`.
- **wdm 0.2.6** — `find_package(wdm 0.2.6 QUIET)` with a **FetchContent
  fallback** that clones `tnagler/wdm`. Installed under
  `<prefix>/include/vinecopulib/wdm/`.
- **Threads** — required.
- **GoogleTest 1.14** — FetchContent, only when `BUILD_TESTING`.
- **Rscript** (optional) — enables the R parity tests; see
  [cmake/findR.cmake](cmake/findR.cmake).
- **Doxygen** — only when `BUILD_DOC`.

Global compile definitions set in
[cmake/compilerDefOpt.cmake](cmake/compilerDefOpt.cmake) (and mirrored in
the umbrella header): `BOOST_NO_AUTO_PTR`,
`BOOST_ALLOW_DEPRECATED_HEADERS`, `BOOST_MATH_PROMOTE_DOUBLE_POLICY=false`,
`BOOST_ALL_NO_LIB`, `USE_BOOST`. Release flags `-O3 -march=native -DNDEBUG`
(an `-mcpu=apple-m1` variant on Apple silicon); always-on warnings
`-Wall -Wextra -Werror=return-type`.

### Build / test / install

The authoritative recipes are in CI
([.github/workflows/continuous_integration.yml](.github/workflows/continuous_integration.yml)):

```bash
# Debug build (header-only) + run tests
mkdir debug && cd debug
cmake .. -DCMAKE_BUILD_TYPE=Debug
make
bin/test_all            # or: ctest

# Release build, precompiled, install, then test
mkdir release && cd release
cmake .. -DCMAKE_BUILD_TYPE=Release -DVINECOPULIB_PRECOMPILED=ON
make && sudo make install
bin/test_all
```

`make` produces the test executables under `<build>/bin/`. ctest is wired
up via `gtest_discover_tests(test_all)` on non-Windows and
`add_test(NAME vinecopulib_tests COMMAND test_all)` on Windows
([test/CMakeLists.txt](test/CMakeLists.txt)).

### Formatting, linting, docs

- **clang-format** — [.clang-format](.clang-format) is exactly
  `BasedOnStyle: Mozilla`. CI's `style` job runs
  `clang-format-14 --dry-run --Werror` over `*.{h,hpp,ipp,cpp,cc}`,
  excluding the vendored
  [misc/nlohmann_json.hpp](include/vinecopulib/misc/nlohmann_json.hpp).
  Format with clang-format 14 before committing.
- **clang-tidy** — [.clang-tidy](.clang-tidy) enables a broad check set
  (`bugprone-*`, `performance-*`, `modernize-*`, `cppcoreguidelines-*`,
  …) but is **advisory**: the CI clang-tidy job is commented out and not
  enforced.
- **docs** — `cmake .. -DBUILD_DOC=ON && make doc` runs Doxygen. The
  website is built with [m.css](https://github.com/mosra/m.css) via
  `doxygen.py docs/Doxyfile-mcss`; see
  [docs/create-docs.md](docs/create-docs.md).

## Working on this repo

### Inspection order

Before changing code, read in this order:

1. `AGENTS.md` (this file) — invariants and boundaries.
2. `docs/*.dox` — the prose primers (`overview.dox`,
   `overview-bicop.dox`, `overview-vinecop.dox`, `setup.dox`) and the
   rendered website.
3. The public header you're about to touch (the `.hpp`), then its
   `implementation/*.ipp`.
4. The matching `test/test_<topic>.cpp` (and its logic in
   `test/src_test/`) for expected behaviour.

Match existing local patterns rather than introducing new ones.

### Definition of done

For any behaviour change:

- Diffs are scoped to the task; no opportunistic refactors that span
  unrelated files.
- Tests added or extended (GoogleTest). Prefer extending an existing
  `test/test_<topic>.cpp` over duplicating logic.
- Public-API or behaviour changes are recorded in [NEWS.md](NEWS.md)
  under the appropriate `### …` heading with a PR number; breaking
  changes use the `DEPRECATED` macro where a soft migration is feasible.
- Doxygen `//!` comments on changed declarations are kept accurate
  (downstream `pyvinecopulib` extracts them via libclang to build its
  Python docstrings, so keep them clean and self-contained).
- Run the validation sequence: `clang-format-14`, a Debug build +
  `bin/test_all`, and — for anything touching the install/precompiled
  path — a `VINECOPULIB_PRECOMPILED=ON` release build.

## Coding conventions

- **Style: `BasedOnStyle: Mozilla`** (2-space indent, ~80-column),
  enforced by clang-format 14.
- **Declaration / implementation split.** Every public `foo.hpp`
  declares the interface and ends with
  `#include <vinecopulib/<module>/implementation/foo.ipp>`. The `.ipp`
  holds the definitions, each marked `inline` (this is what makes
  header-only correct and what the precompiled mode strips). See
  [bicop/class.hpp](include/vinecopulib/bicop/class.hpp) /
  [bicop/clayton.hpp](include/vinecopulib/bicop/clayton.hpp) for the
  pattern.
- **Naming.** Classes and enums are `PascalCase` (`Bicop`, `Vinecop`,
  `RVineStructure`, `FitControlsBicop`, `BicopFamily`); enum *values*,
  methods, and free functions are `snake_case` (`get_family`, `select`,
  `to_pseudo_obs`); private member variables carry a trailing underscore
  (`family_`, `rotation_`, `var_types_`). Tool headers are
  `tools_<domain>.hpp`.
- **Namespaces.** Everything public lives in `vinecopulib`. Infrastructure
  is grouped under `vinecopulib::tools_{stats,eigen,stl,select,…}`; family
  group lists under `vinecopulib::bicop_families`.
- **Doxygen.** Comments use `//! @brief` / `@details` / `@param` /
  `@return` / `@see` / `@literature`, with backticked inline code and
  fenced code blocks — **not** numpydoc. Keep declarations documented;
  see the per-family docs on the `BicopFamily` enum in
  [bicop/family.hpp](include/vinecopulib/bicop/family.hpp) for the house
  style.
- **C++14 + Eigen.** `Eigen::MatrixXd` / `VectorXd` are the workhorse
  types; prefer `.array()`, `.unaryExpr`, and the helpers in
  [misc/tools_eigen.hpp](include/vinecopulib/misc/tools_eigen.hpp)
  (`trim`, `remove_nans`, `unaryExpr_or_nan`, …) over hand-rolled loops.
  Templates are used sparingly; lambdas + `std::function` are the norm in
  implementations.

## Module boundaries

### `bicop/`

- **`Bicop`** ([class.hpp](include/vinecopulib/bicop/class.hpp)) is the
  public facade: a model is `(family, rotation ∈ {0,90,180,270},
  parameter matrix, var_types)`. It exposes `fit` / `select`, the
  evaluation methods (`pdf`, `cdf`, `hfunc1/2`, `hinv1/2`, `simulate`),
  fit statistics (`loglik`, `aic`, `bic`, `mbic`), the parameter↔τ maps
  (`parameters_to_tau`, `tau_to_parameters`), and JSON I/O (`to_json`,
  `to_file`, and constructors from a file / `nlohmann::json`).
- **`AbstractBicop`** ([abstract.hpp](include/vinecopulib/bicop/abstract.hpp))
  is the internal polymorphic base, reached only through `Bicop` (never
  use it directly). The hierarchy is
  `AbstractBicop → {ParBicop → {EllipticalBicop, ArchimedeanBicop,
  ExtremeValueBicop, …}, KernelBicop}`, with one class per family
  (`GaussianBicop`, `ClaytonBicop`, `TllBicop`, …). Subclasses implement
  the `*_raw` primitives (`pdf_raw`, `hfunc1_raw`, …); the base adds
  rotation, discrete handling (`pdf_c_d`, `pdf_d_d`), and numeric
  inverses.
- **`FitControlsBicop`**
  ([fit_controls.hpp](include/vinecopulib/bicop/fit_controls.hpp)) — the
  fit knobs: `family_set` (defaults to `bicop_families::all`),
  `parametric_method` (`"mle"`/`"itau"`), `nonparametric_method`
  (`"constant"`/`"linear"`/`"quadratic"`) + `mult` + `grid_size`,
  `selection_criterion` (`"loglik"`/`"aic"`/`"bic"`/`"mbic"`), `weights`,
  `psi0`, `preselect_families`, `allow_rotations`, `num_threads`.
- **`BicopFamily` + `bicop_families::`**
  ([family.hpp](include/vinecopulib/bicop/family.hpp)) — the 13 family
  values and 14 group lists (`all`, `parametric`, `nonparametric`,
  `one_par`, `two_par`, `three_par`, `elliptical`, `archimedean`,
  `extreme_value`, `bb`, `rotationless`, `lt`, `ut`, `itau`). The group
  lists are the canonical way to constrain a fit's search space.
- **`tools_select`** — bicop candidate generation and quick
  pre-selection used during fitting.

### `vinecop/`

- **`Vinecop`** ([class.hpp](include/vinecopulib/vinecop/class.hpp)) — a
  vine is `(RVineStructure, pair copulas, var_types)`. It exposes `select`
  (structure + families) and `fit` (refit families/params on a fixed
  structure), evaluation (`pdf`, `pdf_full`, `cdf` (Monte-Carlo),
  `simulate`, `rosenblatt`, `inverse_rosenblatt`), per-edge and bulk
  accessors (`get_pair_copula`, `get_all_families`, …), structure
  accessors (`get_matrix`, `get_struct_array`, `get_order`), fit
  statistics (`loglik`/`aic`/`bic`/`mbicv`), `truncate`, the asymptotic
  helpers (`scores`, `hessian`, `scores_cov`), and JSON I/O.
- **`RVineStructure` / `CVineStructure` / `DVineStructure`**
  ([rvine_structure.hpp](include/vinecopulib/vinecop/rvine_structure.hpp))
  — the tree structure as an R-vine matrix / triangular array, with
  truncation, the precomputed `needed_hfunc1/2` masks, validity checks,
  and `RVineStructure::simulate` for random structures. C- and D-vines are
  order-determined specialisations.
- **`FitControlsVinecop`**
  ([fit_controls.hpp](include/vinecopulib/vinecop/fit_controls.hpp))
  inherits `FitControlsBicop` and adds the structure knobs: `trunc_lvl`,
  `tree_criterion` (`"tau"`/`"rho"`/`"hoeffd"`/`"mcor"`), `threshold`,
  `tree_algorithm` (`"mst_prim"` default, `"mst_kruskal"`,
  `"random_weighted"`, `"random_unweighted"`), `select_trunc_lvl` /
  `select_threshold` / `select_families`, `show_trace`, and `seeds`.
- **`VinecopSelector`**
  ([tools_select.hpp](include/vinecopulib/vinecop/tools_select.hpp)) — the
  internal Dissmann selection engine: builds trees level-by-level
  (`select_all_trees`, `sparse_select_all_trees`), weighting edges by the
  chosen dependence criterion and choosing a spanning tree per the
  `tree_algorithm`.

### `misc/`

Infrastructure shared across both modules:

- **`tools_stats`** ([tools_stats.hpp](include/vinecopulib/misc/tools_stats.hpp))
  — `to_pseudo_obs`, `simulate_uniform` (dispatching to `ghalton` for
  `d ≤ 300` and `sobol` otherwise when `qrng=true`), normal/Student-t
  densities·CDFs·quantiles, bivariate normal/t CDFs, dependence-measure
  matrices (via `wdm`), and the `BoxCovering` machinery for discrete-ties
  latent-sample recovery.
- **`tools_optimization` + `tools_bobyqa`** — Powell's BOBYQA
  bound-constrained optimiser behind an `Optimizer` wrapper; the engine
  for parametric MLE.
- **`tools_interpolation`** — `InterpolationGrid`, the bilinear unit-square
  grid backing the `tll` family, with Sinkhorn margin normalisation and
  `integrate_1d` / `integrate_2d`.
- **`tools_serialization` + `nlohmann_json.hpp`** — JSON (de)serialisation
  of matrices, vectors, and triangular arrays; `nlohmann_json.hpp` is
  vendored verbatim (excluded from formatting and coverage).
- **`tools_thread` + `tools_interface` + `tools_batch`** — a `std::thread`
  `ThreadPool`, batched work partitioning, and the `INTERFACED_FROM_R`
  switch that swaps in RcppThread + R printing/interrupt handling when the
  library is compiled inside R.
- **`triangular_array.hpp`** — `TriangularArray<T>`, the ragged container
  that mirrors a (possibly truncated) vine.
- **`tools_eigen`, `tools_stl`, `tools_integration`, `tools_constants`,
  `tools_optional`** — Eigen/STL helpers, 1-d numerical integration
  (Boost.odeint), constants, and a C++14/17 `optional` shim.

## Public API

The Doxygen [website](https://vinecopulib.github.io/vinecopulib/) and the
`docs/*.dox` overviews are the source of truth for signatures and
parameter docs. Quick orientation:

- **Core classes** — `Bicop`, `Vinecop`, `RVineStructure`,
  `CVineStructure`, `DVineStructure`, `FitControlsBicop`,
  `FitControlsVinecop`.
- **Families** — the `BicopFamily` enum and the `bicop_families::` group
  lists.
- **Free functions** (`vinecopulib::tools_stats`) — `to_pseudo_obs`,
  `simulate_uniform`, `sobol`, `ghalton`, plus the distribution helpers
  (`dnorm`/`pnorm`/`qnorm`, `dt`/`pt`/`qt`, …).

All of the above are reachable through the single umbrella include
[`vinecopulib.hpp`](include/vinecopulib.hpp), or via the narrower module
headers (`vinecopulib/bicop/class.hpp`, etc.).

## Tests

GoogleTest, under [test/](test/). The layout is split: thin
`test/test_<topic>.cpp` translation units register the suites, while the
actual logic lives in `test/src_test/*.cpp` with fixtures in
`test/src_test/include/*.hpp`. The executable list (12 binaries, including
the aggregate `test_all`) is in
[cmake/buildTargets.cmake](cmake/buildTargets.cmake); each links
`gtest vinecopulib … src_test`.

Conventions:

- **Run** `bin/test_all` (or `ctest`) from the build directory after
  `make`.
- **Coverage** is collected on the Ubuntu/GNU Debug CI build
  (`-DCODE_COVERAGE=ON`, target `vinecopulib_coverage`) and uploaded to
  Codecov; [codecov.yml](codecov.yml) ignores `test/` and the vendored
  `nlohmann_json.hpp`.
- **R parity tests** are optional: with Rscript available, the
  `cmake/templates/*.R` scripts cross-check parametric bicop / vinecop
  results against the VineCopula R package.
- **Round-trip / parity invariants to preserve** when touching numerics:
  JSON `to_file` → constructor round-trips for `Bicop`, `Vinecop`, and
  `RVineStructure`; copula identities (h-function ↔ inverse, Rosenblatt ↔
  inverse-Rosenblatt, parameter ↔ τ); and agreement with the downstream R
  / Python wrappers, which test against this library's output.

## Extension points

- **New copula family.** Add `bicop/<name>.hpp` (declaration, subclassing
  the right base — `EllipticalBicop`, `ArchimedeanBicop`,
  `ExtremeValueBicop`, or `KernelBicop`) and
  `bicop/implementation/<name>.ipp` (the `inline` definitions of the
  `*_raw` primitives and the τ↔parameter maps). Register the tag in the
  `BicopFamily` enum and the relevant `bicop_families::` group lists in
  [bicop/family.hpp](include/vinecopulib/bicop/family.hpp), wire it into
  candidate generation in `bicop/tools_select`, add a `test/test_*.cpp`
  case, and document it in [NEWS.md](NEWS.md). Downstream wrappers then
  bump their pin.
- **New structure-selection algorithm.** Extend the `tree_algorithm`
  dispatch in `vinecop/tools_select` and accept the new value in
  `FitControlsVinecop`.
- **New dependence / tree criterion.** Add it alongside the wdm-backed
  criteria in `tools_stats` and accept it as a `tree_criterion` value.
- **New quasi-random method.** Slot it in beside `sobol` / `ghalton` in
  `tools_stats`, behind `simulate_uniform`.

## Maintaining this file

If a coding agent or reviewer keeps repeating the same correction — about
a convention this repo enforces — update `AGENTS.md` rather than relying
on tribal knowledge. Do not add ephemeral, user-specific, or
machine-local preferences here. [NEWS.md](NEWS.md) is the place for
release-by-release context; this file is for invariants.
