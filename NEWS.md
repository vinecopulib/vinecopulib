## vinecopulib 1.0.0 (unreleased)

The first stable release. It collects a large amount of work: analytic
derivatives and asymptotic-inference tooling, conditional simulation, discrete
and per-observation-parameter surfaces, substantially faster evaluation and
fitting, a modernized C++17 build, and documentation whose every example is
compiled and run in CI.

Several changes require action. See
[docs/migrating-to-1.0.md](docs/migrating-to-1.0.md) — in particular the
`FitControlsBicop` / `FitControlsVinecop` constructor argument that was inserted
in the middle, which compiles silently against old positional calls and does the
wrong thing.

### BREAKING API CHANGES

* `nonparametric_grid_size` is a new `FitControlsBicop` and `FitControlsVinecop`
  constructor argument, inserted after `nonparametric_mult`. Positional calls
  written against earlier releases still compile and silently rebind every later
  argument; see the migration guide (#654)

* Remove `FitControlsVinecop::get_truncation_level`,
  `get_select_truncation_level`, `set_truncation_level` and
  `set_select_truncation_level`, deprecated since 0.3.1. Use the `trunc_lvl`
  spellings. The four had no definitions, so calling one was a link error (#718)

* Remove `Vinecop::select_all` and `Vinecop::select_families`, deprecated since
  0.3.1. Use `select`, which selects the structure unless the object already
  carries one (#718)

* Rename `Vinecop::hessian_avg` to `hessian` and the former `hessian` to
  `hessian_full`, matching the `pdf` / `pdf_full` pair (#679)

* Replace `tools_serialization::triangular_array_to_json` and
  `json_to_triangular_array` with `TriangularArray::to_json()`, a JSON
  constructor, and `to_list()`. The on-disk format is unchanged (#680)

* Delete `misc/tools_bobyqa.hpp`. No public signature changes, but the header is
  gone for anyone including it directly (#685)

* `Vinecop`'s data constructor no longer defaults its `matrix` argument, which
  made `Vinecop(data)` ambiguous (#709)

* The public umbrella header no longer disables Boost's concept assertions.
  Consumers reaching Boost through `vinecopulib.hpp` may see concept diagnostics
  that were previously suppressed in their own translation units (#714)

* `VINECOPULIB_VERSION` is now a plain decimal literal. It was written with
  leading zeros, making it an octal constant, so any consumer comparing against
  it was comparing a different number than it looked like

* Require C++17, CMake 3.14 and Boost 1.75. The installed target is
  `vinecopulib::vinecopulib`; the unqualified name remains as an alias (#711)

* Rename the CMake options to a `VINECOPULIB_` prefix and remove the unprefixed
  `OPT_ASAN`, `CODE_COVERAGE`, `STRICT_COMPILER` and `BUILD_DOC`. Sanitizers now
  default to off, and `-march=native` is behind `VINECOPULIB_NATIVE_ARCH` so the
  default release build is redistributable (#711)

* Remove the `BUILD_SHARED_LIBS` option; use `VINECOPULIB_PRECOMPILED` (#662)

* Eigen is found in CONFIG mode, so a module-mode-only Eigen install no longer
  configures (#658)

### BEHAVIOR CHANGES

* Every discrete or mixed fit with a nonparametric (`tll`) pair copula changes:
  the fitted grid, density, `loglik`, `aic`, `bic` and `mbicv` all move, so
  family and structure selection can differ. `tll` is the default family, so
  this is the default path for discrete data. Continuous fits are unaffected;
  see BUG FIXES (#739)

* TLL fits change near the boundary of the unit square: the endpoint clamping of
  the interpolation grid was removed, so every `tll` fit differs slightly from
  0.7.3 (#654)

* Kendall's tau of `bb6`, `bb7`, `bb8` and `tawn` changes. Four numerical
  defects in `parameters_to_tau` are fixed, the worst of which returned about
  1e-11 where the true value was 0.33; see BUG FIXES (#713)

* Maximum-likelihood estimates shift in the low digits: BOBYQA is replaced by
  Brent and BFGS (#685)

* `file_to_json` and `json_to_file` now throw on I/O failure instead of failing
  silently, and a `.cbor` extension selects the binary encoding (#684)

* The R-vine structure is stored with the conditioned variable on the diagonal,
  so `get_matrix()`, `get_order()`, `get_struct_array()` and edge orientation
  differ for the same model. Densities and log-likelihoods do not (#702)

### NEW FEATURES

* Add scores, gradient, Hessian and score covariance to `Bicop` and `Vinecop`,
  with `scores_full` and `hessian_full` exposing the per-edge intermediates for
  callers that need several of them (#645, #683, #699)

* Add analytic derivatives of the density and h-functions for every parametric
  family, including the BB families, Tawn and TLL; other families fall back to
  finite differences (#683, #687, #694)

* Add conditional simulation: `Vinecop::simulate_conditional`,
  `FitControlsVinecop::set_conditioning_set` for conditioning-aware structure
  selection, and `reorient` to relabel an already-fitted vine (#696, #697)

* Add conditioning-set overloads to `rosenblatt` and `inverse_rosenblatt`, which
  reorient the vine through non-owning views rather than copying pair copulas
  (#715)

* Accept and document both the expanded `2d` and compact `d + k` data layouts
  consistently across bivariate and vine copula operations, including the
  Rosenblatt transforms and conditional simulation (#729)

* Evaluate copulas with one parameter set per observation: `pdf`, `cdf`,
  h-functions, `loglik`, `simulate` and the derivative methods on `Bicop`, and
  the corresponding `Vinecop` surfaces (#675, #699, #719)

* Add `Vinecop::pdf_full`, returning the density together with the per-edge
  densities and h-functions computed along the way (#669, #699)

* Add `Bicop::get_taildep` and `Bicop::get_beta` for tail-dependence
  coefficients and Blomqvist's beta (#682)

* Allow a user-supplied `tree_criterion` function for structure selection, via
  `set_tree_criterion_function()` and `tree_criterion = "custom"` (#674)

* Add `RVineTrees`, an explicit edge-list view of a structure, with
  `Vinecop::get_trees()` and an `RVineStructure` constructor from it (#698)

* Add a `nonparametric_grid_size` control for the TLL interpolation grid (#654)

* Persist models as CBOR as well as JSON, selected by the file extension (#684)

### PERFORMANCE

* Speed up the bivariate evaluation engine and tighten allocation in the
  derivative cascade (#681)

* Speed up shared Eigen, thread and integration primitives (#689)

* Speed up `tools_stats`: SIMD `qnorm`, the bivariate normal and t kernels,
  pseudo-observations and `BoxCovering` (#690)

* Speed up TLL fitting and evaluation through fused interpolation (#691)

* Speed up vine evaluation and structure selection (#692)

* Compute the Student t score's degrees-of-freedom terms once per call (#693)

* Speed up the `RVineTrees` peel; results are byte-identical (#701)

* Exit structure selection early when the graph is already a tree (#661)

### BUG FIXES

* Validate the pair-copula store in the `RVineTrees` constructor rather than
  indexing it blindly: an empty store means independence on every edge, and one
  that does not cover the structure array is an error instead of a read past its
  end. `Vinecop::fit` no longer returns silently when the model has no pair
  copulas to fit; a vine whose structure is truncated at zero now records the
  independence fit (`loglik` 0) rather than staying unfitted (#748)

* Read an empty pair-copula store as independence in the paths that still
  assumed a full one, completing (#729). `Vinecop::reorient` and `get_trees`
  indexed the store out of bounds, `truncate` gave it tree slots holding no pair
  copula at all, which `pdf` then dereferenced, and deserializing a vine written
  without pair copulas threw a JSON type error instead of round-tripping (#743)

* Compute a real discrete density for kernel pair copulas. `KernelBicop`
  overrode `pdf`, `hfunc1` and `hfunc2` to return the continuous quantity at the
  atom midpoint, so the difference quotients in `AbstractBicop::pdf_c_d` /
  `pdf_d_d` were unreachable. The midpoint rule violates
  `sum_atoms c * (u1 - u1^-) = 1` by up to 10% without converging, reached 40%
  per cell, and shifted a three-dimensional discrete vine's `loglik` by 6.5
  (#739)

* Use the latent sample when fitting a `tll` pair copula to discrete data.
  `TllBicop::fit` assigned `find_latent_sample`'s result to `psobs` but fitted on
  `z_data`, which was already built from the jittered ranks, so declaring a
  variable discrete had no effect on the fitted grid at all (#736, #739)

* Report `select_threshold` rather than `select_trunc_lvl` on the
  "Select threshold" line of `FitControlsVinecop::str()` (#735)

* Remove the declaration of `tools_stats::dependence_matrix`, which had no
  definition anywhere in the library, so any call to it failed to link (#735)

* Treat omitted pair-copulas as implicit independence during evaluation of a
  `Vinecop` constructed from a full structure, preventing conditional
  simulation from accessing an empty pair-copula store (#729)

* Fix four numerical defects in `parameters_to_tau`: Tawn's `A''` evaluated to
  `inf * 0` for small psi and cancelled to a negative value elsewhere, its
  integrand peaks away from 1/2 for asymmetric parameters, and the BB7 and BB8
  integrands cancelled to exactly zero near their boundaries (#713)

* Fix the TLL CDF and the two-dimensional integration it relies on (#659, #667)

* Make `inverse_rosenblatt` thread-safe (#712)

* Evaluate a custom tree criterion serially, since a user-supplied callable is
  not required to be thread-safe (#722)

* Fix Wilson's algorithm for random spanning trees (#653)

* Fix a typo in `pdf_d_d` and adjust the threshold for analytical derivatives
  (#664)

* Fix starting parameters for discrete models (#677)

* Handle edge-case parameter shapes in per-row evaluation (#700)

* Include the thread header explicitly where it was relied on transitively
  (#676)

### BUILD SYSTEM AND DEPENDENCIES

* Modernize the CMake project: `GNUInstallDirs`, a namespaced export,
  `target_compile_features`, an architecture-independent version file,
  `CMakePresets.json`, and appended rather than assigned compiler flags, so
  distribution hardening flags and user `CXXFLAGS` survive (#711)

* Replace the byte-offset arithmetic in the `.ipp` to `.cpp` code generation
  with a regex, and add real dependency edges, so an edited `.ipp` can no longer
  leave a stale generated source behind (#711)

* Reduce Boost to Graph, Math and Random. Odeint, Bimap, Assign, Unordered and
  Concept are gone; `std::optional` replaces the `boost::optional` shim (#711,
  #713, #714)

* Build one test binary plus one per area, make the R parity tests
  parallel-safe, and make R optional behind `VINECOPULIB_R_PARITY_TESTS` (#717)

* Rebuild the CI, release and analysis workflows, with version-consistency,
  header-only-install, benchmark, documentation and clang-tidy jobs, plus
  release automation and a weekly check for a release that was written up but
  never tagged

* Safeguard the vendored `nlohmann_json.hpp` against a newer Xcode (#678)

* Fix several Windows CI issues (#647, #655) and drop the deprecated macos-13
  runner (#663)

### DOCUMENTATION AND TOOLING

* Rewrite the user documentation. Every example is now a `\snippet` of a source
  under `docs/snippets/`, compiled and run in CI, and new pages cover discrete
  variables, conditional simulation, and scores and Hessians. Doxygen runs with
  `WARN_AS_ERROR`

* Publish the website from CI. Every pull request renders the m.css site and
  uploads it as an artifact, and a `v*` tag commits it to `gh-pages`, so the
  published documentation can no longer lag the code (#731)

* Document where the Kendall's tau maps apply. `tau_to_parameters` is available
  only for the one-parameter families in `bicop_families::itau`, notably not for
  `student` (#723)

* Add a benchmark suite, a numerical parity harness and golden-value tests, with
  `benchmarks/README.md` and `scripts/README.md` describing both (#688)

* Add `CONTRIBUTING.md`, `CITATION.cff`, `SECURITY.md`, `CODE_OF_CONDUCT.md`,
  `docs/releasing.md`, and issue and pull-request templates

* Clear the clang-tidy backlog and widen the enabled checks (#710)

* Add `AGENTS.md`, a normative engineering spec for the repository (#672)

* Document per-family parameters, rotations and tail dependence on the
  `BicopFamily` enum, and the property getters and setters that the Python
  bindings extract (#668, #670)

* Seed the RNG-dependent unit tests to remove flakiness (#686)

* Decompose `VinecopSelector` for readability (#695)

* Enforce the clang-format style in CI (#646, #649)

* Fix documentation typos and the mBIC formula (#660, #665, #703, #716)

## vinecopulib 0.7.3 (April 23, 2025)

### BREAKING API CHANGES

* The `mst_algorithm` option to `FitControlsVinecop` has been renamed to `tree_algorithm` to
  allow for alternative spanning tree algorithms (#637).
* `tree_algorithm`'s default value is now `"mst_prim"` instead of `"prim"`, and `"mst_kruskal"`
  replaces `"kruskal"` (#637).
* The CMake option `VINECOPULIB_BUILD_SHARED_LIBS` has been changed to `VINECOPULIB_PRECOMPILED`
  to better reflect its purpose (#641).

### NEW FEATURES

* Allow for random spanning trees as alternatives to the MST-based structure selection using
  `tree_algorithm` in `FitControlsVinecop` with `"random_weighted"` or `"random_unweighted"`
  (#637).

### BUG FIXES

* Decouple edge insertion from criterion computation in `VinecopSelector` to fix randomness
  issues in structure selection when using multiple threads (#640)

## vinecopulib 0.7.2 (March 7, 2025)

### BUG FIXES

* More build system updates by @jschueller (#633)

* Fix deprecation warning in json header (#634)

* Fix TLL speed issues related to FFT (#635)

## vinecopulib 0.7.1 (January 15, 2025)

### NEW FEATURES

* Add `allow_rotation` option to `FitControlsBicop` and `FitControlsVinecop`
  to allow for the rotation of the pair copulas (#628).

* Add a `FitControlsConfig` struct to create flexible and yet safe constructors
  for `FitControlsBicop` and `FitControlsVinecop` (#629).

### BUG FIXES

* Restrict parameter range for fitting Tawn copulas; fix handling of their 
  shape/argument order (#620).

* Compute and save loglik/nobs in `Vinecop::fit()` (#623)

* Disable unwanted compiler output related to BOOST_CONCEPT checks (#624)


## vinecopulib 0.7.0 (January 2, 2025)

### NEW FEATURES

* Use analytical derivatives in discrete pdf/hfuncs (#572)
* Allow for alternative for `"prim"` vs `"kruskal"` in MST-based model selection (#577)
* Improve the dependencies install script to use it in other projects (#576)
* Add tawn copula (#579)
* Improve doc (#580, #585, #607)
* Allow for the discrete Rosenblatt transform (#581)
* Add `Vinecop::fit()` (#584)
* Improve `Bicop::str()` (#588) and `Vinecop::str()` (#589)
* Properly handle discrete variables for the TLL family (#597)
* Weighted pseudo-observations (#602)
* Cross-platform random numbers and add seeds options to `to_pseudo_obs` (#603)
* Improve performance by
    * aligning with the `R` defaults (e.g., `BOOST_NO_AUTO_PTR`, `BOOST_ALLOW_DEPRECATED_HEADERS`, `BOOST_MATH_PROMOTE_DOUBLE_POLICY=false`, `std::string nonparametric_method = "constant"` for the TLL instead of `"quadratic"`, `-O3 -march=native` compiler flags) and add benchmarking example (#592, #611, #613),
    * using `Eigen` element-wise operations instead of `boost` whenever possible (#598, #612),
    * using binary search in the TLL for `get_indices` (#613).

### BUG FIXES

* Improve stability in BB7 PDF (#573)
* Revamped CI/CD pipeline, tests discoverable by CTest, boost version on windows (66cf8b0)
* Fix ASAN issues (#583)
* Fix interface includes and other CMake issue (#586, #599, #601, #608, by @jschueller)

## vinecopulib 0.6.3 (Februrary 20, 2023)

### BUG FIXES

* Replace bitwise ops on boolean operands (#563)

* Handle `NaN`s in `to_pseudo_obs()` (#566)

* Replace calls to `sprintf` in boost libraries (#565)

## vinecopulib 0.6.2 (August 24, 2022)

* updated cmake setup (#549)

* improved documentation (missing data and Rosenblatt transforms)

* better parallelization when there is a small number of edges (#555)

## vinecopulib 0.6.1 (July 13, 2021)

### BREAKING API CHANGES

* refactored serialization using `nlohmann::json` instead of `boost::property_tree` (#539):

    * `to_json` and `to_ptree` methods of both `Bicop` and `Vinecop` objects
    are respectively renamed as `to_file` and `to_json`

    * internal structure of the serialized objects is changed (i.e., not
    possible to read the old files with the new functions)

### BUG FIXES

* use `num_threads` in recursive calls to the inverse Rosenblatt (#535)

* force TLL to be nonnegative (#532)

* fix number of parameters for TLL (#530)

* fix bugs in `Vinecop` print methods (#544)

## vinecopulib 0.5.5 (November 23, 2020)

### BUG FIXES

* fix little bug in copula selection based on mBIC (#527)

* stabilize BB7 copula pdf (#526)

* fix threshold selection for (near-)independent data (#523)

* fix `Vinecop::select()` for 1-dimensional models with discrete variables (#521)

* fix user-visible variable types messed up in `Bicop::flip()` (#519)

## vinecopulib 0.5.4 (September 30, 2020)

### BUG FIXES

* fix uninitialized number of parameters for TLL family (#515)

* fix Kendall's tau of Frank copula for par \<= 3 (#513)

* fix `Vinecop::pdf()` when discrete variables are present (#514)

## vinecopulib 0.5.3 (August 11, 2020)

### NEW FEATURES

* allow 1-dimensional models (#499)

* make AIC default selection criterion (#502)

### BUG FIXES

* make Bicop/Vinecop objects indepent of copied-from-objects (#503)

* enforce parameters bounds in tau_to_parameters for Archimedean families (#507)

## vinecopulib 0.5.2 (May 7, 2020)

### NEW FEATURES

* `str()` methods for `FitControlsBicop` and
  `FitControlsVinecop` (#494)

### BUG FIXES

* fix documentation (#482, #493)

* fix bug in `RVineStructure::simulate()` (#492)

* fix tll family with comonotonic data (#491)

* fix weights handling in family preselection (#490)

* fix archimidean h-functions near independence (#488)

* safety net for NA structure weights (#487)

* fix code qix quality issues (#486)

* fix `Vinecop::str()` (#484)

## vinecopulib 0.5.1 (November 25, 2019)

### BUG FIXES

* fix out of range bug for weighted TLL influence when sample size is small
  (#479).

## vinecopulib 0.5.0 (November 25, 2019)

### NEW FEATURES

* modelling discrete variables with bivariate or vine copulas. (#454)

* selection of partially specified R-vine structures. (#457)

* convenience classes `DVineStructure()`/`CVineStructure()` for D- and
  C-vine structures. (#464)

* new criterion for tree selection: `"joe"` corresponds to -log(1-r^2),
  where r is the pairwise partial correlation. (#426)

* random sampling of R-vine structures. (#429)

* (de)serialization methods `to/from_ptree/json()` for `RVineStructure`
  objects. (#435)

* some improvements in memory efficiency. (#460)

### BREAKING API CHANGES

* `Vinecop` ctors: interchanged order of `structure` and `pair_copulas`,
  removed unpopular versions, new argument `var_types`. (#465)

* removed `tools_stats::simulate_uniform(n, d, seeds)` to avoid implicit
  conversion. (#430)

* `calculate_npars()` becomes `get_npars()`. (#431)

* by default, `RVineStructure::get_struct_array()` and
  `RVineStructure::struct_array()` objects return the structure array in the
  original order (as opposed to natural order). An additional argument is
  available to ask for the old behavior. (#437, #439)

* removed `TriangularArray<T>::operator[]` to access columns.
  `TriangularArray`s are now stored row-major and provide a new constructor
  `TriangularArray<T>(std::vector<std::vector<T>> rows)`. (#433)

### BUG FIXES AND OTHER IMPROVEMENTS

* better support for 0-truncated structures. (#462)

* prevent 0/0 in normalization of `BicopFamily::tll` fits. (#463)

* use `std::shufle()` instead of `std::random_shuffle()` to remain\
  C+17-compliant. (#451)

* ensure consistency in manually created `Bicop(BicopFamily::tll, ...)`,
  with fitted versions. (#446)

* fixed order of ranks in `to_pseudo_obs(.., "first")`. (#469)

* safer computation of multivariate normal cdf. (#475)

## vinecopulib 0.3.2 (July 3, 2019)

### NEW FEATURES

* new `Vinecop::str()` method (#420)

* enhanced extensibility of `RVineStructure`, `Vinecop`, and
  `VinecopSelector` classes and  (#419)

### BUG FIXES

* fix interval adjustment for Brent parameter optimization (#414)

* clean up includes to improve build times (#412)

* better printing for tll family (#415)

* fix batches when `num_threads = 0` (#418)

## vinecopulib 0.3.1 (April 5, 2019)

### NEW FEATURES

* refactoring for enhanced extensibility of the class `Vinecop` (#407)

* add an `str` method to `RVineStructure` (#406)

* simplify algorithms by reversing definition of natural order (#387)

* improve selection of truncation level (#373)

* add truncate methods for `TriangularArray`, `RVineStructure` and `Vinecop`
  (#372)

### DEPRECATED FUNCTIONS

* `FitControlsVinecop` methods like `get_truncation_level()` are now
  deprecated in favor of the shorter `get_trunc_lvl()` versions.

### BUG FIXES

* fix triangular array print method (#405)

* fix potential nan when using Kendall's tau inversion to fit (#403)

* stabilize clayton pdf close to independence (#402)

* fix warning-generating typos in the tests (#391)

* fix deprecated warnings from gtest (#390)

* initializer list constructor for `VinecopSelector` (#384)

* fix baseline criterion in `Bicop::select()` (#382)

* use the `trunc_lvl` of `vine_struct` in `select_families` (#380)

* fix incompatible size (#378)

* remove definitions of `pairwise_rho` and `pairwise_hoeffd` (#375)

* fix windows warnings (#371)

* use `std::log()` (#370)

## vinecopulib 0.3.0 (August 9, 2018)

### NEW FEATURES

* new function `Vinecop::rosenblatt()` for computing the Rosenblatt
  transformation (#367).

* faster algorithms for nonparametric copulas based on bilinear interpolation
  (#357).

* refactor vine structures and related algorithms with triangular arrays to
  improve efficience of truncated models (#347, #354, #365).

* improve random number generation: allow for seeds and quasi-random
  numbers (#342, #356).

* improved parallelization for fitting vine copula models (#338, #344).

* parallelized versions of many algorithms including pdf, cdf and simulation
  (#339, #363).

* allow weights for observations (#336).

### BUG FIXES

* fix cdf of StudentBicop (#353).

* improve numerical stability (#345, #350).

* fix gcc-8 warning (#340).

## vinecopulib 0.2.8 (May 4, 2018)

### NEW FEATURES

* new getters for Kendall's tau (#319).

* log-likelihood as a new data member for AbstractBicop and Vinecop
  (`loglik_`), initialized as NAN and set by fitting
  routines like `fit` and `select` (#327).

* new getters for the log-likelihood (#327).

### BUG FIXES

* fix truncation of pdf values (#320).

* use increased search interval for parameter estimation when initial fit is
  unreasonable (#322).

* ensure that boundaries are respected for Joe's `hinv` methods (#323).

* improve numerical stability by more restrictive parameter bounds for Joe
  and BB7 copulas (#324, #325).

* bandwidth adjustment for family `"tll"` (#326).

## vinecopulib 0.2.7 (March 1, 2018)

### NEW FEATURES

* new criterion for tree selection `"mcor"` (#309).

### BUG FIXES

* fix bandwidth scaling for family `"tll"`(#309).

## vinecopulib 0.2.6 (February 22, 2018)

### NEW FEATURES

* add checks for data in (0, 1) (#305).

* add scaling by `mcor` to estimate the bandwidth of `TllBicop` (#302) .

* add mBICV to select the truncation level and threshold (#304).

### BUG FIXES

* improve Windows build (#301, #302).

* fix `hoeffd` in `calculate_criterion` (#297).

## vinecopulib 0.2.5 (January 14, 2018)

### NEW FEATURES

* speed up vine copula algorithms by pre-computing information related to
  the vine structure (#292).

* the selected threshold parameter can be returned from the `Vinecop`
  object (#290).

### BUG FIXES

* fix storage order of pair copulas when structure is fixed (#289).

* fixed selection algorithm for threshold and truncation level (#290, #294,
  \#295).

## vinecopulib 0.2.4 (December 29, 2017)

### BUG FIXES

* adapt Vinecop's `simulate()` and `pdf()` to truncated vines (#279)

* make bb8 lower bound ensure feasible computations in `parameters_to_tau()`
  (#278 and #280)

* default initialize Rcout (#277).

## vinecopulib 0.2.3 (November 18, 2017)

### NEW FEATURES

* faster implementation of Archimedean pdfs (#274).

### BUG FIXES

* add safeguards for estimation of `Bicop`/`Vinecop` objects with
  insufficient data (#273).

* fix segfault issue in completing a truncated vine fit (#272).

* make `par_method = "itau"` respect the parameter bounds (#271).

## vinecopulib 0.2.2 (November 9, 2017)

### NEW FEATURES

* allow `"loglik"` as selection criterion (#267).

### BUG FIXES

* make interpolation grid symmetric around (0.5, 0.5) again (#268).

## vinecopulib 0.2.1 (November 7, 2017)

### NEW FEATURES

* faster vine copula estimation and selection by parallelizing further
  sub-routines (#259).

* enhanced cross-platform compatibility and addition of a `STRICT_COMPILER`
  option for gcc (#261).

* increased precision of maximum-likelihood estimators (#264).

### BUG FIXES

* made arguments of pairwise dependence measures consistent (#258).

* fixed `itau` estimation method for Frank copulas (only allowed for positive
  parameters) (#263).

## vinecopulib 0.2.0 (October 30, 2017)

### LIBRARY TYPE

* library is now header only by default (#246), with an option to compile it
  as a shared library (#249).

### DEPENDENCIES

* removed dependency on `NLopt` (#239).

### NEW FEATURES

* NA handling (#237, #238).

* parallelized selection/estimation of (pair-) copulas (#240).

* efficient storage and fitting of truncated vines (#248).

* Brent line search for (profile-) maximum-likelihood estimation of one-
  parameter families (#255).

* more restrictive parameter bounds for Archimedean families, ensuring their
  numerical stability (#256).

### BUG FIXES

* error thrown whenever `Vinecop` or `Bicop` constructors are called with
  datasets containing a single row (#251).

* ensure `const` correctness of `Vinecop` and `Bicop` member functions (#225).

* made order of inverse Rosenblatt consistent for d = 2 and d > 2 (#232).

* fixed bug in interpolation near upper right corner (#233).

* interpolation grid is now symmetric around (0.5, 0.5) (#234).

* stabilized quadratic tll estimator near zero (#235).

* stabilized Archimedean pdfs (#256).

## vinecopulib 0.1.0 (August 23, 2017)

### NEW FEATURES

* read/write `Bicop` and `Vinecop` objects (#205) using
  `boost::property_tree::ptree` with `to_ptree()`, `to_json()`, and
  constructors taking `const char   *filename` or a
  `boost::property_tree::ptree` for both classes.

* sparse selection of vine copulas (#206) using new data members in
  `FitControlsVinecop`:
      _`bool select_truncation_level` whether the truncation is selected
       automatically.
      _ `bool select_threshold` whether the threshold parameter is selected
       automatically.
      \* `double threshold` sets a fixed threshold parameter.

* local likelihood estimators (#216) have been implemented by refactoring the
  `tll0` family into a more general `tll` family, where approximations of
  degrees zero, one and two can be fitted by setting the new
  `std::string nonparametric_method` data member of `FitControlsBicop`
  respectively as `constant`, `linear` and `quadratic` (default).

* Kendall's tau (#211) and normalization (#215) for kernel estimators.

* support for clang compiler on linux (#201, #202, #203).

* option to omit R-vine matrix check in `Vinecop` constructors (#198).

### BUG FIXES

* replacing throw `std::string` with throw `std::runtime_error` in
  `tools_opimization.cpp` (#204).

* ensure valid starting parameters in `Bicop::fit()` (#209, #210).

* fix appveyor and travis problems (#208, #212, #213).

## vinecopulib 0.0.3 (June 7, 2017)

### NEW FEATURES

* new functions `Bicop::cdf()` and `Vinecop::cdf()` for evaluating the
  cumulative distribution function of bivariate and vine copulas (#177, #189).

* the constructor of the `RVineMatrix` class now checks whether it is
  provided with a valid R-vine matrix (#192).

* extended documentation to build the library under Windows (#188).

* extended continuous integration tests for Windows (#150, #169).

### BUG FIXES

* vinecopulib.dll is installed to `lib/` instead of `bin/` (#149).

* more pleasing and portable formatting of error messages (#147, #156, #159,
  \#165).

* fixed bugs in `Bicop::select()` caused by `0`s and `1`s or unsufficient
  data (#173, #180).

* fixed compatibility issue with CMake 3.8 (#167).

* fixed uninitialized memory issues on Windows (#169).

## vinecopulib 0.0.2 (March 31, 2017)

### MAJOR CHANGES

* all `tools_xxx` namespaces are no sub-namespaces of `vinecopulib` (#130).

* header files are encapsulate in an addtional `vinecopulib/` folder, i.e.,
  `include/vinecopulib/subdir/file.hpp` (#126).

* removed abitility to extract the git revision (#124).

* new header `misc/tools_interface.hpp` where interface-specific behavior
  can be defined (for example, a custom version of `std::cout`) (#136).

### BUG FIXES

* fix `mat.array() = 0` error on some compilers (#131).

* add missing `<exception>` header (caused errors on not fully C++11
  compliant compilers, #139).

## vinecopulib 0.0.1 (March 29, 2017)

Initial release.
