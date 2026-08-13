# Migrating to vinecopulib 1.0.0 {#migrating-to-1-0}

@tableofcontents

## Do I need to change anything?

Work through this list; if none of it applies, a recompile is all that is
needed.

- [ ] Do you construct `FitControlsBicop` or `FitControlsVinecop` with
      **positional** arguments? → [Constructor arguments](#migrate-ctor).
      **This one compiles silently and does the wrong thing.**
- [ ] Do you call `hessian_avg`, `select_all`, `select_families`, or any
      `*_truncation_level` setter or getter? →
      [Renamed and removed API](#migrate-removed).
- [ ] Do you include `misc/tools_bobyqa.hpp` or call
      `tools_serialization::triangular_array_to_json`? →
      [Renamed and removed API](#migrate-removed).
- [ ] Do you compare against `VINECOPULIB_VERSION`? →
      [Renamed and removed API](#migrate-removed).
- [ ] Do you have regression tests pinning fitted parameters, `tll` densities,
      or Kendall's tau of the BB and Tawn families? →
      [Numerical differences](#migrate-numbers).
- [ ] Do you read `get_matrix()`, `get_order()` or `get_struct_array()`, or
      index pair copulas by position? →
      [R-vine matrix representation](#migrate-structure).
- [ ] Do you build with CMake, or reach Boost through `vinecopulib.hpp`? →
      [Build, CMake, and C++17](#migrate-build).

@section migrate-ctor Constructor arguments

`nonparametric_grid_size` was inserted into both fit-control constructors, after
`nonparametric_mult`. Positional calls written for 0.7.3 still compile, because
the neighbouring parameters have compatible types, and every argument after the
insertion point silently shifts by one.

Before, where the sixth argument was the selection criterion:

```cpp
// 0.7.3
FitControlsBicop controls(bicop_families::all, "mle", "constant", 1.0, "bic");
```

In 1.0.0 that fifth argument is `nonparametric_grid_size`, and `"bic"` no longer
lands where it did. Pass the arguments by name through the setters, which is
immune to further insertions:

```cpp
// 1.0.0
FitControlsBicop controls;
controls.set_family_set(bicop_families::all);
controls.set_parametric_method("mle");
controls.set_nonparametric_method("constant");
controls.set_nonparametric_mult(1.0);
controls.set_selection_criterion("bic");
```

`FitControlsConfig` offers the same by-name construction in one expression.

Note also that `selection_criterion` **defaults to `"aic"`**, not `"bic"`; the
0.7.3 documentation said otherwise, but the code has always defaulted to `"aic"`.

@section migrate-removed Renamed and removed API

| removed | replacement |
| --- | --- |
| `Vinecop::hessian_avg` | `Vinecop::hessian` |
| `Vinecop::hessian` (per-observation) | `Vinecop::hessian_full` |
| `Vinecop::select_all` | `Vinecop::select` |
| `Vinecop::select_families` | `Vinecop::select` |
| `FitControlsVinecop::{get,set}_truncation_level` | `{get,set}_trunc_lvl` |
| `FitControlsVinecop::{get,set}_select_truncation_level` | `{get,set}_select_trunc_lvl` |
| `tools_serialization::triangular_array_to_json` | `TriangularArray::to_json()` |
| `tools_serialization::json_to_triangular_array` | the `TriangularArray` JSON constructor |
| `misc/tools_bobyqa.hpp` | none; the optimizer is internal |

The `*_truncation_level` members were declared but never defined, so a call was
already a link error rather than a deprecation warning.

`Vinecop`'s data constructor no longer defaults its `matrix` argument, since that
made `Vinecop(data)` ambiguous. Pass a structure explicitly, or use the
`RVineStructure` overload.

**`VINECOPULIB_VERSION` was octal.** It was written `000703`, which C++ reads as
octal, so
```cpp
#if VINECOPULIB_VERSION >= 703   // never true: the macro expanded to 451
```
compared against the wrong number. It is now a plain decimal literal, `100000`
for 1.0.0, computed as `major * 100000 + minor * 100 + patch`.

@section migrate-numbers Numerical differences and how to re-baseline

Results change in three places. All three are fixes or deliberate improvements,
but a snapshot test will notice.

- **`tll` fits** differ near the boundary of the unit square: the interpolation
  grid's endpoint clamping was removed.
- **Maximum-likelihood estimates** shift in the low digits: BOBYQA was replaced
  by Brent and BFGS.
- **Kendall's tau of `bb6`, `bb7`, `bb8` and `tawn`** changes, in places by a
  lot. Four numerical defects were fixed; the worst returned about `1e-11` where
  the true value is `0.33`, for Tawn with theta above roughly 50. If you have
  fitted Tawn models with a large theta, their tau was wrong.

To re-baseline, regenerate the expected values rather than widening tolerances,
and check the direction of change against the corrected values. The parity
harness described in `scripts/README.md` compares two builds surface by surface
and is the tool we use for exactly this.

@section migrate-structure R-vine matrix representation

The structure is now stored with the conditioned variable on the diagonal. For
the same fitted model, `get_matrix()`, `get_order()` and `get_struct_array()`
return a different — equivalent — representation than in 0.7.3, and individual
pair copulas may be stored with the opposite orientation.

Densities, log-likelihoods, simulation and the Rosenblatt transform are
unaffected. What breaks is code that compares a matrix against a stored literal,
or that assumes a particular `(tree, edge)` holds a particular pair. Compare
models by their values, or rebuild expectations from `get_trees()`, which gives
each edge as a conditioned pair plus a conditioning set and is representation
independent.

@section migrate-build Build, CMake, and C++17

- **C++17** is required, along with CMake 3.14 and Boost 1.75.
- Link the namespaced target:
  ```cmake
  find_package(vinecopulib REQUIRED)
  target_link_libraries(my_target PRIVATE vinecopulib::vinecopulib)
  ```
  The unqualified `vinecopulib` name still works as an alias.
- **Eigen must provide its CMake config**; module-mode-only installs no longer
  configure.
- The options gained a `VINECOPULIB_` prefix. `OPT_ASAN`, `CODE_COVERAGE`,
  `STRICT_COMPILER` and `BUILD_DOC` are gone — passing one is inert, and CMake
  reports it as an unused variable. `BUILD_SHARED_LIBS` was replaced by
  `VINECOPULIB_PRECOMPILED`.
- **Sanitizers now default to off** (`VINECOPULIB_SANITIZERS`), and
  `-march=native` is behind `VINECOPULIB_NATIVE_ARCH`, so the default release
  build is redistributable. Turn the latter on for benchmarking.
- Compiler flags are appended rather than assigned, so distribution hardening
  flags and a user's `CXXFLAGS` survive.
- **Boost's concept assertions are no longer disabled.** The umbrella header
  used to `#undef BOOST_CONCEPT_ASSERT` in every translation unit that included
  it, including yours. If your own Boost usage now emits concept diagnostics,
  they were always there.
