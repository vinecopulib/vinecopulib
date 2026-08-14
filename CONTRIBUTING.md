# Contributing to vinecopulib

Thanks for taking the time. This file covers the mechanics; the engineering
invariants — module boundaries, API-stability policy, conventions — are in
[AGENTS.md](AGENTS.md), which is worth reading before a first change.

## Pull requests

Pull requests go to `main`. Releases are tags on `main`; there is no long-lived
development branch.

**One pull request, one commit on `main`.** Pull requests are squash-merged, so
the squashed message — not the intermediate commits — is what explains the change
to a future reader. Write it as a summary of the whole pull request: what
changed, why, and anything someone will need later (behavior changes,
follow-ups). Iterations taken to reach a working state are review history, not
project history.

Commit subjects follow the convention already visible in `git log`:

```
feat(vinecop): support conditioning sets in Rosenblatt transforms
refactor(bicop)!: remove members deprecated since 0.3.1
fix(vinecop): evaluate a custom tree criterion serially
docs: rewrite the user documentation for 1.0.0
```

A `!` marks a breaking change. Scopes are the module directories: `bicop`,
`vinecop`, `misc`.

For work that builds on work still in review, prefer
[stacked pull requests](https://docs.github.com/en/pull-requests/how-tos/stacked-pull-requests)
over one large branch: each gets its own review and its own CI run.

## Before you push

```bash
# format (exactly clang-format 14; other versions produce a different result)
pip install clang-format==14.0.6
clang-format --version    # must print 14.0.6
git ls-files '*.hpp' '*.ipp' '*.cpp' | grep -v nlohmann_json | xargs clang-format -i

# build and test
cmake --preset dev        # Debug + sanitizers + strict warnings
cmake --build build-dev -j
cd build-dev && ctest -j$(nproc)
```

While iterating on one area, `make test_<area>_only` builds and runs just that
area, and `--gtest_filter` selects a subset without rebuilding.

The R parity tests need R with **VineCopula >= 2.6.2 from GitHub, not CRAN**.
Without R they skip: configure with `-DVINECOPULIB_R_PARITY_TESTS=OFF` to make
that explicit.

Anything touching numerics should also be checked with the parity harness — see
[scripts/README.md](scripts/README.md).

## Changelog

Every user-visible change adds a bullet to [NEWS.md](NEWS.md) under the
unreleased heading, in the section that fits: `BREAKING API CHANGES`,
`BEHAVIOR CHANGES`, `NEW FEATURES`, `PERFORMANCE`, `BUG FIXES`,
`BUILD SYSTEM AND DEPENDENCIES`, or `DOCUMENTATION AND TOOLING`.

Style, which the existing entries follow:

- one `* ` bullet per change, one blank line between bullets;
- wrapped at 80 columns, continuations indented two spaces;
- **1-3 lines, hard cap 4** — anything longer belongs in the migration guide,
  which the bullet links;
- imperative present tense, identifiers in backticks, no bold;
- a trailing `(#NNN)` with no period after it; several as `(#659, #667)`.

## Documentation

Code examples on the `docs/*.dox` pages are `\snippet`s of compiled sources under
`docs/snippets/`, and CI compiles *and runs* them. Add an example by adding it
there and referencing it — do not paste code into the prose, or it will rot.

Doxygen runs with `WARN_AS_ERROR`, so an undocumented parameter or a broken
reference fails the build. See [docs/create-docs.md](docs/create-docs.md).

## Releasing

See [docs/releasing.md](docs/releasing.md).
