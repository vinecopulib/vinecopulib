# Releasing

Releases are tags on `main`. The steady state is a `NEWS.md` whose top heading
says `(unreleased)`; a dated heading means a release is in flight, and a dated
heading with no tag is a bug — which is what `release-drift.yml` looks for
weekly. That is how 0.8.0 ended up fully written up and never released.

## Steps

1. Confirm `main` is green and every intended pull request has merged.
2. Open a release pull request that retitles the top `NEWS.md` heading from
   `(unreleased)` to the release date and — if the version was not already
   bumped along with the changelog — sets it in the places that must agree:
   - `CMakeLists.txt` — `project(vinecopulib VERSION X.Y.Z ...)`
   - `include/vinecopulib/version.hpp` — `VINECOPULIB_VERSION` (a **plain
     decimal** literal: `000800` is octal and does not mean what it looks like)
     and `VINECOPULIB_LIB_VERSION`
   - `CITATION.cff` — `version` and `date-released`
   - `.zenodo.json` — `version`

   Bumping the version with the changelog and dating the heading here is the
   usual split: `check_version.py` compares the `NEWS.md` heading against the
   project version, so doing both at once keeps every intermediate pull
   request green.
3. Wait for the `version` CI job. `scripts/check_version.py` fails if any of
   those disagree with each other or with `NEWS.md`.
4. Merge, then tag:
   ```bash
   git tag -a vX.Y.Z -m "vinecopulib X.Y.Z"
   git push origin vX.Y.Z
   ```
5. Confirm `release.yml` created the GitHub release from the `NEWS.md` section,
   and `docs.yml` deployed the website.
6. Confirm the Zenodo deposition; update the DOI badge if the concept DOI
   changed.
7. **Open the next cycle immediately**: add a `## vinecopulib X.Y+1.0
   (unreleased)` heading to `NEWS.md`. This is the step that prevents a repeat —
   without it, a dated heading is ambiguous.
8. File pin-bump issues in
   [rvinecopulib](https://github.com/vinecopulib/rvinecopulib) and
   [pyvinecopulib](https://github.com/vinecopulib/pyvinecopulib). Both pin a tag
   of this repository, so every public-API or behavior change here is a real
   change for their users.

## What the automation covers

| workflow | trigger | does |
| --- | --- | --- |
| `continuous_integration.yml` (`version` job) | every pull request | `check_version.py`, i.e. internal consistency |
| `release.yml` | `v*` tag | `check_version.py --released --tag`, re-runs the matrix, extracts the `NEWS.md` section, creates the release |
| `docs.yml` | `v*` tag | builds and publishes the m.css website |
| `release-drift.yml` | weekly | fails if `CMakeLists.txt`'s version has no matching tag |
