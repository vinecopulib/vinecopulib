# Circular copulas and mixed vines: implementation plan

Status: stages 1 to 6 implemented on the integration branch `feat/circulas`;
stage 7 (release validation) is next.
Updated September 22, 2026.

This checklist tracks the first feature release for circular variables. The
release includes circular-circular and circular-linear pair copulas,
nonparametric estimation, and integration with `Vinecop`. Individual stages
are reviewable implementation steps, not separate definitions of a complete
feature release. Follow [AGENTS.md](../../AGENTS.md) throughout.

Files under `docs/plans/` are excluded from the Doxygen input and never
appear on the website.

## Release scope

- Fit, select, evaluate, simulate, and serialize circular-circular and
  circular-linear `Bicop` models.
- Support mixtures of circular and linear variables in `Vinecop`, including
  supplied structures, automatic selection, and Rosenblatt transforms.
- Provide parametric families and a nonparametric estimator that respects
  the geometry of each axis.
- Preserve existing linear models, defaults, serialized models, and the
  continuous/discrete meaning of `"c"` and `"d"` in `var_types`.
- Continue to accept copula-scale data. The library sees a circular variable
  as `u` in `[0, 1]` with `0` identified with `1`. Units, angular direction,
  the cut, and the angle-to-CDF mapping are downstream concerns (marginal
  estimation and the R/Python bindings); the library documents them but does
  not implement them.

The initial implementation targets continuous circular and continuous linear
variables. New combinations involving circular and discrete variables need a
separate mathematical and API decision; reject unsupported combinations
explicitly. Existing discrete functionality must continue to work.

### Proposed family set

These are mathematical family names; enum names and parameter conventions
will be settled in stage 1.

| Pair geometry | Parametric candidates | Nonparametric candidate |
| --- | --- | --- |
| Linear-linear | Existing families | Existing TLL |
| Circular-circular | Cardioid, wrapped Cauchy, and von Mises binding densities | Local likelihood on two periodic axes |
| Circular-linear, either order | Cubic sections with phase (quadratic sections are the case $a = b$); reuse the binding-density construction where appropriate | Local likelihood on one periodic axis and one probit-transformed axis |
| Any supported pair | Existing independence family | — |

The binding-density construction is

$$
c(u,v)=2\pi g\{2\pi(v-qu)-\mu\},\qquad q\in\{-1,1\},
$$

where `g` is a centered circular density, `mu` is a phase, and `q` selects
the association orientation. It is periodic in both coordinates. It is also
valid for circular-linear pairs, but imposes equality at the two ends of
the linear coordinate; it cannot replace the cylindrical sections families.
See [Jones, Pewsey, and Kato](#references) and [Hodel and Fieberg](#references).

Facts about this construction that later stages rely on (all three proposed
`g` are symmetric about zero):

- The h-functions are differences of the circular CDF `G`, and the inverse
  h-functions reduce to `G^{-1}`. Cardioid and wrapped Cauchy have closed-form
  `G`; von Mises needs a series or quadrature for `G` and a root solve for
  `G^{-1}`. Boost.Math, already a dependency, provides the Bessel functions.
- Rotating by 180 degrees maps `mu` to `-mu`; rotating by 90 degrees is
  exactly `q = -1`. The families therefore have two distinct rotations,
  `{0, 90}`, not one and not four.
- `flip()` negates the phase when `q = 1` and is the identity when `q = -1`
  (the density is then exchangeable).
- With `q = 1` and a half-turn phase, a perfectly dependent pair has linear
  Kendall's tau exactly zero. Linear tau must not drive candidate
  generation, preselection, or tree weights for circular pairs.

Multimodal parametric mixtures, general Fourier families, and a separate
Bernstein estimator are later extensions. The nonparametric estimator should
already accommodate patterns beyond a single circular diagonal.

## Stages and dependencies

| Stage | Depends on | Reviewable outcome |
| --- | --- | --- |
| 1. Mathematical and API contract | — | Documented conventions, family formulas, and acceptance cases |
| 2. Geometry and compatibility | 1 | Geometry token, family eligibility, and serialization contracts |
| 3. Parametric pairs | 2 | Complete circular-circular and circular-linear `Bicop` support |
| 4. Nonparametric pairs | 2, 3 | Periodic and mixed-axis fits, usable wherever a parametric circular pair is |
| 5. Supplied vine structures | 3 | Mixed-vine fitting, evaluation, transforms, and simulation; `select` rejects circular input |
| 6. Automatic selection | 5 | Geometry-aware family and structure selection |
| 7. Release validation | 1–6 | Validated mixed-vine feature with documentation and examples |

Stage 4 is a parallel track. It is exploratory by nature and must not sit on
the critical path: stages 5 and 6 complete with parametric circular families
only, and stage 4 plugs in through the eligibility filtering from stage 2.

All stages are collected on the integration branch `feat/circulas`. Each
stage is one pull request against that branch, stacked on its predecessor
where the dependency table requires it, and made of meaningful intermediate
commits rather than a single squashed step. `feat/circulas` merges to `main`
after stage 7. Collect entries in the unreleased section of `NEWS.md` and cut
the release after the merge. Record actual PR numbers beside completed tasks.
Do not merge, tag, or publish without express authorization.

## 1. Mathematical and API contract

- [x] (#785) Fix the spelling of the new `var_types` value; the token approach in
  stage 2 is decided. The spelling reaches the R and Python signatures, so it
  is settled here, not later. Settled: `"a"` (angular); see the decision
  record.
- [x] (#785) Document that the library takes `u` in `[0, 1]` with `0` identified
  with `1`, and that the cut is wherever the caller's marginal CDF puts it.
  State that ranks (`to_pseudo_obs`) place the cut at the caller's zero angle. Done in [circular-copulas-contract.md](circular-copulas-contract.md), *Copula-scale conventions*.
- [x] (#785) Define endpoint behavior: densities join at `0` and `1` on circular
  axes; CDFs and h-functions retain their anchored probability semantics.
  Inverses must return the correct branch on the chosen interval. Done, *Endpoint behavior*.
- [x] (#785) Derive density, CDF, both h-functions, and both inverses for each
  proposed parametric family. Identify numerical integrations and root solves
  (von Mises `G` and `G^{-1}`; cubic sections constraints). Done, *Binding-density circulas* and *Cylindrical sections copulas*;
  verified by `tools/circulas/check_contract.py`.
- [x] (#785) Fix concentration bounds, the independence limit, parameter
  identifiability, and the number of fitted parameters. The wrapped Cauchy
  concentration needs an upper bound strictly below its singular limit.
  Derive parameter constraints for the phased cubic-sections family explicitly. Done, *Parameter domains* and *Cubic sections*: the cubic feasible region
  is exactly the box, so only box constraints are needed.
- [x] (#785) Fix the orientation and phase conventions: `q` is represented by the
  copula rotation in `{0, 90}`; `mu` is a periodic, unbounded parameter in
  the optimizer (no bound to hit, no finite-difference problem at a cut).
  Specify how `allow_rotations` maps onto the two distinct rotations. Done, *Phase parameters* and *Orientation and rotation*.
- [x] (#785) Record the `flip()` transformations for every family, including the
  phased cylindrical families, which need not be exchangeable. Done, *Flip* and *Symmetries of the cylindrical families*.
- [x] (#785) Define ordinary Kendall's tau and Blomqvist's beta relative to the
  chosen cut. Keep their current meanings; expose circular dependence
  summaries separately. Specify tail-dependence behavior. `itau` stays
  unavailable for circular families; `parameters_to_tau` returns `NaN` and
  its inverse is not available. Done, *Dependence measures*.
- [x] (#785) State the periodicity condition for a mixed vine: the conditional CDF
  of a circular variable wraps from `0` to `1`, so every pair copula in which
  that variable (or its h-transform) is a conditioned argument must have
  equal density at both ends of that coordinate. A circular variable that
  appears only in a conditioning set imposes nothing on the pair copula. Done, *Mixed vines: the periodicity condition*, with a proof sketch.
- [x] (#785) Distinguish density periodicity from changing the cut and refitting.
  Investigate how changing cuts affects the simplifying assumption: a shift
  of a conditional circular CDF can depend on the conditioning values.
  Do not promise arbitrary-cut invariance of simplified vines without proof. Done, *Changing the cut*: no invariance is promised beyond tree 1.
- [x] (#785) Specify tolerances and independent numerical reference calculations
  for the identities used in the later test stages. Done, *Acceptance cases and tolerances*; references in `tools/circulas/`.

## 2. Geometry, family eligibility, and compatibility

Proposal: a third `var_types` value marking a continuous circular variable.
`"c"` and `"d"` keep their meaning; the new value is a kind of continuous
variable. This reuses the existing propagation through trees
([vinecop/class.ipp](../../include/vinecopulib/vinecop/implementation/class.ipp),
`set_var_types_internal`, and the selector's per-edge `var_types`), the
constructors' default arguments, the `"vt"` JSON field, the views, and the
binding signatures. Legacy JSON reads work unchanged. A separate per-variable
attribute would duplicate all of that; the circular-discrete combination,
when it comes, can be a further value.

- [x] (#790) Accept the new value in `Bicop::check_var_types` and
  `Vinecop::check_var_types`.  Audit every `== "d"` and `get_n_discrete()`
  branch to confirm it stays correct with a third value; `as_continuous()`
  and the `BicopView` / `VinecopView` continuous paths must preserve it.
  Audit result (September 21, 2026): the `== "d"` and `get_n_discrete()`
  branches are safe; the hazards are the twelve tests that check continuity
  positively (`== "c"`, `== {"c", "c"}` and their negations) in
  `abstract.ipp`, `parametric.ipp`, `Bicop::check_deriv_preconditions`, and
  the vine score path, plus the literal `{"c", "c"}` returned by
  `Bicop::as_continuous` and `BicopView::get_var_types`. Replace all literal
  comparisons by named type predicates before any circular code lands.
- [x] (#790) Decide whether circular families infer geometry when constructed
  explicitly, and validate contradictory family/geometry combinations. Decided: no inference; a circular family constructed with linear
  types throws, as does a linear family with a circular type.
- [x] (#790) Define family capabilities for both argument orders, including the
  geometry-dependent nonparametric estimator. Include independence in every
  supported candidate set. `family_accepts_var_types()` / `eligible_families()`.
- [x] (#790) Update family registration, name conversion, convenience groups, and
  factory dispatch. Add the new families to `bicop_families::all`; the
  eligibility filter keeps the effective default search unchanged for
  linear callers. Add a group for the two-rotation families. Enum values, names, `two_rotations`, `circular`, and `cylindrical`
  groups added; the families join `all`, `parametric`, `two_par`, and
  `three_par` in stage 3, when they become constructible.
- [x] (#790) Thread geometry into `tools_select::create_candidate_bicops`, which
  today receives only data and controls, chooses rotations from the sign of
  linear tau, and would drop one orientation of a circular family. 
- [x] (#790) Specify filtering for explicit and empty family sets and behavior when
  no compatible candidate remains. Apply the same rules to fitting and selection. 
- [x] (#791) Preserve geometry through copying, views, flipping, rotations,
  resetting, truncation, structure conversions, and conditional/reoriented
  interfaces. Done at the `Bicop` level (copy, `BicopView::as_continuous`,
  `flip()` with the two-rotation canonicalization, `set_rotation`, JSON);
  the vine-level interfaces (truncation, relabeling, `VinecopView`,
  conditional simulation and reorientation) are verified by the stage 5
  tasks once a circular vine can be built.
- [x] (#791) Define the JSON field for nonparametric grid geometry. The
  `Bicop` JSON carries only family, rotation, parameters, `var_types`, and fit
  statistics, and `KernelBicop::set_parameters` rebuilds knots from the row
  count alone, so a circular grid cannot reload without it. Test legacy reads
  and complete round-trips; reject unsupported formats clearly. Defined in the
  contract (*Nonparametric grid serialization*); stage 4 implements and tests
  it.
- [x] Audit all special cases keyed on `BicopFamily::tll` before choosing
  between geometry-dependent `tll` and separate family identifiers. Four
  sites: the factory in `abstract.ipp`, the two parameter-skipping branches
  in `Bicop::hfunc1_continuous` / `hinv2_continuous` (which should test
  `bicop_families::nonparametric`), and the grid line of both `str()`
  methods. Decision: one geometry-dependent `tll`; see the decision record.

Main code: [bicop/class.hpp](../../include/vinecopulib/bicop/class.hpp),
[family.hpp](../../include/vinecopulib/bicop/family.hpp),
[abstract.hpp](../../include/vinecopulib/bicop/abstract.hpp),
[vinecop/class.hpp](../../include/vinecopulib/vinecop/class.hpp), their `.ipp`
implementations, and fit controls.

## 3. Parametric pair copulas

- [x] (#791) Implement cardioid, wrapped Cauchy, and von Mises binding families in
  the existing `.hpp` / inline `.ipp` pattern, sharing the circular primitives
  (`g`, `G`, `G^{-1}`) that have identical semantics.
- [x] (#791) Implement cubic sections (quadratic sections are the case $a = b$),
  phase handling, and both axis orders. Support circular-linear use of binding families without duplicating
  their mathematical implementations.
- [x] (#791) Implement stable lifted circular CDFs and inverses, retaining the number
  of completed turns. Audit cancellation and concentration limits, including
  the von Mises normalizing constant and wrapped Cauchy near-singular limit.
- [x] (#791) Implement likelihood fitting with circular initialization and periodic
  phase optimization. The existing tau-based starting values and search
  bounds must not constrain these fits incorrectly.
- [x] (#791) Implement or validate derivative support. Support broadcast and
  per-observation parameter matrices under the existing contract.
- [x] (#791) Honor observation weights, missing-value handling, parameter validation,
  seeded simulation, and fit statistics.
- [x] (#791) Update candidate construction and preselection: produce exactly the two
  rotations for circular families and bypass the tau-sign and tail
  (`lt` / `ut`) heuristics for them.
- [x] (#791) Test normalization, uniform marginals, CDF boundary values, CDF/density
  derivatives, h/inverse identities, independence limits, flip identities,
  simulation moments, and recovery of phases near the cut. Include the
  zero-tau half-turn case as a fit-recovery test. `BindingBicop` holds the leaves in terms of the lifted CDF; the
  families supply `g`, the lifted CDF, its inverse, and the moment map. `SectionsBicop` reads the argument order from `var_types`. Bracketed Newton for the cardioid and von Mises inverses; closed
  form for the wrapped Cauchy. `CircularBicop::fit` starts from moment estimates and calls the
  extracted `ParBicop::fit_mle` with an unbounded phase. Finite-difference fallback of `ParBicop`; per-row parameters via
  `binaryExpr_or_nan`. Done in stage 2b. `test_circular`, with golden values from `tools/circulas`.

Main code: [parametric.ipp](../../include/vinecopulib/bicop/implementation/parametric.ipp),
[tools_select.ipp](../../include/vinecopulib/bicop/implementation/tools_select.ipp),
new family files, and numerical helpers under `misc/`.

## 4. Nonparametric estimation (parallel track)

The current TLL fit transforms both coordinates with `qnorm`. Its local
likelihood formulas, bandwidth selection, and influence calculation assume
that transformed geometry. A circular extension needs new fitting formulas;
matching the two boundary rows after an ordinary TLL fit is insufficient.
A circular axis has no boundary, so it needs no transform; the probit
transform and its Jacobian apply to linear axes only.
- [x] (#795) Prototype periodic smoothing in circular coordinates and probit
  smoothing in linear coordinates. Compare wrapped Gaussian and von Mises
  kernels.
  Chosen: the von Mises kernel. Its local log-linear fit in
  $(\cos\theta, \sin\theta)$ has a closed form through the Bessel
  ratio $A = I_1 / I_0$ and its inverse, the circular analog of the
  Gaussian closed form on a linear axis; a wrapped Gaussian kernel has no
  such form. The two kernels are close for every concentration in use.
- [x] (#795) Derive the periodic/mixed local-likelihood estimator and its supported
  polynomial orders. Use a positive periodic kernel estimate as a numerical
  reference. Specify unsupported method choices rather than silently
  substituting a different estimator.
  Product kernel (von Mises on circular axes, Gaussian on the probit scale
  of linear axes) with a local log-linear model per axis and no interaction
  term, so the local likelihood equations separate and each axis contributes
  a closed-form correction. Orders `constant` and `linear` have closed
  forms; `quadratic` (the wrappers' default) adds the second harmonic on a
  circular axis, whose normalizer has no closed form and is handled by
  trapezoid quadrature and a damped Newton solve of the moment equations
  (#795).
- [x] (#795) Define bandwidth selection for both geometries, including observation
  weights. Evaluate held-out likelihood, cut sensitivity, and strong or
  multimodal dependence before choosing defaults.
  The kernel variance is the same fraction of the transformed margin's
  variance as on a linear axis ($n^{-1/3}$, times 1.5 for the linear
  fit), the angle of a uniform variable having variance $(2\pi)^2/12$;
  a factor $1 - R$ with the circular association $R$ narrows it under
  strong dependence, and `nonparametric_mult` scales it. Weights enter the
  kernel sums. Held-out log-likelihood minus the truth, n = 500, 4
  replications, multipliers 0.25 / 0.5 / 1 / 2 / 4: von Mises kappa = 2,
  constant: -0.08 / -0.03 / -0.03 / -0.07 / -0.20, linear: -0.18 / -0.06 /
  -0.02 / -0.02 / -0.06, quadratic: -0.07 / -0.03 / -0.03 / -0.09 / -0.21;
  quadratic sections a = 0.8, constant: -0.05 / -0.02 / -0.01 / -0.02 /
  -0.03, quadratic: -0.06 / -0.02 / -0.01 / -0.01 / -0.02; a bimodal
  mixture of two circulas prefers
  0.25 to 0.5; the half-turn wrapped Cauchy with rho = 0.95 loses about 0.7
  at every multiplier, a ridge too narrow for a product kernel. The default
  multiplier 1 is at or next to the optimum in every other case, at n = 500
  and n = 2000.
- [x] (#795) Give each axis its own knot vector. `InterpolationGrid` takes one
  vector for both axes and `make_normal_grid` concentrates knots in the
  tails; a circular axis needs uniform knots on `[0, 1]`.
  `InterpolationGrid` holds one knot vector per axis (PR #795, first
  commit); `KernelBicop::make_grid_points` gives uniform knots on `[0, 1]`
  for a circular axis and the normal grid otherwise.
- [x] (#795) Represent a circular grid boundary consistently, either with shared
  endpoint values or an explicit wraparound cell. Keep density nonnegative
  and enforce uniform marginal integrals with the appropriate grid weights.
  Shared endpoint values: the knot at 1 repeats the knot at 0, the
  estimator copies the value, and the trapezoid weights of the two ends add
  up to one interior weight. Values are positive by construction.
- [x] (#795) Ensure margin normalization preserves periodic endpoint equality.
  Require both constraints to meet tolerance before accepting the fit.
  The Sinkhorn passes rescale rows and columns; identical end rows receive
  identical factors, so equality survives (tested).
- [x] (#795) Adapt interpolation, one- and two-dimensional integration, and direct
  conditional inversion. Preserve anchored CDF values at `0` and `1`.
  Nothing to adapt: the anchored integrals, conditional masses, and direct
  inversion work on any ascending knot vector; the periodic ends only need
  equal values.
- [x] (#795) Preserve knots, axis geometry, density values, and effective degrees of
  freedom through `get_parameters`, `set_parameters`, flipping, and JSON.
  Values alone must not reload onto a different grid.
  `get_parameters` returns the values; `set_parameters` places them on the
  default knots of the variable types (rebuilding the grid when the shape
  changes); `set_var_types` rebuilds the knots when the geometry of an axis
  changes; `flip` exchanges knots with values; JSON records `"grid"` for a
  circular pair and validates it on reading.
- [x] (#795) Derive or validate the influence/effective-degrees-of-freedom calculation
  before using AIC, BIC, or mBIC penalties. Reusing the Gaussian TLL influence
  formula for a different smoother needs justification.
  Derived: the influence is $K(0)\,[M^{-1}]_{00} / n$ with the local
  information $M = f_0\, E[(1, \psi)(1, \psi)^\top]$ under the fitted
  local model, whose moments are Gaussian on a linear axis and Bessel ratios
  ($I_1/I_0$, $I_2/I_0$) on a circular one; it reduces to the
  existing formula when both axes are linear.
- [x] (#795) Test periodic densities and arbitrary asymmetric/multimodal patterns,
  h/inverse identities, grid refinement, bandwidth extremes, fit statistics,
  axis swaps, and compatibility with existing linear TLL reference fits.
- [x] (#795) Register the estimator in the eligibility tables from stage 2 and add
  it to the mixed-vine tests of stages 5 and 6 once those have landed.

Main code: [tll.ipp](../../include/vinecopulib/bicop/implementation/tll.ipp),
[kernel.ipp](../../include/vinecopulib/bicop/implementation/kernel.ipp), and
[tools_interpolation.ipp](../../include/vinecopulib/misc/implementation/tools_interpolation.ipp).
  `family_accepts_var_types` already admitted `tll` everywhere; the vine
  test fits a mixed D-vine with `tll` on every edge through `select`.

## 5. Vine models with supplied structures

- [x] (#793) Derive pair geometry from the conditioned variable identities at every
  tree level, following the existing `var_types` propagation. A circular
  variable retains its geometry after its conditional probability transform;
  a circular variable only in the conditioning set does not make two linear
  conditioned variables circular. Implemented as `Vinecop::edge_var_types`,
  read directly from the structure, and used both for the propagation and
  for omitted pair copulas.
- [x] (#790) Make `Vinecop::select` (and the threshold and truncation searches) throw
  on circular input until stage 6 lands. With the default tau criterion a
  perfectly dependent half-turn pair has weight zero and any positive
  `threshold` sets it to independence.
- [x] (#793) Validate supplied pair families and propagate geometry through order
  changes, edge flips, truncation, and omitted independence edges.
- [x] (#793) Support sequential refitting (`fit` on a fixed structure) and verify
  `pdf`, `pdf_full`, `loglik`, `cdf`, `rosenblatt`, `inverse_rosenblatt`,
  and simulation.
- [x] (#793) Cover conditional simulation, supported reorientations, and the
  associated views without losing axis geometry or nonparametric grids.
- [x] (#793) Verify scores and Hessians for supported parametric models; retain
  explicit unsupported-operation behavior for nonparametric derivatives. The
  circular families use the finite-difference derivative leaves; the full
  (`step_wise = false`) vine scores agree with finite differences of the
  joint log-likelihood.
- [x] (#793) Test small circular-linear-linear and circular-circular-linear models
  against an independently assembled pair-density product and numerical
  marginalization. Include higher-tree circular edges and a larger mixed vine.
- [x] (#793) Verify joint density equality at each circular boundary, Rosenblatt
  round-trips, simulated uniform marginals, and serial/threaded consistency.

Main code: [vinecop/class.hpp](../../include/vinecopulib/vinecop/class.hpp),
[vinecop/tools_select.hpp](../../include/vinecopulib/vinecop/tools_select.hpp),
their implementations, and structure/view conversion code.

## 6. Automatic family and structure selection

- [x] (#794) Apply compatible candidate sets at every tree level, including the
  nonparametric estimator when available and both circular-linear argument
  orders.
- [x] (#794) Benchmark the existing bounded criteria first: `cxi` (Chatterjee's xi,
  detects any functional relationship and is nearly cut-invariant) and
  `hoeffd`, on the zero-tau half-turn case, reflective association, and
  multimodal patterns. Compare against circular correlation statistics and
  likelihood gain over independence; the latter is unbounded and costs one
  full candidate fit per edge per tree, so it is an option, not the baseline.
  Benchmarked in the PR (n = 500, 20 replications): |tau| and rho_S
  vanish on the half-turn pair (0.06, 0.48 -> 0.02 under a shifted cut) and
  change with the cut for every family; `cxi` is nearly cut-invariant and
  strong for functional relationships (0.81 on the half-turn pair) but weak
  for diffuse dependence (0.05 for von Mises kappa = 0.5, 0.10 for the
  quadratic sections with a = 1, against 0.02 under independence); `hoeffd`
  is weaker still. The moment-based circular measure is exactly
  cut-invariant and separates every case (0.25, 0.29, 0.95 against 0.05).
- [x] (#794) Choose the default criterion for pairs involving a circular variable
  and specify weighting, missing observations, small samples, and
  independence behavior.
  Chosen: `tools_stats::pairwise_circular()` for every pair with a
  circular conditioned variable, under every built-in criterion; `"custom"`
  receives the pair data unchanged. Weights enter as weighted means; missing
  observations are removed by `calculate_criterion` as for the linear
  criteria; below 11 observations the weight is 0 as for the linear criteria;
  independence gives a value of order n^{-1/2}.
- [x] (#794) Integrate the criterion with existing spanning-tree algorithms and
  custom criteria. Cache fitted edge candidates if likelihood-based weights
  are offered.
  No caching needed: the measure costs one pass over the pair.
- [x] (#794) Audit thresholds, automatic threshold search, truncation selection,
  and mBICv against the chosen criterion's scale. Remove the stage 5 rejection
  of circular input.
  The measure is bounded by one like tau, so thresholds keep their
  meaning; the threshold search and truncation selection compare it against
  the same thresholds; mBICv does not depend on the criterion.
- [x] (#794) Test selection on known mixed vines and nonlinear circular associations,
  including the zero-tau half-turn case. Verify reproducibility and behavior
  for explicit family restrictions and incompatible candidate sets.
- [x] (#794) Benchmark selection cost and confirm existing linear defaults retain
  their current effective candidate sets and numerical behavior.
  Selection on the three-dimensional test vines takes well under a
  second; the linear criteria are untouched, since the circular branch is
  entered only when a pair has a circular variable.

Main code: [vinecop/tools_select.ipp](../../include/vinecopulib/vinecop/implementation/tools_select.ipp),
bivariate candidate selection, and both fit-control classes.

## 7. Release checklist

- [ ] All preceding stages have implementation PRs and their acceptance tests.
- [ ] Extend existing GoogleTest areas where suitable; keep circular cases
  independent of R parity fixtures that only recognize existing families.
- [ ] Add numerical reference tests for formulas not supported by those R
  fixtures. Use `cylcop` comparisons where conventions match, with core
  correctness tests that do not require that package.
- [ ] Run clang-format 14 and the required lint/spelling checks, the Debug
  build and `bin/test_all`, and the precompiled Release build and tests.
  Verify an installed consumer can include and use the new families.
- [ ] Benchmark pair fitting/evaluation, grid resolution, and mixed-vine
  fitting/simulation. Record accuracy and runtime together.
- [ ] Add Doxygen documentation and compiled examples for circular-circular
  pairs, circular-linear pairs in both orders, and a mixed parametric/TLL vine.
- [ ] Document marginal transforms, cuts, phase versus copula rotation,
  dependence summaries, simplifying-assumption limitations, supported
  variable types, and nonparametric method choices.
- [ ] Document new API and JSON fields and coordinate the downstream R/Python
  work. Consolidate the unreleased `NEWS.md` entries for the release.
- [ ] Confirm the release delivers mixed-vine selection and nonparametric
  fits as well as the bivariate parametric families.

## References

These sources motivate the design. Application papers provide examples of
mixed vines; the mathematical contracts and numerical tests above remain
required for a general library implementation.

- [Jones, Pewsey, and Kato (2015), *On a class of circulas: copulas for circular
  distributions*](https://doi.org/10.1007/s10463-014-0493-6): binding-density
  construction, phase and orientation, CDFs, and circular dependence measures.
- [Hodel and Fieberg (2022), *Circular-linear copulae for animal movement
  data*](https://doi.org/10.1111/2041-210X.13821): cylindrical quadratic/cubic
  sections, reflection constructions, and the accompanying `cylcop` software.
- [Wang et al. (2021), *Circular-linear-linear probabilistic model based on
  vine copulas*](https://doi.org/10.1016/j.jweia.2021.104704): combining
  cylindrical and ordinary pair copulas in a trivariate model.
- [Nagar et al., *A dependent circular-linear model for multivariate
  biomechanical data: Ilizarov ring fixator study*](https://arxiv.org/abs/2312.10159):
  a six-variable application combining all three pair geometries and truncation.
- [Carnicero, Ausin, and Wiper (2013), *Non-parametric copulas for
  circular-linear and circular-circular data: an application to wind
  directions*](https://doi.org/10.1007/s00477-013-0733-y): continuity constraints
  for circular Bernstein copulas; a comparison for periodic nonparametric fits.
- [Garcia-Portugues, Crujeiras, and Gonzalez-Manteiga, *Kernel density
  estimation for directional-linear data*](https://arxiv.org/abs/1210.3214):
  directional-linear kernel estimation as background for mixed-axis smoothing.

## Decision record

Resolved design choices: decision, rationale, validation evidence, and
implementing PR. Entries without a PR were settled during planning.

| Decision | Rationale | Evidence / PR |
| --- | --- | --- |
| Pair copulas use the binding-density construction for circular-circular pairs and cylindrical sections for circular-linear pairs. | Closed-form or one-dimensional numerics throughout; covers the published mixed-vine applications. | Jones, Pewsey, and Kato; Hodel and Fieberg |
| Continuous circular and continuous linear variables only; circular-discrete is rejected explicitly. | Mixed circular/discrete needs its own mathematical and API decision. | planning |
| Marginal transforms, cuts, and angular conventions stay downstream; the library takes `u` in `[0, 1]` with `0` identified with `1`. | Matches the existing copula-scale contract and the `kde1d` exclusion in AGENTS.md. | planning |
| Geometry is a third `var_types` value. | Reuses propagation, JSON, views, and binding signatures; a parallel attribute would duplicate them. `Vinecop::set_var_types_internal` and the selector's edge inheritance already propagate any token without literal comparisons. | planning (confirmed by maintainer, September 15, 2026; propagation verified September 21, 2026) |
| The token is `"a"` (angular). | Short and parallel to `"c"` / `"d"`. The literal is confined to the type predicates, so the spelling can change later at one site. | maintainer, September 21, 2026 |
| Orientation `q` is the copula rotation in `{0, 90}`; phase `mu` is periodic and unbounded. | Symmetric `g` makes 180 and 270 redundant; a periodic phase has no bound to hit. | planning (confirmed by maintainer, September 15, 2026) |
| `itau` stays unavailable for circular families. | No valid identification result from linear tau. | planning |
| Stage 4 (nonparametric) is a parallel track; stages 5 and 6 ship with parametric circular families. | Exploratory work must not gate the release. | planning |
| One geometry-dependent `tll`; the estimator reads the axis geometry from `var_types`. | The `"vt"` JSON field already records the geometry; a separate identifier would be a second source of truth that can disagree with it. `family_set = {tll}` keeps meaning "nonparametric" for every geometry, as `{gaussian}` needs no circular twin. | audit, September 21, 2026 |
| Family names: `cardioid`, `wrapped_cauchy`, `von_mises` (binding), `cubic_sections` (cylindrical); new group `bicop_families::two_rotations`. | Bare snake_case values match the flat existing enum; `create_candidate_bicops` branches on rotation arity, so the two-rotation families need their own group beside `rotationless`. | maintainer, September 21, 2026 |
| No `quad_sections` family: quadratic sections are `cubic_sections` with $a = b$. | One parameter saved under BIC did not justify a fifth family in every list, docstring, and downstream wrapper. | maintainer, September 22, 2026 (#791) |
| `parameters_to_tau` returns `NaN` for the circular families. | Kendall's tau depends on the cut, so it misreports circular dependence; its nested quadrature also dominated the cost of printing a vine. | maintainer, September 22, 2026 (#791) |
| The CDF of a binding family is closed form through the Fourier coefficients of `g`. | Replaces per-point quadrature of $h_1$ at a thousandth of the cost, with no accuracy loss. | #791 |
| Stage pull requests target the integration branch `feat/circulas`, stacked where a stage depends on its predecessor; `feat/circulas` merges to `main` after stage 7. | Keeps `main` free of a half-finished feature while each stage still gets its own review. | maintainer, September 21, 2026 |

Still open: parameter conventions, nonparametric fitting formulas, and the
default tree criterion for circular pairs.
