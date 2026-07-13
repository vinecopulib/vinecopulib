# Analytic derivative code generation (BB and Tawn families)

The analytic derivative leaves (`pdf_deriv_raw`, `pdf_deriv2_raw`,
`hfunc1_deriv_raw`, `hfunc1_deriv2_raw`, and for Tawn also `hfunc2_deriv_raw` /
`hfunc2_deriv2_raw`) of the BB1/BB6/BB7/BB8 and Tawn copulas are **machine
generated**. The generated C++ lives between the

```
// [BEGIN generated derivative leaves - regenerate with tools/derivs/generate.py]
...
// [END generated derivative leaves]
```

markers in `include/vinecopulib/bicop/implementation/{bb1,bb6,bb7,bb8,tawn}.ipp`.
**Do not edit that code by hand** — edit the generator and regenerate.

## Why generated

VineCopula has no closed-form derivatives for these families, so there is no C
source to port. Differentiating the closed-form copula directly is intractable
(4th-order mixed derivatives of the BB6/7/8 inverses). Instead we differentiate
the *generator* / *Pickands* function — cheap — and feed it through the
family-independent "derivative calculus" of the derivation note
(`bb_tawn_analytic_derivatives.tex`), expanded per selector and CSE'd. This
matches the speed of the hand-written families (Clayton, Joe, ...).

## Files

- `assembly.py` — the symbolic derivative calculus (Archimedean and
  extreme-value), i.e. the `T_i`/`T_ij`, implicit `w_i`/`w_ij`, coordinate
  pullback, and `D f = f · L`, `D² f = f · (L'' + L'⊗L')` assembly. Family
  independent.
- `generate.py` — per-family specs (generator `phi` / Pickands `A`, inverse
  preamble, power-rewrite base) + the emitter that CSE's each selector and
  splices the C++ into the `.ipp` files.

## Regenerate

From the repo root:

```bash
python3 tools/derivs/generate.py            # all families
python3 tools/derivs/generate.py bb7 tawn   # a subset
```

Requires `sympy`. Re-running is idempotent (it replaces the marked region).

## Correctness

`generate.py` asserts that its power-rewrite is an exact symbolic identity. The
generated leaves are then gated in CI by
`ParBicopTest.derivatives_match_finite_differences` (analytic vs. independent
central finite differences, all families and rotations). The derivation itself
was verified by hand, symbolically (sympy), and against 45-digit numerical
partials of the closed-form copula.
