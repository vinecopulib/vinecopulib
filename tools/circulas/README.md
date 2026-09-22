# Circular copulas: reference calculations

Independent numerical references for the formulas fixed in
[docs/plans/circular-copulas-contract.md](../../docs/plans/circular-copulas-contract.md).
Nothing here ships in the library; the scripts exist so that the C++
implementation is tested against numbers it did not produce.

## Files

- `circulas.py` — the formulas of the contract in NumPy: the three circular
  densities with their lifted CDFs and inverses, the binding-density circula
  (density, CDF, h-functions, inverses), and the two cylindrical
  sections families. Each function follows the contract's notation.
- `check_contract.py` — asserts every identity of the contract's acceptance
  table against numerical integration and differentiation of the density:
  normalization, margins, periodicity, `h = dC`, inverse round-trips, the
  lifted-CDF properties, the rotation and flip tables, the identifiability
  normalizations, the closed-form dependence measures, and the half-turn case.
  Exits nonzero on the first failure.
- `dump_golden.py` — writes reference values (density, CDF, h-functions, and
  inverses on a fixed lattice of arguments and parameters, plus Kendall's tau)
  to a JSON file, or with `--cpp` regenerates
  `test/src_test/include/circular_golden.hpp`, the header the C++ golden tests
  read them from.

## Run

From the repo root, with NumPy, SciPy, and SymPy available:

```bash
python3 tools/circulas/check_contract.py
python3 tools/circulas/dump_golden.py --cpp   # regenerates the test header
```

The checks take under a minute. Any change to a formula in the contract must
be mirrored in `circulas.py` and pass `check_contract.py` before the C++ side
is touched.
