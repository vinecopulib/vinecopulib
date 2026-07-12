#!/usr/bin/env python3
"""Compare two parity_dump JSON files with per-key tolerance classes.

Generate the two dumps with parity_dump into bench_results/ (gitignored), e.g.
    bin/parity_dump bench_results/parity_master.json   # on the baseline
    bin/parity_dump bench_results/parity_branch.json   # on the branch

Usage:
    compare_parity.py master.json branch.json [--strict] [--tol-config cfg.json]

Tolerance classes (first matching prefix wins; longest-prefix order below):
    A = 0.0        bit-identical (default for everything)
    B = 1e-12      relative; floating-point reassociation
    C = 1e-8       relative; quantified tolerance/algorithm changes
Keys are slash-joined JSON paths, e.g. "bicop_eval/clayton/rot0/pdf".

--strict forces class A everywhere (use during Phase 1).
--tol-config points to a JSON file of {"prefix": tolerance} overrides that
replace the built-in table (use as gates are relaxed per phase).
"""

import argparse
import json
import math
import sys

# Built-in gate table; update as optimization phases land (see plan).
# Ordered by specificity: first match wins.
DEFAULT_GATES = [
    # Phase 2.6: integration tolerance 1e-12 -> 1e-9 on tau
    # ("tau/", 1e-7),
    # Phase 2: log-pdf pathway + vectorized leaves (fp reassociation)
    # ("bicop_eval/", 1e-12),
    # ("bicop_fit/", 1e-8),
    # Phase 3: tll overhaul
    # ("tll/", 1e-6),
    # ("vinecop/", 1e-9),
]


def flatten(obj, prefix=""):
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.update(flatten(v, f"{prefix}{k}/"))
    elif isinstance(obj, list):
        if all(isinstance(x, (int, float)) for x in obj):
            out[prefix.rstrip("/")] = obj
        else:
            for i, v in enumerate(obj):
                out.update(flatten(v, f"{prefix}{i}/"))
    else:
        out[prefix.rstrip("/")] = [obj]
    return out


def tolerance_for(key, gates):
    for prefix, tol in gates:
        if key.startswith(prefix):
            return tol
    return 0.0


def rel_diff(a, b):
    # nlohmann serializes NaN/Inf as JSON null -> Python None; coerce to NaN so a
    # value going non-finite is reported as a FAIL instead of raising TypeError
    # (which would abort the whole comparison) in the isnan checks below.
    if a is None:
        a = float("nan")
    if b is None:
        b = float("nan")
    if a == b:
        return 0.0
    if math.isnan(a) and math.isnan(b):
        return 0.0
    if math.isnan(a) != math.isnan(b):
        return math.inf
    denom = max(abs(a), abs(b), 1.0)
    return abs(a - b) / denom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("master")
    ap.add_argument("branch")
    ap.add_argument("--strict", action="store_true",
                    help="require bit-identical everywhere")
    ap.add_argument("--tol-config", default=None,
                    help="JSON file: {prefix: tol} overriding the gate table")
    args = ap.parse_args()

    with open(args.master) as f:
        master = flatten(json.load(f))
    with open(args.branch) as f:
        branch = flatten(json.load(f))

    gates = [] if args.strict else list(DEFAULT_GATES)
    if args.tol_config and not args.strict:
        with open(args.tol_config) as f:
            gates = sorted(json.load(f).items(),
                           key=lambda kv: -len(kv[0]))

    missing = sorted(set(master) - set(branch))
    added = sorted(set(branch) - set(master))
    for k in missing:
        print(f"MISSING in branch: {k}")
    for k in added:
        print(f"NEW in branch:     {k}")

    n_fail = 0
    worst = {}
    for key in sorted(set(master) & set(branch)):
        a, b = master[key], branch[key]
        if len(a) != len(b):
            print(f"FAIL {key}: length {len(a)} vs {len(b)}")
            n_fail += 1
            continue
        tol = tolerance_for(key, gates)
        max_d = max((rel_diff(x, y) for x, y in zip(a, b)), default=0.0)
        group = key.split("/")[0]
        worst[group] = max(worst.get(group, 0.0), max_d)
        if max_d > tol:
            print(f"FAIL {key}: max rel diff {max_d:.3e} > tol {tol:.1e}")
            n_fail += 1

    print("\nworst relative difference per group:")
    for group in sorted(worst):
        print(f"  {group:<14} {worst[group]:.3e}")

    if missing or n_fail:
        print(f"\n{n_fail} comparison failure(s), {len(missing)} missing key(s)")
        return 1
    print("\nall parity gates passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
