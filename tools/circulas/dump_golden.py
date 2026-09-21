"""Dumps reference values of the circular families for the C++ golden tests.

Run from the repo root: ``python3 tools/circulas/dump_golden.py OUT.json``.
For each family, parameter point, and rotation, the file records the density,
CDF, both h-functions, and both inverse h-functions on a fixed lattice of
arguments, plus Kendall's tau by numerical integration. The parameter lattice
covers interiors; near-boundary points carry an ``_edge`` suffix, following
the convention of ``scripts/README.md``.
"""
import json
import sys
import numpy as np
from scipy import integrate

import circulas as cc

ARGS = np.array([0.02, 0.2, 0.5, 0.8, 0.98])
U, V = [x.ravel() for x in np.meshgrid(ARGS, ARGS, indexing="ij")]
W = np.array([0.05, 0.3, 0.5, 0.7, 0.95])

BINDING_PARAMS = {
    "cardioid": {"interior": [(0.2, 0.0), (0.4, 1.0), (0.3, -2.5)], "edge": [(0.5, 3.0), (0.0, 0.0)]},
    "wrapped_cauchy": {"interior": [(0.3, 0.0), (0.7, 1.0), (0.5, np.pi)], "edge": [(0.99, -3.0), (0.0, 0.0)]},
    "von_mises": {"interior": [(1.0, 0.0), (4.0, 1.0), (2.0, np.pi)], "edge": [(100.0, 0.5), (0.0, 0.0)]},
}
SECTIONS_PARAMS = {
    "quad_sections": {"interior": [(0.5, 0.5, 0.0), (0.8, 0.8, 1.0), (0.3, 0.3, -2.0)], "edge": [(1.0, 1.0, 3.0), (0.0, 0.0, 0.0)]},
    "cubic_sections": {"interior": [(0.5, -0.3, 0.0), (0.8, 0.6, 1.0), (0.3, -0.9, -2.0)], "edge": [(1.0, -1.0, 3.0), (1.0, 1.0, -1.0), (0.0, 0.0, 0.0)]},
}


def kendall_tau_numeric(h1, h2):
    """tau = 1 - 4 int int h1(v | u) h2(u | v) du dv."""
    val = integrate.dblquad(lambda v, u: h1(u, v) * h2(u, v), 0, 1, 0, 1, epsabs=1e-9)[0]
    return 1.0 - 4.0 * val


def binding_entry(Fam, par, mu, rotation):
    G = Fam(par)
    q, phase = (1, mu) if rotation == 0 else (-1, -mu)  # rotation 90 is q = -1 with phase -mu
    entry = {
        "parameters": [par, mu], "rotation": rotation, "u": U.tolist(), "v": V.tolist(), "w": W.tolist(),
        "pdf": cc.circula_pdf(G, q, phase, U, V).tolist(),
        "cdf": cc.circula_cdf(G, q, phase, U, V).tolist(),
        "hfunc1": cc.circula_hfunc1(G, q, phase, U, V).tolist(),
        "hfunc2": cc.circula_hfunc2(G, q, phase, U, V).tolist(),
        "hinv1": [cc.circula_hinv1(G, q, phase, U, w).tolist() for w in W],
        "hinv2": [cc.circula_hinv2(G, q, phase, V, w).tolist() for w in W],
        "tau": kendall_tau_numeric(lambda u, v: cc.circula_hfunc1(G, q, phase, u, v),
                                   lambda u, v: cc.circula_hfunc2(G, q, phase, u, v)),
    }
    return entry


def sections_entry(a, b, mu):
    return {
        "parameters": [a, mu] if a == b else [a, b, mu], "rotation": 0,
        "u": U.tolist(), "v": V.tolist(), "w": W.tolist(),
        "pdf": cc.sections_pdf(a, b, mu, U, V).tolist(),
        "cdf": cc.sections_cdf(a, b, mu, U, V).tolist(),
        "hfunc1": cc.sections_hfunc1(a, b, mu, U, V).tolist(),
        "hfunc2": cc.sections_hfunc2(a, b, mu, U, V).tolist(),
        "hinv1": [cc.sections_hinv1(a, b, mu, U, w).tolist() for w in W],
        "hinv2": [cc.sections_hinv2(a, b, mu, V, w).tolist() for w in W],
        "tau": cc.sections_kendall_tau(a, b, mu),
        "spearman": cc.sections_spearman_rho(a, b, mu),
    }


def main(path):
    out = {"lattice_note": "u and v are the flattened 5 x 5 product of ARGS; hinv rows correspond to w"}
    for name, groups in BINDING_PARAMS.items():
        Fam = cc.FAMILIES[name]
        for suffix, points in (("", groups["interior"]), ("_edge", groups["edge"])):
            for rotation in (0, 90):
                out[f"{name}{suffix}_rot{rotation}"] = [binding_entry(Fam, par, mu, rotation) for par, mu in points]
    for name, groups in SECTIONS_PARAMS.items():
        for suffix, points in (("", groups["interior"]), ("_edge", groups["edge"])):
            out[f"{name}{suffix}"] = [sections_entry(a, b, mu) for a, b, mu in points]
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {path}: {sum(len(v) for v in out.values() if isinstance(v, list))} parameter points")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    main(sys.argv[1])
