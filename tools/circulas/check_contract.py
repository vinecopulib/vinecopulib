"""Checks every identity of the contract's acceptance table.

Run from the repo root: ``python3 tools/circulas/check_contract.py``. Each
formula in ``circulas.py`` is compared with numerical integration and
differentiation of the density, and the symmetry, identifiability, and
dependence-measure claims of the contract are asserted. The script stops at
the first failure with a nonzero exit code.
"""
import sys
import numpy as np
import sympy as sp
from scipy import integrate, stats

import circulas as cc

TWO_PI = cc.TWO_PI
rng = np.random.default_rng(20260921)


def dblquad(f, u_hi=1.0, v_hi=1.0):
    return integrate.dblquad(lambda t, s: f(s, t), 0.0, u_hi, 0.0, v_hi, epsabs=1e-11)[0]


# --------------------------------------------------------------------------
def check_lifted_cdf(G, name):
    t = np.linspace(-7.0, 7.0, 41)
    quad = np.array([integrate.quad(G.g, 0.0, ti, epsabs=1e-12)[0] for ti in t])
    assert np.allclose(G.Gt(t), quad, atol=1e-9), f"{name}: Gt' != g"
    assert np.allclose(G.Gt(t + TWO_PI), G.Gt(t) + 1.0, atol=1e-12), f"{name}: lift"
    assert np.allclose(G.Gt(-t), -G.Gt(t), atol=1e-12), f"{name}: oddness"
    w = np.linspace(-2.3, 2.3, 47)
    assert np.allclose(G.Gt(G.Gtinv(w)), w, atol=1e-9), f"{name}: Gt(Gtinv) over several turns"
    # Ghat' = Gt
    d = 1e-5
    assert np.allclose((G.Ghat(t + d) - G.Ghat(t - d)) / (2 * d), G.Gt(t), atol=1e-7), f"{name}: Ghat' != Gt"
    print(f"  lifted CDF {name}: OK")


def check_circula(G, name):
    for q in (1, -1):
        for mu in (0.0, 0.7, np.pi, -2.5):
            pdf = lambda u, v: cc.circula_pdf(G, q, mu, u, v)
            assert abs(dblquad(pdf) - 1.0) < 1e-8, f"{name}: normalization"
            for x in (0.13, 0.5, 0.91):
                assert abs(integrate.quad(lambda t: pdf(x, t), 0, 1)[0] - 1) < 1e-9
                assert abs(integrate.quad(lambda s: pdf(s, x), 0, 1)[0] - 1) < 1e-9
            grid = np.linspace(0, 1, 9)
            assert np.allclose(pdf(0.0, grid), pdf(1.0, grid), atol=1e-12), f"{name}: periodic in u"
            assert np.allclose(pdf(grid, 0.0), pdf(grid, 1.0), atol=1e-12), f"{name}: periodic in v"
            for (u0, v0) in ((0.2, 0.35), (0.6, 0.8), (0.95, 0.1), (0.5, 0.5)):
                C = lambda a, b: dblquad(pdf, a, b)
                d = 1e-5
                assert abs((C(u0 + d, v0) - C(u0 - d, v0)) / (2 * d) - cc.circula_hfunc1(G, q, mu, u0, v0)) < 1e-6
                assert abs((C(u0, v0 + d) - C(u0, v0 - d)) / (2 * d) - cc.circula_hfunc2(G, q, mu, u0, v0)) < 1e-6
                assert abs(C(u0, v0) - cc.circula_cdf(G, q, mu, u0, v0)) < 1e-8, f"{name}: closed-form CDF"
            u = rng.uniform(size=40)
            w = rng.uniform(size=40)
            vv = cc.circula_hinv1(G, q, mu, u, w)
            assert np.all((vv >= -1e-12) & (vv <= 1 + 1e-12)), f"{name}: hinv1 in [0, 1]"
            assert np.allclose(cc.circula_hfunc1(G, q, mu, u, vv), w, atol=1e-8), f"{name}: h1(hinv1)"
            uu = cc.circula_hinv2(G, q, mu, u, w)
            assert np.all((uu >= -1e-12) & (uu <= 1 + 1e-12)), f"{name}: hinv2 in [0, 1]"
            assert np.allclose(cc.circula_hfunc2(G, q, mu, uu, u), w, atol=1e-8), f"{name}: h2(hinv2)"
            # endpoints
            assert abs(cc.circula_hfunc1(G, q, mu, 0.3, 0.0)) < 1e-12
            assert abs(cc.circula_hfunc1(G, q, mu, 0.3, 1.0) - 1) < 1e-12
            assert abs(cc.circula_hinv1(G, q, mu, np.array([0.3]), np.array([0.0]))[0]) < 1e-9
            assert abs(cc.circula_hinv1(G, q, mu, np.array([0.3]), np.array([1.0]))[0] - 1) < 1e-9
            assert abs(cc.circula_cdf(G, q, mu, 1.0, 0.4) - 0.4) < 1e-9 and abs(cc.circula_cdf(G, q, mu, 0.4, 1.0) - 0.4) < 1e-9
    print(f"  circula {name}: normalization, margins, periodicity, C, h = dC, inverses, endpoints: OK")


def check_circula_symmetries(G, name):
    U, V = rng.uniform(size=200), rng.uniform(size=200)
    mu = 0.9
    base = lambda u, v: cc.circula_pdf(G, 1, mu, u, v)
    # rotation table of the contract, in the library's convention
    assert np.allclose(base(*cc.rotate_args(90, U, V)), cc.circula_pdf(G, -1, -mu, U, V))
    assert np.allclose(base(*cc.rotate_args(180, U, V)), cc.circula_pdf(G, 1, -mu, U, V))
    assert np.allclose(base(*cc.rotate_args(270, U, V)), cc.circula_pdf(G, -1, mu, U, V))
    # flip table
    assert np.allclose(base(V, U), cc.circula_pdf(G, 1, -mu, U, V))
    assert np.allclose(cc.circula_pdf(G, -1, mu, V, U), cc.circula_pdf(G, -1, mu, U, V))
    # generic library flip (90 -> 270, mu -> -mu) equals the rotation-90 copula
    assert np.allclose(cc.circula_pdf(G, 1, -mu, *cc.rotate_args(270, U, V)), base(*cc.rotate_args(90, U, V)))
    print(f"  symmetries {name}: rotations, flip, canonicalization: OK")


def check_identifiability_and_independence():
    U, V = rng.uniform(size=20), rng.uniform(size=20)
    for Fam, par in ((cc.Cardioid, 0.3), (cc.WrappedCauchy, 0.5), (cc.VonMises, 2.0)):
        assert np.allclose(cc.circula_pdf(Fam(-par), 1, 0.4, U, V), cc.circula_pdf(Fam(par), 1, 0.4 + np.pi, U, V))
        assert np.allclose(cc.circula_pdf(Fam(0.0), 1, 0.4, U, V), 1.0)
        assert np.allclose(cc.circula_pdf(Fam(par), 1, 0.4, U, V), cc.circula_pdf(Fam(par), 1, 0.4 + TWO_PI, U, V))
    assert np.allclose(cc.sections_pdf(0.4, -0.7, 1.0, U, V), cc.sections_pdf(-0.4, 0.7, 1.0 + np.pi, U, V))
    assert np.allclose(cc.sections_pdf(0.0, 0.0, 1.0, U, V), 1.0)
    print("  identifiability (negative concentration == phase + pi), independence, 2 pi periodicity of mu: OK")


def check_half_turn_and_moments():
    G = cc.WrappedCauchy(0.95)
    u, v = cc.circula_simulate(G, 1, np.pi, 200000, rng)
    tau, rho = stats.kendalltau(u, v)[0], stats.spearmanr(u, v)[0]
    beta = 4 * np.mean((u <= 0.5) & (v <= 0.5)) - 1
    # values at rho = 0.95; the limits as rho -> 1 are 0, -1/2, and -1
    assert abs(tau + 0.054) < 0.03 and abs(rho + 0.479) < 0.03 and abs(beta + 0.902) < 0.03, (tau, rho, beta)
    print(f"  half-turn case (rho = 0.95, mu = pi): tau = {tau:+.3f}, rho_S = {rho:+.3f}, beta = {beta:+.3f}: OK")
    for Fam, par in ((cc.Cardioid, 0.4), (cc.WrappedCauchy, 0.7), (cc.VonMises, 3.0)):
        G = Fam(par)
        for q, mu in ((1, 2.0), (-1, -2.9)):
            u, v = cc.circula_simulate(G, q, mu, 100000, rng)
            mu_hat, rbar = cc.circula_moment_start(u, v, q)
            assert abs(np.angle(np.exp(1j * (mu_hat - mu)))) < 0.03, (Fam.__name__, mu_hat, mu)
            assert abs(rbar - G.mean_resultant_length()) < 0.01, (Fam.__name__, rbar)
    kappa = cc.von_mises_kappa_from_resultant(cc.VonMises(3.0).mean_resultant_length())
    assert abs(kappa - 3.0) < 1e-6
    print("  moment initialization of the binding families (phase, mean resultant length, A^-1): OK")


def check_sections():
    for (a, b, mu) in ((0.7, 0.7, 0.4), (1.0, -1.0, 2.0), (0.3, -0.9, -1.1), (0.0, 1.0, np.pi)):
        pdf = lambda u, v: cc.sections_pdf(a, b, mu, u, v)
        assert abs(dblquad(pdf) - 1.0) < 1e-9
        for x in (0.1, 0.5, 0.93):
            assert abs(integrate.quad(lambda t: pdf(x, t), 0, 1)[0] - 1) < 1e-10
            assert abs(integrate.quad(lambda s: pdf(s, x), 0, 1)[0] - 1) < 1e-10
        grid = np.linspace(0, 1, 11)
        assert np.allclose(pdf(0.0, grid), pdf(1.0, grid))
        if a != -b:
            assert not np.allclose(pdf(grid, 0.0), pdf(grid, 1.0)), "should not be periodic in v"
        U, V, W = rng.uniform(size=30), rng.uniform(size=30), rng.uniform(size=30)
        Cn = np.array([dblquad(pdf, u0, v0) for u0, v0 in zip(U, V)])
        assert np.allclose(cc.sections_cdf(a, b, mu, U, V), Cn, atol=1e-9)
        d = 1e-6
        C = lambda u, v: cc.sections_cdf(a, b, mu, u, v)
        assert np.allclose((C(U + d, V) - C(U - d, V)) / (2 * d), cc.sections_hfunc1(a, b, mu, U, V), atol=1e-7)
        assert np.allclose((C(U, V + d) - C(U, V - d)) / (2 * d), cc.sections_hfunc2(a, b, mu, U, V), atol=1e-7)
        assert np.all(pdf(rng.uniform(size=5000), rng.uniform(size=5000)) >= -1e-12)
        assert np.allclose(C(1.0, V), V) and np.allclose(C(U, 1.0), U) and np.allclose(C(0.0, V), 0) and np.allclose(C(U, 0.0), 0)
        vv = cc.sections_hinv1(a, b, mu, U, W)
        assert np.allclose(cc.sections_hfunc1(a, b, mu, U, vv), W, atol=1e-10)
        uu = cc.sections_hinv2(a, b, mu, V, W)
        assert np.allclose(cc.sections_hfunc2(a, b, mu, uu, V), W, atol=1e-10)
    print("  sections: normalization, margins, periodicity in u only, C, h = dC, inverses, endpoints: OK")

    # feasible region is the box
    v = np.linspace(0, 1, 20001)
    for a, b in ((1, 1), (1, -1), (-1, 1), (-1, -1), (1, 0), (0, 1), (0.5, -1), (1, 0.3)):
        assert np.max(np.abs(cc.sections_p(a, b, v))) <= 1 + 1e-12
    for a, b in ((1.001, 0), (0, -1.001), (1.01, 1.01)):
        assert np.max(np.abs(cc.sections_p(a, b, v))) > 1
    print("  sections: feasible region is exactly |a| <= 1, |b| <= 1: OK")

    # symmetry table
    U, V = rng.uniform(size=50), rng.uniform(size=50)
    a, b, mu = 0.5, -0.2, 0.8
    f = lambda a_, b_, mu_, u, v: cc.sections_pdf(a_, b_, mu_, u, v)
    assert np.allclose(f(a, b, mu, 1 - U, V), f(a, b, -mu, U, V))
    assert np.allclose(f(a, b, mu, U, 1 - V), f(-b, -a, mu, U, V))
    assert np.allclose(f(a, b, mu, U, 1 - V), f(b, a, mu + np.pi, U, V))
    assert np.allclose(f(a, b, mu, 1 - U, 1 - V), f(b, a, np.pi - mu, U, V))
    assert np.allclose(f(a, a, mu, 1 - U, 1 - V), f(a, a, np.pi - mu, U, V))
    assert np.allclose(f(a, a, mu, U, 1 - V), f(a, a, mu + np.pi, U, V))
    assert np.allclose(f(0.6, 0.6, mu, U, V), 1 + 0.6 * np.cos(TWO_PI * U - mu) * (1 - 2 * V))
    print("  sections: symmetry table, quadratic sections at a == b: OK")

    # closed-form dependence measures (symbolic)
    u_, v_, a_, b_, mu_ = sp.symbols("u v a b mu", real=True)
    P = a_ * v_ * (1 - v_) ** 2 + b_ * v_ ** 2 * (1 - v_)
    C = u_ * v_ + (sp.sin(2 * sp.pi * u_ - mu_) + sp.sin(mu_)) / (2 * sp.pi) * P
    c = sp.diff(C, u_, v_)
    rho_s = sp.simplify(12 * sp.integrate(sp.integrate(C, (u_, 0, 1)), (v_, 0, 1)) - 3)
    tau = sp.simplify(4 * sp.integrate(sp.integrate(sp.expand(C * c), (v_, 0, 1)), (u_, 0, 1)) - 1)
    assert sp.simplify(rho_s - (a_ + b_) * sp.sin(mu_) / (2 * sp.pi)) == 0
    assert sp.simplify(tau - (a_ + b_) * sp.sin(mu_) / (3 * sp.pi)) == 0
    print(f"  sections: Kendall tau = {tau}, Spearman rho = {rho_s}: OK")

    # moment initialization
    for a, b, mu, cubic in ((0.8, 0.8, 1.2, False), (0.9, -0.5, -2.0, True), (0.2, 0.9, 0.3, True)):
        U = rng.uniform(size=200000)
        V = cc.sections_hinv1(a, b, mu, U[:20000], rng.uniform(size=20000))  # inverse is slow; 2e4 suffices
        a_hat, b_hat, mu_hat = cc.sections_moment_start(U[:20000], V, cubic=cubic)
        assert abs(a_hat - a) < 0.08 and abs(b_hat - b) < 0.08, (a_hat, b_hat)
        assert abs(np.angle(np.exp(1j * (mu_hat - mu)))) < 0.1, (mu_hat, mu)
    print("  sections: moment initialization recovers (a, b, mu): OK")


def main():
    print("binding-density circulas")
    for G, name in ((cc.Cardioid(0.45), "cardioid rho=0.45"),
                    (cc.WrappedCauchy(0.8), "wrapped Cauchy rho=0.8"),
                    (cc.VonMises(6.0), "von Mises kappa=6")):
        check_lifted_cdf(G, name)
        check_circula(G, name)
        check_circula_symmetries(G, name)
    check_lifted_cdf(cc.VonMises(100.0), "von Mises kappa=100 (upper bound)")
    check_lifted_cdf(cc.WrappedCauchy(0.99), "wrapped Cauchy rho=0.99 (upper bound)")
    check_lifted_cdf(cc.Cardioid(0.5), "cardioid rho=0.5 (upper bound)")
    check_identifiability_and_independence()
    check_half_turn_and_moments()
    print("cylindrical sections copulas")
    check_sections()
    print("all contract identities hold")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as err:
        print(f"FAILED: {err}", file=sys.stderr)
        sys.exit(1)
