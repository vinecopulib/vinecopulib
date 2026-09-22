"""Reference implementation of the circular pair-copula formulas.

This module mirrors ``docs/plans/circular-copulas-contract.md`` one function
per formula, in NumPy, with no attempt at speed. It exists so that the C++
implementation can be tested against numbers that were not produced by it.
Nothing here ships in the library.

Notation follows the contract: ``q`` is the orientation (the copula rotation
in the library, ``q = 1`` for rotation 0), ``mu`` the phase in radians, and
``Gt`` the lifted distribution function of the circular density ``g``.
"""
import numpy as np
from scipy import optimize, special

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# circular densities g, symmetric about zero, with lifted CDFs Gt (Gt(0) = 0)
# ---------------------------------------------------------------------------
class Cardioid:
    """g(t) = (1 + 2 rho cos t) / (2 pi), rho in [0, 1/2]."""

    lower, upper = 0.0, 0.5

    def __init__(self, rho):
        self.rho = rho

    def g(self, t):
        return (1.0 + 2.0 * self.rho * np.cos(t)) / TWO_PI

    def Gt(self, t):
        return (t + 2.0 * self.rho * np.sin(t)) / TWO_PI

    def Ghat(self, t):
        """An antiderivative of Gt; only the differences in C(u, v) matter."""
        return (t ** 2 / 2.0 - 2.0 * self.rho * np.cos(t)) / TWO_PI

    def Gtinv(self, w):
        w = np.atleast_1d(np.asarray(w, float))
        out = np.empty_like(w)
        for i, wi in enumerate(w):
            lo, hi = TWO_PI * wi - 2.0 * self.rho - 1e-12, TWO_PI * wi + 2.0 * self.rho + 1e-12
            out[i] = optimize.brentq(lambda t: self.Gt(t) - wi, lo, hi, xtol=1e-14)
        return out

    def mean_resultant_length(self):
        return self.rho


class WrappedCauchy:
    """g(t) = (1 - rho^2) / (2 pi (1 + rho^2 - 2 rho cos t)), rho in [0, 0.99]."""

    lower, upper = 0.0, 0.99

    def __init__(self, rho):
        self.rho = rho

    def g(self, t):
        r = self.rho
        return (1.0 - r ** 2) / (TWO_PI * (1.0 + r ** 2 - 2.0 * r * np.cos(t)))

    def Gt(self, t):
        r = self.rho
        return (t + 2.0 * np.arctan(r * np.sin(t) / (1.0 - r * np.cos(t)))) / TWO_PI

    def Ghat(self, t, nterms=4000):
        r = self.rho
        k = np.arange(1, nterms + 1)
        t = np.asarray(t, float)
        series = np.sum((r ** k)[:, None] * np.cos(np.outer(k, t.ravel())) / (k ** 2)[:, None], axis=0)
        return (t ** 2 / 2.0 - 2.0 * series.reshape(t.shape)) / TWO_PI

    def Gtinv(self, w):
        w = np.asarray(w, float)
        k = np.round(w)
        return TWO_PI * k + 2.0 * np.arctan((1.0 - self.rho) / (1.0 + self.rho) * np.tan(np.pi * (w - k)))

    def mean_resultant_length(self):
        return self.rho


class VonMises:
    """g(t) = exp(kappa cos t) / (2 pi I_0(kappa)), kappa in [0, 100]."""

    lower, upper = 0.0, 100.0

    def __init__(self, kappa, nterms=None):
        self.kappa = kappa
        # terms decay like exp(-j^2 / (2 kappa)); this reaches 1e-12 with margin
        self.nterms = nterms or int(np.ceil(np.sqrt(2.0 * max(kappa, 1.0) * 30.0))) + 10

    def g(self, t):
        return np.exp(self.kappa * np.cos(t)) / (TWO_PI * special.i0(self.kappa))

    def _ratios(self):
        j = np.arange(1, self.nterms + 1)
        # I_j / I_0 through the exponentially scaled Bessel functions
        return j, special.ive(j, self.kappa) / special.ive(0, self.kappa)

    def Gt(self, t):
        j, ratio = self._ratios()
        t = np.asarray(t, float)
        series = np.sum(ratio[:, None] * np.sin(np.outer(j, t.ravel())) / j[:, None], axis=0)
        return (t + 2.0 * series.reshape(t.shape)) / TWO_PI

    def Ghat(self, t):
        j, ratio = self._ratios()
        t = np.asarray(t, float)
        series = np.sum(ratio[:, None] * (1.0 - np.cos(np.outer(j, t.ravel()))) / (j ** 2)[:, None], axis=0)
        return (t ** 2 / 2.0 + 2.0 * series.reshape(t.shape)) / TWO_PI

    def Gtinv(self, w):
        w = np.atleast_1d(np.asarray(w, float))
        out = np.empty_like(w)
        for i, wi in enumerate(w):
            lo, hi = TWO_PI * wi - np.pi, TWO_PI * wi + np.pi
            out[i] = optimize.brentq(lambda t: self.Gt(t) - wi, lo, hi, xtol=1e-14)
        return out

    def mean_resultant_length(self):
        return special.i1e(self.kappa) / special.i0e(self.kappa)


def von_mises_kappa_from_resultant(rbar):
    """Inverse of A(kappa) = I_1(kappa) / I_0(kappa), clipped to the domain."""
    rbar = min(max(rbar, 0.0), 0.999999)
    if rbar == 0.0:
        return 0.0
    f = lambda k: special.i1e(k) / special.i0e(k) - rbar
    if f(VonMises.upper) < 0:
        return VonMises.upper
    return optimize.brentq(f, 1e-12, VonMises.upper)


FAMILIES = {"cardioid": Cardioid, "wrapped_cauchy": WrappedCauchy, "von_mises": VonMises}


# ---------------------------------------------------------------------------
# binding-density circula
# ---------------------------------------------------------------------------
def circula_pdf(G, q, mu, u, v):
    return TWO_PI * G.g(TWO_PI * (v - q * u) - mu)


def circula_cdf(G, q, mu, u, v):
    gh = G.Ghat
    return (q / TWO_PI) * (gh(TWO_PI * v - mu) - gh(TWO_PI * (v - q * u) - mu)
                           - gh(-mu) + gh(-TWO_PI * q * u - mu))


def circula_hfunc1(G, q, mu, u, v):
    """h_1(v | u) = dC/du."""
    return G.Gt(TWO_PI * (v - q * u) - mu) - G.Gt(-TWO_PI * q * u - mu)


def circula_hfunc2(G, q, mu, u, v):
    """h_2(u | v) = dC/dv."""
    return q * (G.Gt(TWO_PI * v - mu) - G.Gt(TWO_PI * (v - q * u) - mu))


def circula_hinv1(G, q, mu, u, w):
    return q * u + (mu + G.Gtinv(w + G.Gt(-TWO_PI * q * u - mu))) / TWO_PI


def circula_hinv2(G, q, mu, v, w):
    return q * (v - (mu + G.Gtinv(G.Gt(TWO_PI * v - mu) - q * w)) / TWO_PI)


def circula_simulate(G, q, mu, n, rng):
    u = rng.uniform(size=n)
    v = circula_hinv1(G, q, mu, u, rng.uniform(size=n))
    return u, v


def circula_moment_start(u, v, q, weights=None):
    """Phase and mean resultant length of the residual angles 2 pi (v - q u)."""
    theta = TWO_PI * (v - q * u)
    w = np.ones_like(u) if weights is None else weights
    z = np.sum(w * np.exp(1j * theta)) / np.sum(w)
    return np.angle(z), np.abs(z)


# rotations as the library applies them: c_rot(u, v) = c_0(rotated arguments)
def rotate_args(rotation, u, v):
    if rotation == 0:
        return u, v
    if rotation == 90:
        return v, 1.0 - u
    if rotation == 180:
        return 1.0 - u, 1.0 - v
    if rotation == 270:
        return 1.0 - v, u
    raise ValueError(rotation)


# ---------------------------------------------------------------------------
# cylindrical sections copulas: circular u (first argument), linear v
# ---------------------------------------------------------------------------
def sections_p(a, b, v):
    """p(v) = a (1 - v)(1 - 3v) + b v (2 - 3v)."""
    return a * (1.0 - v) * (1.0 - 3.0 * v) + b * v * (2.0 - 3.0 * v)


def sections_P(a, b, v):
    """Antiderivative of p with P(0) = P(1) = 0."""
    return a * v * (1.0 - v) ** 2 + b * v ** 2 * (1.0 - v)


def sections_pdf(a, b, mu, u, v):
    return 1.0 + np.cos(TWO_PI * u - mu) * sections_p(a, b, v)


def sections_cdf(a, b, mu, u, v):
    return u * v + (np.sin(TWO_PI * u - mu) + np.sin(mu)) / TWO_PI * sections_P(a, b, v)


def sections_hfunc1(a, b, mu, u, v):
    """h_1(v | u) = dC/du."""
    return v + np.cos(TWO_PI * u - mu) * sections_P(a, b, v)


def sections_hfunc2(a, b, mu, u, v):
    """h_2(u | v) = dC/dv."""
    return u + (np.sin(TWO_PI * u - mu) + np.sin(mu)) / TWO_PI * sections_p(a, b, v)


def _bracketed_root(f, lo=0.0, hi=1.0):
    flo, fhi = f(lo), f(hi)
    if flo >= 0:
        return lo
    if fhi <= 0:
        return hi
    return optimize.brentq(f, lo, hi, xtol=1e-14)


def sections_hinv1(a, b, mu, u, w):
    u, w = np.broadcast_arrays(np.asarray(u, float), np.asarray(w, float))
    return np.array([_bracketed_root(lambda t: sections_hfunc1(a, b, mu, ui, t) - wi)
                     for ui, wi in zip(u.ravel(), w.ravel())]).reshape(u.shape)


def sections_hinv2(a, b, mu, v, w):
    v, w = np.broadcast_arrays(np.asarray(v, float), np.asarray(w, float))
    return np.array([_bracketed_root(lambda s: sections_hfunc2(a, b, mu, s, vi) - wi)
                     for vi, wi in zip(v.ravel(), w.ravel())]).reshape(v.shape)


def sections_kendall_tau(a, b, mu):
    return (a + b) * np.sin(mu) / (3.0 * np.pi)


def sections_spearman_rho(a, b, mu):
    return (a + b) * np.sin(mu) / TWO_PI


def sections_moment_start(u, v, weights=None, cubic=True):
    """Moment estimates (a, b, mu) from E[exp(i 2 pi U) p_k(V)]."""
    w = np.ones_like(u) if weights is None else weights
    z = np.exp(1j * TWO_PI * u)
    if not cubic:
        m = np.sum(w * z * (1.0 - 2.0 * v)) / np.sum(w)
        a = min(6.0 * np.abs(m), 1.0)
        return a, a, np.angle(m)
    m1 = np.sum(w * z * (1.0 - v) * (1.0 - 3.0 * v)) / np.sum(w)
    m2 = np.sum(w * z * v * (2.0 - 3.0 * v)) / np.sum(w)
    mu = np.angle(m1 + m2)
    rhs = np.array([np.real(m1 * np.exp(-1j * mu)), np.real(m2 * np.exp(-1j * mu))])
    gram = np.array([[1.0 / 15.0, 1.0 / 60.0], [1.0 / 60.0, 1.0 / 15.0]])
    a, b = np.linalg.solve(gram, rhs)
    if a < 0:  # (a, b, mu) ~ (-a, -b, mu + pi)
        a, b, mu = -a, -b, mu + np.pi
    return min(a, 1.0), min(max(b, -1.0), 1.0), np.angle(np.exp(1j * mu))
