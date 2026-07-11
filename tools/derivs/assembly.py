"""Symbolic derivative "assembly" for the parametric bicop families.

This is the family-independent derivative calculus of the BB/Tawn derivation
note (``bb_tawn_analytic_derivatives.tex``), ported to sympy. It is used only at
code-generation time (see ``generate.py``): each requested derivative selector is
expanded symbolically in terms of the generator/Pickands derivatives, with the
copula value ``w = C(u, v)`` kept as a free runtime symbol, then CSE'd and
emitted as a fast C++ closed form. Nothing here ships in the library.

Coordinate indices ``xi``: Archimedean 0=u, 1=v, 2=theta, 3=delta; Tawn
0=u, 1=v, 2=psi1, 3=psi2, 4=theta.
"""
import sympy as sp

th, de = sp.symbols('theta delta', positive=True)
u, v, w, z = sp.symbols('u v w z', positive=True)


def archimedean(phi_of_z):
    """Assembly for a strict Archimedean copula with generator ``phi_of_z``.

    Returns a dict of callables ``c_d1(xi)``, ``c_d2(xi, xj)``, ``h1_d1(xi)``,
    ``h1_d2(xi, xj)``, ``logc_d1(xi)``, ``logc_d2(xi, xj)`` giving the requested
    derivative as a sympy expression in (u, v, w, theta, delta).
    """
    dcache = {}

    def dz(a, b, c):
        key = (a, b, c)
        if key not in dcache:
            dcache[key] = sp.diff(phi_of_z, z, a, th, b, de, c)
        return dcache[key]

    def P(pt, a, b, c):
        return dz(a, b, c).subs(z, pt)

    def phi_field(pt):
        return dict(fz=P(pt, 1, 0, 0), fzz=P(pt, 2, 0, 0), fth=P(pt, 0, 1, 0),
                    fde=P(pt, 0, 0, 1), fzth=P(pt, 1, 1, 0), fzde=P(pt, 1, 0, 1),
                    fthth=P(pt, 0, 2, 0), fthde=P(pt, 0, 1, 1), fdede=P(pt, 0, 0, 2))

    def G_field(pt):
        g, r, zzz = P(pt, 1, 0, 0), P(pt, 2, 0, 0), P(pt, 3, 0, 0)
        zth, zde = P(pt, 1, 1, 0), P(pt, 1, 0, 1)
        zzth, zzde = P(pt, 2, 1, 0), P(pt, 2, 0, 1)
        zthth, zthde, zdede = P(pt, 1, 2, 0), P(pt, 1, 1, 1), P(pt, 1, 0, 2)
        return dict(fz=r / g, fzz=zzz / g - (r / g)**2, fth=zth / g, fde=zde / g,
                    fzth=zzth / g - r * zth / g**2, fzde=zzde / g - r * zde / g**2,
                    fthth=zthth / g - zth**2 / g**2, fthde=zthde / g - zth * zde / g**2,
                    fdede=zdede / g - zde**2 / g**2)

    def R_field(pt):
        r, zzz, zzzz = P(pt, 2, 0, 0), P(pt, 3, 0, 0), P(pt, 4, 0, 0)
        zzth, zzde = P(pt, 2, 1, 0), P(pt, 2, 0, 1)
        zzzth, zzzde = P(pt, 3, 1, 0), P(pt, 3, 0, 1)
        zzthth, zzthde, zzdede = P(pt, 2, 2, 0), P(pt, 2, 1, 1), P(pt, 2, 0, 2)
        return dict(fz=zzz / r, fzz=zzzz / r - (zzz / r)**2, fth=zzth / r, fde=zzde / r,
                    fzth=zzzth / r - zzz * zzth / r**2, fzde=zzzde / r - zzz * zzde / r**2,
                    fthth=zzthth / r - zzth**2 / r**2, fthde=zzthde / r - zzth * zzde / r**2,
                    fdede=zzdede / r - zzde**2 / r**2)

    def eth(xi): return 1 if xi == 2 else 0
    def ede(xi): return 1 if xi == 3 else 0
    def Ti(F, zi, xi): return F['fz'] * zi + F['fth'] * eth(xi) + F['fde'] * ede(xi)

    def Tij(F, zi, zj, zij, xi, xj):
        return (F['fzz'] * zi * zj + F['fz'] * zij
                + F['fzth'] * (zi * eth(xj) + zj * eth(xi))
                + F['fzde'] * (zi * ede(xj) + zj * ede(xi)) + F['fthth'] * eth(xi) * eth(xj)
                + F['fthde'] * (eth(xi) * ede(xj) + eth(xj) * ede(xi)) + F['fdede'] * ede(xi) * ede(xj))

    PHu, PHv, PHw = phi_field(u), phi_field(v), phi_field(w)
    Gu, Gv, Gw, Rw = G_field(u), G_field(v), G_field(w), R_field(w)
    gw = P(w, 1, 0, 0)

    def zu(xi): return 1 if xi == 0 else 0
    def zv(xi): return 1 if xi == 1 else 0

    def w_i(xi):
        num = Ti(PHu, zu(xi), xi) + Ti(PHv, zv(xi), xi) - (P(w, 0, 1, 0) * eth(xi) + P(w, 0, 0, 1) * ede(xi))
        return num / gw

    def w_ij(xi, xj, wi, wj):
        Tu = Tij(PHu, zu(xi), zu(xj), 0, xi, xj)
        Tv = Tij(PHv, zv(xi), zv(xj), 0, xi, xj)
        Rij = (PHw['fzz'] * wi * wj + PHw['fzth'] * (wi * eth(xj) + wj * eth(xi))
               + PHw['fzde'] * (wi * ede(xj) + wj * ede(xi)) + PHw['fthth'] * eth(xi) * eth(xj)
               + PHw['fthde'] * (eth(xi) * ede(xj) + eth(xj) * ede(xi)) + PHw['fdede'] * ede(xi) * ede(xj))
        return (Tu + Tv - Rij) / gw

    def Li(xi):
        wi = w_i(xi)
        return Ti(Gu, zu(xi), xi) + Ti(Gv, zv(xi), xi) + Ti(Rw, wi, xi) - 3 * Ti(Gw, wi, xi)

    def H1i(xi):
        wi = w_i(xi)
        return Ti(Gu, zu(xi), xi) - Ti(Gw, wi, xi)

    def Lij(xi, xj):
        wi, wj = w_i(xi), w_i(xj)
        wij = w_ij(xi, xj, wi, wj)
        return (Tij(Gu, zu(xi), zu(xj), 0, xi, xj) + Tij(Gv, zv(xi), zv(xj), 0, xi, xj)
                + Tij(Rw, wi, wj, wij, xi, xj) - 3 * Tij(Gw, wi, wj, wij, xi, xj))

    def H1ij(xi, xj):
        wi, wj = w_i(xi), w_i(xj)
        wij = w_ij(xi, xj, wi, wj)
        return Tij(Gu, zu(xi), zu(xj), 0, xi, xj) - Tij(Gw, wi, wj, wij, xi, xj)

    cval = -PHu['fz'] * PHv['fz'] * PHw['fzz'] / gw**3
    h1val = PHu['fz'] / gw
    return dict(
        vars=(u, v, w, th, de),
        c_d1=lambda xi: cval * Li(xi),
        c_d2=lambda xi, xj: cval * (Lij(xi, xj) + Li(xi) * Li(xj)),
        h1_d1=lambda xi: h1val * H1i(xi),
        h1_d2=lambda xi, xj: h1val * (H1ij(xi, xj) + H1i(xi) * H1i(xj)),
        logc_d1=lambda xi: Li(xi),
        logc_d2=lambda xi, xj: Lij(xi, xj),
    )


# ---------------------------------------------------------------------------
# Extreme-value assembly (Tawn). Coordinate indices: 0=u, 1=v, 2=psi1, 3=psi2,
# 4=theta. The Pickands argument t and r = -log(uv) are kept as runtime symbols
# (computed in the leaf preamble); A_of_t is the Pickands function A(t; params).
# ---------------------------------------------------------------------------
r_s, t_s, p1, p2, pth = sp.symbols('r_s t_s psi1 psi2 theta_ev', positive=True)


def extreme_value(A_of_t):
    """Assembly for an extreme-value copula C = exp(-r A(t)).

    ``A_of_t`` is A(t_s; p1, p2, pth). Returns callables c_d1/c_d2/h1_d1/h1_d2/
    h2_d1/h2_d2/logc_d1/logc_d2 as sympy expressions in (u, v, r_s, t_s, p1, p2,
    pth).
    """
    PIDX = [2, 3, 4]
    PU = {2: p1, 3: p2, 4: pth}

    def A(k):
        return sp.diff(A_of_t, t_s, k)

    def Aa(k, a):
        return sp.diff(sp.diff(A_of_t, t_s, k), PU[a])

    def Aab(k, a, b):
        return sp.diff(sp.diff(sp.diff(A_of_t, t_s, k), PU[a]), PU[b])

    d_ = t_s * (1 - t_s)
    dp = 1 - 2 * t_s
    dpp = sp.Integer(-2)

    E = {'val': -r_s * A(0), 'r': -A(0), 't': -r_s * A(1), 'rr': sp.Integer(0),
         'rt': -A(1), 'tt': -r_s * A(2), 'a': lambda a: -r_s * Aa(0, a),
         'ra': lambda a: -Aa(0, a), 'ta': lambda a: -r_s * Aa(1, a),
         'ab': lambda a, b: -r_s * Aab(0, a, b)}
    P = {'val': A(0) - t_s * A(1), 'r': sp.Integer(0), 't': -t_s * A(2),
         'rr': sp.Integer(0), 'rt': sp.Integer(0), 'tt': -A(2) - t_s * A(3),
         'a': lambda a: Aa(0, a) - t_s * Aa(1, a), 'ra': lambda a: sp.Integer(0),
         'ta': lambda a: -t_s * Aa(2, a), 'ab': lambda a, b: Aab(0, a, b) - t_s * Aab(1, a, b)}
    Q = {'val': A(0) + (1 - t_s) * A(1), 'r': sp.Integer(0), 't': (1 - t_s) * A(2),
         'rr': sp.Integer(0), 'rt': sp.Integer(0), 'tt': -A(2) + (1 - t_s) * A(3),
         'a': lambda a: Aa(0, a) + (1 - t_s) * Aa(1, a), 'ra': lambda a: sp.Integer(0),
         'ta': lambda a: (1 - t_s) * Aa(2, a), 'ab': lambda a, b: Aab(0, a, b) + (1 - t_s) * Aab(1, a, b)}
    Pv, Qv = P['val'], Q['val']
    B = {'val': Pv * Qv + d_ * A(2) / r_s, 'r': -d_ * A(2) / r_s**2,
         't': P['t'] * Qv + Pv * Q['t'] + (dp * A(2) + d_ * A(3)) / r_s,
         'rr': 2 * d_ * A(2) / r_s**3, 'rt': -(dp * A(2) + d_ * A(3)) / r_s**2,
         'tt': P['tt'] * Qv + 2 * P['t'] * Q['t'] + Pv * Q['tt'] + (dpp * A(2) + 2 * dp * A(3) + d_ * A(4)) / r_s,
         'a': lambda a: P['a'](a) * Qv + Pv * Q['a'](a) + d_ * Aa(2, a) / r_s,
         'ra': lambda a: -d_ * Aa(2, a) / r_s**2,
         'ta': lambda a: P['ta'](a) * Qv + P['t'] * Q['a'](a) + P['a'](a) * Q['t'] + Pv * Q['ta'](a) + (dp * Aa(2, a) + d_ * Aa(3, a)) / r_s,
         'ab': lambda a, b: P['ab'](a, b) * Qv + P['a'](a) * Q['a'](b) + P['a'](b) * Q['a'](a) + Pv * Q['ab'](a, b) + d_ * Aab(2, a, b) / r_s}

    def a_i(xi): return -(1 if xi == 0 else 0) / u
    def b_i(xi): return -(1 if xi == 1 else 0) / v
    def a_ij(xi, xj): return (1 if (xi == 0 and xj == 0) else 0) / u**2
    def b_ij(xi, xj): return (1 if (xi == 1 and xj == 1) else 0) / v**2
    def r_i(xi): return a_i(xi) + b_i(xi)
    def r_ij(xi, xj): return a_ij(xi, xj) + b_ij(xi, xj)
    def t_i(xi): return (b_i(xi) - t_s * r_i(xi)) / r_s
    def t_ij(xi, xj):
        return (b_ij(xi, xj) - t_i(xj) * r_i(xi) - t_i(xi) * r_i(xj) - t_s * r_ij(xi, xj)) / r_s
    def e(xi, a): return 1 if xi == a else 0

    def DiF(F, xi):
        return F['r'] * r_i(xi) + F['t'] * t_i(xi) + sum(F['a'](a) * e(xi, a) for a in PIDX)

    def DijF(F, xi, xj):
        out = (F['rr'] * r_i(xi) * r_i(xj) + F['rt'] * (r_i(xi) * t_i(xj) + t_i(xi) * r_i(xj))
               + F['tt'] * t_i(xi) * t_i(xj) + F['r'] * r_ij(xi, xj) + F['t'] * t_ij(xi, xj))
        out += sum(F['ra'](a) * (r_i(xi) * e(xj, a) + r_i(xj) * e(xi, a)) for a in PIDX)
        out += sum(F['ta'](a) * (t_i(xi) * e(xj, a) + t_i(xj) * e(xi, a)) for a in PIDX)
        out += sum(F['ab'](a, b) * e(xi, a) * e(xj, b) for a in PIDX for b in PIDX)
        return out

    def LiF(F, xi): return DiF(F, xi) / F['val']
    def LijF(F, xi, xj): return DijF(F, xi, xj) / F['val'] - DiF(F, xi) * DiF(F, xj) / F['val']**2

    def Lam1i(xi): return DiF(E, xi) + LiF(P, xi) + a_i(xi)
    def Lam2i(xi): return DiF(E, xi) + LiF(Q, xi) + b_i(xi)
    def Lci(xi): return DiF(E, xi) + LiF(B, xi) + a_i(xi) + b_i(xi)
    def Lam1ij(xi, xj): return DijF(E, xi, xj) + LijF(P, xi, xj) + a_ij(xi, xj)
    def Lam2ij(xi, xj): return DijF(E, xi, xj) + LijF(Q, xi, xj) + b_ij(xi, xj)
    def Lcij(xi, xj): return DijF(E, xi, xj) + LijF(B, xi, xj) + a_ij(xi, xj) + b_ij(xi, xj)

    C = sp.exp(E['val'])
    h1v, h2v, cv = C * Pv / u, C * Qv / v, C * B['val'] / (u * v)
    return dict(
        vars=(u, v, r_s, t_s, p1, p2, pth),
        c_d1=lambda xi: cv * Lci(xi),
        c_d2=lambda xi, xj: cv * (Lcij(xi, xj) + Lci(xi) * Lci(xj)),
        h1_d1=lambda xi: h1v * Lam1i(xi),
        h1_d2=lambda xi, xj: h1v * (Lam1ij(xi, xj) + Lam1i(xi) * Lam1i(xj)),
        h2_d1=lambda xi: h2v * Lam2i(xi),
        h2_d2=lambda xi, xj: h2v * (Lam2ij(xi, xj) + Lam2i(xi) * Lam2i(xj)),
        logc_d1=lambda xi: Lci(xi),
        logc_d2=lambda xi, xj: Lcij(xi, xj),
    )
