# Circular copulas: mathematical and API contract

Status: stage 1 of [circular-copulas.md](circular-copulas.md). Updated
September 21, 2026.

This document fixes the conventions, formulas, parameter domains, and
acceptance cases that stages 2 to 7 implement. It is normative for those
stages: an implementation that disagrees with a formula here is wrong until
this document is changed. Every formula below has been checked numerically by
the reference scripts in [tools/circulas/](../../tools/circulas/), which also
generate the golden values for the later test stages.

Throughout, a *circular* variable is one with `var_types` entry `"a"`, a
*linear* variable one with entry `"c"`. Discrete variables (`"d"`) are outside
the scope of the circular feature; see [Variable types](#variable-types).

## Copula-scale conventions

### The circle on the unit interval

The library sees a circular variable as a value $u \in [0, 1]$ with $0$
identified with $1$. The caller's marginal CDF decides where the cut sits: if
the angle $\phi$ has CDF $F$ on $[\phi_0, \phi_0 + 2\pi)$, then $u = F(\phi)$
and the cut is at $\phi_0$. Pseudo-observations computed by `to_pseudo_obs`
rank the angles as numbers, so they place the cut at the caller's zero angle.
Units, the direction of rotation, and the choice of $\phi_0$ are downstream
concerns. The library documents them and takes copula-scale input as it does
for linear variables.

A copula density $c$ is *periodic in its first coordinate* if
$c(0, v) = c(1, v)$ for all $v$, and analogously for the second coordinate.
Every family below that accepts a circular variable in a given argument
position is periodic in that coordinate. This is the property that makes the
joint density of a mixed vine continuous across each cut; see
[Mixed vines](#mixed-vines-the-periodicity-condition).

### Endpoint behavior

- **Densities** are periodic in each circular coordinate. They need not be
  periodic in a linear coordinate, and in general are not.
- **CDFs and h-functions** keep their anchored meaning: $C(0, v) = 0$,
  $C(u, 0) = 0$, $C(1, v) = v$, $C(u, 1) = u$, and
  $h_1(v \mid u) = \partial C / \partial u$ runs from $0$ at $v = 0$ to $1$ at
  $v = 1$. A conditional CDF of a circular variable therefore wraps from $0$ to
  $1$ across the cut, as the marginal one does.
- **Inverse h-functions** return the branch in $[0, 1]$: for
  $w \in [0, 1]$, $h_1^{-1}(w \mid u) \in [0, 1]$ with the value $0$ at
  $w = 0$ and $1$ at $w = 1$. No modular reduction happens after the fact; the
  formulas below produce the correct branch directly.
- **Parameters** that are phases are periodic and accepted on all of
  $\mathbb{R}$; see [Phase parameters](#phase-parameters).

### Variable types

| `var_types` | Meaning | Status |
| --- | --- | --- |
| `"c"` | continuous linear | unchanged |
| `"d"` | discrete | unchanged |
| `"a"` | continuous circular (angular) | new |

A circular variable is a kind of continuous variable: it takes one data
column, needs no left-limit column, and every code path that asks "is this
variable discrete?" answers no. Code paths that ask "is this variable linear
continuous?" must not treat `"a"` as `"c"`; the audit in the plan lists the
sites. Both questions are answered through named predicates, never by
comparing against the literal.

Supported pair geometries are `{"a", "a"}` (circular-circular),
`{"a", "c"}` and `{"c", "a"}` (circular-linear in either order), and the
existing linear and discrete combinations. The pairs `{"a", "d"}` and
`{"d", "a"}` are rejected by `check_var_types` with a message that names the
combination; supporting them needs a separate mathematical decision.

## Binding-density circulas (circular-circular)

### Construction

Let $g$ be a circular density on $\mathbb{R}$, that is, a $2\pi$-periodic
nonnegative function with $\int_{-\pi}^{\pi} g = 1$, and assume $g$ is
symmetric about zero. For an orientation $q \in \{-1, 1\}$ and a phase
$\mu \in \mathbb{R}$ the binding density

$$
c(u, v) = 2\pi\, g\bigl(2\pi (v - q u) - \mu\bigr), \qquad u, v \in [0, 1],
$$

is a copula density that is periodic in both coordinates. Uniform margins
follow from $\int_0^1 2\pi g(2\pi t - \phi)\, dt = \int_{-\pi}^{\pi} g = 1$
for every $\phi$. The density is constant along the lines $v - qu = $ const,
so $q = 1$ concentrates mass around $v \equiv u + \mu / 2\pi$ and $q = -1$
around $v \equiv -u + \mu / 2\pi$ (mod 1).

The three proposed choices of $g$ are the cardioid, wrapped Cauchy, and von
Mises densities; all three are symmetric about zero, which the rotation and
flip rules below rely on.

### The lifted CDF

The h-functions are differences of the *lifted* distribution function

$$
\tilde G(\theta) = \int_0^{\theta} g(t)\, dt, \qquad \theta \in \mathbb{R},
$$

which satisfies $\tilde G(0) = 0$, $\tilde G(-\theta) = -\tilde G(\theta)$,
$\tilde G(\theta + 2\pi) = \tilde G(\theta) + 1$, and $\tilde G' = g$. When
$g > 0$ everywhere, $\tilde G$ is a strictly increasing bijection of
$\mathbb{R}$, and $\tilde G^{-1}$ retains the number of completed turns: the
integer part of $w$ is the number of turns and the fractional part the
position within the turn. Implementations must compute $\tilde G$ and
$\tilde G^{-1}$ on all of $\mathbb{R}$, not on one period, because the
h-function formulas below evaluate them at arguments outside $[-\pi, \pi)$.

### Density, CDF, h-functions, and inverses

For $q \in \{-1, 1\}$ and $\mu \in \mathbb{R}$, with $\hat G$ any
antiderivative of $\tilde G$,

$$
\begin{aligned}
c(u, v) &= 2\pi\, g\bigl(2\pi(v - qu) - \mu\bigr), \\
C(u, v) &= \frac{q}{2\pi}\Bigl[\hat G(2\pi v - \mu) - \hat G\bigl(2\pi(v - qu) - \mu\bigr)
          - \hat G(-\mu) + \hat G(-2\pi q u - \mu)\Bigr], \\
h_1(v \mid u) &= \tilde G\bigl(2\pi(v - qu) - \mu\bigr) - \tilde G(-2\pi q u - \mu), \\
h_2(u \mid v) &= q\,\Bigl[\tilde G(2\pi v - \mu) - \tilde G\bigl(2\pi(v - qu) - \mu\bigr)\Bigr], \\
h_1^{-1}(w \mid u) &= q u + \frac{1}{2\pi}\Bigl[\mu + \tilde G^{-1}\bigl(w + \tilde G(-2\pi q u - \mu)\bigr)\Bigr], \\
h_2^{-1}(w \mid v) &= q\,\Bigl[v - \frac{1}{2\pi}\Bigl(\mu + \tilde G^{-1}\bigl(\tilde G(2\pi v - \mu) - q w\bigr)\Bigr)\Bigr].
\end{aligned}
$$

Here $h_1(v \mid u) = \partial C / \partial u$ is the conditional CDF of the
second argument given the first, and $h_2(u \mid v) = \partial C / \partial v$
the reverse, matching the library's convention. The inverse formulas produce
values in $[0, 1]$ for $w \in [0, 1]$ without any wrapping: at $w = 0$ the
argument of $\tilde G^{-1}$ collapses to the point where the two $\tilde G$
terms cancel, giving $0$, and at $w = 1$ the lifted inverse advances by
exactly one turn, giving $1$. Nothing beyond $g$, $\tilde G$, $\tilde G^{-1}$
(and $\hat G$ for the CDF) is family specific, so one base class implements
the primitives for all three families.

### The three families

| Family | $g(\theta)$ | $\tilde G(\theta)$ | $\hat G(\theta)$ | $\tilde G^{-1}$ |
| --- | --- | --- | --- | --- |
| `cardioid` | $\dfrac{1 + 2\rho\cos\theta}{2\pi}$ | $\dfrac{\theta + 2\rho\sin\theta}{2\pi}$ | $\dfrac{\theta^2/2 - 2\rho\cos\theta}{2\pi}$ | root of a Kepler-type equation |
| `wrapped_cauchy` | $\dfrac{1 - \rho^2}{2\pi(1 + \rho^2 - 2\rho\cos\theta)}$ | $\dfrac{1}{2\pi}\Bigl[\theta + 2\arctan\dfrac{\rho\sin\theta}{1 - \rho\cos\theta}\Bigr]$ | series, see below | closed form, see below |
| `von_mises` | $\dfrac{e^{\kappa\cos\theta}}{2\pi I_0(\kappa)}$ | $\dfrac{1}{2\pi}\Bigl[\theta + \dfrac{2}{I_0(\kappa)}\displaystyle\sum_{j \ge 1}\dfrac{I_j(\kappa)}{j}\sin(j\theta)\Bigr]$ | series | root solve |

The wrapped Cauchy expression for $\tilde G$ is a valid lift: the arctangent
term is continuous and bounded on $\mathbb{R}$ because
$1 - \rho\cos\theta > 0$ for $\rho < 1$, so the formula inherits
$\tilde G(\theta + 2\pi) = \tilde G(\theta) + 1$ from the linear term.
Its inverse is closed form by reducing to one turn: with $k$ the nearest
integer to $w$,

$$
\tilde G^{-1}(w) = 2\pi k + 2\arctan\Bigl(\frac{1 - \rho}{1 + \rho}\tan\bigl(\pi (w - k)\bigr)\Bigr).
$$

The antiderivative is
$\hat G(\theta) = \frac{1}{2\pi}\bigl[\theta^2 / 2 - 2\sum_{k \ge 1} \rho^k \cos(k\theta) / k^2\bigr]$,
a geometrically convergent series (the real part of a dilogarithm). In general,
if $g(\theta) = \bigl(1 + 2\sum_{j \ge 1} \rho_j \cos j\theta\bigr) / 2\pi$, then
$\hat G(\theta) = \theta^2 / 4\pi - \sum_{j \ge 1} \rho_j \cos(j\theta) / (\pi j^2)$,
so every binding family with known Fourier coefficients ($\rho$ for the
cardioid, $\rho^j$ for the wrapped Cauchy, $I_j(\kappa) / I_0(\kappa)$ for the
von Mises) has a closed-form CDF.

For the von Mises family the series for $\tilde G$ uses the Bessel ratios
$I_j(\kappa) / I_0(\kappa)$, which decay like $\exp(-j^2 / 2\kappa)$; the
number of terms for a tolerance $\varepsilon$ is about
$\sqrt{2\kappa \ln(1/\varepsilon)}$, below 70 for $\kappa \le 100$ and
$\varepsilon = 10^{-10}$. The ratios must be computed without forming
$I_0(\kappa)$ itself when $\kappa$ is large ($I_0$ overflows a double near
$\kappa = 713$); Boost.Math provides `cyl_bessel_i`, and the backward
recurrence of Hill (1977) is the standard alternative. Within the bound
$\kappa \le 100$ fixed below, direct evaluation is safe:
$I_0(100) \approx 10^{42}$.

The cardioid inverse solves $\theta + 2\rho\sin\theta = 2\pi w$. The left side
is increasing with derivative $1 + 2\rho\cos\theta \ge 0$, which vanishes at
one point when $\rho = 1/2$, so the solver must be bracketed (the root lies in
$[2\pi w - 2\rho, 2\pi w + 2\rho]$) rather than a bare Newton iteration. The
von Mises inverse is likewise a bracketed root solve of $\tilde G(\theta) = w$
on $[2\pi w - \pi, 2\pi w + \pi]$, with Newton steps using $g$ as the
derivative once inside the bracket.

### Parameter domains, independence, and identifiability

| Family | Parameters, in order | Domain | Independence |
| --- | --- | --- | --- |
| `cardioid` | $(\rho, \mu)$ | $\rho \in [0, 1/2]$, $\mu \in \mathbb{R}$ | $\rho = 0$ |
| `wrapped_cauchy` | $(\rho, \mu)$ | $\rho \in [0, 0.99]$, $\mu \in \mathbb{R}$ | $\rho = 0$ |
| `von_mises` | $(\kappa, \mu)$ | $\kappa \in [0, 100]$, $\mu \in \mathbb{R}$ | $\kappa = 0$ |

A negative concentration is not a separate model: for each of the three
densities, $g(\theta; -\rho) = g(\theta + \pi; \rho)$, so
$(-\rho, \mu)$ and $(\rho, \mu + \pi)$ give the same copula. Restricting the
concentration to be nonnegative is therefore a normalization, not a loss of
generality. At concentration zero the phase is not identified; fits report
whatever value the optimizer ends at, and the number of parameters counted by
`get_npars()` is $2$ regardless.

The upper bounds are chosen so that the density stays representable and the
optimizer well conditioned. The cardioid bound $1/2$ is the largest value for
which $g \ge 0$. For the wrapped Cauchy, $\rho = 0.99$ gives a ratio of largest to smallest
density value of $((1 + \rho)/(1 - \rho))^2 \approx 4 \cdot 10^4$ and a
Kendall's $\tau$ of about $0.96$ at $\mu = 0$, comparable to the caps of the
existing one-parameter families (Clayton at $\theta = 28$ gives $0.93$). For
the von Mises family, $\kappa = 100$ gives a peak density of about $25$, a
minimum of about $10^{-86}$, and a circular standard deviation of $0.1$ radians;
larger values buy no useful additional dependence. Stage 3 may tighten these
bounds after numerical experiments and must record any change here.

### Phase parameters

Phases are stored in radians, following Jones, Pewsey, and Kato (2015). The
optimizer treats a phase as unbounded (lower bound $-\infty$, upper bound
$+\infty$), which `tools_transforms` maps to the identity transform, so there
is no bound to hit and no finite-difference step crosses a cut in parameter
space. `set_parameters` accepts any real phase and stores it as given. `fit`
reports the phase reduced to $[-\pi, \pi)$. Two parameter vectors whose phases
differ by a multiple of $2\pi$ define the same copula, and tests compare
copulas, not parameter vectors, whenever a phase is involved.

### Orientation and rotation

The orientation $q$ is not a parameter; it is the copula rotation. The base
(rotation $0$) copula is the $q = 1$ density. The library rotates by evaluating
the base copula on transformed arguments, with $c_{90}(u, v) = c_0(v, 1 - u)$,
$c_{180}(u, v) = c_0(1 - u, 1 - v)$, and $c_{270}(u, v) = c_0(1 - v, u)$. For
a base density with parameters $(\rho, \mu)$, symmetry and periodicity of $g$
give

$$
\begin{aligned}
c_{90}(u, v) &= 2\pi\, g\bigl(2\pi(v + u) + \mu\bigr), &\quad&\text{the } q = -1 \text{ density with phase } -\mu, \\
c_{180}(u, v) &= 2\pi\, g\bigl(2\pi(v - u) + \mu\bigr), &&\text{the } q = 1 \text{ density with phase } -\mu, \\
c_{270}(u, v) &= 2\pi\, g\bigl(2\pi(v + u) - \mu\bigr), &&\text{the } q = -1 \text{ density with phase } \mu.
\end{aligned}
$$

Because the phase is a free parameter, rotations $180$ and $270$ add nothing:
the family has exactly two distinct rotations, $0$ and $90$. The three
families form the group `bicop_families::two_rotations`. `check_rotation`
rejects $180$ and $270$ for them with a message naming the allowed set, as it
rejects nonzero rotations for `rotationless` families. With
`allow_rotations = true`, candidate generation produces both rotations; with
`allow_rotations = false`, only rotation $0$. There is no sign of a linear
dependence measure that could choose between them, so both are always fitted.

### Flip

`Bicop::flip()` exchanges the two arguments. For the base copula,
$c_0(v, u) = 2\pi g(2\pi(u - v) - \mu) = 2\pi g(2\pi(v - u) + \mu)$, so the
flip of rotation $0$ negates the phase. The rotation-$90$ copula is
exchangeable, so its flip is the identity. The generic flip in `Bicop`
switches $90 \leftrightarrow 270$ and then flips the base copula; for a
`two_rotations` family this yields rotation $270$ with phase $-\mu$, which is
the same density as rotation $90$ with phase $\mu$. The implementation must
canonicalize that result back to rotation $90$ with the original phase, so
that a flipped model never carries a rotation outside $\{0, 90\}$.

| Rotation | `flip()` | Density unchanged? |
| --- | --- | --- |
| $0$ | $\mu \to -\mu$ | no (unless $\mu \in \pi\mathbb{Z}$) |
| $90$ | identity | yes |

### Dependence measures

Kendall's $\tau$, Spearman's $\rho$, and Blomqvist's $\beta$ keep their
ordinary definitions. They depend on the cut and are misleading for circular
dependence: for $q = 1$ and $\mu = \pi$ the copula concentrates on
$v \equiv u + 1/2$, and in the perfectly dependent limit $\tau = 0$,
$\rho_S = -1/2$, and $\beta = -1$. Consequently:

- `parameters_to_tau` returns `NaN` for every circular family, so that a
  cut-dependent number is never reported as a dependence measure.
- `tau_to_parameters` is not available (`no_tau_to_parameters`), the
  `bicop_families::itau` list is unchanged, and `parametric_method = "itau"`
  excludes the circular families as it excludes the two-parameter ones.
- Linear $\tau$ must not drive the choice of rotation, family preselection, or
  tree weights for pairs involving a circular variable. The `lt` / `ut`
  tail-dependence heuristics do not apply either: all three densities are
  bounded, so both tail-dependence coefficients are zero.

Circular dependence summaries (the circular correlation coefficients of
Fisher and Lee or of Jammalamadaka and Sarma) are a stage 6 concern and are
exposed separately when the tree criterion is chosen.

### Fitting

`ParBicop::fit` initializes from Kendall's $\tau$, which is useless here. The
binding families use a moment initialization instead. Under the model, the
residual angles $\theta_i = 2\pi(v_i - q u_i) \bmod 2\pi$ have density
$g(\cdot - \mu)$, so with the (weighted) mean resultant
$\bar R e^{i\bar\mu} = \sum_i w_i e^{i\theta_i} / \sum_i w_i$,

- $\hat\mu = \bar\mu$ for all three families;
- `cardioid`: $\hat\rho = \min(\bar R, 1/2)$, since the mean resultant length
  of the cardioid is $\rho$;
- `wrapped_cauchy`: $\hat\rho = \bar R$, since its mean resultant length is
  $\rho$;
- `von_mises`: $\hat\kappa = A^{-1}(\bar R)$ with $A = I_1 / I_0$, solved
  numerically and clipped to the domain.

Each rotation is a separate candidate, so $q$ is fixed within a fit. The
maximum-likelihood fit then starts from these values with the phase
unbounded. The `adjust_parameters_bounds` step, which narrows the search box
around a $\tau$-implied value, is skipped for these families.

## Cylindrical sections copulas (circular-linear)

### Convention

A cylindrical family has one circular and one linear argument. Internally the
circular variable is the first argument. When `var_types` is `{"c", "a"}` the
primitives are evaluated on swapped columns and $h_1 \leftrightarrow h_2$,
$h_1^{-1} \leftrightarrow h_2^{-1}$ are exchanged, exactly as `flip()` would
do; the parameter vector is the same in both orders. The family is
`rotationless`: its parameter set is closed under every reflection of either
axis, as shown below, so a rotation would only relabel parameters.

### Cubic sections

The family `cubic_sections` with parameters $(a, b, \mu)$,
$a \in [0, 1]$, $b \in [-1, 1]$, $\mu \in \mathbb{R}$, is

$$
\begin{aligned}
c(u, v) &= 1 + \cos(2\pi u - \mu)\, p(v), &\quad p(v) &= a\,(1 - v)(1 - 3v) + b\, v\,(2 - 3v), \\
C(u, v) &= uv + \frac{1}{2\pi}\bigl[\sin(2\pi u - \mu) + \sin\mu\bigr]\, P(v), & P(v) &= a\, v(1 - v)^2 + b\, v^2(1 - v), \\
h_1(v \mid u) &= v + \cos(2\pi u - \mu)\, P(v), \\
h_2(u \mid v) &= u + \frac{1}{2\pi}\bigl[\sin(2\pi u - \mu) + \sin\mu\bigr]\, p(v),
\end{aligned}
$$

with $P' = p$ and $P(0) = P(1) = 0$. It has cubic sections in $v$ in the
sense of Nelsen, Quesada-Molina, and Rodríguez-Lallena (1997). The parameter
$a$ controls the direction preference near $v = 0$ and $b$ the one near
$v = 1$; the two ends may prefer different directions with different
strengths. Setting $a = b$ gives quadratic sections, $p(v) = a(1 - 2v)$, in
the sense of Quesada-Molina and Rodríguez-Lallena (1995): the circular-linear
analog of the Farlie-Gumbel-Morgenstern copula, with mass shifted toward
$u \equiv \mu / 2\pi$ for small $v$ and toward the opposite direction for
large $v$. The term $\sin\mu$ in $C$ enforces $C(0, v) = 0$ and $C(1, v) = v$
for every phase.

**Constraint region.** Because $\cos(2\pi u - \mu)$ attains every value in
$[-1, 1]$ as $u$ ranges over $[0, 1]$, nonnegativity of $c$ is equivalent to
$|p(v)| \le 1$ for all $v \in [0, 1]$. This region is exactly the box
$|a| \le 1$, $|b| \le 1$. The endpoint values $p(0) = a$ and $p(1) = -b$ show
that the box is necessary. For sufficiency, $\max_v |p(v)|$ is a convex
function of $(a, b)$ (a maximum of absolute values of linear functions), so
its maximum over the box is attained at a corner; at the four corners
$p(v) = \pm(1 - 2v)$ or $\pm(1 - 6v + 6v^2)$, both bounded by $1$ on
$[0, 1]$. The fit therefore needs only box constraints. The pair
$(a, b, \mu)$ and $(-a, -b, \mu + \pi)$ define the same copula, so $a \ge 0$
is a normalization, with the residual ambiguity $(0, b, \mu) \sim (0, -b, \mu + \pi)$
on the boundary $a = 0$ only. Independence is $a = b = 0$.

$h_1^{-1}(w \mid u)$ is the root in $[0, 1]$ of a cubic in $v$ and
$h_2^{-1}(w \mid v)$ a bracketed root solve in $u$; both are bracketed on
$[0, 1]$ and monotone because the derivative is $c \ge 0$.

In closed form, $\tau = (a + b)\sin\mu / (3\pi)$ and
$\rho_S = (a + b)\sin\mu / (2\pi)$; both vanish at $\mu = 0$ for every $a, b$,
so the same warning as for the circulas applies and `parameters_to_tau`
returns `NaN` here as well.

### Symmetries of the cylindrical family

Under the reflections of the two axes and the $180$-degree rotation,

| Operation | `cubic_sections` $(a, b, \mu)$ |
| --- | --- |
| $u \to 1 - u$ (reflect the circular axis) | $(a, b, -\mu)$ |
| $v \to 1 - v$ (reflect the linear axis) | $(-b, -a, \mu)$, i.e. $(b, a, \mu + \pi)$ |
| rotation $180$ | $(b, a, \pi - \mu)$ |
| `flip()` | swap `var_types`; parameters unchanged |

Each row maps the parameter domain onto itself (after the sign normalization
of $a$), which is why no rotation is needed. `flip()` exchanges the axis
order, and the primitives read the order from `var_types`, so the parameters
do not change.

### Binding families on the cylinder

A binding density is also a valid copula density for a circular-linear pair:
it is periodic in the circular coordinate, and periodicity in the linear one
does no harm. It is admitted for `{"a", "c"}` and `{"c", "a"}` with the same
parameters, rotations, and flip rules as for `{"a", "a"}`. It imposes equal
density at the two ends of the linear coordinate, which the cylindrical
family does not, so it complements rather than replaces it.

### Fitting

The moment initialization uses the identities

$$
\begin{aligned}
\mathbb{E}\bigl[e^{i 2\pi U} p_1(V)\bigr] &= \bigl(\tfrac{a}{15} + \tfrac{b}{60}\bigr) e^{i\mu}, \quad
\mathbb{E}\bigl[e^{i 2\pi U} p_2(V)\bigr] = \bigl(\tfrac{a}{60} + \tfrac{b}{15}\bigr) e^{i\mu},
\end{aligned}
$$

with $p_1(v) = (1 - v)(1 - 3v)$ and $p_2(v) = v(2 - 3v)$. Replacing the
expectations by weighted sample means gives $\hat\mu$ as the argument of the
sum of the two moments and the amplitudes by a $2 \times 2$ linear solve; the
results are clipped to
the parameter domain and the sign of $a$ normalized. The maximum-likelihood
fit starts from these values.

## Mixed vines: the periodicity condition

### Statement

Let $x_i$ be a circular variable of a vine with marginal density $f_i$ that
is periodic across the cut. The joint density of the vine, viewed as a
function of $x_i$, is continuous across the cut if every pair copula in which
$x_i$, or a conditional CDF of $x_i$, is a *conditioned* argument is periodic
in that argument. A pair copula in which $x_i$ appears only in the
conditioning set is not constrained by $x_i$.

### Why the condition suffices

The vine density is $\prod_k f_k(x_k) \prod_e c_e\bigl(F_{a_e \mid D_e}, F_{b_e \mid D_e}\bigr)$
with conditioned variables $a_e, b_e$ and conditioning set $D_e$. Consider
the factors as functions of $x_i$.

First, the marginal $f_i$ is periodic by assumption.

Second, in tree $1$ the conditional CDFs are the marginals, and the argument
$u_i = F_i(x_i)$ of a pair copula with $a_e = i$ wraps from $0$ to $1$ across
the cut. The factor $c_e(u_i, u_j)$ is continuous across the cut if and only if
$c_e(0, \cdot) = c_e(1, \cdot)$, which is periodicity in the first argument.

Third, the arguments passed to tree $t + 1$ are h-functions of tree $t$. For
an edge with $a_e = i$, the next-tree argument $h_{i \mid j D}(u_i \mid u_j) = \int_0^{u_i} c_e(s, u_j)\, ds$
is a CDF in $u_i$ and again wraps from $0$ to $1$; the same periodicity
requirement therefore applies to every later pair copula in which this
transform is a conditioned argument, and the argument repeats up the trees.
The *other* h-function of the same edge, $h_{j \mid i D}(u_j \mid u_i) = \int_0^{u_j} c_e(u_i, t)\, dt$,
is continuous in $u_i$ across the cut exactly when $c_e(0, \cdot) = c_e(1, \cdot)$,
which is the same condition. So once every edge that has $x_i$ (or a
transform of it) as a conditioned argument is periodic in that argument, all
conditional CDFs computed through the vine are continuous in $x_i$ across the
cut, including those of other variables given $x_i$.

Finally, an edge with $i \in D_e$ but $i \notin \{a_e, b_e\}$ receives $x_i$
only through the conditional CDFs $F_{a_e \mid D_e}$ and $F_{b_e \mid D_e}$,
which the previous step showed to be continuous in $x_i$. Composing with any
continuous $c_e$ preserves continuity, so this edge imposes nothing on
$c_e$.

### Consequences for eligibility

For a pair copula with conditioned variables $(a_e, b_e)$ the relevant
geometry is that of $a_e$ and $b_e$ alone. The library already propagates
`var_types` this way: `Vinecop::set_var_types_internal` assigns each pair's
types from the conditioned variables of the previous tree without any literal
comparison, and the selector's edge inheritance does the same. The eligible
families for an edge are those of [Family eligibility](#family-eligibility)
for the pair's geometry. The independence copula is periodic in both
coordinates and is eligible everywhere.

### Changing the cut

Shifting the cut of $x_i$ by a fraction $\delta$ of a turn maps
$u_i \to (u_i + \delta) \bmod 1$. In tree $1$ every family here absorbs the
shift into its phase: $\mu \to \mu + 2\pi q\delta$ for a binding density and
$\mu \to \mu + 2\pi\delta$ for a cylindrical one. The parametric families are
therefore closed under cut changes at the pair level, and a nonparametric
estimator on a uniform periodic grid is closed up to interpolation error.

In higher trees the picture changes. Under the new cut, the conditional CDF
becomes
$F^{\text{new}}_{i \mid D}(x_i \mid x_D) = \bigl(F_{i \mid D}(x_i \mid x_D) - F_{i \mid D}(\phi_0 + 2\pi\delta \mid x_D)\bigr) \bmod 1$,
and the subtracted term depends on $x_D$ unless $x_i$ is independent of $x_D$.
The pair copula of $(i, j)$ given $D$ under the new cut is thus the old one
composed with a shift in its first argument that varies with the conditioning
values. It is a pair copula of a *non-simplified* vine in general, and a
simplified vine fitted under one cut is not the reparameterization of a
simplified vine fitted under another. The documentation must state this and
must not promise invariance of a simplified mixed vine under cut changes.
Cut sensitivity of higher-tree fits is a stage 6 validation item.

## API contract

### Family eligibility

| Pair geometry | Eligible families |
| --- | --- |
| `{"c", "c"}`, `{"c", "d"}`, `{"d", "c"}`, `{"d", "d"}` | all existing families; unchanged |
| `{"a", "a"}` | `indep`, `cardioid`, `wrapped_cauchy`, `von_mises`, `tll` |
| `{"a", "c"}`, `{"c", "a"}` | `indep`, `cardioid`, `wrapped_cauchy`, `von_mises`, `cubic_sections`, `tll` |
| `{"a", "d"}`, `{"d", "a"}` | rejected |

Eligibility is a function of `(family, var_types)` exposed in one place and
used by construction, `set_var_types`, fitting, selection, and the vine. The
new families are added to `bicop_families::all`; because eligibility filters
the candidate set, the effective search for a linear pair is unchanged.
Explicit family sets are intersected with the eligible set. If the
intersection is empty, `select` throws a message listing the eligible
families for the pair's geometry rather than fitting nothing or silently
substituting; `indep` is always eligible, so a family set containing it never
triggers this.

### Construction and validation

- The `Bicop` constructor keeps its default `var_types = {"c", "c"}`.
  Constructing a circular family with an ineligible geometry, including the
  default, throws with a message naming the family and its eligible
  geometries. Geometry is never inferred from the family.
- `set_var_types` applies the same check, so a circular pair copula inside a
  `Vinecop` must be constructed with its geometry, and `Vinecop`'s propagation
  of `var_types` to its pair copulas validates each edge.
- A linear family with an `"a"` in its `var_types` is rejected for the same
  reason: its density is not periodic, and the vine condition would be
  violated silently.
- `check_var_types` in both classes accepts `"a"` and rejects the
  circular-discrete pairs listed above.
- `Vinecop::select` and the threshold and truncation searches throw when any
  variable is circular until stage 6 replaces the tree criterion; `Vinecop::fit`
  on a supplied structure is supported from stage 5.

### Parameters

Parameter order and domains are as in the tables above. The optimizer bounds
are the domain bounds, with $\pm\infty$ for phases. `get_npars()` returns the
parameter count. JSON serialization of the parametric families is unchanged:
`"fam"`, `"rot"`, `"par"`, `"vt"`, and the fit statistics suffice, and the new
family names round-trip through `get_family_name` / `get_family_enum`.

### Nonparametric grid serialization

A `tll` model whose grid is not the default normal grid on two linear axes
carries an additional JSON object `"grid"` with two fields: `"knots"`, a list
of two vectors holding the knot positions of the first and second axis, and
`"types"`, the two axis types (`"c"` or `"a"`) the knots were built for. The
`"par"` field keeps the density values on the knot lattice, row `i` and
column `j` belonging to the `i`-th knot of the first axis and the `j`-th of
the second. On reading, a model with a `"grid"` field rebuilds its
interpolation grid from those knots; a model without one rebuilds the normal
grid from the row count of `"par"`, as today, so every existing file reloads
unchanged. A `"grid"` whose knot counts disagree with the shape of `"par"`,
whose types disagree with `"vt"`, or whose knots are not increasing in
`[0, 1]` is rejected with a message naming the field. Writing a linear `tll`
fit does not add the field, so files written by a linear-only library are
read unchanged.

### Candidate generation and preselection

`create_candidate_bicops` receives the pair's geometry. For a
`two_rotations` family it emits rotations $\{0, 90\}$ when
`allow_rotations` is set and $\{0\}$ otherwise, ignoring the sign of
Kendall's $\tau$. `preselect_candidates` does not apply the `lt` / `ut`
heuristics to circular families. Linear families are never candidates for a
pair with a circular variable.

### Tree criterion

An edge whose pair has a circular conditioned variable is weighted by
`tools_stats::pairwise_circular` under every built-in `tree_criterion`. For
two circular variables the measure is
$\max_{q = \pm 1} \bigl| E\, e^{i 2\pi (V - qU)} \bigr|$, the larger
mean resultant length of the angle differences and sums; it is one exactly
for a rotation or reflection and does not depend on the cut of either
variable. For a circular $U$ and a linear $V$ it is
$\bigl(|E\, e^{i 2\pi U} f_1(V)|^2 + |E\, e^{i 2\pi U} f_2(V)|^2\bigr)^{1/2}$
with the orthonormal Legendre polynomials $f_1(v) = \sqrt{3}(2v - 1)$ and
$f_2(v) = \sqrt{5}(6v^2 - 6v + 1)$, the moments that identify the sections
copulas; it does not depend on the cut of $U$ and is invariant under
$V \to 1 - V$. Both are bounded by one and vanish under independence. The
`"custom"` criterion receives the pair data unchanged. Rank-based criteria are
not used for such pairs because they depend on the cut: Kendall's $\tau$ of
the half-turn pair is zero.

### Derivatives and views

`Bicop::check_deriv_preconditions` accepts circular variable types; the
parametric circular families obtain their parameter derivatives from the
existing finite-difference fallback unless stage 3 adds analytic leaves. The
`as_continuous()` operations drop discreteness element-wise and preserve
`"a"`.

## Acceptance cases and tolerances

Reference values come from the scripts in
[tools/circulas/](../../tools/circulas/), which evaluate every formula with
NumPy and SciPy against numerical integration and differentiation of the
density, independently of the C++ implementation. Golden values for the C++
tests are generated by those scripts and stored with the tests, following the
existing golden-value convention in `scripts/README.md`.

| Identity | Families | Tolerance | Reference |
| --- | --- | --- | --- |
| $\int_0^1\!\int_0^1 c = 1$; uniform margins | all | $10^{-8}$ | adaptive quadrature |
| $c(0, v) = c(1, v)$ on circular axes; $c(u, 0) \ne c(u, 1)$ for cylindrical | all | $10^{-12}$ | exact |
| $h_1 = \partial C / \partial u$, $h_2 = \partial C / \partial v$ | all | $10^{-6}$ | central differences of the integrated $C$ |
| $h_k(h_k^{-1}(w)) = w$, $h_k^{-1} \in [0, 1]$, endpoint values $0$ and $1$ | all | $10^{-10}$ closed form; $10^{-8}$ root-solved | exact |
| $\tilde G' = g$, $\tilde G(\theta + 2\pi) = \tilde G(\theta) + 1$, oddness, $\tilde G(\tilde G^{-1}(w)) = w$ on several turns | binding | $10^{-9}$ | quadrature of $g$ |
| von Mises series against quadrature up to $\kappa = 100$ | `von_mises` | $10^{-9}$ | quadrature |
| closed-form $C$ via $\hat G$ against double quadrature | `cardioid` | $10^{-9}$ | quadrature |
| rotation and flip identities of the tables above | all | $10^{-12}$ | exact |
| independence at zero concentration or amplitude | all | exact | — |
| $(-\rho, \mu) = (\rho, \mu + \pi)$; $(-a, -b, \mu) = (a, b, \mu + \pi)$ | all | $10^{-12}$ | exact |
| closed-form $\tau$, $\rho_S$ of the cubic sections family (reference only; the library reports `NaN`) | `cubic_sections` | exact | symbolic |
| the half-turn case at wrapped Cauchy $\rho = 0.95$: $\tau = -0.054$, $\rho_S = -0.479$, $\beta = -0.902$ (limits $0$, $-1/2$, $-1$) | binding | $\pm 0.03$ | Monte Carlo, $n = 2 \cdot 10^5$ |
| fit recovery of $(\text{concentration}, \mu)$ at $n = 2000$, including $\mu$ near $\pm\pi$ and the half-turn case | all | $\pm 0.05$ / $\pm 0.1$ rad | simulation |
| Rosenblatt round-trip of a mixed vine; joint density equal at each circular cut | vine | $10^{-10}$ | exact |
| mixed-vine density against the product of pair densities assembled by hand | vine | $10^{-10}$ | exact |

Tests compare copulas, not parameter vectors, whenever phases are involved,
by evaluating densities on a fixed grid or by reducing phases to
$[-\pi, \pi)$ first.

## References

- Jones, M. C., Pewsey, A., and Kato, S. (2015). On a class of circulas:
  copulas for circular distributions. *Annals of the Institute of Statistical
  Mathematics* 67, 843–862. The binding-density construction, its CDF and
  conditionals, and circular dependence measures.
- Hodel, F. H. and Fieberg, J. R. (2022). Circular-linear copulae for animal
  movement data. *Methods in Ecology and Evolution* 13, 1001–1013. Cylindrical
  copulas with quadratic and cubic sections and the `cylcop` package. The
  parameterizations above were derived independently; constants may differ
  from `cylcop`'s and must be reconciled before any stage 7 comparison.
- Quesada-Molina, J. J. and Rodríguez-Lallena, J. A. (1995). Bivariate copulas
  with quadratic sections. *Journal of Nonparametric Statistics* 5, 323–337.
- Nelsen, R. B., Quesada-Molina, J. J., and Rodríguez-Lallena, J. A. (1997).
  Bivariate copulas with cubic sections. *Journal of Nonparametric Statistics*
  7, 205–220.
- Hill, G. W. (1977). Algorithm 518: Incomplete Bessel function $I_0$: the von
  Mises distribution. *ACM Transactions on Mathematical Software* 3, 279–284.
  The backward recurrence for the von Mises CDF.
