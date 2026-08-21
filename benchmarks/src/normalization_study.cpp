// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Accuracy and cost of the schemes that could renormalize a TLL grid to
// uniform margins, and of the knot clamp in the conditional cdf.
//
//   build:  cmake --preset bench && cmake --build build/bench --target
//           normalization_study
//   run:    build/bench/bin/normalization_study
//
// It is a study, not a regression test: run it when changing the
// normalization scheme, its iteration bound, or the knot clamp.
//
// The arms are measured from the grid a `tll` fit stores, pushed back off the
// uniform-margin manifold by a row and column rescaling: the same scaling
// class, hence the same limit, but the same amount of work to reach it as the
// raw local-likelihood surface. To measure that surface itself, pass
// `norm_maxiter = 0` where `TllBicop::fit` builds the grid
// (`bicop/implementation/tll.ipp`), rebuild, and run again.
//
// Every arm is a positive row/column rescaling, so all of them stay in one
// scaling class and share the same doubly stochastic limit; they differ only
// in how fast they get there, how the residual is split between the two
// margins, and whether they commute with transposition. Section 0 checks the
// study's restatement of the shipped sweep against the library first, so a
// mismatch there means the study is wrong rather than the library.

#include "helpers.hpp"
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <vinecopulib.hpp>

using namespace vinecopulib;
using tools_interpolation::InterpolationGrid;

namespace {

constexpr double kMinMass = 1e-20; // the guard normalize_margins applies

// --------------------------------------------------------------- the grid

//! The grid `KernelBicop::make_normal_grid` builds, with the endpoint
//! snapping `InterpolationGrid`'s constructor applies on top. Restated here
//! rather than reached through the protected member; section 0 checks it.
Eigen::VectorXd
grid_points(size_t m)
{
  Eigen::VectorXd z(m);
  for (size_t i = 0; i < m; ++i)
    z(i) = -3.25 + static_cast<double>(i) * (6.5 / static_cast<double>(m - 1));
  Eigen::VectorXd g = tools_stats::pnorm(z);
  g(0) = 0.0;
  g(m - 1) = 1.0;
  return g;
}

Eigen::VectorXd
trapezoid_weights(const Eigen::VectorXd& g)
{
  const ptrdiff_t m = g.size();
  Eigen::VectorXd w(m);
  w(0) = (g(1) - g(0)) / 2.0;
  w.segment(1, m - 2) = (g.tail(m - 2) - g.head(m - 2)) / 2.0;
  w(m - 1) = (g(m - 1) - g(m - 2)) / 2.0;
  return w;
}

// ------------------------------------------------------------- the metrics

struct Residual
{
  double r1; // max_i |int c(., u_i) - 1|, the first margin (rows)
  double r2; // max_j |int c(u_j, .) - 1|, the second margin (columns)
  double worst() const { return std::max(r1, r2); }
};

Residual
residual(const Eigen::MatrixXd& v, const Eigen::VectorXd& w)
{
  return { (v * w).array().abs().maxCoeff() == 0.0
             ? 1.0
             : ((v * w).array() - 1.0).abs().maxCoeff(),
           ((v.transpose() * w).array() - 1.0).abs().maxCoeff() };
}

// ---------------------------------------------------------------- the arms

Eigen::MatrixXd
row_scaled(const Eigen::MatrixXd& v, const Eigen::VectorXd& w)
{
  Eigen::VectorXd s = (v * w).cwiseMax(kMinMass).cwiseInverse();
  Eigen::MatrixXd out = v;
  out.array().colwise() *= s.array();
  return out;
}

Eigen::MatrixXd
col_scaled(const Eigen::MatrixXd& v, const Eigen::VectorXd& w)
{
  Eigen::VectorXd s = (v.transpose() * w).cwiseMax(kMinMass).cwiseInverse();
  Eigen::MatrixXd out = v;
  out.array().rowwise() *= s.transpose().array();
  return out;
}

//! rows then columns -- one shipped sweep
void
step_alt(Eigen::MatrixXd& v, const Eigen::VectorXd& w)
{
  v = col_scaled(row_scaled(v, w), w);
}

//! the geometric mean of the two sweep orderings. Transposing swaps them, so
//! the mean is exactly transposition equivariant at any iteration count.
void
step_altsym(Eigen::MatrixXd& v, const Eigen::VectorXd& w)
{
  const Eigen::MatrixXd a = col_scaled(row_scaled(v, w), w);
  const Eigen::MatrixXd b = row_scaled(col_scaled(v, w), w);
  v = (a.array() * b.array()).sqrt();
}

//! the geometric mean of the row-only and column-only rescalings, optionally
//! over-relaxed. alpha = 1 is the plain simultaneous (Jacobi) update; without
//! the square root (alpha = 2) the total-mass mode oscillates forever.
void
step_sym(Eigen::MatrixXd& v, const Eigen::VectorXd& w, double alpha)
{
  const Eigen::VectorXd r = (v * w).cwiseMax(kMinMass);
  const Eigen::VectorXd c = (v.transpose() * w).cwiseMax(kMinMass);
  const Eigen::VectorXd sr = r.array().pow(-alpha / 2.0);
  const Eigen::VectorXd sc = c.array().pow(-alpha / 2.0);
  v.array().colwise() *= sr.array();
  v.array().rowwise() *= sc.transpose().array();
}

// --- the dual (log-domain) formulation, shared by `logdomain` and `anderson`

//! log(sum_j exp(a_j)), returning -inf if every term is -inf
double
logsumexp(const Eigen::VectorXd& a)
{
  const double hi = a.maxCoeff();
  if (!std::isfinite(hi))
    return hi;
  return hi + std::log((a.array() - hi).exp().sum());
}

//! log w_j + log V_ij, the matrix the dual iteration works on
Eigen::MatrixXd
log_kernel(const Eigen::MatrixXd& v, const Eigen::VectorXd& w)
{
  Eigen::MatrixXd lk = v.array().log();
  lk.array().rowwise() += w.array().log().transpose();
  return lk;
}

//! one alternating sweep on the dual potentials: the same iteration as
//! `step_alt`, computed through logsumexp instead of products and quotients
void
dual_step(const Eigen::MatrixXd& lk,
          const Eigen::MatrixXd& lkt,
          Eigen::VectorXd& x,
          Eigen::VectorXd& y)
{
  const ptrdiff_t m = x.size();
  for (ptrdiff_t i = 0; i < m; ++i) {
    const double s = logsumexp((lk.row(i).transpose() + y).eval());
    x(i) = std::isfinite(s) ? -s : 0.0;
  }
  for (ptrdiff_t j = 0; j < m; ++j) {
    const double s = logsumexp((lkt.row(j).transpose() + x).eval());
    y(j) = std::isfinite(s) ? -s : 0.0;
  }
}

Eigen::MatrixXd
from_potentials(const Eigen::MatrixXd& v,
                const Eigen::VectorXd& x,
                const Eigen::VectorXd& y)
{
  Eigen::MatrixXd out = v;
  out.array().colwise() *= x.array().exp();
  out.array().rowwise() *= y.array().exp().transpose();
  return out;
}

//! Anderson-accelerated dual iteration (type II, depth `depth`): the next
//! iterate is the least-squares combination of the last `depth` fixed-point
//! residuals rather than the last one alone.
int
run_anderson(Eigen::MatrixXd& v,
             const Eigen::VectorXd& w,
             int max_iter,
             double tol,
             int depth)
{
  const ptrdiff_t m = v.rows();
  const Eigen::MatrixXd lk = log_kernel(v, w);
  const Eigen::MatrixXd lkt = log_kernel(v.transpose(), w);
  Eigen::VectorXd z = Eigen::VectorXd::Zero(2 * m), gz(2 * m), f(2 * m);
  std::vector<Eigen::VectorXd> dz, df;
  Eigen::VectorXd f_prev, z_prev;
  int used = 0;

  for (int k = 0; k < max_iter; ++k) {
    Eigen::VectorXd x = z.head(m), y = z.tail(m);
    dual_step(lk, lkt, x, y);
    gz.head(m) = x;
    gz.tail(m) = y;
    f = gz - z;
    used = k + 1;

    if (k > 0) {
      dz.push_back(z - z_prev);
      df.push_back(f - f_prev);
      if (static_cast<int>(dz.size()) > depth) {
        dz.erase(dz.begin());
        df.erase(df.begin());
      }
    }
    f_prev = f;
    z_prev = z;

    if (dz.empty()) {
      z = gz;
    } else {
      const ptrdiff_t d = static_cast<ptrdiff_t>(dz.size());
      Eigen::MatrixXd dF(2 * m, d), dZ(2 * m, d);
      for (ptrdiff_t i = 0; i < d; ++i) {
        dF.col(i) = df[static_cast<size_t>(i)];
        dZ.col(i) = dz[static_cast<size_t>(i)];
      }
      Eigen::VectorXd gamma = dF.colPivHouseholderQr().solve(f);
      if (!gamma.array().isFinite().all())
        gamma.setZero();
      z = gz - (dZ + dF) * gamma;
    }

    if (tol > 0.0) {
      const Eigen::MatrixXd cur = from_potentials(v, z.head(m), z.tail(m));
      if (residual(cur, w).worst() < tol)
        break;
    }
  }

  v = from_potentials(v, z.head(m), z.tail(m));
  return used;
}

//! One arm: `run` applies at most `max_iter` iterations, stopping early when
//! both margins are within `tol` (pass tol <= 0 to disable), and returns the
//! number of iterations it used.
struct Arm
{
  std::string name;
  std::function<int(Eigen::MatrixXd&, const Eigen::VectorXd&, int, double)> run;
};

Arm
simple_arm(
  const std::string& name,
  const std::function<void(Eigen::MatrixXd&, const Eigen::VectorXd&)>& step)
{
  return { name,
           [step](Eigen::MatrixXd& v,
                  const Eigen::VectorXd& w,
                  int max_iter,
                  double tol) {
             for (int k = 0; k < max_iter; ++k) {
               step(v, w);
               if ((tol > 0.0) && (residual(v, w).worst() < tol))
                 return k + 1;
             }
             return max_iter;
           } };
}

std::vector<Arm>
all_arms()
{
  std::vector<Arm> arms;
  arms.push_back(simple_arm("alt", step_alt));
  arms.push_back(simple_arm("altsym", step_altsym));
  arms.push_back(
    simple_arm("sym", [](Eigen::MatrixXd& v, const Eigen::VectorXd& w) {
      step_sym(v, w, 1.0);
    }));
  arms.push_back(
    simple_arm("sym(1.3)", [](Eigen::MatrixXd& v, const Eigen::VectorXd& w) {
      step_sym(v, w, 1.3);
    }));
  arms.push_back(
    simple_arm("sym(1.8)", [](Eigen::MatrixXd& v, const Eigen::VectorXd& w) {
      step_sym(v, w, 1.8);
    }));
  arms.push_back({ "logdomain",
                   [](Eigen::MatrixXd& v,
                      const Eigen::VectorXd& w,
                      int max_iter,
                      double tol) {
                     const ptrdiff_t m = v.rows();
                     const Eigen::MatrixXd lk = log_kernel(v, w);
                     const Eigen::MatrixXd lkt = log_kernel(v.transpose(), w);
                     Eigen::VectorXd x = Eigen::VectorXd::Zero(m);
                     Eigen::VectorXd y = Eigen::VectorXd::Zero(m);
                     int used = max_iter;
                     for (int k = 0; k < max_iter; ++k) {
                       dual_step(lk, lkt, x, y);
                       if (tol > 0.0) {
                         if (residual(from_potentials(v, x, y), w).worst() <
                             tol) {
                           used = k + 1;
                           break;
                         }
                       }
                     }
                     v = from_potentials(v, x, y);
                     return used;
                   } });
  arms.push_back(
    { "anderson(3)",
      [](Eigen::MatrixXd& v,
         const Eigen::VectorXd& w,
         int max_iter,
         double tol) { return run_anderson(v, w, max_iter, tol, 3); } });
  return arms;
}

double
asymmetry(const Arm& arm,
          const Eigen::MatrixXd& v,
          const Eigen::VectorXd& w,
          int k)
{
  Eigen::MatrixXd a = v, b = v.transpose();
  arm.run(a, w, k, -1.0);
  arm.run(b, w, k, -1.0);
  const double scale = a.array().abs().maxCoeff();
  return (a - b.transpose()).array().abs().maxCoeff() / scale;
}

// ------------------------------------------------------------- the fixtures

struct Case
{
  std::string label;
  BicopFamily family;
  double par;
  size_t n;
  size_t m;
  std::string method;
};

Eigen::MatrixXd
sim(BicopFamily family, double par, size_t n, int seed)
{
  Bicop bc(family, 0, Eigen::MatrixXd::Constant(1, 1, par));
  return bc.simulate(n, false, { seed });
}

//! Rescales a grid off the uniform-margin manifold, so the arms face the work
//! the raw local-likelihood surface would give them. Any positive row/column
//! rescaling leaves the doubly stochastic limit unchanged.
Eigen::MatrixXd
unnormalize(Eigen::MatrixXd v)
{
  const ptrdiff_t m = v.rows();
  Eigen::VectorXd s(m);
  for (ptrdiff_t i = 0; i < m; ++i) {
    s(i) = std::exp(std::sin(3.0 * static_cast<double>(i)));
  }
  v.array().colwise() *= s.array();
  v.array().rowwise() *= s.transpose().array();
  return v;
}

//! the grid a `tll` fit stores, i.e. after the shipped normalization
Eigen::MatrixXd
fitted_grid(const Case& cs, int seed, bool swap_args)
{
  Eigen::MatrixXd u = sim(cs.family, cs.par, cs.n, seed);
  if (swap_args)
    u.col(0).swap(u.col(1));
  FitControlsBicop controls({ BicopFamily::tll });
  controls.set_nonparametric_method(cs.method);
  controls.set_nonparametric_grid_size(cs.m);
  Bicop bc(BicopFamily::tll);
  bc.fit(u, controls);
  return bc.get_parameters();
}

std::vector<Case>
study_cases()
{
  std::vector<Case> cs;
  for (double rho : { 0.0, 0.3, 0.5, 0.7, 0.9, 0.95 })
    cs.push_back({ "gauss/rho=" + std::to_string(rho).substr(0, 4),
                   BicopFamily::gaussian,
                   rho,
                   1000,
                   30,
                   "constant" });
  cs.push_back(
    { "clayton/th=5", BicopFamily::clayton, 5.0, 1000, 30, "constant" });
  cs.push_back(
    { "gumbel/th=3", BicopFamily::gumbel, 3.0, 1000, 30, "constant" });
  cs.push_back(
    { "gauss/rho=0.9/n=200", BicopFamily::gaussian, 0.9, 200, 30, "constant" });
  cs.push_back({ "gauss/rho=0.9/n=5000",
                 BicopFamily::gaussian,
                 0.9,
                 5000,
                 30,
                 "constant" });
  cs.push_back({ "gauss/rho=0.9/quadratic",
                 BicopFamily::gaussian,
                 0.9,
                 1000,
                 30,
                 "quadratic" });
  cs.push_back(
    { "gauss/rho=0.9/m=50", BicopFamily::gaussian, 0.9, 1000, 50, "constant" });
  return cs;
}

} // namespace

// ============================================================== section 0

namespace {

bool
validate()
{
  std::cout << "== 0. validation (a failure here means the study is wrong, "
               "not the library)\n\n";
  bool ok = true;

  const size_t m = 30;
  const Eigen::VectorXd g = grid_points(m);
  const Eigen::VectorXd w = trapezoid_weights(g);
  const Case cs{ "gauss/rho=0.7", BicopFamily::gaussian, 0.7, 1000, m,
                 "constant" };
  const Eigen::MatrixXd v3 = fitted_grid(cs, 11, false);

  // the grid points, via the density the library reports for that grid
  {
    const Eigen::MatrixXd u =
      tools_stats::simulate_uniform(500, 2, true, { 3 });
    InterpolationGrid ig(g, v3, 0);
    Bicop bc(BicopFamily::tll, 0, v3);
    const double d = (ig.interpolate(u) - bc.pdf(u)).array().abs().maxCoeff();
    std::cout << "   grid points reproduce Bicop::pdf   max |diff| = " << d
              << (d < 1e-15 ? "   ok\n" : "   MISMATCH\n");
    ok = ok && (d < 1e-15);
  }

  // the shipped iteration, exactly: this validates the weights, the trapezoid
  // rule, the averaging and the guard at once
  {
    const Eigen::MatrixXd raw = unnormalize(v3);
    double worst = 0.0;
    for (int k : { 1, 2, 3, 5 }) {
      InterpolationGrid ig(g, raw, 0);
      ig.set_values(raw, k);
      Eigen::MatrixXd mine = raw;
      for (int i = 0; i < k; ++i)
        step_altsym(mine, w);
      worst = std::max(worst,
                       ((mine - ig.get_values()).array().abs() /
                        ig.get_values().array().abs().max(1e-300))
                         .maxCoeff());
    }
    std::cout << "   `altsym` reproduces set_values(V, k) max rel = " << worst
              << (worst < 1e-13 ? "   ok\n" : "   MISMATCH\n");
    ok = ok && (worst < 1e-13);
  }

  // the sweep map is memoryless, which licenses starting every arm from the
  // grid the library already normalized three times
  {
    const Eigen::MatrixXd raw = unnormalize(v3);
    double worst = 0.0;
    for (int k : { 2, 7, 22 }) {
      InterpolationGrid a(g, raw, 0), b(g, raw, 0);
      a.set_values(raw, 3 + k);
      b.set_values(raw, 3);
      const Eigen::MatrixXd mid = b.get_values();
      b.set_values(mid, k);
      worst = std::max(
        worst, (a.get_values() - b.get_values()).array().abs().maxCoeff());
    }
    std::cout << "   set_values(V,3+k) == (V,3) then (V,k)  max |diff| = "
              << worst << (worst < 1e-14 ? "   ok\n" : "   MISMATCH\n");
    ok = ok && (worst < 1e-14);
  }

  // every arm reaches the same doubly stochastic limit; if this fails, the
  // issue's diagnosis (the order dependence is a normalization artifact) is
  // wrong
  {
    Eigen::MatrixXd ref = v3;
    for (int i = 0; i < 2000; ++i)
      step_alt(ref, w);
    double worst = 0.0;
    for (const auto& arm : all_arms()) {
      Eigen::MatrixXd cur = v3;
      arm.run(cur, w, 2000, 1e-14);
      worst = std::max(worst, (cur - ref).array().abs().maxCoeff());
    }
    std::cout << "   every arm reaches the same limit   max |diff| = " << worst
              << (worst < 1e-9 ? "   ok\n" : "   MISMATCH\n");
    ok = ok && (worst < 1e-9);
  }

  // hygiene: the guard must never fire on a real fit, or every residual
  // below is meaningless
  {
    double lo = 1e300, hi = 0.0;
    size_t subnormal = 0, guarded = 0;
    for (const auto& cs2 : study_cases()) {
      const Eigen::MatrixXd v = fitted_grid(cs2, 11, false);
      const Eigen::VectorXd w2 = trapezoid_weights(grid_points(cs2.m));
      lo = std::min(lo, v.minCoeff());
      hi = std::max(hi, v.maxCoeff());
      for (ptrdiff_t i = 0; i < v.size(); ++i)
        subnormal +=
          (v(i) != 0.0) && (std::fabs(v(i)) < 2.2250738585072014e-308);
      guarded += ((v * w2).array() <= kMinMass).count();
      guarded += ((v.transpose() * w2).array() <= kMinMass).count();
    }
    std::cout << "   fitted values in [" << lo << ", " << hi << "], "
              << subnormal << " subnormal, guard fires " << guarded << " times"
              << (guarded == 0 ? "   ok\n" : "   WARNING\n");
  }

  std::cout << "\n";
  if (!ok)
    std::cout << "   *** MISMATCH(ES); the results below are not "
                 "trustworthy ***\n\n";
  return ok;
}

// ============================================================== section A

void
section_a()
{
  std::cout << "== A. residual and asymmetry after three passes\n\n"
            << "   R1 = max_i |int c(., u_i) du - 1| (first argument), R2 "
               "the same for\n"
            << "   the second, both after three `alt` passes -- which is what "
               "used to\n"
            << "   ship, and why R2 is exact. asym is\n"
            << "   max |N(V)' - N(V')| / max |N(V)| after three passes.\n\n";
  std::cout << std::left << std::setw(26) << "case" << std::right
            << std::setw(12) << "R1" << std::setw(12) << "R2" << std::setw(12)
            << "asym(alt)" << std::setw(12) << "asym(altsym)"
            << "\n";
  std::cout << std::string(74, '-') << "\n"
            << std::scientific << std::setprecision(2);

  const auto arms = all_arms();
  for (const auto& cs : study_cases()) {
    const Eigen::VectorXd g = grid_points(cs.m);
    const Eigen::VectorXd w = trapezoid_weights(g);
    const Eigen::MatrixXd v3 = unnormalize(fitted_grid(cs, 11, false));
    Eigen::MatrixXd after_alt = v3;
    for (int k = 0; k < 3; ++k) {
      step_alt(after_alt, w);
    }
    const Residual r = residual(after_alt, w);
    const double a_alt = asymmetry(arms[0], v3, w, 3);
    const double a_as = asymmetry(arms[1], v3, w, 3);
    std::cout << std::left << std::setw(26) << cs.label << std::right
              << std::setw(12) << r.r1 << std::setw(12) << r.r2 << std::setw(12)
              << a_alt << std::setw(12) << a_as << "\n";
  }
  std::cout << "\n";
}

// ============================================================== section B

void
section_b()
{
  std::cout << "== B. worst-margin residual after k further iterations, from "
               "the fitted grid\n\n";
  const std::vector<int> ks = { 1, 2, 3, 5, 10, 25, 50 };
  const auto arms = all_arms();

  for (const auto& cs : study_cases()) {
    if (cs.label.find("n=") != std::string::npos ||
        cs.label.find("quadratic") != std::string::npos ||
        cs.label.find("m=50") != std::string::npos)
      continue;
    const Eigen::VectorXd w = trapezoid_weights(grid_points(cs.m));
    const Eigen::MatrixXd v3 = unnormalize(fitted_grid(cs, 11, false));
    std::cout << cs.label << "  (residual at k=0 is " << residual(v3, w).worst()
              << ")\n";
    std::cout << std::left << std::setw(14) << "  arm" << std::right;
    for (int k : ks)
      std::cout << std::setw(11) << ("k=" + std::to_string(k));
    std::cout << "\n";
    for (const auto& arm : arms) {
      std::cout << std::left << "  " << std::setw(12) << arm.name << std::right;
      for (int k : ks) {
        Eigen::MatrixXd v = v3;
        arm.run(v, w, k, -1.0);
        std::cout << std::setw(11) << residual(v, w).worst();
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
}

// ============================================================== section C

void
section_c()
{
  std::cout << "== C. iterations to tolerance, and time, from the fitted "
               "grid\n\n"
            << "   ns/iter is measured on the m of the case; the cap is "
               "250. These are\n"
            << "   reference implementations, so they allocate where the "
               "shipped one does\n"
            << "   not -- compare iteration counts, and take timings as "
               "relative.\n\n";
  const auto arms = all_arms();
  std::cout << std::left << std::setw(26) << "case" << std::setw(14) << "arm"
            << std::right << std::setw(9) << "k(1e-8)" << std::setw(9)
            << "k(1e-10)" << std::setw(9) << "k(1e-12)" << std::setw(12)
            << "ns/iter"
            << "\n";
  std::cout << std::string(79, '-') << "\n";

  for (const auto& cs : study_cases()) {
    const Eigen::VectorXd w = trapezoid_weights(grid_points(cs.m));
    const Eigen::MatrixXd v3 = unnormalize(fitted_grid(cs, 11, false));
    for (const auto& arm : arms) {
      int k8 = 0, k10 = 0, k12 = 0;
      {
        Eigen::MatrixXd v = v3;
        k8 = arm.run(v, w, 250, 1e-8);
      }
      {
        Eigen::MatrixXd v = v3;
        k10 = arm.run(v, w, 250, 1e-10);
      }
      {
        Eigen::MatrixXd v = v3;
        k12 = arm.run(v, w, 250, 1e-12);
      }
      const int reps = 200;
      const auto t0 = std::chrono::steady_clock::now();
      for (int r = 0; r < reps; ++r) {
        Eigen::MatrixXd v = v3;
        arm.run(v, w, 10, -1.0);
      }
      const auto t1 = std::chrono::steady_clock::now();
      const double ns =
        std::chrono::duration<double, std::nano>(t1 - t0).count() /
        (reps * 10.0);
      std::cout << std::left << std::setw(26) << cs.label << std::setw(14)
                << arm.name << std::right << std::fixed << std::setprecision(0)
                << std::setw(9) << k8 << std::setw(9) << k10 << std::setw(9)
                << k12 << std::setw(12) << ns << std::scientific
                << std::setprecision(2) << "\n";
    }
  }
  std::cout << "\n";
}

// ============================================================== section D

//! `cond_cdf` with and without the 1e-4 knot clamp. The integrand is
//! piecewise linear, so the unclamped trapezoid is exact.
double
cond_cdf_ref(const Eigen::MatrixXd& v,
             const Eigen::VectorXd& g,
             double u_cond,
             double u,
             bool clamp)
{
  const ptrdiff_t m = g.size();
  ptrdiff_t i = 0;
  while ((i < m - 2) && (g(i + 1) <= u_cond))
    ++i;
  const double x2x = g(i + 1) - u_cond, xx1 = u_cond - g(i);
  const double x2x1 = g(i + 1) - g(i);
  auto knot = [&](ptrdiff_t j) {
    const double val = (v(i, j) * x2x + v(i + 1, j) * xx1) / x2x1;
    return clamp ? std::max(val, 1e-4) : std::max(val, 0.0);
  };
  double tmpint = 0.0, int1 = 0.0, v_k = knot(0);
  for (ptrdiff_t k = 0; k < m - 1; ++k) {
    const double v_k1 = knot(k + 1);
    int1 += (v_k1 + v_k) * (g(k + 1) - g(k)) / 2.0;
    if (!(u < g(k))) {
      if (u < g(k + 1))
        tmpint += (2 * v_k + (v_k1 - v_k) * (u - g(k)) / (g(k + 1) - g(k))) *
                  (u - g(k)) / 2.0;
      else
        tmpint += (v_k1 + v_k) * (g(k + 1) - g(k)) / 2.0;
    }
    v_k = v_k1;
  }
  return tmpint / std::max(int1, 1e-20);
}

void
section_d()
{
  std::cout << "== D. the 1e-4 knot clamp in cond_cdf / cond_quantile\n\n"
            << "   hfunc is clamped, pdf trims at 1e-20 and integrate_2d "
               "does not clamp\n"
            << "   at all, so hfunc is not the integral of the pdf the same "
               "object reports.\n"
            << "   err is max |hfunc1_clamped - hfunc1_exact| over a 21x21 "
               "interior lattice.\n\n";
  std::cout << std::left << std::setw(26) << "case" << std::right
            << std::setw(14) << "frac<1e-4" << std::setw(14) << "err"
            << std::setw(14) << "check vs lib"
            << "\n";
  std::cout << std::string(68, '-') << "\n";

  Eigen::MatrixXd lattice(21 * 21, 2);
  for (int a = 0; a < 21; ++a)
    for (int b = 0; b < 21; ++b) {
      lattice(a * 21 + b, 0) = 0.025 + 0.0475 * a;
      lattice(a * 21 + b, 1) = 0.025 + 0.0475 * b;
    }

  for (const auto& cs : study_cases()) {
    const Eigen::VectorXd g = grid_points(cs.m);
    const Eigen::MatrixXd v3 = fitted_grid(cs, 11, false);
    Bicop bc(BicopFamily::tll, 0, v3);
    const Eigen::VectorXd h = bc.hfunc1(lattice);

    double err = 0.0, check = 0.0;
    for (ptrdiff_t r = 0; r < lattice.rows(); ++r) {
      const double u1 = lattice(r, 0), u2 = lattice(r, 1);
      const double clamped = cond_cdf_ref(v3, g, u1, u2, true);
      const double exact = cond_cdf_ref(v3, g, u1, u2, false);
      err = std::max(err, std::fabs(clamped - exact));
      check = std::max(check, std::fabs(clamped - h(r)));
    }
    const double frac = static_cast<double>((v3.array() < 1e-4).count()) /
                        static_cast<double>(v3.size());
    std::cout << std::left << std::setw(26) << cs.label << std::right
              << std::fixed << std::setprecision(1) << std::setw(13)
              << 100.0 * frac << "%" << std::scientific << std::setprecision(2)
              << std::setw(14) << err << std::setw(14) << check << "\n";
  }
  std::cout << "\n   (check vs lib is the study's clamped restatement against "
               "Bicop::hfunc1;\n    it must be at machine precision or the "
               "err column means nothing)\n\n";
}

// ============================================================== section E

void
section_e()
{
  std::cout << "== E. does converging change the estimate? ISE against the "
               "truth\n\n"
            << "   30 replicates, gaussian, n=1000, m=30, constant. ISE = "
               "int (chat - c)^2,\n"
            << "   trapezoid on the grid. Every arm shares one limit "
               "(section 0), so this\n"
            << "   is about the iteration count, not the scheme.\n\n";
  std::cout << std::left << std::setw(12) << "rho" << std::right
            << std::setw(14) << "ISE(x3)" << std::setw(12) << "x5"
            << std::setw(12) << "x10" << std::setw(12) << "x25" << std::setw(12)
            << "converged"
            << "\n";
  std::cout << std::string(74, '-') << "\n";

  const size_t m = 30, n = 1000;
  const Eigen::VectorXd g = grid_points(m);
  const Eigen::VectorXd w = trapezoid_weights(g);

  for (double rho : { 0.3, 0.5, 0.7, 0.9, 0.95 }) {
    // the truth on the grid, normalized the same way, so the comparison is
    // between two doubly stochastic surfaces
    Eigen::MatrixXd truth(m, m);
    Eigen::VectorXd z(m);
    for (size_t i = 0; i < m; ++i)
      z(i) =
        -3.25 + static_cast<double>(i) * (6.5 / static_cast<double>(m - 1));
    for (size_t i = 0; i < m; ++i)
      for (size_t j = 0; j < m; ++j)
        truth(i, j) = std::exp(-(rho * rho * (z(i) * z(i) + z(j) * z(j)) -
                                 2 * rho * z(i) * z(j)) /
                               (2 * (1 - rho * rho))) /
                      std::sqrt(1 - rho * rho);
    for (int k = 0; k < 500; ++k)
      step_alt(truth, w);

    std::vector<double> ise(5, 0.0);
    const int reps = 30;
    for (int r = 0; r < reps; ++r) {
      Eigen::MatrixXd u = sim(BicopFamily::gaussian, rho, n, 100 + r);
      FitControlsBicop controls({ BicopFamily::tll });
      controls.set_nonparametric_grid_size(m);
      Bicop bc(BicopFamily::tll);
      bc.fit(u, controls);
      const Eigen::MatrixXd v3 = bc.get_parameters();
      const Eigen::MatrixXd raw = unnormalize(v3);
      const std::vector<int> extra = { 3, 5, 10, 25, 500 };
      for (size_t s = 0; s < extra.size(); ++s) {
        Eigen::MatrixXd v = raw;
        for (int k = 0; k < extra[s]; ++k)
          step_altsym(v, w);
        const Eigen::MatrixXd d = (v - truth).array().square();
        ise[s] += w.transpose() * d * w;
      }
    }
    std::cout << std::left << std::fixed << std::setprecision(2)
              << std::setw(12) << rho << std::right << std::scientific
              << std::setprecision(4) << std::setw(14) << ise[0] / reps;
    for (size_t s = 1; s < ise.size(); ++s)
      std::cout << std::fixed << std::setprecision(1) << std::setw(11)
                << 100.0 * (ise[s] / ise[0] - 1.0) << "%";
    std::cout << "\n";
  }
  std::cout << "\n";
}

} // namespace

int
main()
{
  std::cout << "\nInterpolationGrid margin normalization: accuracy and cost\n"
            << std::string(74, '=') << "\n\n";
  const bool ok = validate();
  section_a();
  section_b();
  section_c();
  section_d();
  section_e();
  if (!ok)
    std::cout << "*** section 0 reported a MISMATCH; results above are not "
                 "trustworthy ***\n";
  return 0;
}
