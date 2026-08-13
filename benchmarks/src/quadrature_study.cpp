// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

// Accuracy and speed of candidate quadrature rules for the four
// parameters_to_tau integrals (BB6, BB7, BB8, Tawn).
//
// The integrands are restated here, templated on the scalar type, so that
// references can be computed in 50-digit arithmetic. Keep them in sync with
// bb6.ipp, bb7.ipp, bb8.ipp and extreme_value.ipp/tawn.ipp.
//
//   build:  cmake --preset bench && cmake --build --preset bench --target
//   quadrature_study run:    build/bench/bin/quadrature_study

#include <vinecopulib/bicop/class.hpp>

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/numeric/odeint.hpp>

#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using Real50 = boost::multiprecision::cpp_bin_float_50;

// The bounds the library integrates between: the integrands are singular at one
// or both endpoints, so the interval is clipped.
template<typename Real>
Real
lb()
{
  return Real(1) / 1000000000000LL; // 1e-12
}

// ---------------------------------------------------------------- integrands

template<typename Real>
struct Bb6
{
  Real theta, delta;
  Real operator()(Real v) const
  {
    using std::log1p;
    using std::pow;
    Real omv = 1 - v;
    Real res = -4 * (omv - pow(omv, -theta) + pow(omv, -theta) * v);
    return 1 / (delta * theta) * log1p(-pow(omv, theta)) * res;
  }
  static constexpr const char* name = "bb6";
  Real split() const { return Real(1) / 2; }
};

template<typename Real>
struct Bb7
{
  Real theta, delta;
  Real operator()(Real v) const
  {
    using std::exp;
    using std::expm1;
    using std::log1p;
    Real log1mv = log1p(-v);
    Real tmp = exp(theta * log1mv);
    Real log1mtmp = log1p(-tmp);
    Real res = -4 * expm1(-delta * log1mtmp) / (theta * delta);
    return res / exp((theta - 1) * log1mv + (-delta - 1) * log1mtmp);
  }
  static constexpr const char* name = "bb7";
  Real split() const { return Real(1) / 2; }
};

template<typename Real>
struct Bb8
{
  Real theta, delta;
  Real operator()(Real t) const
  {
    using std::exp;
    using std::log1p;
    Real log1mtd = log1p(-t * delta);
    Real num = log1p(-exp(theta * log1mtd));
    Real den = log1p(-exp(theta * log1p(-delta)));
    return (num - den) * (1 - t * delta - exp((1 - theta) * log1mtd));
  }
  static constexpr const char* name = "bb8";
  Real split() const { return Real(1) / 2; }
};

// t (1 - t) A''(t) / A(t) with Tawn's Pickands dependence function.
template<typename Real>
struct Tawn
{
  Real psi1, psi2, theta;
  Real operator()(Real t) const
  {
    using std::exp;
    using std::log;
    using std::log1p;
    using std::max;
    using std::min;
    using std::pow;
    Real A = (1 - psi1) * (1 - t) + (1 - psi2) * t +
             pow(pow(psi2 * t, theta) + pow(psi1 * (1 - t), theta), 1 / theta);
    // A'' = (theta - 1) (x + y)^(1/theta - 2) x y / (t (1-t))^2 with
    // x = (psi2 t)^theta, y = (psi1 (1-t))^theta; see tawn.ipp.
    Real lx = theta * log(psi2 * t), ly = theta * log(psi1 * (1 - t));
    Real hi = max(lx, ly);
    Real ls = hi + log1p(exp(min(lx, ly) - hi));
    Real App = exp(log(theta - 1) + (1 / theta - 2) * ls + lx + ly -
                   2 * log(t) - 2 * log1p(-t));
    return t * (1 - t) * App / A;
  }
  static constexpr const char* name = "tawn";
  // A'' peaks where (psi2 t)^theta and (psi1 (1-t))^theta balance.
  Real split() const { return psi1 / (psi1 + psi2); }
};

// ---------------------------------------------------------------------- arms

// odeint over [a, b], so the split variant can reuse it.
template<typename F>
double
arm_odeint_on(const F& f, double a, double b)
{
  boost::numeric::odeint::runge_kutta_dopri5<double> stepper;
  double x = 0.0;
  auto ifunc = [&f](const double, double& dxdt, const double t) {
    dxdt = f(t);
  };
  integrate_adaptive(
    boost::numeric::odeint::make_controlled(1e-9, 1e-9, stepper),
    ifunc,
    x,
    a,
    b,
    1e-12);
  return x;
}

template<typename F>
double
arm_odeint_split(const F& f)
{
  const double s = static_cast<double>(f.split());
  return arm_odeint_on(f, 1e-12, s) + arm_odeint_on(f, s, 1.0 - 1e-12);
}

template<typename F>
double
arm_tanh_sinh_split(const F& f)
{
  const double s = static_cast<double>(f.split());
  boost::math::quadrature::tanh_sinh<double> integrator;
  return integrator.integrate(f, 1e-12, s) +
         integrator.integrate(f, s, 1.0 - 1e-12);
}

// The rule currently in tools_integration.hpp.
template<typename F>
double
arm_odeint(const F& f)
{
  boost::numeric::odeint::runge_kutta_dopri5<double> stepper;
  const double a = 1e-12, b = 1.0 - 1e-12, tol = 1e-9;
  double x = 0.0;
  auto ifunc = [&f](const double, double& dxdt, const double t) {
    dxdt = f(t);
  };
  integrate_adaptive(boost::numeric::odeint::make_controlled(tol, tol, stepper),
                     ifunc,
                     x,
                     a,
                     b,
                     a);
  return x;
}

// Tolerance 1e-10, not 1e-12: two of the four integrands cannot reach 1e-12 in
// double, and asking for it makes the adaptive rule subdivide to max_depth.
template<unsigned N, typename F>
double
arm_gauss_kronrod(const F& f)
{
  using boost::math::quadrature::gauss_kronrod;
  double err = 0;
  return gauss_kronrod<double, N>::integrate(
    f, 1e-12, 1.0 - 1e-12, 12, 1e-10, &err);
}

template<typename F>
double
arm_tanh_sinh(const F& f)
{
  boost::math::quadrature::tanh_sinh<double> integrator;
  return integrator.integrate(f, 1e-12, 1.0 - 1e-12);
}

// The 2019 attempt (archive/add-quadrature): a fixed-order Gauss-Legendre rule.
template<typename F>
double
arm_gauss_legendre20(const F& f)
{
  return boost::math::quadrature::gauss<double, 20>::integrate(
    f, 1e-12, 1.0 - 1e-12);
}

// Reference: adaptive Gauss-Kronrod at 50 digits, deep enough that the endpoint
// singularities are resolved far below double precision. Cross-checked against
// a 50-digit tanh_sinh wherever that converges -- two independent
// high-precision methods agreeing is the evidence; either alone is not.
std::string
fmt_sci(double x)
{
  std::ostringstream os;
  os << std::scientific << std::setprecision(2) << x;
  return os.str();
}

template<typename F>
Real50
reference(const F& f, double* cross_check_rel)
{
  using boost::math::quadrature::gauss_kronrod;
  const Real50 a = lb<Real50>(), b = Real50(1) - a;
  // Split at the peak here as well: an unsplit reference misses a narrow
  // interior peak just as badly as the arms being judged do.
  const Real50 s = f.split();
  const bool do_split = s > a && s < b;
  Real50 err = 0;
  auto gk = [&](const Real50& x, const Real50& y) {
    return gauss_kronrod<Real50, 61>::integrate(
      f, x, y, 25, Real50(1e-30), &err);
  };
  const Real50 value = do_split ? gk(a, s) + gk(s, b) : gk(a, b);

  // Cross-check the 50-digit GK61 value against 50-digit tanh-sinh, and report
  // how far apart they are: the reference is only as good as that agreement.
  *cross_check_rel = std::numeric_limits<double>::quiet_NaN();
  try {
    boost::math::quadrature::tanh_sinh<Real50> ts(25);
    const Real50 other = do_split
                           ? ts.integrate(f, a, s) + ts.integrate(f, s, b)
                           : ts.integrate(f, a, b);
    const Real50 scale = abs(value) > Real50(1e-300) ? abs(value) : Real50(1);
    *cross_check_rel = static_cast<double>(abs(other - value) / scale);
  } catch (const std::exception&) {
    // tanh_sinh declines on some endpoint singularities; the GK value stands.
  }
  return value;
}

// ----------------------------------------------------------------- validation

// The restated integrand must reproduce the library. `transform` applies the
// family's outer expression, e.g. 1 + I for BB6/BB7. Without this the study
// would happily measure a mis-transcribed integrand.
template<typename F, typename Transform>
bool
matches_library(vinecopulib::BicopFamily family,
                const Eigen::VectorXd& par,
                const F& f,
                const Transform& transform,
                bool split = false)
{
  vinecopulib::Bicop bicop(family, 0, par);
  const double library = bicop.get_tau();
  // Use the rule the library uses, so a mismatch means the restated integrand
  // disagrees with the shipped one rather than the two rules disagreeing.
  const double restated =
    transform(split ? arm_tanh_sinh_split(f) : arm_tanh_sinh(f));
  const double scale = std::max(1.0, std::abs(library));
  const bool ok = std::abs(library - restated) / scale < 1e-6;
  if (!ok) {
    std::cout << "  MISMATCH vs library: parameters_to_tau=" << library
              << " restated=" << restated << " (par " << par.transpose()
              << ")\n";
  }
  return ok;
}

// --------------------------------------------------------------------- driver

struct Stats
{
  double max_rel = 0.0;
  double total_ns = 0.0;
  size_t n = 0;
  std::string worst_case;
};

template<typename MakeDouble, typename MakeRef>
void
evaluate(const std::string& family,
         const std::vector<std::string>& labels,
         const MakeDouble& make_double,
         const MakeRef& make_ref,
         std::vector<Stats>& stats,
         const std::vector<std::string>& arm_names)
{
  const size_t n_cases = labels.size();
  for (size_t c = 0; c < n_cases; ++c) {
    std::cout << "  [" << family << " " << (c + 1) << "/" << n_cases << "] "
              << labels[c] << std::flush;
    double cross_check_rel = 0.0;
    const Real50 ref = reference(make_ref(c), &cross_check_rel);
    std::cout << (std::isnan(cross_check_rel)
                    ? std::string("  (ref cross-check unavailable)\n")
                    : "  (ref cross-check " + fmt_sci(cross_check_rel) + ")\n")
              << std::flush;
    auto fd = make_double(c);

    std::vector<std::function<double()>> arms = {
      [&] { return arm_odeint(fd); },
      [&] { return arm_odeint_split(fd); },
      [&] { return arm_tanh_sinh(fd); },
      [&] { return arm_tanh_sinh_split(fd); },
      [&] { return arm_gauss_kronrod<15>(fd); },
      [&] { return arm_gauss_legendre20(fd); },
    };

    for (size_t a = 0; a < arms.size(); ++a) {
      // one untimed call so any lazily built tables are warm
      double value = 0.0;
      try {
        value = arms[a]();
      } catch (const std::exception&) {
        stats[a].max_rel = std::numeric_limits<double>::infinity();
        stats[a].worst_case = labels[c] + " (threw)";
        continue;
      }
      const int reps = 10;
      auto t0 = std::chrono::steady_clock::now();
      for (int r = 0; r < reps; ++r)
        value = arms[a]();
      auto t1 = std::chrono::steady_clock::now();

      const Real50 denom = abs(ref) > Real50(1e-300) ? abs(ref) : Real50(1);
      const double rel = static_cast<double>(abs(Real50(value) - ref) / denom);
      if (rel > stats[a].max_rel) {
        stats[a].max_rel = rel;
        stats[a].worst_case = labels[c];
      }
      stats[a].total_ns +=
        std::chrono::duration<double, std::nano>(t1 - t0).count() / reps;
      stats[a].n += 1;
    }
  }

  std::cout << "\n=== " << family << " (" << n_cases
            << " parameter sets) ===\n";
  std::cout << std::left << std::setw(16) << "arm" << std::right
            << std::setw(14) << "max rel err" << std::setw(12) << "ns/call"
            << "   worst case\n";
  for (size_t a = 0; a < arm_names.size(); ++a) {
    std::cout << std::left << std::setw(16) << arm_names[a] << std::right
              << std::setw(14) << std::scientific << std::setprecision(3)
              << stats[a].max_rel << std::setw(12) << std::fixed
              << std::setprecision(0)
              << (stats[a].n ? stats[a].total_ns / stats[a].n : 0.0) << "   "
              << stats[a].worst_case << "\n";
  }
}

int
main()
{
  const std::vector<std::string> arm_names = { "odeint",    "odeint split",
                                               "tanh_sinh", "tanh_sinh split",
                                               "GK15",      "GL20 (2019)" };

  // Parameter grids inside each family's admissible region, including the
  // boundary corners where the integrands are most singular.
  // bb6: theta in [1, 6], delta in [1, 8]
  const std::vector<std::pair<double, double>> bb6 = {
    { 1.0, 1.0 }, { 1.0001, 1.0001 }, { 1.5, 1.5 }, { 3.0, 4.0 },
    { 6.0, 8.0 }, { 1.0, 8.0 },       { 6.0, 1.0 }
  };
  // bb7: theta in [1, 6], delta in [0.01, 25]
  const std::vector<std::pair<double, double>> bb7 = {
    { 1.0, 0.01 }, { 1.0001, 0.01 }, { 1.5, 1.0 }, { 3.0, 5.0 },
    { 6.0, 25.0 }, { 1.0, 25.0 },    { 6.0, 0.01 }
  };
  // bb8: theta in [1, 8], delta in [1e-4, 1]
  const std::vector<std::pair<double, double>> bb8 = {
    { 1.0, 1e-4 }, { 1.0001, 0.9999 }, { 3.0, 0.8 }, { 8.0, 1.0 },
    { 8.0, 1e-4 }, { 1.0, 1.0 },       { 4.0, 0.5 }
  };
  // tawn: psi1, psi2 in (0, 1], theta in [1, 60]
  const std::vector<std::array<double, 3>> tawn = {
    { 0.8, 0.5, 2.0 },
    { 0.5, 0.5, 1.0001 },
    { 0.999, 0.999, 60.0 },
    { 0.01, 0.99, 5.0 },
    { 1.0, 1.0, 10.0 },
    { 0.3, 0.7, 30.0 },
    { 0.999, 0.001, 1.5 },
    { 0.5, 0.5, 60.0 },
    // strongly asymmetric psi: A'' peaks at psi1 / (psi1 + psi2), decades away
    // from 1/2, and narrows as theta grows
    { 1e-4, 0.5, 60.0 },
    { 1e-3, 0.5, 10.0 },
    { 0.01, 0.5, 60.0 },
    { 1e-6, 0.5, 60.0 }
  };

  auto label2 = [](const char* a, const char* b, double x, double y) {
    std::ostringstream os;
    os << a << "=" << x << ", " << b << "=" << y;
    return os.str();
  };

  std::cout << "Validating the restated integrands against Bicop::get_tau()\n";
  size_t mismatches = 0;
  {
    for (auto& p : bb6) {
      Eigen::VectorXd par(2);
      par << p.first, p.second;
      mismatches += !matches_library(vinecopulib::BicopFamily::bb6,
                                     par,
                                     Bb6<double>{ p.first, p.second },
                                     [](double i) { return 1 + i; });
    }
    for (auto& p : bb7) {
      Eigen::VectorXd par(2);
      par << p.first, p.second;
      mismatches += !matches_library(vinecopulib::BicopFamily::bb7,
                                     par,
                                     Bb7<double>{ p.first, p.second },
                                     [](double i) { return 1 + i; });
    }
    for (auto& p : bb8) {
      Eigen::VectorXd par(2);
      par << p.first, p.second;
      const double th = p.first, de = p.second;
      mismatches +=
        !matches_library(vinecopulib::BicopFamily::bb8,
                         par,
                         Bb8<double>{ th, de },
                         [th, de](double i) { return 1 - 4 / (de * th) * i; });
    }
    for (auto& p : tawn) {
      Eigen::VectorXd par(3);
      par << p[0], p[1], p[2];
      mismatches += !matches_library(
        vinecopulib::BicopFamily::tawn,
        par,
        Tawn<double>{ p[0], p[1], p[2] },
        [](double i) { return i; },
        /* split = */ true);
    }
    std::cout << (mismatches ? "  -> " : "  -> no mismatches; ")
              << (mismatches
                    ? std::to_string(mismatches) +
                        " MISMATCH(ES); results below are not trustworthy\n"
                    : std::string("integrands agree with the library\n"));
  }

  {
    std::vector<std::string> labels;
    for (auto& p : bb6)
      labels.push_back(label2("theta", "delta", p.first, p.second));
    std::vector<Stats> stats(arm_names.size());
    evaluate(
      "BB6",
      labels,
      [&](size_t c) {
        return Bb6<double>{ bb6[c].first, bb6[c].second };
      },
      [&](size_t c) {
        return Bb6<Real50>{ Real50(bb6[c].first), Real50(bb6[c].second) };
      },
      stats,
      arm_names);
  }
  {
    std::vector<std::string> labels;
    for (auto& p : bb7)
      labels.push_back(label2("theta", "delta", p.first, p.second));
    std::vector<Stats> stats(arm_names.size());
    evaluate(
      "BB7",
      labels,
      [&](size_t c) {
        return Bb7<double>{ bb7[c].first, bb7[c].second };
      },
      [&](size_t c) {
        return Bb7<Real50>{ Real50(bb7[c].first), Real50(bb7[c].second) };
      },
      stats,
      arm_names);
  }
  {
    std::vector<std::string> labels;
    for (auto& p : bb8)
      labels.push_back(label2("theta", "delta", p.first, p.second));
    std::vector<Stats> stats(arm_names.size());
    evaluate(
      "BB8",
      labels,
      [&](size_t c) {
        return Bb8<double>{ bb8[c].first, bb8[c].second };
      },
      [&](size_t c) {
        return Bb8<Real50>{ Real50(bb8[c].first), Real50(bb8[c].second) };
      },
      stats,
      arm_names);
  }
  {
    std::vector<std::string> labels;
    for (auto& p : tawn) {
      std::ostringstream os;
      os << "psi1=" << p[0] << ", psi2=" << p[1] << ", theta=" << p[2];
      labels.push_back(os.str());
    }
    std::vector<Stats> stats(arm_names.size());
    evaluate(
      "Tawn",
      labels,
      [&](size_t c) {
        return Tawn<double>{ tawn[c][0], tawn[c][1], tawn[c][2] };
      },
      [&](size_t c) {
        return Tawn<Real50>{ Real50(tawn[c][0]),
                             Real50(tawn[c][1]),
                             Real50(tawn[c][2]) };
      },
      stats,
      arm_names);
  }

  return 0;
}
