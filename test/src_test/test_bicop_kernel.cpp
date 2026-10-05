// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "include/kernel_test.hpp"
#include "include/r_parity.hpp"
#include "include/test_utils.hpp"
#include <limits>

namespace test_bicop_kernel {
using namespace vinecopulib;
using test_utils::all_close;

TEST_P(TrafokernelTest, sanity_checks)
{
  auto values = bicop_.get_parameters();
  EXPECT_ANY_THROW(bicop_.set_parameters(values.block(0, 0, 30, 1)));
  EXPECT_ANY_THROW(bicop_.set_parameters(values.block(0, 0, 1, 30)));
  EXPECT_NO_THROW(bicop_.set_parameters(values.block(0, 0, 10, 10)));
  EXPECT_ANY_THROW(bicop_.set_parameters(values.block(0, 0, 2, 2)));
  EXPECT_ANY_THROW(bicop_.set_parameters(-1 * values));
}

TEST_P(TrafokernelTest, fit)
{
  bicop_.fit(u, controls);
  EXPECT_EQ(bicop_.get_parameters().rows(), 30);
  EXPECT_EQ(bicop_.get_parameters().cols(), 30);

  // make sure that npars are copied
  auto bicop_cpy = bicop_;
  EXPECT_GT(bicop_cpy.get_npars(), 1.0);

  // catches bugs when n < (grid size)^2
  controls.set_weights(Eigen::VectorXd::Constant(20, 1.0));
  bicop_.fit(u.topRows(20), controls);

  controls.set_nonparametric_grid_size(10);
  bicop_.fit(u, controls);
  EXPECT_EQ(bicop_.get_parameters().rows(), 10);
  EXPECT_EQ(bicop_.get_parameters().cols(), 10);
}

TEST_P(TrafokernelTest, serialization)
{
  test_r_parity::RTempDir dir;
  bicop_.to_file(dir.file("bicop.json"));
  Bicop pc(dir.file("bicop.json"));

  EXPECT_EQ(bicop_.get_rotation(), pc.get_rotation());
  EXPECT_EQ(bicop_.get_family_name(), pc.get_family_name());
  EXPECT_EQ(bicop_.get_var_types(), pc.get_var_types());
  EXPECT_EQ(bicop_.get_npars(), pc.get_npars());
  ASSERT_TRUE(
    all_close(bicop_.get_parameters(), pc.get_parameters(), 1e-4, 1e-4));
}

TEST_P(TrafokernelTest, eval_funcs)
{
  bicop_.fit(u, controls);

  EXPECT_GE(bicop_.pdf(u).minCoeff(), 0.0);
  EXPECT_GE(bicop_.cdf(u).minCoeff(), 0.0);
  EXPECT_GE(bicop_.hfunc1(u).minCoeff(), 0.0);
  EXPECT_GE(bicop_.hfunc2(u).minCoeff(), 0.0);
  EXPECT_GE(bicop_.hinv1(u).minCoeff(), 0.0);
  EXPECT_GE(bicop_.hinv2(u).minCoeff(), 0.0);
  EXPECT_LE(bicop_.cdf(u).maxCoeff(), 1.0);
  EXPECT_LE(bicop_.hfunc1(u).maxCoeff(), 1.0);
  EXPECT_LE(bicop_.hfunc2(u).maxCoeff(), 1.0);
  EXPECT_LE(bicop_.hinv1(u).maxCoeff(), 1.0);
  EXPECT_LE(bicop_.hinv2(u).maxCoeff(), 1.0);
  EXPECT_GE(bicop_.get_npars(), 0.0);
  EXPECT_LE(bicop_.get_npars(), 100.0);
  EXPECT_NEAR(bicop_.get_loglik(), bicop_.loglik(u), 1e-5);

  u(0, 0) = std::numeric_limits<double>::quiet_NaN();
  u(1, 1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_NO_THROW(bicop_.pdf(u.block(0, 0, 10, 2)));
  EXPECT_TRUE(bicop_.pdf(u.block(0, 0, 1, 2)).array().isNaN()(0));
  EXPECT_NO_THROW(bicop_.cdf(u.block(0, 0, 10, 2)));
  EXPECT_TRUE(bicop_.cdf(u.block(0, 0, 1, 2)).array().isNaN()(0));
  EXPECT_NO_THROW(bicop_.hfunc1(u.block(0, 0, 10, 2)));
  EXPECT_TRUE(bicop_.hfunc1(u.block(0, 0, 1, 2)).array().isNaN()(0));
  EXPECT_NO_THROW(bicop_.hinv1(u.block(0, 0, 10, 2)));
  EXPECT_TRUE(bicop_.hinv1(u.block(0, 0, 1, 2)).array().isNaN()(0));
  EXPECT_NO_THROW(bicop_.hfunc2(u.block(0, 0, 10, 2)));
  EXPECT_TRUE(bicop_.hfunc2(u.block(0, 0, 1, 2)).array().isNaN()(0));
  EXPECT_NO_THROW(bicop_.hinv2(u.block(0, 0, 10, 2)));
  EXPECT_TRUE(bicop_.hinv2(u.block(0, 0, 1, 2)).array().isNaN()(0));
  EXPECT_NO_THROW(bicop_.loglik(u.block(0, 0, 10, 2)));
}

TEST_P(TrafokernelTest, select)
{
  auto newcop = Bicop(u, controls);
  EXPECT_NEAR(newcop.loglik(u), newcop.get_loglik(), 1e-5);
  EXPECT_EQ(newcop.get_family(), BicopFamily::tll);
}

TEST_P(TrafokernelTest, flip)
{
  auto pdf = bicop_.pdf(u);
  u.col(0).swap(u.col(1));
  bicop_.flip();
  auto pdf_flipped = bicop_.pdf(u);
  EXPECT_TRUE(all_close(pdf, pdf_flipped, 1e-10, 1e-10));
}

TEST_P(TrafokernelTest, tau)
{
  double tau = bicop_.parameters_to_tau(bicop_.get_parameters());
  EXPECT_GE(tau, -1.0);
  EXPECT_LE(tau, 1.0);
}

TEST_P(TrafokernelTest, reset)
{
  // this is essentially what we do when converting between C++ and R objects
  auto cop = Bicop(u, controls);
  auto cop_new = Bicop(BicopFamily::tll, 0, cop.get_parameters());
  EXPECT_EQ(cop.get_parameters(), cop_new.get_parameters());
  EXPECT_EQ(cop.get_family(), cop_new.get_family());
  EXPECT_EQ(cop.get_rotation(), cop_new.get_rotation());
  EXPECT_NEAR(cop.loglik(u), cop_new.loglik(u), 1e-8);
}

INSTANTIATE_TEST_SUITE_P(TrafokernelTest,
                         TrafokernelTest,
                         ::testing::Values("constant", "linear", "quadratic"));

// Golden-value regression guards: fixed-seed fits and evaluations pinned to
// the current implementation. Floating-point reassociation stays well within
// the 1e-8 tolerances; anything larger signals an unintended behavior change.
// Re-baselined when the margin normalization and the conditional cdf changed.
namespace {

//! The grid a TLL fit stores its density on: `KernelBicop::make_normal_grid`
//! with the endpoint snapping `InterpolationGrid`'s constructor applies.
Eigen::VectorXd
tll_grid_points(int m)
{
  Eigen::VectorXd z(m);
  for (int i = 0; i < m; ++i) {
    z(i) = -3.25 + i * (6.5 / (m - 1));
  }
  Eigen::VectorXd g = tools_stats::pnorm(z);
  g(0) = 0.0;
  g(m - 1) = 1.0;
  return g;
}

//! Weights such that `v.dot(w)` is the trapezoid integral over [0, 1] of the
//! piecewise linear function through `(g, v)` -- exact, since that is what the
//! grid interpolates.
Eigen::VectorXd
tll_grid_weights(const Eigen::VectorXd& g)
{
  const int m = static_cast<int>(g.size());
  Eigen::VectorXd w(m);
  w(0) = (g(1) - g(0)) / 2.0;
  w.segment(1, m - 2) = (g.tail(m - 2) - g.head(m - 2)) / 2.0;
  w(m - 1) = (g(m - 1) - g(m - 2)) / 2.0;
  return w;
}

Bicop
fit_tll(double rho, size_t n, int seed = 5)
{
  Bicop gauss(BicopFamily::gaussian, 0, Eigen::MatrixXd::Constant(1, 1, rho));
  FitControlsBicop controls({ BicopFamily::tll });
  Bicop tll(BicopFamily::tll);
  tll.fit(gauss.simulate(n, false, { seed }), controls);
  return tll;
}

} // namespace

// A fitted grid must be a copula density: both margins integrate to one. The
// margins are checked separately and against each other, because rescaling
// rows and then columns satisfies one of them exactly and leaves the other
// carrying the entire residual.
TEST(test_bicop_kernel_margins, both_margins_are_uniform)
{
  const Eigen::VectorXd g = tll_grid_points(30);
  const Eigen::VectorXd w = tll_grid_weights(g);
  for (double rho : { 0.3, 0.7, 0.9 }) {
    const Eigen::MatrixXd v = fit_tll(rho, 1000).get_parameters();
    const double r1 = ((v * w).array() - 1.0).abs().maxCoeff();
    const double r2 = ((v.transpose() * w).array() - 1.0).abs().maxCoeff();
    EXPECT_LT(r1, 1e-5) << "rho = " << rho;
    EXPECT_LT(r2, 1e-5) << "rho = " << rho;
    EXPECT_LT(std::max(r1, r2), 10.0 * std::min(r1, r2))
      << "margins are not equally accurate at rho = " << rho << ": " << r1
      << " vs " << r2;
  }
}

// Fitting (u1, u2) and fitting (u2, u1) must give transposed grids, or a vine
// that reuses a pair copula in the opposite orientation gets a different model
// from one that refits it.
TEST(test_bicop_kernel_margins, fit_does_not_depend_on_argument_order)
{
  for (double rho : { 0.3, 0.7, 0.9 }) {
    Bicop gauss(BicopFamily::gaussian, 0, Eigen::MatrixXd::Constant(1, 1, rho));
    Eigen::MatrixXd u = gauss.simulate(1000, false, { 5 });
    Eigen::MatrixXd u_swapped = u;
    u_swapped.col(0).swap(u_swapped.col(1));

    FitControlsBicop controls({ BicopFamily::tll });
    Bicop a(BicopFamily::tll);
    a.fit(u, controls);
    Bicop b(BicopFamily::tll);
    b.fit(u_swapped, controls);
    b.flip();

    EXPECT_TRUE(all_close(a.get_parameters(), b.get_parameters(), 1e-12, 1e-12))
      << "rho = " << rho;
  }
}

// `flip` transposes without renormalizing, which is only correct if the two
// margins are interchangeable.
TEST(test_bicop_kernel_margins, flip_is_an_involution)
{
  Bicop tll = fit_tll(0.7, 1000);
  const Eigen::MatrixXd before = tll.get_parameters();
  tll.flip();
  tll.flip();
  EXPECT_EQ(before, tll.get_parameters());
}

// The default grid is already uniform, so normalizing it must be a no-op.
TEST(test_bicop_kernel_margins, independence_grid_is_exactly_uniform)
{
  const Eigen::MatrixXd v = Bicop(BicopFamily::tll).get_parameters();
  EXPECT_EQ(v, Eigen::MatrixXd::Constant(v.rows(), v.cols(), 1.0));
}

// `hfunc1` must be the conditional cdf of the density the same object reports.
// It used to floor the interpolated density at 1e-4 while `pdf` floored it at
// 1e-20, which on a strongly dependent fit moves the h-function by ~1e-4.
TEST(test_bicop_kernel_margins, hfunc1_is_the_conditional_cdf_of_the_pdf)
{
  const int m = 30;
  const Eigen::VectorXd g = tll_grid_points(m);
  for (double rho : { 0.3, 0.7, 0.9 }) {
    const Bicop tll = fit_tll(rho, 1000);
    for (double u1 : { 0.1, 0.5, 0.9 }) {
      // hfunc1 conditions on the first argument and integrates the second;
      // the density is piecewise linear in u2 along a fixed u1, so trapezoid
      // integration over the grid knots is exact
      Eigen::MatrixXd line(m, 2);
      line.col(0).setConstant(u1);
      line.col(1) = g;
      const Eigen::VectorXd c = tll.pdf(line);

      Eigen::VectorXd cum(m);
      cum(0) = 0.0;
      for (int k = 0; k < m - 1; ++k) {
        cum(k + 1) = cum(k) + (c(k + 1) + c(k)) * (g(k + 1) - g(k)) / 2.0;
      }
      for (int k = 1; k < m - 1; ++k) {
        Eigen::MatrixXd probe(1, 2);
        probe << u1, g(k);
        EXPECT_NEAR(tll.hfunc1(probe)(0), cum(k) / cum(m - 1), 1e-8)
          << "rho = " << rho << ", u1 = " << u1 << ", knot " << k;
      }
    }
  }
}

namespace {

// The probability of a rectangle, from `(grid_points, values)` alone: the
// four-corner definition `C(u1, u2) = M(u1, u2) * u2 / M(1, u2)` for the
// interpolant's mass `M`, evaluated in `long double`.
//
// Differencing four corners of an order-one function costs an absolute `eps`
// whatever the summation accuracy, so a `double` evaluation of this cannot
// serve as truth for a route that avoids that difference -- it would carry the
// very error being measured, and correlate with the route it is meant to judge.
// The caller therefore requires a `long double` wider than `double`, which
// x86-64 has (64-bit mantissa) and arm64 and MSVC do not.
long double
rect_prob_reference(const Eigen::VectorXd& g,
                    const Eigen::MatrixXd& v,
                    double a1,
                    double b1,
                    double a2,
                    double b2)
{
  using ld = long double;
  const int m = static_cast<int>(g.size());
  auto mass = [&](ld x, ld y) {
    ld total = 0.0L;
    for (int i = 0; i + 1 < m && x > static_cast<ld>(g(i)); ++i) {
      const ld hx = static_cast<ld>(g(i + 1)) - static_cast<ld>(g(i));
      const ld sx = (std::min<ld>(x, g(i + 1)) - g(i)) / hx;
      for (int j = 0; j + 1 < m && y > static_cast<ld>(g(j)); ++j) {
        const ld hy = static_cast<ld>(g(j + 1)) - static_cast<ld>(g(j));
        const ld sy = (std::min<ld>(y, g(j + 1)) - g(j)) / hy;
        total += hx * hy *
                 ((sx - sx * sx / 2) * (sy - sy * sy / 2) * v(i, j) +
                  (sx * sx / 2) * (sy - sy * sy / 2) * v(i + 1, j) +
                  (sx - sx * sx / 2) * (sy * sy / 2) * v(i, j + 1) +
                  (sx * sx / 2) * (sy * sy / 2) * v(i + 1, j + 1));
      }
    }
    return total;
  };
  auto cdf = [&](ld x, ld y) {
    return (x <= 0.0L || y <= 0.0L) ? 0.0L : mass(x, y) * y / mass(1.0L, y);
  };
  return (cdf(b1, b2) + cdf(a1, a2)) - (cdf(a1, b2) + cdf(b1, a2));
}

} // namespace

// A mixed-discrete density is a rectangle probability divided by the atom's
// area, so it amplifies any absolute error in the corners by the reciprocal of
// that area. For a kernel pair the probability is a sum of nonnegative hat
// weights against a nonnegative grid, which costs one power of the atom width
// instead of two -- and the density is what the fit and every vine evaluation
// consume, so the difference is not academic.
TEST(test_bicop_kernel_accuracy, discrete_density_beats_a_cdf_difference)
{
  if (std::numeric_limits<long double>::digits <=
      std::numeric_limits<double>::digits) {
    GTEST_SKIP() << "needs a long double wider than double for the reference";
  }

  const Bicop tll = fit_tll(0.7, 1000);
  const Eigen::MatrixXd v = tll.get_parameters();
  const Eigen::VectorXd g = tll_grid_points(static_cast<int>(v.rows()));

  Bicop dd = tll;
  dd.set_var_types({ "d", "d" });

  // Dyadic atoms, so every bound is exactly representable and the reference
  // sees the same rectangle the library does. Widths start at 1/64: at 1/8 both
  // routes sit on the rounding floor of the corners themselves and the ratio
  // between them is noise, so there is nothing there to assert.
  for (int k : { 64, 512, 4096 }) {
    const double w = 1.0 / k;
    const int step = std::max(1, k / 12);
    double worst_exact = 0.0;
    double worst_corners = 0.0;
    for (int i = step; i < k; i += step) {
      for (int j = step; j < k; j += step) {
        const double b1 = i * w, a1 = b1 - w, b2 = j * w, a2 = b2 - w;
        // away from the boundary, where the atom stops being the small quantity
        if (a1 < 0.1 || b1 > 0.9 || a2 < 0.1 || b2 > 0.9) {
          continue;
        }
        // the density is the absolute value of the rectangle probability, so
        // that is what both routes are measured against
        const long double truth =
          std::abs(rect_prob_reference(g, v, a1, b1, a2, b2));
        ASSERT_GT(truth, 0.0L);

        Eigen::MatrixXd atom(1, 4);
        atom << b1, b2, a1, a2;
        const double exact = dd.pdf(atom)(0) * w * w;

        Eigen::MatrixXd corners(4, 2);
        corners << b1, b2, a1, b2, b1, a2, a1, a2;
        const Eigen::VectorXd c = tll.cdf(corners);
        const double differenced = std::abs((c(0) + c(3)) - (c(1) + c(2)));

        worst_exact = std::max(
          worst_exact, std::abs(static_cast<double>((exact - truth) / truth)));
        worst_corners = std::max(
          worst_corners,
          std::abs(static_cast<double>((differenced - truth) / truth)));
      }
    }
    // one power of the atom width rather than two. The margin is an order of
    // magnitude and more on every width tested, but the fitted grid differs
    // between platforms, so the assertion is deliberately not on the ratio's
    // size
    EXPECT_GT(worst_corners, 3.0 * worst_exact)
      << "atom width 1/" << k << ": exact " << worst_exact << ", four-corner "
      << worst_corners;
    EXPECT_LT(worst_exact, 1e-9) << "atom width 1/" << k << ": " << worst_exact;
  }
}

TEST(test_bicop_kernel_golden, golden_fits)
{
  auto bc =
    Bicop(BicopFamily::gaussian, 0, Eigen::MatrixXd::Constant(1, 1, 0.5));
  const auto data = bc.simulate(500, false, { 5 });
  Eigen::MatrixXd probe(6, 2);
  // clang-format off
  probe << 0.1, 0.1,  0.25, 0.75,  0.5, 0.5,
           0.75, 0.25,  0.9, 0.9,  0.05, 0.95;
  // clang-format on

  struct GoldenCase
  {
    std::string method;
    double npars;
    double loglik;
    Eigen::VectorXd pdf, hfunc1, hinv1, cdf, grid_probes;
  };
  std::vector<GoldenCase> cases(3);
  cases[0].method = "constant";
  cases[0].npars = 37.569790843197538;
  cases[0].loglik = 83.802106687208749;
  cases[0].pdf = Eigen::VectorXd(6);
  cases[0].pdf << 1.8940729582511884, 0.70935284843621338, 1.1148071649101758,
    0.77757379162372142, 2.2676315165456322, 0.17741750671483902;
  cases[0].hfunc1 = Eigen::VectorXd(6);
  cases[0].hfunc1 << 0.21899579248568388, 0.87479181886029278,
    0.52400104653037682, 0.13936087025008306, 0.80142252950786008,
    0.9957097179806309;
  cases[0].hinv1 = Eigen::VectorXd(6);
  cases[0].hinv1 << 0.042889211995842169, 0.6046520170087788,
    0.47846417922982398, 0.38336414037387623, 0.94387210622534923,
    0.79087117919954331;
  cases[0].cdf = Eigen::VectorXd(6);
  cases[0].cdf << 0.025367960010206622, 0.22763666123632589,
    0.32481401349537559, 0.23074282202648569, 0.83110985565914741,
    0.049805007976438359;
  cases[0].grid_probes = Eigen::VectorXd(5);
  cases[0].grid_probes << 15.827682893775718, 2.3984959626577322,
    1.1010338412868619, 2.8078410066920165, 5.790622493853256e-05;

  cases[1].method = "linear";
  cases[1].npars = 38.518228881716112;
  cases[1].loglik = 82.474046581040042;
  cases[1].pdf = Eigen::VectorXd(6);
  cases[1].pdf << 1.9188084368755076, 0.69510840161844223, 1.137204918323985,
    0.75619683069637178, 2.2307219461331593, 0.15415538783697091;
  cases[1].hfunc1 = Eigen::VectorXd(6);
  cases[1].hfunc1 << 0.22103482105128525, 0.88281900042762274,
    0.5218681353794673, 0.12728269288068772, 0.79688343566644859,
    0.99685079203148297;
  cases[1].hinv1 = Eigen::VectorXd(6);
  cases[1].hinv1 << 0.042836645988880998, 0.59177073332846897,
    0.48076405869388239, 0.39623272868149678, 0.94683523543492631,
    0.78562273644698866;
  cases[1].cdf = Eigen::VectorXd(6);
  cases[1].cdf << 0.026760106976800743, 0.22963293048386801,
    0.33036297224802852, 0.23276056961803065, 0.83294864128931778,
    0.04988387277881965;
  cases[1].grid_probes = Eigen::VectorXd(5);
  cases[1].grid_probes << 18.699710006479634, 2.4188195976400464,
    1.1236231416223819, 2.7977314845160701, 2.1478682342732693e-06;

  cases[2].method = "quadratic";
  cases[2].npars = 33.370124518148273;
  cases[2].loglik = 76.876425793478063;
  cases[2].pdf = Eigen::VectorXd(6);
  cases[2].pdf << 1.8251813429745112, 0.70329649116715842, 1.1450619754169749,
    0.74608219048986524, 2.235954375660707, 0.12960910101829565;
  cases[2].hfunc1 = Eigen::VectorXd(6);
  cases[2].hfunc1 << 0.2153720489450402, 0.88690062222533828,
    0.52421540595933014, 0.12748621382581807, 0.78167702122552829,
    0.99737633410425253;
  cases[2].hinv1 = Eigen::VectorXd(6);
  cases[2].hinv1 << 0.042851904268459945, 0.58854252057678169,
    0.47885517777009479, 0.39640746891648543, 0.95272760367705434,
    0.78121889618964591;
  cases[2].cdf = Eigen::VectorXd(6);
  cases[2].cdf << 0.027415045547700051, 0.23036402155514504,
    0.33225543730670176, 0.23248452233614378, 0.83480392402231929,
    0.049904458110113306;
  cases[2].grid_probes = Eigen::VectorXd(5);
  cases[2].grid_probes << 27.471550075175461, 2.4263410494972639,
    1.1361441410553463, 3.1363504740526791, 2.4694365156066189e-11;

  for (const auto& gc : cases) {
    FitControlsBicop controls({ BicopFamily::tll });
    controls.set_nonparametric_method(gc.method);
    Bicop tll(BicopFamily::tll);
    tll.fit(data, controls);
    EXPECT_NEAR(tll.get_npars(), gc.npars, 1e-8) << gc.method;
    EXPECT_NEAR(tll.get_loglik(), gc.loglik, 1e-8) << gc.method;
    EXPECT_TRUE(all_close(tll.pdf(probe), gc.pdf, 1e-8, 1e-10)) << gc.method;
    EXPECT_TRUE(all_close(tll.hfunc1(probe), gc.hfunc1, 1e-8, 1e-10))
      << gc.method;
    EXPECT_TRUE(all_close(tll.hinv1(probe), gc.hinv1, 1e-8, 1e-8)) << gc.method;
    EXPECT_TRUE(all_close(tll.cdf(probe), gc.cdf, 1e-8, 1e-10)) << gc.method;
    const Eigen::MatrixXd grid = tll.get_parameters();
    Eigen::VectorXd gp(5);
    gp << grid(0), grid(217), grid(435), grid(653), grid(880);
    EXPECT_TRUE(all_close(gp, gc.grid_probes, 1e-8, 1e-10)) << gc.method;
  }
}

namespace {

class TllLocalLikelihood : public TllBicop
{
public:
  using TllBicop::fit_local_likelihood;
};

//! The local likelihood fit at `x` by Newton's method, with the integral of
//! the kernel computed by the trapezoid rule in the coordinates whitened by
//! the bandwidth. Returns the estimate and the influence
//! `W_B(0) / n * (M^{-1})_{00}`.
std::pair<double, double>
brute_force_local_likelihood(const Eigen::Vector2d& x,
                             const Eigen::MatrixXd& data,
                             const Eigen::Matrix2d& B,
                             const std::string& method,
                             const Eigen::VectorXd& weights)
{
  const Eigen::Index p =
    (method == "constant") ? 1 : ((method == "linear") ? 3 : 6);
  const auto basis = [p](const Eigen::Vector2d& v) {
    Eigen::VectorXd a(p);
    a(0) = 1.0;
    if (p > 1) {
      a(1) = v(0);
      a(2) = v(1);
    }
    if (p > 3) {
      a(3) = v(0) * v(0) / 2.0;
      a(4) = v(0) * v(1);
      a(5) = v(1) * v(1) / 2.0;
    }
    return a;
  };
  const double two_pi = 2.0 * 3.14159265358979323846;
  const double kernel0 = 1.0 / (two_pi * std::sqrt(B.determinant()));
  const Eigen::Matrix2d Bi = B.inverse();
  const auto n = static_cast<double>(data.rows());

  // the data side of the score: sum_i w_i W_B(Z_i - x) A(Z_i - x)
  Eigen::VectorXd moments = Eigen::VectorXd::Zero(p);
  for (Eigen::Index i = 0; i < data.rows(); ++i) {
    const Eigen::Vector2d v = data.row(i).transpose() - x;
    const double w = (weights.size() > 0) ? weights(i) : 1.0;
    moments += w * kernel0 * std::exp(-0.5 * v.dot(Bi * v)) * basis(v);
  }

  // int W_B(v) g(v) dv = E[g(L z)] with z ~ N(0, I) and B = L L'
  const Eigen::Matrix2d L = B.llt().matrixL();
  const double h = 0.1;
  const int half = 100;
  std::vector<Eigen::VectorXd> nodes;
  std::vector<double> node_weights;
  for (int i1 = -half; i1 <= half; ++i1) {
    for (int i2 = -half; i2 <= half; ++i2) {
      const Eigen::Vector2d z(h * i1, h * i2);
      nodes.push_back(basis(L * z));
      node_weights.push_back(h * h * std::exp(-0.5 * z.squaredNorm()) / two_pi);
    }
  }

  Eigen::VectorXd a = Eigen::VectorXd::Zero(p);
  a(0) = std::log(moments(0) / n);
  Eigen::MatrixXd M(p, p);
  for (int it = 0; it < 100; ++it) {
    Eigen::VectorXd score = moments;
    M.setZero();
    for (size_t q = 0; q < nodes.size(); ++q) {
      const double e = node_weights[q] * std::exp(a.dot(nodes[q]));
      score -= n * e * nodes[q];
      M += e * nodes[q] * nodes[q].transpose();
    }
    const Eigen::VectorXd step = (n * M).ldlt().solve(score);
    a += step;
    if (step.lpNorm<Eigen::Infinity>() < 1e-14) {
      break;
    }
  }
  return { std::exp(a(0)), kernel0 / n * M.inverse()(0, 0) };
}

} // namespace

TEST(test_bicop_kernel, local_likelihood_maximizes_the_local_likelihood)
{
  // a correlated bandwidth, on the normal scale
  auto bc =
    Bicop(BicopFamily::gaussian, 0, Eigen::MatrixXd::Constant(1, 1, 0.6));
  const Eigen::MatrixXd data =
    tools_stats::qnorm(bc.simulate(200, false, { 3 }));
  Eigen::Matrix2d B;
  B << 0.3, 0.18, 0.18, 0.3;
  Eigen::MatrixXd x(4, 2);
  x << 0.4, -0.2, 0.0, 0.0, 1.5, 1.2, -1.0, 1.0;
  const Eigen::VectorXd weights =
    1.0 + tools_stats::simulate_uniform(200, 1, false, { 4 }).array();

  TllLocalLikelihood tll;
  for (const std::string method : { "constant", "linear", "quadratic" }) {
    const Eigen::MatrixXd fit =
      tll.fit_local_likelihood(x, data, B, method, Eigen::VectorXd());
    const Eigen::MatrixXd weighted =
      tll.fit_local_likelihood(x, data, B, method, weights);
    for (Eigen::Index k = 0; k < x.rows(); ++k) {
      const auto [estimate, influence] = brute_force_local_likelihood(
        x.row(k).transpose(), data, B, method, Eigen::VectorXd());
      EXPECT_NEAR(fit(k, 0) / estimate, 1.0, 1e-9) << method << " " << k;
      EXPECT_NEAR(fit(k, 1) / influence, 1.0, 1e-9) << method << " " << k;
      const double weighted_estimate =
        brute_force_local_likelihood(
          x.row(k).transpose(), data, B, method, weights)
          .first;
      EXPECT_NEAR(weighted(k, 0) / weighted_estimate, 1.0, 1e-9)
        << method << " " << k;
    }
  }
}

TEST(test_bicop_kernel_golden, hfunc_hinv_identity)
{
  auto bc =
    Bicop(BicopFamily::gaussian, 0, Eigen::MatrixXd::Constant(1, 1, 0.5));
  const auto data = bc.simulate(500, false, { 5 });

  for (const auto& method : { std::string("constant"),
                              std::string("linear"),
                              std::string("quadratic") }) {
    FitControlsBicop controls({ BicopFamily::tll });
    controls.set_nonparametric_method(method);
    Bicop tll(BicopFamily::tll);
    tll.fit(data, controls);

    // round trip on a deterministic interior grid
    Eigen::MatrixXd u(9 * 9, 2);
    Eigen::VectorXd g = Eigen::VectorXd::LinSpaced(9, 0.1, 0.9);
    size_t k = 0;
    for (long i = 0; i < 9; ++i) {
      for (long j = 0; j < 9; ++j) {
        u(k, 0) = g(i);
        u(k, 1) = g(j);
        ++k;
      }
    }
    Eigen::MatrixXd u_inv = u;
    u_inv.col(1) = tll.hinv1(u);
    Eigen::VectorXd q = tll.hfunc1(u_inv);
    EXPECT_TRUE(all_close(q, u.col(1), 0.0, 1e-8)) << method;

    // monotonicity of hinv1 in the target quantile
    Eigen::MatrixXd u_mono(50, 2);
    u_mono.col(0) = Eigen::VectorXd::Constant(50, 0.4);
    u_mono.col(1) = Eigen::VectorXd::LinSpaced(50, 0.01, 0.99);
    Eigen::VectorXd v = tll.hinv1(u_mono);
    for (long i = 1; i < v.size(); ++i) {
      EXPECT_GE(v(i), v(i - 1) - 1e-12) << method;
    }
  }
}
}
