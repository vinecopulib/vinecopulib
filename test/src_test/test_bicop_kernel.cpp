// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "include/kernel_test.hpp"
#include "include/r_parity.hpp"
#include "include/test_utils.hpp"

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
  cases[0].npars = 37.000327042690529;
  cases[0].loglik = 88.538467766120675;
  cases[0].pdf = Eigen::VectorXd(6);
  cases[0].pdf << 1.9919574139287595, 0.68216086502430573, 1.1300814135201216,
    0.67960363298559812, 1.9930902734906932, 0.18239227425062596;
  cases[0].hfunc1 = Eigen::VectorXd(6);
  cases[0].hfunc1 << 0.23659470476088243, 0.89667232461375157,
    0.49423919599774013, 0.11674185749009439, 0.79123264641260593,
    0.99491617243809083;
  cases[0].hinv1 = Eigen::VectorXd(6);
  cases[0].hinv1 << 0.039281387072030501, 0.57152308973806543,
    0.50509203639013178, 0.4213165287079067, 0.95234734906032714,
    0.72588964669310063;
  cases[0].cdf = Eigen::VectorXd(6);
  cases[0].cdf << 0.026001734522761211, 0.23436511876766064,
    0.33879498981043821, 0.2328428085812691, 0.83052099047340289,
    0.04978895320650701;
  cases[0].grid_probes = Eigen::VectorXd(5);
  cases[0].grid_probes << 0.91541618290656168, 2.8004708797060593,
    1.079243827522306, 2.6823808892119843, 0.001656625154939985;

  cases[1].method = "linear";
  cases[1].npars = 37.734218596805036;
  cases[1].loglik = 87.89070543544986;
  cases[1].pdf = Eigen::VectorXd(6);
  cases[1].pdf << 2.1147484360269111, 0.62557840649005281, 1.1605887631415703,
    0.66223063946627514, 2.1311487707942747, 0.16367327862097206;
  cases[1].hfunc1 = Eigen::VectorXd(6);
  cases[1].hfunc1 << 0.24426649002304113, 0.91078930227766286,
    0.49656848950533444, 0.10322468716386059, 0.77883392509654781,
    0.99639672660723977;
  cases[1].hinv1 = Eigen::VectorXd(6);
  cases[1].hinv1 << 0.03900409965324849, 0.55224589559366488,
    0.50295513860106622, 0.43633650157871079, 0.9543796762149559,
    0.71189513687558392;
  cases[1].cdf = Eigen::VectorXd(6);
  cases[1].cdf << 0.027257703816971766, 0.23694283302649019,
    0.34670484688762238, 0.23658250630015026, 0.8337177927126862,
    0.049884234805975383;
  cases[1].grid_probes = Eigen::VectorXd(5);
  cases[1].grid_probes << 0.68899768218664004, 2.9389943791369473,
    1.1122401226577394, 3.0236851112882945, 3.6668397922157365e-05;

  cases[2].method = "quadratic";
  cases[2].npars = 19.315912883957658;
  cases[2].loglik = 84.142992201245306;
  cases[2].pdf = Eigen::VectorXd(6);
  cases[2].pdf << 2.0414794193343084, 0.61291098953624357, 1.1570502697684848,
    0.65845304894037671, 2.1608602349656523, 0.11875599689850495;
  cases[2].hfunc1 = Eigen::VectorXd(6);
  cases[2].hfunc1 << 0.23462644843576383, 0.9082303915564236,
    0.4980031311021596, 0.10433205987706211, 0.77281318369365759,
    0.99677732517569217;
  cases[2].hinv1 = Eigen::VectorXd(6);
  cases[2].hinv1 << 0.040079953487945796, 0.55232703807108607,
    0.50172540287002709, 0.43568361010984546, 0.95598455230608914,
    0.72272856136941177;
  cases[2].cdf = Eigen::VectorXd(6);
  cases[2].cdf << 0.025255464924933916, 0.23595672916972246,
    0.34513102474705393, 0.23551736010615565, 0.83292765435424554,
    0.049892244841626389;
  cases[2].grid_probes = Eigen::VectorXd(5);
  cases[2].grid_probes << 0.048036494984137437, 2.6465744594545155,
    1.1181605901828306, 2.9472169281733644, 2.6790149427227697e-07;

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
