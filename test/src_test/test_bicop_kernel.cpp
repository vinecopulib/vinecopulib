// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
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

// Golden-value regression guards for the performance work: fixed-seed fits
// and evaluations pinned to the pre-optimization implementation. Class-B
// optimizations (floating-point reassociation) stay well within the 1e-8
// tolerances; anything larger signals an unintended behavior change.
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
  cases[0].loglik = 88.537611609339933;
  cases[0].pdf = Eigen::VectorXd(6);
  cases[0].pdf << 1.9920693244173888, 0.68223264104149317, 1.1300894677022009,
    0.67953894472340226, 1.992993165308208, 0.18242727854749052;
  cases[0].hfunc1 = Eigen::VectorXd(6);
  cases[0].hfunc1 << 0.23658607047125324, 0.89666511364276091,
    0.49422315377271525, 0.11673440978420759, 0.79122136895234918,
    0.99491543977136443;
  cases[0].hinv1 = Eigen::VectorXd(6);
  cases[0].hinv1 << 0.039283122081542388, 0.57153584840125404,
    0.50510621539433487, 0.42133118640049361, 0.95235049651819281,
    0.72590232422226109;
  cases[0].cdf = Eigen::VectorXd(6);
  cases[0].cdf << 0.026003049149675449, 0.23438161606545282, 0.3388084584458943,
    0.23284538470946511, 0.83053173864168195, 0.049793141885435759;
  cases[0].grid_probes = Eigen::VectorXd(5);
  cases[0].grid_probes << 0.91548412857751971, 2.8006107481433404,
    1.0792311697691204, 2.6820004917790174, 0.001657074800914232;

  cases[1].method = "linear";
  cases[1].npars = 37.734218596805036;
  cases[1].loglik = 87.892204122740054;
  cases[1].pdf = Eigen::VectorXd(6);
  cases[1].pdf << 2.1144069092419406, 0.62538328502263485, 1.1605849912836177,
    0.66242802925973743, 2.1314471201078322, 0.16357640564084397;
  cases[1].hfunc1 = Eigen::VectorXd(6);
  cases[1].hfunc1 << 0.24429353041499491, 0.91080867527660725,
    0.49661789283091529, 0.10324553681290974, 0.7788659362205681,
    0.99639765654827461;
  cases[1].hinv1 = Eigen::VectorXd(6);
  cases[1].hinv1 << 0.038999194483039901, 0.55220823522540741,
    0.50291259869118221, 0.43629357797908597, 0.9543718233180698,
    0.71185221741325222;
  cases[1].cdf = Eigen::VectorXd(6);
  cases[1].cdf << 0.027253298336657676, 0.23689322100612506, 0.3466670852339287,
    0.23657671687757489, 0.83369279538193652, 0.04987082389178156;
  cases[1].grid_probes = Eigen::VectorXd(5);
  cases[1].grid_probes << 0.68879142973403928, 2.9385513732236999,
    1.1122992213940723, 3.0246339629663934, 3.6640428424338637e-05;

  cases[2].method = "quadratic";
  cases[2].npars = 19.315912883957651;
  cases[2].loglik = 84.141138454688559;
  cases[2].pdf = Eigen::VectorXd(6);
  cases[2].pdf << 2.0419253004032352, 0.61317262541702111, 1.1570504390386878,
    0.65818314435011982, 2.1604228094451736, 0.11884962833032209;
  cases[2].hfunc1 = Eigen::VectorXd(6);
  cases[2].hfunc1 << 0.23459154280716302, 0.90820398859517293,
    0.49793765968004844, 0.10430422980907414, 0.77277085813672075,
    0.99677548140291572;
  cases[2].hinv1 = Eigen::VectorXd(6);
  cases[2].hinv1 << 0.04008680660626851, 0.55237759070587344,
    0.50178197349305265, 0.43574022615212016, 0.95599419446079992,
    0.72277853023842908;
  cases[2].cdf = Eigen::VectorXd(6);
  cases[2].cdf << 0.025261223375066241, 0.23602522243110072,
    0.34518361647512608, 0.2355258899124702, 0.83296077233960586,
    0.049910964418499251;
  cases[2].grid_probes = Eigen::VectorXd(5);
  cases[2].grid_probes << 0.048059981571689885, 2.6471395430672491,
    1.1180743763936611, 2.9461094294257286, 2.6815474736691195e-07;

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
