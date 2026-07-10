// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "test_utils.hpp"
#include "vinecop_test.hpp"
#include <string>
#include <vinecopulib.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace test_vinecop_class {
using namespace vinecopulib;
using test_utils::all_close;

TEST_F(VinecopTest, constructors_without_error)
{
  Vinecop vinecop(5);
  Vinecop vinecop_indep(model_matrix);

  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 90);
    }
  }

  Vinecop vinecop_parametrized(model_matrix, pair_copulas);
}

TEST_F(VinecopTest, copy)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 90);
    }
  }

  Vinecop vinecop1(model_matrix, pair_copulas);
  auto pc = vinecop1.get_pair_copula(0, 0);
  pc.set_parameters(pc.get_parameters().array() + 1);
  EXPECT_EQ(vinecop1.get_parameters(0, 0), vinecop1.get_parameters(0, 1));

  Vinecop vinecop2 = vinecop1;
  pair_copulas[0][0].set_parameters(pc.get_parameters().array() + 1);
  vinecop2.set_all_pair_copulas(pair_copulas);
  EXPECT_EQ(vinecop1.get_parameters(0, 0), vinecop1.get_parameters(0, 1));
}

TEST_F(VinecopTest, print)
{
  auto cvine = CVineStructure(std::vector<size_t>({ 5, 4, 3, 2, 1 }));
  auto vc1 = Vinecop(cvine);

  // check if first, second and last line are correct
  std::string expected_first_line = "Vinecop model with 5 variables";
  std::string expected_second_line =
    "tree edge conditioned variables conditioning variables var_types       "
    "family rotation parameters df tau ";
  std::string expected_last_line =
    "   4    1                  5, 4                3, 2, 1      c, c "
    "Independence                        0.0 ";

  std::istringstream input;
  input.str(vc1.str());

  std::string line;
  // get first, second and last line
  std::getline(input, line);
  EXPECT_EQ(line, expected_first_line);
  std::getline(input, line);
  EXPECT_EQ(line, expected_second_line);
  std::string last_line;
  while (std::getline(input, line)) {
    last_line = line;
  }
  EXPECT_EQ(last_line, expected_last_line);

  // create vine with 7 variables
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::tawn, 270);
    }
  }

  Vinecop vc2(model_matrix, pair_copulas);

  // check if first, second and last line are correct
  expected_first_line = "Vinecop model with 7 variables";
  expected_second_line =
    "tree edge conditioned variables conditioning variables var_types family "
    "rotation       parameters  df   tau ";
  expected_last_line = "   6    1                  4, 7          3, 1, 2, 6, 5 "
                       "     c, c   Tawn      270 1.00, 1.00, 1.00 3.0 -0.00 ";

  input.clear();
  input.str(vc2.str());

  // get first, second and last line
  std::getline(input, line);
  EXPECT_EQ(line, expected_first_line);
  std::getline(input, line);
  EXPECT_EQ(line, expected_second_line);
  while (std::getline(input, line)) {
    last_line = line;
  }
  EXPECT_EQ(last_line, expected_last_line);

  auto data = tools_stats::simulate_uniform(100, 5);
  auto controls = FitControlsVinecop({ BicopFamily::tll });
  vc1.select(data, controls);

  // check if first and second are correct
  // we don't check the last line as it is random
  // but at least we see it doesn't crash
  expected_first_line = "Vinecop model with 5 variables";
  expected_second_line =
    "tree edge conditioned variables conditioning variables var_types family "
    "rotation   parameters   df   tau ";

  input.clear();
  input.str(vc1.str());

  // get first, second and last line
  std::getline(input, line);
  EXPECT_EQ(line, expected_first_line) << vc1.str();
  std::getline(input, line);
  EXPECT_EQ(line, expected_second_line) << vc1.str();
}

TEST_F(VinecopTest, serialization)
{

  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(7, 7);
  mat << 5, 2, 6, 6, 6, 6, 6, 6, 6, 1, 2, 5, 5, 0, 2, 5, 2, 5, 2, 0, 0, 1, 1, 5,
    1, 0, 0, 0, 3, 7, 7, 0, 0, 0, 0, 7, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0;

  // create vine with 7 variables, 2-truncated
  size_t d = 7;
  auto pc_store = Vinecop::make_pair_copula_store(d, 5);
  for (auto& tree : pc_store) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::bb1, 90);
    }
  }

  auto vc = Vinecop(mat, pc_store);
  vc.truncate(3);

  // serialize
  vc.to_file(std::string("temp"));

  // unserialize
  auto vc2 = Vinecop(std::string("temp"));

  // Remove temp file
  std::string cmd = rm + "temp";
  int sys_exit_code = system(cmd.c_str());
  if (sys_exit_code != 0) {
    throw std::runtime_error("error in system call");
  }

  EXPECT_EQ(vc.get_all_rotations(), vc2.get_all_rotations());
  EXPECT_EQ(vc.get_all_families(), vc2.get_all_families());
  EXPECT_EQ(vc.get_var_types(), vc2.get_var_types());
  EXPECT_EQ(vc.get_matrix(), vc2.get_matrix());
}

TEST_F(VinecopTest, getters_are_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 90);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  for (auto& tree : vinecop.get_all_families()) {
    for (auto& fam : tree) {
      EXPECT_EQ(fam, BicopFamily::clayton);
    }
  }

  for (auto& tree : vinecop.get_all_pair_copulas()) {
    for (auto& pc : tree) {
      EXPECT_EQ(pc.get_family(), BicopFamily::clayton);
      EXPECT_EQ(pc.get_rotation(), 90);
    }
  }

  for (auto& tree : vinecop.get_all_parameters()) {
    for (auto& par : tree) {
      EXPECT_EQ(par.size(), 1);
      EXPECT_EQ(par(0), 1e-10);
    }
  }

  for (auto& tree : vinecop.get_all_rotations()) {
    for (auto& rot : tree) {
      EXPECT_EQ(rot, 90);
    }
  }

  for (auto& tree : vinecop.get_all_taus()) {
    for (auto& tau : tree) {
      ASSERT_TRUE(fabs(tau) < 1e-4);
    }
  }

  EXPECT_NO_THROW(vinecop.get_dim());
  EXPECT_NO_THROW(vinecop.get_rvine_structure());
  EXPECT_ANY_THROW(vinecop.get_loglik());
  EXPECT_ANY_THROW(vinecop.get_nobs());
  EXPECT_ANY_THROW(vinecop.get_aic());
  EXPECT_ANY_THROW(vinecop.get_bic());
  EXPECT_ANY_THROW(vinecop.get_mbicv());
  EXPECT_ANY_THROW(vinecop.aic());
  EXPECT_ANY_THROW(vinecop.bic());
  EXPECT_ANY_THROW(vinecop.mbicv());
}

TEST_F(VinecopTest, 1dim)
{
  auto data = tools_stats::simulate_uniform(3, 1);
  auto vc = Vinecop(1);
  vc.select(data);
  EXPECT_TRUE((vc.pdf(data).array() == 1).all());
  EXPECT_TRUE((vc.rosenblatt(data).array() == data.array()).all());
  EXPECT_TRUE((vc.inverse_rosenblatt(data).array() == data.array()).all());
  vc.loglik();
  vc.aic();
  vc.simulate(3);
  Vinecop::make_pair_copula_store(1, 2);
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(1, 1);
  mat(0, 0) = 1;
  RVineStructure rvine_structure(mat);
}

TEST_F(VinecopTest, fit_statistics_getters_are_correct)
{
  auto data = tools_stats::simulate_uniform(100, 3);
  auto vc = Vinecop(
    data, RVineStructure(), {}, FitControlsVinecop({ BicopFamily::clayton }));
  EXPECT_NEAR(vc.get_loglik(), vc.loglik(data), 1e-10);
  EXPECT_NEAR(static_cast<double>(vc.get_nobs()), 100, 1e-10);
  EXPECT_NEAR(vc.get_aic(), vc.aic(data), 1e-10);
  EXPECT_NEAR(vc.get_bic(), vc.bic(data), 1e-10);
  EXPECT_NEAR(vc.get_mbicv(0.6), vc.mbicv(data, 0.6), 1e-10);
}

TEST_F(VinecopTest, truncate_methods_works)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);
  vinecop.truncate(2);
  EXPECT_EQ(vinecop.get_all_pair_copulas().size(), 2);
  EXPECT_EQ(vinecop.get_rvine_structure().get_trunc_lvl(), 2);
  vinecop.truncate(0);
  EXPECT_EQ(vinecop.get_all_pair_copulas().size(), 0);
  EXPECT_EQ(vinecop.get_rvine_structure().get_trunc_lvl(), 0);
}

TEST_F(VinecopTest, pdf_is_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  ASSERT_TRUE(all_close(vinecop.pdf(u), f, 1e-4, 1e-4));
}

TEST_F(VinecopTest, hfuncs_is_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  auto r = vinecop.pdf_full(u, 1, true);

  auto rvine_structure = vinecop.get_rvine_structure();
  size_t d = vinecop.get_dim();
  size_t trunc_lvl = rvine_structure.get_trunc_lvl();

  ASSERT_EQ(r.hfunc1.get_trunc_lvl(), trunc_lvl);
  ASSERT_EQ(r.hfunc2.get_trunc_lvl(), trunc_lvl);
  ASSERT_EQ(r.hfunc1_sub.get_trunc_lvl(), trunc_lvl);
  ASSERT_EQ(r.hfunc2_sub.get_trunc_lvl(), trunc_lvl);
  ASSERT_EQ(r.pdf_edges.get_trunc_lvl(), trunc_lvl);

  for (size_t tree = 0; tree < trunc_lvl; ++tree) {
    for (size_t edge = 0; edge < d - tree - 1; ++edge) {
      EXPECT_EQ(r.pdf_edges(tree, edge).size(), u.rows());
      if (rvine_structure.needed_hfunc1(tree, edge)) {
        EXPECT_EQ(r.hfunc1(tree, edge).size(), u.rows());
        EXPECT_EQ(r.hfunc1_sub(tree, edge).size(), u.rows());
      } else {
        EXPECT_EQ(r.hfunc1(tree, edge).size(), 0);
        EXPECT_EQ(r.hfunc1_sub(tree, edge).size(), 0);
      }
      if (rvine_structure.needed_hfunc2(tree, edge)) {
        EXPECT_EQ(r.hfunc2(tree, edge).size(), u.rows());
        EXPECT_EQ(r.hfunc2_sub(tree, edge).size(), u.rows());
      } else {
        EXPECT_EQ(r.hfunc2(tree, edge).size(), 0);
        EXPECT_EQ(r.hfunc2_sub(tree, edge).size(), 0);
      }
    }
  }
}

TEST_F(VinecopTest, cdf_is_correct)
{
  // Create a bivariate copula and a corresponding vine with two variables
  auto pair_copulas = Vinecop::make_pair_copula_store(2);
  auto par = Eigen::VectorXd::Constant(1, 0.5);
  auto bicop = Bicop(BicopFamily::gaussian, 0, par);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::gaussian, 0, par);
    }
  }
  Eigen::Matrix<size_t, 2, 2> matrix;
  matrix << 1, 1, 2, 0;
  Vinecop vinecop(matrix, pair_copulas);

  // Test whether the analytic and simulated versions are "close" enough
  auto u2 = vinecop.simulate(10);
  ASSERT_TRUE(all_close(vinecop.cdf(u2, 10000), bicop.cdf(u2), 1e-2, 1e-2));

  // verify that qrng stuff works
  Vinecop vinecop2(301);
  vinecop.simulate(10, true);
}

TEST_F(VinecopTest, simulate_is_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  // only check if it works
  vinecop.simulate(10);
  // check the underlying transformation from independent samples
  ASSERT_TRUE(all_close(vinecop.inverse_rosenblatt(u), sim, 1e-4, 1e-4));

  // verify that qrng stuff works
  vinecop.simulate(10, true);
  Vinecop vinecop2(301);
  vinecop.simulate(10, true);
}

TEST_F(VinecopTest, rosenblatt_is_correct)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);
  auto u2 = vinecop.simulate(5);
  ASSERT_TRUE(all_close(
    vinecop.rosenblatt(vinecop.inverse_rosenblatt(u2)), u2, 1e-6, 1e-6));

  // truncated multivariate
  pair_copulas = Vinecop::make_pair_copula_store(7, 2);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  vinecop = Vinecop(model_matrix, pair_copulas);
  ASSERT_TRUE(all_close(
    vinecop.rosenblatt(vinecop.inverse_rosenblatt(u)), u, 1e-6, 1e-6));

  // bivariate case
  pair_copulas = Vinecop::make_pair_copula_store(2);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Eigen::Matrix<size_t, 2, 2> mat;
  mat << 1, 1, 2, 0;
  vinecop = Vinecop(mat, pair_copulas);
  u = vinecop.simulate(5);
  ASSERT_TRUE(all_close(
    vinecop.rosenblatt(vinecop.inverse_rosenblatt(u)), u, 1e-6, 1e-6));
}

TEST_F(VinecopTest, scores_stepwise)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  auto uu = vinecop.simulate(100, false, 1, { 1 });

  auto J = vinecop.hessian(uu, true);
  EXPECT_TRUE(J.isUpperTriangular());

  // hessian() is the observation-average of hessian_full()
  auto H = vinecop.hessian_full(uu, true);
  Eigen::MatrixXd J_manual(J.rows(), J.cols());
  size_t d = vinecop.get_dim();
  size_t trunc_lvl = vinecop.get_trunc_lvl();
  size_t ipar = 0;
  for (size_t t = 0; t < trunc_lvl; t++) {
    for (size_t e = 0; e < d - 1 - t; e++) {
      size_t np = static_cast<size_t>(
        vinecop.get_pair_copula(t, e).get_parameters().size());
      for (size_t p = 0; p < np; p++) {
        J_manual.row(ipar++) = H(t, e)[p].colwise().mean();
      }
    }
  }
  EXPECT_TRUE(all_close(J, J_manual, 1e-10, 1e-10));

  // scores_cov() is the mean-centered covariance of scores()
  auto s = vinecop.scores(uu, true);
  Eigen::MatrixXd sc = s.rowwise() - s.colwise().mean();
  Eigen::MatrixXd I_manual =
    (sc.adjoint() * sc) / static_cast<double>(s.rows());
  EXPECT_TRUE(all_close(vinecop.scores_cov(uu, true), I_manual, 1e-10, 1e-10));
}

TEST_F(VinecopTest, scores_joint)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  auto uu = vinecop.simulate(100, false, 1, { 1 });
  auto J = vinecop.hessian(uu, false);
  auto I = vinecop.scores_cov(uu, false);

  EXPECT_FALSE(J.isUpperTriangular());

  Eigen::MatrixXd Jinv = J.triangularView<Eigen::Upper>()
                           .solve(Eigen::MatrixXd::Identity(J.cols(), J.cols()))
                           .triangularView<Eigen::Upper>();
  // std::cout << J << std::endl << std::endl;
  // std::cout << Jinv * I << std::endl;
  I = vinecop.scores_cov(uu, true);
  // std::cout << (Jinv * I * Jinv.transpose() /
  // u.rows()).diagonal().cwiseSqrt()
  //           << std::endl
  //           << std::endl;
  // std::cout << vinecop.str() << std::endl;
}

// The analytic full gradient (step_wise = false) must match VineCopula's
// RVineGrad and the analytic joint Hessian must match RVineHessian, both
// element-wise (analytic vs analytic). Requires VineCopula >= 2.6.2, whose
// #101 fix makes RVineHessian depend on the data so it can be used directly
// (not finite-differenced). The R oracle (test_vinecop_derivatives.R)
// hardcodes the same model and reorders RVineGrad and RVineHessian into
// vinecopulib's parameter order.
//
// The model uses only 0/180 families: vinecopulib and VineCopula do not share
// the same structure / swapped-rotation conventions, so a 90/270 vine cannot
// be built identically in both (see the R oracle's header). 0/180 families
// are symmetric under the argument swap that distinguishes the conventions.
// Vine-level rotations are validated against brute force instead
// (hessian_matches_brute_force, below), and rotated bicop derivatives against
// BiCopDeriv* (ParBicopTest).
TEST(VinecopDerivatives, full_scores_match_RVineGrad_and_RVineHessian)
{
  std::string cmd =
    std::string(RSCRIPT) + std::string(TEST_VINECOP_DERIV) + " 300";
  int sys_exit_code = system(cmd.c_str());
  if (sys_exit_code != 0) {
    throw std::runtime_error("error in system call");
  }
  auto mat =
    tools_eigen::read_matxs("temp_vderiv_matrix").colwise().reverse().eval();
  auto u = tools_eigen::read_matxd("temp_vderiv_data");
  auto grads = tools_eigen::read_matxd("temp_vderiv_grad");
  auto hess_r = tools_eigen::read_matxd("temp_vderiv_hess");
  cmd = rm + "temp_vderiv_data temp_vderiv_matrix temp_vderiv_grad "
             "temp_vderiv_hess";
  sys_exit_code += system(cmd.c_str());
  if (sys_exit_code != 0) {
    throw std::runtime_error("error in system call");
  }

  // must stay in sync with the R script (0/180 families only; pair copula
  // (t, e) <-> VineCopula matrix cell [d - t, e + 1])
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 0, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 0, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  pcs[3][0] = Bicop(BicopFamily::clayton, 0, P1(2.5));
  Vinecop vc(RVineStructure(mat), pcs);
  ASSERT_EQ(static_cast<size_t>(vc.get_npars()), grads.rows());

  // full gradient vs RVineGrad: sum over observations, and a single row
  Eigen::MatrixXd s = vc.scores(u, false, 1);
  EXPECT_TRUE(all_close(
    s.colwise().sum().transpose().eval(), grads.col(0).eval(), 1e-4, 1e-4));
  EXPECT_TRUE(
    all_close(s.row(0).transpose().eval(), grads.col(1).eval(), 1e-4, 1e-4));

  // gradient() is the observation-average of the scores
  EXPECT_TRUE(all_close(vc.gradient(u, false, 1),
                        s.colwise().mean().transpose().eval(),
                        1e-12,
                        1e-12));

  // scores_full(keep_all) returns the same scores plus non-empty caches
  auto full = vc.scores_full(u, false, 1, true);
  EXPECT_TRUE(all_close(full.scores, s, 1e-12, 1e-12));
  EXPECT_EQ(full.pdf_edges.get_dim(), vc.get_dim());
  EXPECT_EQ(full.logpdf_deriv_pars(0, 0).size(), 1u);
  EXPECT_EQ(full.logpdf_deriv_pars(0, 0)[0].size(), u.rows());

  // deterministic across the number of threads
  EXPECT_TRUE(all_close(vc.scores(u, false, 3), s, 1e-12, 1e-12));

  // joint Hessian vs RVineHessian, element-wise. n * hessian() is the summed
  // observed information, as RVineHessian returns.
  Eigen::MatrixXd H = vc.hessian(u, false, 1) * static_cast<double>(u.rows());
  EXPECT_TRUE(all_close(H, hess_r, 1e-4, 1e-4));
  EXPECT_TRUE(all_close(H, H.transpose().eval(), 1e-10, 1e-10)); // symmetric
}

// The analytic joint Hessian must match brute-force central finite
// differences of the whole-vine log-likelihood — the unforgiving ground
// truth, and (unlike the R oracle) valid for 90/270 rotations and
// truncation, which VineCopula's RVineHessian cannot provide.
TEST(VinecopDerivatives, hessian_matches_brute_force)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 90, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 270, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(-0.3));
  pcs[3][0] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  Vinecop vc(DVineStructure({ 1, 2, 3, 4, 5 }), pcs);
  auto u = vc.simulate(150, false, 1, { 42, 1, 2 });

  size_t npars = static_cast<size_t>(vc.get_npars());
  Eigen::MatrixXd H = vc.hessian(u, false, 1) * static_cast<double>(u.rows());

  // flat (t, e, p) list in hessian()'s column order
  std::vector<std::array<size_t, 3>> pl;
  for (size_t t = 0; t < 4; ++t) {
    for (size_t e = 0; e < 4 - t; ++e) {
      size_t np = static_cast<size_t>(pcs[t][e].get_parameters().size());
      for (size_t p = 0; p < np; ++p) {
        pl.push_back({ t, e, p });
      }
    }
  }
  auto structure = vc.get_rvine_structure();
  auto loglik = [&](const std::vector<std::vector<Bicop>>& pp) {
    return Vinecop(structure, pp).pdf(u).array().max(1e-300).log().sum();
  };
  auto bump = [&](std::vector<std::vector<Bicop>> pp, size_t a, double h) {
    auto pr = pp[pl[a][0]][pl[a][1]].get_parameters();
    pr(pl[a][2]) += h;
    pp[pl[a][0]][pl[a][1]].set_parameters(pr);
    return pp;
  };
  double base = loglik(pcs);
  double h = 1e-4;
  for (size_t a = 0; a < npars; ++a) {
    for (size_t bb = a; bb < npars; ++bb) {
      double fd;
      if (a == bb) {
        fd = (loglik(bump(pcs, a, h)) - 2 * base + loglik(bump(pcs, a, -h))) /
             (h * h);
      } else {
        fd = (loglik(bump(bump(pcs, a, h), bb, h)) -
              loglik(bump(bump(pcs, a, h), bb, -h)) -
              loglik(bump(bump(pcs, a, -h), bb, h)) +
              loglik(bump(bump(pcs, a, -h), bb, -h))) /
             (4 * h * h);
      }
      EXPECT_NEAR(H(a, bb), fd, 1e-3 * (1.0 + std::abs(fd)))
        << "H(" << a << ", " << bb << ")";
    }
  }
  EXPECT_TRUE(all_close(H, H.transpose().eval(), 1e-10, 1e-10));
  EXPECT_TRUE(all_close(vc.hessian(u, false, 3) * static_cast<double>(u.rows()),
                        H,
                        1e-12,
                        1e-12)); // threading determinism
}

// The analytic step-wise Hessian must match central finite differences of
// the step-wise scores (the quantity the finite-difference path produced),
// including a two-parameter (Student) edge's own-edge cross block.
TEST(VinecopDerivatives, stepwise_hessian_matches_fd_of_scores)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 90, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 270, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(-0.3));
  pcs[3][0] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  auto structure = DVineStructure({ 1, 2, 3, 4, 5 });
  Vinecop vc(structure, pcs);
  auto u = vc.simulate(120, false, 1, { 5, 6, 7 });

  auto Ha = vc.hessian_full(u, true, 1);

  std::vector<std::array<size_t, 3>> pl;
  for (size_t t = 0; t < 4; ++t) {
    for (size_t e = 0; e < 4 - t; ++e) {
      size_t np = static_cast<size_t>(pcs[t][e].get_parameters().size());
      for (size_t p = 0; p < np; ++p) {
        pl.push_back({ t, e, p });
      }
    }
  }
  auto bump = [&](std::vector<std::vector<Bicop>> pp, size_t a, double h) {
    auto pr = pp[pl[a][0]][pl[a][1]].get_parameters();
    pr(pl[a][2]) += h;
    pp[pl[a][0]][pl[a][1]].set_parameters(pr);
    return pp;
  };
  for (size_t a = 0; a < pl.size(); ++a) {
    double base = std::abs(pcs[pl[a][0]][pl[a][1]].get_parameters()(pl[a][2]));
    double h = 1e-4 * std::max(1.0, base);
    Eigen::MatrixXd sp = Vinecop(structure, bump(pcs, a, h)).scores(u, true, 1);
    Eigen::MatrixXd sm =
      Vinecop(structure, bump(pcs, a, -h)).scores(u, true, 1);
    Eigen::MatrixXd fd = (sp - sm) / (2 * h);
    EXPECT_TRUE(all_close(Ha(pl[a][0], pl[a][1])[pl[a][2]], fd, 1e-3, 1e-3))
      << "param " << a;
  }
  // threading determinism
  auto Ha3 = vc.hessian_full(u, true, 3);
  for (size_t t = 0; t < 4; ++t) {
    for (size_t e = 0; e < 4 - t; ++e) {
      for (size_t p = 0; p < Ha(t, e).size(); ++p) {
        EXPECT_TRUE(all_close(Ha(t, e)[p], Ha3(t, e)[p], 1e-12, 1e-12));
      }
    }
  }
}

// A vine with a family outside the analytic set (bb1) still gets an analytic
// joint Hessian: the second-order cascade calls the pair copula's derivative
// leaves, which fall back to finite differences inside the leaf for bb1.
TEST(VinecopDerivatives, hessian_bb1_edge_matches_brute_force)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(3);
  Eigen::VectorXd bb1_par(2);
  bb1_par << 1.0, 1.5;
  pcs[0][0] = Bicop(BicopFamily::bb1, 0, bb1_par);
  pcs[0][1] = Bicop(BicopFamily::gaussian, 0, P1(0.5));
  pcs[1][0] = Bicop(BicopFamily::gumbel, 0, P1(1.5));
  Vinecop vc(DVineStructure({ 1, 2, 3 }), pcs);
  auto u = vc.simulate(100, false, 1, { 3 });

  size_t npars = static_cast<size_t>(vc.get_npars());
  Eigen::MatrixXd H = vc.hessian(u, false, 1) * static_cast<double>(u.rows());
  EXPECT_TRUE(H.allFinite());

  std::vector<std::array<size_t, 3>> pl;
  for (size_t t = 0; t < 2; ++t) {
    for (size_t e = 0; e < 2 - t; ++e) {
      size_t np = static_cast<size_t>(pcs[t][e].get_parameters().size());
      for (size_t p = 0; p < np; ++p) {
        pl.push_back({ t, e, p });
      }
    }
  }
  auto structure = vc.get_rvine_structure();
  auto loglik = [&](const std::vector<std::vector<Bicop>>& pp) {
    return Vinecop(structure, pp).pdf(u).array().max(1e-300).log().sum();
  };
  auto bump = [&](std::vector<std::vector<Bicop>> pp, size_t a, double hh) {
    auto pr = pp[pl[a][0]][pl[a][1]].get_parameters();
    pr(pl[a][2]) += hh;
    pp[pl[a][0]][pl[a][1]].set_parameters(pr);
    return pp;
  };
  double base = loglik(pcs);
  double h = 1e-4;
  // looser tolerance: bb1's second derivatives come from the leaf-level
  // finite-difference fallback, so the cascade carries that error
  for (size_t a = 0; a < npars; ++a) {
    for (size_t bb = a; bb < npars; ++bb) {
      double fd;
      if (a == bb) {
        fd = (loglik(bump(pcs, a, h)) - 2 * base + loglik(bump(pcs, a, -h))) /
             (h * h);
      } else {
        fd = (loglik(bump(bump(pcs, a, h), bb, h)) -
              loglik(bump(bump(pcs, a, h), bb, -h)) -
              loglik(bump(bump(pcs, a, -h), bb, h)) +
              loglik(bump(bump(pcs, a, -h), bb, -h))) /
             (4 * h * h);
      }
      EXPECT_NEAR(H(a, bb), fd, 5e-2 * (1.0 + std::abs(fd)))
        << "H(" << a << ", " << bb << ")";
    }
  }
}

// scores() and hessian() reject models with nonparametric pair copulas
// (differentiating w.r.t. the interpolation grid is meaningless; mirrors
// the guard added in the perf-overhaul PR).
TEST(VinecopDerivatives, reject_nonparametric)
{
  auto data = tools_stats::simulate_uniform(200, 3);
  Vinecop vc(DVineStructure({ 1, 2, 3 }));
  vc.select(data, FitControlsVinecop({ BicopFamily::tll }));
  EXPECT_THROW(vc.scores(data, true), std::runtime_error);
  EXPECT_THROW(vc.scores(data, false), std::runtime_error);
  EXPECT_THROW(vc.hessian(data, true), std::runtime_error);
  EXPECT_THROW(vc.hessian(data, false), std::runtime_error);
  EXPECT_THROW(vc.gradient(data, false), std::runtime_error);
}

// The step-wise scores must match per-edge finite differences of the edge
// log-densities, with the edge arguments reconstructed from pdf_full().
TEST(VinecopDerivatives, stepwise_scores_match_per_edge_reference)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  auto pcs = Vinecop::make_pair_copula_store(5);
  pcs[0][0] = Bicop(BicopFamily::gaussian, 0, P1(0.6));
  pcs[0][1] = Bicop(BicopFamily::clayton, 90, P1(2.5));
  pcs[0][2] = Bicop(BicopFamily::gumbel, 180, P1(1.8));
  pcs[0][3] = Bicop(BicopFamily::joe, 270, P1(2.2));
  pcs[1][0] = Bicop(BicopFamily::frank, 0, P1(4.0));
  Eigen::VectorXd st_par(2);
  st_par << 0.5, 6.0;
  pcs[1][1] = Bicop(BicopFamily::student, 0, st_par);
  pcs[1][2] = Bicop(BicopFamily::clayton, 180, P1(1.2));
  pcs[2][0] = Bicop(BicopFamily::gumbel, 0, P1(1.4));
  pcs[2][1] = Bicop(BicopFamily::gaussian, 0, P1(0.2));
  pcs[3][0] = Bicop(BicopFamily::gaussian, 0, P1(-0.3));
  Vinecop vc(DVineStructure({ 1, 2, 3, 4, 5 }), pcs);

  auto u = vc.simulate(200, false, 1, { 7, 8, 9 });
  Eigen::MatrixXd s = vc.scores(u, true, 1);

  // reconstruct each edge's arguments from the cached h-functions
  auto full = vc.pdf_full(u, 1, true);
  auto structure = vc.get_rvine_structure();
  auto order = structure.get_order();
  size_t d = vc.get_dim();
  size_t ipar = 0;
  for (size_t t = 0; t < d - 1; ++t) {
    for (size_t e = 0; e < d - 1 - t; ++e) {
      size_t m = structure.min_array(t, e);
      bool direct = (m == structure.struct_array(t, e, true));
      Eigen::MatrixXd u_e(u.rows(), 2);
      if (t == 0) {
        u_e.col(0) = u.col(order[e] - 1);
        u_e.col(1) = u.col(order[m - 1] - 1);
      } else {
        u_e.col(0) = full.hfunc2(t - 1, e);
        u_e.col(1) =
          direct ? full.hfunc2(t - 1, m - 1) : full.hfunc1(t - 1, m - 1);
      }

      // old-style central differences of the edge log-density
      Bicop ec = vc.get_pair_copula(t, e);
      auto pars = ec.get_parameters();
      auto ub = ec.get_parameters_upper_bounds();
      auto lb = ec.get_parameters_lower_bounds();
      for (size_t p = 0; p < static_cast<size_t>(pars.size()); ++p) {
        auto pars_tmp = pars;
        double eps = 0;
        pars_tmp(p) = std::min(pars(p) + 1e-3, ub(p));
        eps += pars_tmp(p) - pars(p);
        ec.set_parameters(pars_tmp);
        Eigen::VectorXd f1 = ec.pdf(u_e).array().max(1e-20).log();
        pars_tmp(p) = std::max(pars(p) - 1e-3, lb(p));
        eps -= pars_tmp(p) - pars(p);
        ec.set_parameters(pars_tmp);
        Eigen::VectorXd f2 = ec.pdf(u_e).array().max(1e-20).log();
        ec.set_parameters(pars);
        EXPECT_TRUE(
          all_close(s.col(ipar).eval(), ((f1 - f2) / eps).eval(), 1e-3, 1e-4))
          << "tree " << t << ", edge " << e << ", parameter " << p;
        ipar++;
      }
    }
  }
  EXPECT_EQ(ipar, static_cast<size_t>(vc.get_npars()));

  EXPECT_TRUE(all_close(vc.scores(u, true, 3), s, 1e-12, 1e-12));
}

// The analytic cascade must agree with brute-force finite differences of
// the whole-vine log-likelihood on a truncated vine with rotated pair
// copulas (covers 90/270 rotations and truncation, which the R oracle
// cannot: VineCopula's RVineGrad is broken for rotated families).
TEST_F(VinecopTest, full_scores_match_brute_force_on_truncated_vine)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7, 3);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);
  auto uu = vinecop.simulate(100, false, 1, { 17 });

  Eigen::VectorXd grad =
    vinecop.scores(uu, false, 1).colwise().sum().transpose();

  auto loglik_of = [&](const std::vector<std::vector<Bicop>>& pcs) {
    Vinecop v(model_matrix, pcs);
    return v.pdf(uu).array().max(1e-300).log().sum();
  };
  size_t ipar = 0;
  for (size_t t = 0; t < 3; ++t) {
    for (size_t e = 0; e < 7 - 1 - t; ++e) {
      double h = 1e-5 * 3.0;
      auto pp = pair_copulas;
      auto pm = pair_copulas;
      pp[t][e].set_parameters(Eigen::VectorXd::Constant(1, 3.0 + h));
      pm[t][e].set_parameters(Eigen::VectorXd::Constant(1, 3.0 - h));
      double fd = (loglik_of(pp) - loglik_of(pm)) / (2 * h);
      EXPECT_NEAR(grad(ipar), fd, 1e-4 * (1.0 + std::abs(fd)))
        << "tree " << t << ", edge " << e;
      ipar++;
    }
  }
  EXPECT_EQ(ipar, static_cast<size_t>(vinecop.get_npars()));
}

// Vines with families outside the analytic set or with discrete variables
// keep working through the fallback paths.
TEST(VinecopDerivatives, fallback_paths_run)
{
  auto P1 = [](double v) { return Eigen::VectorXd::Constant(1, v).eval(); };
  // a bb1 edge is parametric, so scores/hessian still take the analytic
  // cascade — the bb1 pair copula's derivative leaves simply fall back to
  // finite differences internally
  auto pcs = Vinecop::make_pair_copula_store(3);
  Eigen::VectorXd bb1_par(2);
  bb1_par << 1.0, 1.5;
  pcs[0][0] = Bicop(BicopFamily::bb1, 0, bb1_par);
  pcs[0][1] = Bicop(BicopFamily::gaussian, 0, P1(0.5));
  pcs[1][0] = Bicop(BicopFamily::gumbel, 0, P1(1.5));
  Vinecop vc(DVineStructure({ 1, 2, 3 }), pcs);
  auto u = vc.simulate(50, false, 1, { 3 });

  Eigen::MatrixXd s_full = vc.scores(u, false, 1);
  Eigen::MatrixXd s_sw = vc.scores(u, true, 1);
  EXPECT_EQ(s_full.cols(), 4);
  EXPECT_EQ(s_sw.cols(), 4);
  EXPECT_TRUE(s_full.allFinite());
  EXPECT_TRUE(s_sw.allFinite());
  // full and step-wise gradients agree for the last tree's parameters (no
  // deeper edges depend on them), giving the fallback a value check
  EXPECT_TRUE(all_close(s_full.col(3).eval(), s_sw.col(3).eval(), 1e-3, 1e-4));

  // discrete variables keep the finite-difference paths
  Vinecop vc_disc(DVineStructure({ 1, 2, 3 }), pcs);
  vc_disc.set_var_types({ "c", "d", "c" });
  Eigen::MatrixXd u4(u.rows(), 4);
  u4.leftCols(3) = u;
  u4.col(3) = (u.col(1).array() * 0.9).matrix();
  Eigen::MatrixXd s_disc = vc_disc.scores(u4, true, 1);
  EXPECT_EQ(s_disc.cols(), 4);
  EXPECT_TRUE(s_disc.allFinite());
  Eigen::MatrixXd s_disc_full = vc_disc.scores(u4, false, 1);
  EXPECT_TRUE(s_disc_full.allFinite());
}

TEST_F(VinecopTest, aic_bic_are_correct)
{
  int d = 7;
  auto data = tools_stats::simulate_uniform(100, 7);
  Vinecop true_model(d);

  auto pair_copulas = Vinecop::make_pair_copula_store(d);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 0, Eigen::VectorXd::Constant(1, 3.0));
    }
  }
  Vinecop complex_model(model_matrix, pair_copulas);

  ASSERT_TRUE(true_model.aic(data) < complex_model.aic(data));
  ASSERT_TRUE(true_model.bic(data) < complex_model.bic(data));
  true_model.select(data);

  FitControlsVinecop controls({ BicopFamily::gaussian, BicopFamily::tll });
  complex_model.select(data);
  ASSERT_NEAR(complex_model.get_aic(), complex_model.aic(data), 1e-2);
  ASSERT_NEAR(complex_model.get_bic(), complex_model.bic(data), 1e-2);
}

TEST_F(VinecopTest, fit_parameters_is_correct)
{
  u.conservativeResize(50, 7);
  auto controls = FitControlsVinecop({ BicopFamily::clayton }, "itau");
  Vinecop vc(7);
  vc.select(u, controls);
  auto rvine_structure = vc.get_rvine_structure();

  auto pcs = vc.get_all_pair_copulas();
  for (auto& pc : pcs[0])
    pc.set_parameters(Eigen::VectorXd::Constant(1, 1));
  Vinecop vc2(rvine_structure, pcs);
  vc2.fit(u, controls);

  ASSERT_NEAR(vc2.get_loglik(), vc2.loglik(u), 1e-2);
  ASSERT_NEAR(vc2.get_aic(), vc2.aic(u), 1e-2);
  ASSERT_NEAR(vc2.get_bic(), vc2.bic(u), 1e-2);
  ASSERT_TRUE(vc.str() == vc2.str());

  Vinecop vc3(rvine_structure, pcs);
  controls.set_select_families(false);
  vc3.select(u, controls);
  ASSERT_TRUE(vc.str() == vc3.str());
}

TEST_F(VinecopTest, family_select_finds_true_rotations)
{
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);
  auto data = vinecop.simulate(2000);

  auto controls = FitControlsVinecop({ BicopFamily::clayton }, "itau");
  // controls.set_show_trace(true);
  Vinecop fit(data, model_matrix, {}, controls);

  // don't check last two trees to avoid random failures because of
  // estimation uncertainty
  auto true_rots = vinecop.get_all_rotations();
  auto fitd_rots = fit.get_all_rotations();
  true_rots.erase(true_rots.end() - 2, true_rots.end());
  fitd_rots.erase(fitd_rots.end() - 2, fitd_rots.end());
  EXPECT_EQ(true_rots, fitd_rots);
}

TEST_F(VinecopTest, family_select_returns_pcs_in_right_order)
{
  u.conservativeResize(50, 7);
  auto pair_copulas = Vinecop::make_pair_copula_store(7);
  auto par = Eigen::VectorXd::Constant(1, 3.0);
  for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
      pc = Bicop(BicopFamily::clayton, 270, par);
    }
  }
  Vinecop vinecop(model_matrix, pair_copulas);

  auto controls = FitControlsVinecop(bicop_families::itau, "itau");
  // controls.set_show_trace(true);
  Vinecop fit_struct(u, RVineStructure(), {}, controls);
  Vinecop fit_fam(u, fit_struct.get_matrix(), {}, controls);

  EXPECT_EQ(fit_struct.get_all_parameters(), fit_fam.get_all_parameters());
}

TEST_F(VinecopTest, trace_works)
{
  u.conservativeResize(10, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_show_trace(true);
  controls.set_select_threshold(true);
  controls.set_trunc_lvl(3);
  testing::internal::CaptureStdout();
  Vinecop fit(u, RVineStructure(), {}, controls);
  std::string output = testing::internal::GetCapturedStdout();
  EXPECT_TRUE(!output.empty());
}

TEST_F(VinecopTest, works_multi_threaded)
{
  u.conservativeResize(100, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_select_trunc_lvl(true);

  Vinecop fit1(u, RVineStructure(), {}, controls);
  controls.set_num_threads(2);
  Vinecop fit2(u, RVineStructure(), {}, controls);

  auto pcs = fit1.get_all_pair_copulas();
  for (auto& pc : pcs[0])
    pc.set_parameters(Eigen::VectorXd::Constant(1, 1));
  Vinecop fit3(fit1.get_rvine_structure(), pcs);
  fit3.fit(u, controls, 2);

  // check for equality in likelihood, since the pair copulas may be stored
  // in a different order when running in parallel
  EXPECT_NEAR(fit1.loglik(u), fit2.loglik(u), 1e-2);
  EXPECT_NEAR(fit1.loglik(u), fit3.loglik(u), 1e-2);

  // check if parallel evaluators have same output as single threaded ones
  EXPECT_TRUE(all_close(fit2.pdf(u, 2), fit2.pdf(u), 1e-10, 1e-10));
  EXPECT_TRUE(all_close(
    fit2.inverse_rosenblatt(u, 2), fit2.inverse_rosenblatt(u), 1e-10, 1e-10));
  EXPECT_TRUE(
    all_close(fit2.rosenblatt(u, 2), fit2.rosenblatt(u), 1e-10, 1e-10));

  // just check that it works
  fit2.simulate(2, false, 3);
  fit2.cdf(u, 100, 2);
}

// check if the same conditioned sets appear for each tree
inline size_t
get_pairs_unequal(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix1,
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix2,
  size_t trunc_lvl)
{
  std::vector<std::vector<std::vector<size_t>>> vc_sets(trunc_lvl),
    vcl_sets(trunc_lvl);
  size_t pairs_unequal = 0;
  for (size_t tree = 0; tree < trunc_lvl; ++tree) {
    vc_sets[tree].resize(trunc_lvl - tree);
    vcl_sets[tree].resize(trunc_lvl - tree);
    for (size_t edge = 0; edge < trunc_lvl - tree; ++edge) {
      vc_sets[tree][edge].resize(2);
      vc_sets[tree][edge][0] = matrix1(tree, edge);
      vc_sets[tree][edge][1] = matrix1(trunc_lvl - edge, edge);
      vcl_sets[tree][edge].resize(2);
      vcl_sets[tree][edge][0] = matrix2(tree, edge);
      vcl_sets[tree][edge][1] = matrix2(trunc_lvl - edge, edge);
    }
    for (auto s1 : vc_sets[tree]) {
      bool is_in_both = false;
      for (auto s2 : vcl_sets[tree]) {
        if (tools_stl::is_same_set(s1, s2)) {
          is_in_both = true;
        }
      }
      if (!is_in_both) {
        ++pairs_unequal;
      }
    }
  }
  return pairs_unequal;
}

TEST_F(VinecopTest, select_finds_right_structure_prim)
{
  // check whether the same structure appears if we only allow for
  // independence (pair-copula estimates differ otherwise)

  // select structure and get matrix
  Vinecop fit(7);
  fit.select(u, FitControlsVinecop({ BicopFamily::indep }));
  auto vcl_matrix = fit.get_matrix();

  // check if the same conditioned sets appear for each tree
  size_t pairs_unequal = get_pairs_unequal(vc_matrix, vcl_matrix, 6);
  EXPECT_EQ(pairs_unequal, 0);
}

TEST_F(VinecopTest, select_finds_right_structure_kruskal)
{
  // check whether the same structure appears if we only allow for
  // independence (pair-copula estimates differ otherwise)
  FitControlsVinecop controls({ BicopFamily::indep });
  EXPECT_EQ(controls.get_tree_algorithm(), "mst_prim");
  EXPECT_ANY_THROW(controls.set_tree_algorithm("foobar"));
  controls.set_tree_algorithm("mst_kruskal");

  // select structure and get matrix
  Vinecop fit(7);
  fit.select(u, controls);
  auto vcl_matrix = fit.get_matrix();

  // check if the same conditioned sets appear for each tree
  size_t pairs_unequal = get_pairs_unequal(vc_matrix, vcl_matrix, 6);
  EXPECT_EQ(pairs_unequal, 0);
}

TEST_F(VinecopTest, tree_criterion_custom_works)
{
  auto tau_fn = [](const Eigen::MatrixXd& data,
                   const Eigen::VectorXd& weights) {
    return wdm::wdm(data, "tau", weights)(0, 1);
  };

  // setter path: a custom criterion equal to Kendall's tau must recover the
  // same structure as the built-in "tau" (independence families so the
  // structure is determined solely by the criterion).
  Vinecop fit(7);
  FitControlsVinecop controls({ BicopFamily::indep });
  controls.set_tree_criterion("custom");
  controls.set_tree_criterion_function(tau_fn);
  fit.select(u, controls);
  EXPECT_EQ(get_pairs_unequal(vc_matrix, fit.get_matrix(), 6), 0);

  // FitControlsConfig construction path covers the new optional field.
  FitControlsConfig cfg;
  cfg.family_set = std::vector<BicopFamily>{ BicopFamily::indep };
  cfg.tree_criterion = "custom";
  cfg.tree_criterion_function = tau_fn;
  Vinecop fit_cfg(7);
  fit_cfg.select(u, FitControlsVinecop(cfg));
  EXPECT_EQ(get_pairs_unequal(vc_matrix, fit_cfg.get_matrix(), 6), 0);

  // "custom" without a callable must throw during selection.
  FitControlsVinecop bad({ BicopFamily::indep });
  bad.set_tree_criterion("custom");
  Vinecop fit_bad(7);
  EXPECT_ANY_THROW(fit_bad.select(u, bad));
}

TEST_F(VinecopTest, select_finds_different_structures_random)
{
  // Initialize the controls
  FitControlsVinecop controls_weighted({ BicopFamily::tll });
  controls_weighted.set_tree_algorithm("random_weighted");

  FitControlsVinecop controls_unweighted({ BicopFamily::tll });
  controls_unweighted.set_tree_algorithm("random_unweighted");

  // For reseeding the random number generator
  std::random_device rd;
  std::vector<int> seeds(20);

  // To store the unique structures
  std::set<TriangularArray<size_t>> unique_structures_weighted;
  std::set<TriangularArray<size_t>> unique_structures_unweighted;

  // To store the first RNG value after each reseeding
  std::set<uint32_t> first_rng_outputs;

  const size_t num_trials = 10;

  for (size_t i = 0; i < num_trials; ++i) {
    // Seed controls randomly for each test run
    std::generate(
      seeds.begin(), seeds.end(), [&]() { return static_cast<int>(rd()); });
    controls_weighted.set_seeds(seeds);
    controls_unweighted.set_seeds(seeds);

    // Check RNG output changes
    auto rng_sample_weighted =
      controls_weighted.get_rng()(); // Get first sample
    first_rng_outputs.insert(rng_sample_weighted);

    auto rng_sample_unweighted =
      controls_unweighted.get_rng()(); // Get first sample
    first_rng_outputs.insert(rng_sample_unweighted);

    // Select a random structure for the weighted method
    Vinecop fit_weighted(u, RVineStructure(), {}, controls_weighted);
    auto struct_array_weighted = fit_weighted.get_struct_array();
    unique_structures_weighted.insert(struct_array_weighted);

    // Select a random structure for the unweighted method
    Vinecop fit_unweighted(u, RVineStructure(), {}, controls_unweighted);
    auto struct_array_unweighted = fit_unweighted.get_struct_array();
    unique_structures_unweighted.insert(struct_array_unweighted);
  }

  // The probability that any 2 samples are the same by chance is very low
  EXPECT_EQ(first_rng_outputs.size(), num_trials);
  EXPECT_EQ(unique_structures_weighted.size(), num_trials);
  EXPECT_EQ(unique_structures_unweighted.size(), num_trials);
}

TEST_F(VinecopTest, fixed_truncation)
{
  u.conservativeResize(10, 7);
  FitControlsVinecop controls({ BicopFamily::indep });
  controls.set_trunc_lvl(2);
  // controls.set_show_trace(true);
  Vinecop fit(7);
  fit.select(u, controls);
  EXPECT_EQ(fit.get_all_pair_copulas().size(), 2);

  fit.select(u, controls);
  EXPECT_EQ(fit.get_all_pair_copulas().size(), 2);

  Vinecop fit2(u, fit.get_rvine_structure(), {}, controls);
  EXPECT_EQ(fit2.get_all_pair_copulas().size(), 2);

  Vinecop fit3(u, fit.get_rvine_structure());
  EXPECT_EQ(fit3.get_all_pair_copulas().size(), 6);

  fit3.truncate(2);
  size_t pairs_unequal =
    get_pairs_unequal(fit3.get_matrix(), fit.get_matrix(), 2);
  EXPECT_EQ(pairs_unequal, 0);
}

TEST_F(VinecopTest, sparse_threshold_selection)
{
  u.conservativeResize(20, 7);

  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_select_threshold(true);
  controls.set_threshold(NAN);
  // controls.set_show_trace(true);
  controls.set_selection_criterion("mbicv");

  Vinecop fit(7);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);

  u = tools_stats::simulate_uniform(100, 7);
  fit.select(u, controls);
}

TEST_F(VinecopTest, sparse_truncation_selection)
{
  u.conservativeResize(50, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_select_trunc_lvl(true);
  // controls.set_show_trace(true);
  u = tools_stats::simulate_uniform(100, 7);
  Vinecop fit(7);
  fit.select(u, controls);
  EXPECT_LE(fit.get_rvine_structure().get_trunc_lvl(), 6);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);
  fit.select(u, controls);
  EXPECT_LE(fit.get_rvine_structure().get_trunc_lvl(), 6);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);
}

TEST_F(VinecopTest, sparse_both_selection)
{
  u.conservativeResize(20, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  controls.set_select_trunc_lvl(true);
  controls.set_select_threshold(true);
  controls.set_tree_criterion("joe");
  controls.set_selection_criterion("mbicv");
  // controls.set_show_trace(true);
  Vinecop fit(7);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);
  fit.select(u, controls);
  EXPECT_NEAR(fit.get_loglik(), fit.loglik(u), 0.001);

  // 1d models
  Eigen::MatrixXd uu = u.col(0).matrix();
  fit = Vinecop(1);
  fit.select(uu, controls);
}

TEST_F(VinecopTest, partial_selection)
{
  u.conservativeResize(20, 7);
  FitControlsVinecop controls(bicop_families::itau, "itau");
  // controls.set_show_trace(true);
  auto fixed = CVineStructure(std::vector<size_t>{ 5, 4, 7, 1, 3, 6, 2 });
  fixed.truncate(1);
  Vinecop fit(fixed);
  fit.select(u, controls);

  // a C-vine with root node 2 contains 6 edges with vertex 2.
  auto rvm = fit.get_rvine_structure().get_matrix();
  int count2 = 0;
  for (int i = 0; i < 6; i++) {
    //  diagonal element         base element
    if ((rvm(6 - i, i) == 2) || (rvm(0, i) == 2))
      count2++;
  }
  EXPECT_EQ(count2, 6);
}

TEST_F(VinecopTest, tawn_flipping)
{
  FitControlsVinecop controls({ BicopFamily::tawn });
  Vinecop fit1(7);
  fit1.select(u, controls);
  Vinecop fit2(fit1.get_rvine_structure());
  fit2.select(u, controls);

  for (size_t tree = 0; tree < fit1.get_trunc_lvl(); ++tree) {
    for (size_t edge = 0; edge < 6 - tree; ++edge) {
      auto pc1 = fit1.get_pair_copula(tree, edge);
      auto pc2 = fit2.get_pair_copula(tree, edge);
      ASSERT_EQ(pc1.get_family(), pc2.get_family());
      ASSERT_EQ(pc1.get_rotation(), pc2.get_rotation());
      ASSERT_TRUE(
        all_close(pc1.get_parameters(), pc2.get_parameters(), 0.01, 0.01));
    }
  }
}

}
