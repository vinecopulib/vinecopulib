// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <cmath>
#include <vinecopulib.hpp>
#include <vinecopulib/misc/tools_optimization.hpp>
#include <vinecopulib/misc/tools_transforms.hpp>

namespace test_tools_optimization {

using namespace vinecopulib;
using namespace vinecopulib::tools_transforms;
using namespace vinecopulib::tools_optimization;

// ---------------------------------------------------------------------------
// Parameter transforms
// ---------------------------------------------------------------------------

TEST(test_tools_transforms, roundtrip_and_derivative)
{
  const double eta_max = 30.0;
  std::vector<Transform> transforms;
  { // interval-bounded -> logistic
    Transform t;
    t.kind = TransformKind::Logistic;
    t.a = -1.0;
    t.b = 1.0;
    t.eta_max = eta_max;
    transforms.push_back(t);
  }
  { // lower-bounded -> softplus
    Transform t;
    t.kind = TransformKind::SoftplusLower;
    t.a = 1.0;
    t.eta_max = eta_max;
    transforms.push_back(t);
  }
  { // upper-bounded -> mirrored softplus
    Transform t;
    t.kind = TransformKind::SoftplusUpper;
    t.b = 5.0;
    t.eta_max = eta_max;
    transforms.push_back(t);
  }
  { // unbounded -> identity
    Transform t;
    t.kind = TransformKind::Identity;
    t.eta_max = eta_max;
    transforms.push_back(t);
  }

  for (const auto& t : transforms) {
    for (double eta : { -8.0, -3.0, -0.5, 0.0, 0.5, 3.0, 8.0 }) {
      double theta = t.to_theta(eta);
      // theta stays within the bounds
      EXPECT_GE(theta, t.a - 1e-12);
      EXPECT_LE(theta, t.b + 1e-12);
      // eta -> theta -> eta round-trips (within the eta clamp range)
      EXPECT_NEAR(t.to_eta(theta), eta, 1e-8);
      // analytic derivative matches central finite differences
      double fd = (t.to_theta(eta + 1e-6) - t.to_theta(eta - 1e-6)) / 2e-6;
      EXPECT_NEAR(t.dtheta_deta(eta), fd, 1e-4);
    }
  }
}

TEST(test_tools_transforms, boundary_values_are_finite)
{
  Transform t;
  t.kind = TransformKind::Logistic;
  t.a = -1.0;
  t.b = 1.0;
  // inverse-transforming a value on the bound gives a large but finite eta
  EXPECT_TRUE(std::isfinite(t.to_eta(t.a)));
  EXPECT_TRUE(std::isfinite(t.to_eta(t.b)));
  // forward-transforming an extreme eta stays inside the bounds
  EXPECT_LE(t.to_theta(1e6), t.b);
  EXPECT_GE(t.to_theta(-1e6), t.a);
}

TEST(test_tools_transforms, from_bounds_infers_kinds)
{
  Eigen::VectorXd lb(3), ub(3);
  lb << -1.0, 1e-10, -std::numeric_limits<double>::infinity();
  ub << 1.0, 28.0, std::numeric_limits<double>::infinity();
  auto pt = ParameterTransform::from_bounds(lb, ub);
  ASSERT_EQ(pt.size(), 3u);

  // an unbounded coordinate uses the identity map
  Eigen::VectorXd eta(3);
  eta << 0.3, 2.0, 7.5;
  Eigen::VectorXd theta = pt.to_theta(eta);
  EXPECT_NEAR(theta(2), 7.5, 1e-12);

  // interior round-trip on the bounded coordinates
  Eigen::VectorXd theta_in(3);
  theta_in << 0.4, 3.0, -2.0;
  EXPECT_TRUE(
    (pt.to_theta(pt.to_eta(theta_in)) - theta_in).cwiseAbs().maxCoeff() < 1e-8);
}

// ---------------------------------------------------------------------------
// BFGS optimizer
// ---------------------------------------------------------------------------

TEST(test_tools_optimization, constant_function_returns_start)
{
  Objective f = [](const Eigen::VectorXd&, Eigen::VectorXd& grad) {
    if (grad.size() > 0)
      grad.setZero();
    return 0.0;
  };
  Eigen::VectorXd lb(2), ub(2), x0(2);
  lb << -1.0, -1.0;
  ub << 1.0, 1.0;
  x0 << 0.2, -0.3;

  Optimizer opt;
  auto res = opt.optimize(x0, lb, ub, f);
  EXPECT_NEAR(opt.get_objective_max(), 0.0, 1e-10);
  EXPECT_NEAR(res(0), 0.2, 1e-6);
  EXPECT_NEAR(res(1), -0.3, 1e-6);
}

TEST(test_tools_optimization, quadratic_interior_optimum)
{
  for (int n : { 1, 2, 3 }) {
    Eigen::VectorXd c = Eigen::VectorXd::LinSpaced(n, -1.0, 2.0);
    // maximize f(x) = -||x - c||^2, unique maximizer at x = c, f = 0
    Objective f = [&c](const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
      if (grad.size() > 0)
        grad = -2.0 * (x - c);
      return -(x - c).squaredNorm();
    };
    Eigen::VectorXd lb = c.array() - 10.0;
    Eigen::VectorXd ub = c.array() + 10.0;
    Eigen::VectorXd x0 = c.array() + 1.5;

    Optimizer opt;
    auto res = opt.optimize(x0, lb, ub, f);
    EXPECT_TRUE((res - c).cwiseAbs().maxCoeff() < 1e-5);
    EXPECT_NEAR(opt.get_objective_max(), 0.0, 1e-8);
  }
}

TEST(test_tools_optimization, optimum_on_boundary)
{
  // maximize f(x) = -(x - 2)^2 over [-1, 1]: constrained maximizer at x = 1
  Objective f = [](const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    if (grad.size() > 0)
      grad(0) = -2.0 * (x(0) - 2.0);
    return -(x(0) - 2.0) * (x(0) - 2.0);
  };
  Eigen::VectorXd lb(1), ub(1), x0(1);
  lb << -1.0;
  ub << 1.0;
  x0 << 0.0;

  Optimizer opt;
  auto res = opt.optimize(x0, lb, ub, f);
  EXPECT_NEAR(res(0), 1.0, 1e-3);
  EXPECT_NEAR(opt.get_objective_max(), -1.0, 1e-3);
}

TEST(test_tools_optimization, rosenbrock)
{
  // maximize the negative Rosenbrock function; maximizer at (1, 1), f = 0
  Objective f = [](const Eigen::VectorXd& p, Eigen::VectorXd& grad) {
    double x = p(0), y = p(1);
    double val = -(100.0 * (y - x * x) * (y - x * x) + (1.0 - x) * (1.0 - x));
    if (grad.size() > 0) {
      grad(0) = 400.0 * x * (y - x * x) + 2.0 * (1.0 - x);
      grad(1) = -200.0 * (y - x * x);
    }
    return val;
  };
  Eigen::VectorXd lb(2), ub(2), x0(2);
  lb << -5.0, -5.0;
  ub << 5.0, 5.0;
  x0 << -1.2, 1.0;

  Optimizer opt;
  auto res = opt.optimize(x0, lb, ub, f);
  EXPECT_NEAR(res(0), 1.0, 1e-3);
  EXPECT_NEAR(res(1), 1.0, 1e-3);
  EXPECT_GT(opt.get_objective_max(), -1e-6);
}

TEST(test_tools_optimization, size_checks_throw)
{
  Objective f = [](const Eigen::VectorXd&, Eigen::VectorXd&) { return 0.0; };
  Eigen::VectorXd x0(2), lb(2), ub(1);
  x0 << 0.0, 0.0;
  lb << -1.0, -1.0;
  ub << 1.0;
  Optimizer opt;
  EXPECT_ANY_THROW(opt.optimize(x0, lb, ub, f));
}

} // namespace test_tools_optimization
