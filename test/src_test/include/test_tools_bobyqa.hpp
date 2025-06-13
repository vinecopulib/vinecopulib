// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <cmath>
#include <vinecopulib.hpp>
#include <vinecopulib/misc/tools_bobyqa.hpp>

namespace test_tools_bobyqa {

using namespace vinecopulib;

TEST(test_tools_bobyqa, const_function)
{

  auto f = [](long /*n*/, const double* /*x*/) -> double { return 0.0; };

  const long variables_count = 2;
  const long number_of_interpolation_conditions = variables_count + 2;

  Eigen::VectorXd lb(2);
  Eigen::VectorXd ub(2);
  Eigen::VectorXd x(2);

  lb << -1.0, -1.0;
  ub << 1.0, 1.0;
  x << 0.0, 0.0;

  const double initial_trust_region_radius = 1e-3;
  const double final_trust_region_radius = 1e3;
  const long max_function_calls_count = 5;

  auto result = tools_bobyqa::bobyqa(f,
                                     variables_count,
                                     number_of_interpolation_conditions,
                                     x,
                                     lb,
                                     ub,
                                     initial_trust_region_radius,
                                     final_trust_region_radius,
                                     max_function_calls_count);

  ASSERT_TRUE(fabs(result.second) < 1e-5);
  ASSERT_TRUE(fabs(result.first(0)) < 1e-5);
  ASSERT_TRUE(fabs(result.first(1)) < 1e-5);
}

TEST(test_tools_bobyqa, complex_quadratic_function)
{

  auto f = [](long /*n*/, const double* x) -> double {
    return -4 * x[0] * x[1] + 5 * x[0] * x[0] + 8 * x[1] * x[1] +
           16 * sqrt(5.0) * x[0] + 8 * sqrt(5.0) * x[1] - 44.0;
  };

  const long variables_count = 2;
  const long number_of_interpolation_conditions =
    (variables_count + 1) * (variables_count + 2) / 2;

  Eigen::VectorXd lb(2);
  Eigen::VectorXd ub(2);
  Eigen::VectorXd x(2);

  lb << -4.0, -3.0;
  ub << 5.0, 5.0;
  x << 0.0, -sqrt(5.0);

  const double initial_trust_region_radius = 1e-3;
  const double final_trust_region_radius = 1e3;
  const long max_function_calls_count = 22;

  auto result = tools_bobyqa::bobyqa(f,
                                     variables_count,
                                     number_of_interpolation_conditions,
                                     x,
                                     lb,
                                     ub,
                                     initial_trust_region_radius,
                                     final_trust_region_radius,
                                     max_function_calls_count);

  ASSERT_TRUE(fabs(result.second + 142.99689) < 1e-5);
  ASSERT_TRUE(fabs(result.first(0) + 4.0) < 1e-5);
  ASSERT_TRUE(fabs(result.first(1) + 2.11803) < 1e-5);
}

TEST(test_tools_bobyqa, quadratic_function_with_jump)
{

  auto f = [](long /*n*/, const double* x) -> double {
    return x[1] < 0.5 ? x[0] * x[0] + x[1] * x[1]
                      : x[0] * x[0] + x[1] * x[1] + 10;
  };

  const long variables_count = 2;
  const long number_of_interpolation_conditions =
    (variables_count + 1) * (variables_count + 2) / 2;

  Eigen::VectorXd lb(2);
  Eigen::VectorXd ub(2);
  Eigen::VectorXd x(2);

  lb << -1.0, -1.0;
  ub << 1.0, 1.0;
  x << 0.5, 0.5;

  const double initial_trust_region_radius = 1e-3;
  const double final_trust_region_radius = 1e3;
  const long max_function_calls_count = 22;

  auto result = tools_bobyqa::bobyqa(f,
                                     variables_count,
                                     number_of_interpolation_conditions,
                                     x,
                                     lb,
                                     ub,
                                     initial_trust_region_radius,
                                     final_trust_region_radius,
                                     max_function_calls_count);

  ASSERT_TRUE(fabs(result.second - 0.49700) < 1e-5);
  ASSERT_TRUE(fabs(result.first(0) - 0.49899) < 1e-5);
  ASSERT_TRUE(fabs(result.first(1) - 0.49800) < 1e-5);
}
}
