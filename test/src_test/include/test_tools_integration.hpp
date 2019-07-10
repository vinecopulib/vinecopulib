// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/bicop/class.hpp>
#include <vinecopulib/misc/tools_integration.hpp>
#include <vinecopulib/misc/tools_stats.hpp>

namespace test_tools_integration {

using namespace vinecopulib;

TEST(test_tools_integration, integrate_uniform)
{
  auto myfun = [](double x, double i = 1) { return std::pow(x, i); };
  auto moments_uniform = [](double a = 0, double b = 1, double i = 1) {
    double output = 0.0;
    for (size_t j = 0; j <= static_cast<size_t>(i); j++) {
      double j_dbl = static_cast<double>(j);
      output += std::pow(a, j_dbl) * std::pow(b, i - j_dbl);
    }
    return output / (i + 1);
  };
  EXPECT_ANY_THROW(tools_integration::integrate(myfun, 0, 1, 42));

  for (size_t i = 1; i < 5; i++) {
    double i_dbl = static_cast<double>(i);
    auto f = [i_dbl, myfun](double x) { return myfun(x, i_dbl); };
    EXPECT_NEAR(
      tools_integration::integrate(f), moments_uniform(0, 1, i_dbl), 1e-8);
    EXPECT_NEAR(tools_integration::integrate(f, 3, 5) / 2,
                moments_uniform(3, 5, i_dbl),
                1e-8);
  }
}

TEST(test_tools_integration, integrate_hfunc)
{
  auto cop = Bicop(BicopFamily::clayton, 0, Eigen::VectorXd::Constant(1, 10));
  auto myfun = [cop](const double u1, const double u2) {
    Eigen::MatrixXd input(1, 2);
    input << u1, u2;
    auto output = cop.pdf(input);
    return output(0);
  };
  auto u = tools_stats::simulate_uniform(20, 2, { 1 });
  for (size_t i = 0; i < 20; i++) {
    auto f = [i, u, myfun](double x) { return myfun(x, u(i, 1)); };
    auto truth = cop.hfunc2(u.row(i));
    EXPECT_NEAR(tools_integration::integrate(f, 0, 1, 50), 1.0, 1e-5);
    EXPECT_NEAR(
      tools_integration::integrate(f, 0, u(i, 0), 50), truth(0, 0), 1e-5);
  }
}
}
