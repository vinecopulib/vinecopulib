// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "gtest/gtest.h"
#include <vinecopulib/misc/tools_integration.hpp>

namespace test_tools_integration {

using namespace vinecopulib;

TEST(test_tools_integration, integrate_is_correct)
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
    EXPECT_NEAR(tools_integration::integrate(f, 3, 5),
                moments_uniform(3, 5, i_dbl),
                1e-8);
  }
}
}
