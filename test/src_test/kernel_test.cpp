// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include "include/kernel_test.hpp"
#include <vinecopulib.hpp>

TrafokernelTest::TrafokernelTest()
  : bicop_(Bicop(vinecopulib::BicopFamily::tll, 0))
  , controls()
  , u(tools_stats::simulate_uniform(20, 2, true, { 1 }))
{
  controls.family_set = std::vector<BicopFamily>{ vinecopulib::BicopFamily::tll };
  controls.nonparametric_method = GetParam();

}
